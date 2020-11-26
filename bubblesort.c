#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>

#include <unistd.h>

/* Define DEBUG by compiler's command line (-DDEBUG=1) */
#ifndef DEBUG
	#define DEBUG 0
#endif

/* If true, only use bubblesort in first iteration, then use merge */
#ifndef PREFER_MERGE
	#define PREFER_MERGE 1
#endif

/* Try to skip unnecessary converge phase at cost of more broadcasts */
#ifndef SKIP_CONVERGE
	#define SKIP_CONVERGE 1
#endif

/* Size of the converge phase slice. Default is 2 (length will be process' array size divided by 2) */
#ifndef CONV_DIV
	#define CONV_DIV 2
#endif

void bubblesort(int *array, const int SIZE);
void merge_2(int *src_a, int len_a, int *src_b, int len_b, int *dst, int length);
#if PREFER_MERGE == 1
void merge_3(int *src_a, int len_a, int *src_b, int len_b, int *src_c, int len_c, int *dst, int length);
#endif
#if DEBUG == 1
void print_array(int *array, int len);
#endif

int main(int argc, char *argv[])
{
	/* Initialize MPI */
	MPI_Init(&argc, &argv);

	/* Get MPI processor identification */
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	/* Get MPI processor number */
	int mpi_size;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

	/* Default with 1.000.000 entries */
	/* Or configure with compiler's command line (-DN=N) */
#ifndef N
	#if DEBUG == 1
		/* DEBUG array with 40 entries */
		const unsigned ARRAY_LEN = 40;
	#else
		const unsigned ARRAY_LEN = 1000000;
	#endif
#else
	const unsigned ARRAY_LEN = N;
#endif

	/* Size of each slice. Last slice can be bigger due to an odd array length */
	const unsigned SLICE_LEN = ARRAY_LEN / mpi_size + (mpi_rank == mpi_size - 1 ? ARRAY_LEN % mpi_size : 0);

	/* Allocate memory for the slice */
	int *values = malloc(SLICE_LEN*sizeof(int));
	if(values == NULL){
		printf("P%d: Not enough memory to allocate array of length %lu\n", mpi_rank, SLICE_LEN*sizeof(int));
		MPI_Finalize();
		exit(1);
	}

	/* Allocate memory for sorting status of all process */
	unsigned char *sorted = malloc(mpi_size*sizeof(unsigned char));
	if(sorted == NULL){
		printf("P%d: Not enough memory to allocate array of length %lu\n", mpi_rank, mpi_size*sizeof(unsigned char));
		free(values);
		values = NULL;
		MPI_Finalize();
		exit(1);
	}

	/* Alocate memory for combining two sorted vectors. Size of worst case */
	int *combined = malloc((ARRAY_LEN/mpi_size + ARRAY_LEN % mpi_size)*sizeof(int));
	if(combined == NULL){
		printf("P%d: Not enough memory to allocate array of length %lu\n", mpi_rank, (ARRAY_LEN/mpi_size + ARRAY_LEN % mpi_size)*sizeof(int));
		free(values);
		values = NULL;
		free(sorted);
		sorted = NULL;
		MPI_Finalize();
		exit(1);
	}

	/* Communication buffer. Size of the worst case */
	int *right_val = malloc((ARRAY_LEN/mpi_size + ARRAY_LEN % mpi_size)/CONV_DIV*sizeof(int));
	if(right_val == NULL){
		printf("P%d: Not enough memory to allocate array of length %lu\n", mpi_rank, (ARRAY_LEN/mpi_size + ARRAY_LEN%mpi_size)/CONV_DIV*sizeof(int));
		free(values);
		values = NULL;
		free(sorted);
		sorted = NULL;
		free(combined);
		combined = NULL;
		MPI_Finalize();
		exit(1);
	}

#if DEBUG == 1
	printf("P%d: Populating array\n", mpi_rank);
	MPI_Barrier(MPI_COMM_WORLD);
	usleep(100);
#endif

	/* Populate array in decreasing order to execute worst case */
	for(int i = 0; i < SLICE_LEN; i++)
		values[i] = ARRAY_LEN - ARRAY_LEN/mpi_size*mpi_rank - i;

#if DEBUG == 1
	printf("\nP%d: Array is -> ", mpi_rank);
	print_array(values, SLICE_LEN);
	printf("\n");
	MPI_Barrier(MPI_COMM_WORLD);
	usleep(100);
#else
	double then = MPI_Wtime();
#endif

	while(true){
		/* Length of received message. Will be set on first iteration */
		static int mpi_count = 0;

	#if PREFER_MERGE == 1
		static bool first_sort = true;
		if(first_sort){
	#endif
			/* Sort the local array */
			bubblesort(values, SLICE_LEN);
	#if PREFER_MERGE == 1
			first_sort = false;
		} else {
			/* The array is partially sorted. Only combine the values */
			merge_3(values, SLICE_LEN/CONV_DIV, &values[SLICE_LEN/CONV_DIV], SLICE_LEN - SLICE_LEN/CONV_DIV*2, &values[SLICE_LEN - SLICE_LEN/CONV_DIV], SLICE_LEN/CONV_DIV, combined, SLICE_LEN);

			/* Now put the values back to the main array */
			memcpy(values, combined, SLICE_LEN*sizeof(int));
		}
	#endif

	#if DEBUG == 1
		printf("\nP%d: sorted local array -> ", mpi_rank);
		print_array(values, SLICE_LEN);
		printf("\n");
		MPI_Barrier(MPI_COMM_WORLD);
		usleep(100);
	#endif

		/* Send last (biggest) element to right */
		if(mpi_rank != mpi_size - 1){
		#if DEBUG == 1
			printf("P%d: Sending value %d to right\n", mpi_rank, values[SLICE_LEN - 1]);
		#endif
			MPI_Send(&values[SLICE_LEN - 1], 1, MPI_INT, mpi_rank + 1, 0, MPI_COMM_WORLD);
		}
	
	#if DEBUG == 1
		MPI_Barrier(MPI_COMM_WORLD);
		usleep(100);
	#endif

		sorted[mpi_rank] = true;
		if(mpi_rank != 0){
			/* Receive biggest element from left */
			int biggest;
			MPI_Status mpi_status;
			MPI_Recv(&biggest, 1, MPI_INT, mpi_rank - 1, 0, MPI_COMM_WORLD, &mpi_status);

			/* Check if left's biggest number is smaller than this slice smallest number */
			sorted[mpi_rank] = values[0] > biggest;

		#if DEBUG == 1
			printf("P%d: %d > %d = %d\n", mpi_rank, values[0], biggest, sorted[mpi_rank]);
		#endif
		}

	#if DEBUG == 1
		MPI_Barrier(MPI_COMM_WORLD);
		usleep(100);
	#endif

		/* Let every process know */
		bool finished = true;
		for(int i = 1; i < mpi_size; i++){
			/* If the current rank is the root of broadcast, send if it is sorted */
			MPI_Bcast(&sorted[i], 1, MPI_UNSIGNED_CHAR, i, MPI_COMM_WORLD);
			finished = finished && sorted[i];
		#if SKIP_CONVERGE == 0
			/* Don't continue with broadcast if one of the arrays is not sorted */
			if(!finished)
				break;
		#endif
		}

	#if DEBUG == 1
		MPI_Barrier(MPI_COMM_WORLD);
		usleep(100);
	#endif

		/* If everybody is sorted, finish execution. */
		if(finished)
			break;

		/* Send the lesser values to left */
		if(
			mpi_rank != 0 
		#if SKIP_CONVERGE == 1
			/* OPTIMIZE: Only send if not sorted with my left neighbor */
			&& !sorted[mpi_rank]
		#endif
		){
			MPI_Send(values, SLICE_LEN/CONV_DIV, MPI_INT, mpi_rank - 1, 0, MPI_COMM_WORLD);
		#if DEBUG == 1
			printf("P%d: Sending values to left\n", mpi_rank);
		#endif
		}

	#if DEBUG == 1
		MPI_Barrier(MPI_COMM_WORLD);
		usleep(100);
	#endif

		if(
			mpi_rank != mpi_size - 1 
		#if SKIP_CONVERGE == 1
			/* OPTIMIZE: Only converge if not sorted with right neighbor */
			&& !sorted[mpi_rank + 1]
		#endif
		){
			MPI_Status mpi_status;
			/* Only check for message size if not received before */
			if(mpi_count == 0){
				/* Probe incoming message */
				MPI_Probe(mpi_rank + 1, 0, MPI_COMM_WORLD, &mpi_status);

				/* Check for its length */
				MPI_Get_count(&mpi_status, MPI_INT, &mpi_count);
			}

			/* Receive the values */
			MPI_Recv(right_val, mpi_count, MPI_INT, mpi_rank + 1, 0, MPI_COMM_WORLD, &mpi_status);

			/* Combine lesser values from right with greater values from here */
			merge_2(&values[SLICE_LEN - SLICE_LEN/CONV_DIV], SLICE_LEN/CONV_DIV, right_val, mpi_count, combined, SLICE_LEN/2 + SLICE_LEN%2 + mpi_count);

			/* Put the lesser values back to this process */
			memcpy(&values[SLICE_LEN - SLICE_LEN/CONV_DIV], combined, (SLICE_LEN/CONV_DIV)*sizeof(int));

		#if DEBUG == 1
			printf("P%d: Combined array (sent) -> ", mpi_rank);
			print_array(values, SLICE_LEN);
			printf("\n");
		#endif

			/* Send back the values to right */
			MPI_Send(&combined[SLICE_LEN/CONV_DIV], mpi_count, MPI_INT, mpi_rank + 1, 0, MPI_COMM_WORLD);
		}

	#if DEBUG == 1
		MPI_Barrier(MPI_COMM_WORLD);
		usleep(100);
	#endif

		/* Receive back the values from the left */
		if(
			mpi_rank != 0
		#if SKIP_CONVERGE == 1
			/* OPTIMIZE: Only send if not sorted with my left neighbor */
			&& !sorted[mpi_rank]
		#endif
		){
			MPI_Status mpi_status;
			MPI_Recv(values, SLICE_LEN / CONV_DIV, MPI_INT, mpi_rank - 1, 0, MPI_COMM_WORLD, &mpi_status);
		#if DEBUG == 1
			printf("P%d: Combined array (received) -> ", mpi_rank);
			print_array(values, SLICE_LEN);
			printf("\n");
		#endif
		}

	#if DEBUG == 1
		MPI_Barrier(MPI_COMM_WORLD);
		usleep(100);
	#endif

	}

	#if DEBUG == 1
		if(mpi_rank == 0){
			printf("All process sorted\n");
			printf("\nSorted array\n");
		}
		MPI_Barrier(MPI_COMM_WORLD);
		usleep(100);
		for(int i = 0; i < mpi_size; i++){
			if(i == mpi_rank){
				print_array(values, SLICE_LEN);
				fflush(stdout);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			usleep(100);
		}
		if(mpi_rank == 0)
			printf("\n");
	#else
		double now = MPI_Wtime();
		if(mpi_rank == 0)
			printf("Array sorted in %f\n", now - then);
	#endif

	free(values);
	values = NULL;
	free(sorted);
	sorted = NULL;
	free(combined);
	combined = NULL;
	free(right_val);
	right_val = NULL;

	MPI_Finalize();

}

void bubblesort(int *array, const int SIZE)
{
	/* Sort the array: bubblesort */
	bool swapped = true;
	for(int i = 0; swapped && i < SIZE; i++){
		swapped = false;
		for(int j = i + 1; j < SIZE; j++){
			if(array[j] < array[i]){
				int swap = array[i];
				array[i] = array[j];
				array[j] = swap;
				swapped = true;
			}
		}
	}
}

void merge_2(int *src_a, int len_a, int *src_b, int len_b, int *dst, int length)
{
	int ia = 0;
	int ib = 0;
	for(int i = 0; i < length; i++){
		if((ia < len_a && src_a[ia] <= src_b[ib]) || ib == len_b)
			dst[i] = src_a[ia++];
		else
			dst[i] = src_b[ib++];
	}
}

#if PREFER_MERGE == 1
void merge_3(int *src_a, int len_a, int *src_b, int len_b, int *src_c, int len_c, int *dst, int length)
{
	int ia = 0;
	int ib = 0;
	int ic = 0;
	for(int i = 0; i < length; i++){
		if(ia < len_a && (ib == len_b || src_a[ia] <= src_b[ib]) && (ic == len_c || src_a[ia] <= src_c[ic]))
			dst[i] = src_a[ia++];
		else if(ib < len_b && (ia == len_a || src_b[ib] <= src_a[ia]) && (ic == len_c || src_b[ib] <= src_c[ic]))
			dst[i] = src_b[ib++];
		else
			dst[i] = src_c[ic++];
	}
}
#endif

#if DEBUG == 1
void print_array(int *array, int len)
{
	for(int i = 0; i < len; i++)
		printf("%d ", array[i]);
}
#endif


