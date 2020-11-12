#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>

/* Define DEBUG by compiler's command line (-DDEBUG=1) */
#ifndef DEBUG
	#define DEBUG 1
#endif

void bubblesort(int *array, const int SIZE);
void combine(int *src_a, int len_a, int *src_b, int len_b, int *dst, int length);
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
		printf("P%d: Not enough memory to allocate array of length %d\n", mpi_rank, SLICE_LEN);
		MPI_Finalize();
		exit(1);
	}

	/* Communication buffers. Avoid reallocation */
	int *right_val = NULL;
	int *combined = NULL;

#if DEBUG == 1
	printf("P%d: Populating array\n", mpi_rank);
#endif

	/* Populate array in decreasing order to execute worst case */
	for(int i = 0; i < SLICE_LEN; i++)
		values[i] = ARRAY_LEN - ARRAY_LEN/mpi_size*mpi_rank - i;

#if DEBUG == 1
	print_array(values, SLICE_LEN);
#else
	double then = MPI_Wtime();
#endif

	while(true){
		/* Sort the local array */
		bubblesort(values, SLICE_LEN);

	#if DEBUG == 1
		printf("\nP%d: sorted array\n", mpi_rank);
		print_array(values, SLICE_LEN);
	#endif

		/* Send last (biggest) element to right */
		if(mpi_rank != mpi_size - 1){
		#if DEBUG == 1
			printf("P%d: Sending value %d to right\n", mpi_rank, values[SLICE_LEN - 1]);
		#endif
			MPI_Send(&values[SLICE_LEN - 1], 1, MPI_INT, mpi_rank + 1, 0, MPI_COMM_WORLD);
		}

		bool sorted = true;
		if(mpi_rank != 0){
			/* Receive biggest element from left */
			int biggest;
			MPI_Status mpi_status;
			MPI_Recv(&biggest, 1, MPI_INT, mpi_rank - 1, 0, MPI_COMM_WORLD, &mpi_status);

			/* Check if left's biggest number is smaller than this slice smallest number */
			sorted = values[0] > biggest;

		#if DEBUG == 1
			printf("P%d: %d > %d = %d\n", mpi_rank, values[0], biggest, sorted);
		#endif
		}

		/* Let every process know */
		bool finished = true;
		for(int i = 1; i < mpi_size; i++){
			bool aux;
			/* If the current rank is the root of broadcast, send if it is sorted */
			if(i == mpi_rank)
				aux = sorted;
			MPI_Bcast(&aux, 1, MPI_C_BOOL, i, MPI_COMM_WORLD);

			/* If some of the processes is not sorted, don't continue the broadcast */
			if(!aux){
				finished = false;
			#if DEBUG == 1
				printf("P%d: process %d is not sorted\n", mpi_rank, i);
			#endif
				break;
			}
		}

		/* If everybody is sorted, finish execution. */
		if(finished)
			break;

		/* Send the lesser values to left */
		if(mpi_rank != 0)
			MPI_Send(values, SLICE_LEN / 2, MPI_INT, mpi_rank - 1, 0, MPI_COMM_WORLD);

		if(mpi_rank != mpi_size - 1){
			MPI_Status mpi_status;
			/* Only check for message size and allocate buffer if not allocated before */
			static int mpi_count;
			if(right_val == NULL){
				/* Probe incoming message */
				MPI_Probe(mpi_rank + 1, 0, MPI_COMM_WORLD, &mpi_status);

				/* Check for its length */
				MPI_Get_count(&mpi_status, MPI_INT, &mpi_count);

				/* Create an array for lesser values of right */
				right_val = malloc(mpi_count*sizeof(int));
				if(right_val == NULL){
					printf("P%d: Not enough memory to allocate array of length %d\n", mpi_rank, mpi_count);
					MPI_Finalize();
					free(values);
					values = NULL;
					exit(1);
				}
			}

			/* Receive the values */
			MPI_Recv(right_val, mpi_count, MPI_INT, mpi_rank + 1, 0, MPI_COMM_WORLD, &mpi_status);

			/* Only allocate if not allocated before */
			if(combined == NULL){
				/* Create aux array for combined values */
				combined = malloc((SLICE_LEN/2 + SLICE_LEN%2 + mpi_count)*sizeof(int));
				if(combined == NULL){
					printf("P%d: Not enough memory to allocate array of length %d\n", mpi_rank, SLICE_LEN/2 + SLICE_LEN%2 + mpi_count);
					MPI_Finalize();
					free(right_val);
					right_val = NULL;
					free(values);
					values = NULL;
					exit(1);
				}
			}

			/* Combine lesser values from right with greater values from here */
			combine(&values[SLICE_LEN/2], SLICE_LEN/2 + SLICE_LEN%2, right_val, mpi_count, combined, SLICE_LEN/2 + SLICE_LEN%2 + mpi_count);

			/* Put the lesser values back to this process */
			memcpy(&values[SLICE_LEN/2], combined, (SLICE_LEN/2 + SLICE_LEN%2)*sizeof(int));

			/* Put the greater values back to the right values buffer */
			memcpy(right_val, &combined[SLICE_LEN/2 + SLICE_LEN%2], (mpi_count)*sizeof(int));

			/* Send back the values to right */
			MPI_Send(right_val, mpi_count, MPI_INT, mpi_rank + 1, 0, MPI_COMM_WORLD);
		}

		/* Receive back the values from the left */
		if(mpi_rank != 0){
			MPI_Status mpi_status;
			MPI_Recv(values, SLICE_LEN / 2, MPI_INT, mpi_rank - 1, 0, MPI_COMM_WORLD, &mpi_status);
		}

	}

	#if DEBUG == 1
		printf("P%d: All process sorted\n", mpi_rank);
		printf("\nP%d: sorted array\n", mpi_rank);
		print_array(values, SLICE_LEN);
	#else
		double now = MPI_Wtime();
	#endif

	free(values);
	values = NULL;

	if(right_val != NULL){
		free(right_val);
		right_val = NULL;
	}

	if(combined != NULL){
		free(combined);
		combined = NULL;
	}

#if DEBUG == 0
	if(mpi_rank == 0)
		printf("Array sorted in %f\n", now - then);
#endif

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

void combine(int *src_a, int len_a, int *src_b, int len_b, int *dst, int length)
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

#if DEBUG == 1
void print_array(int *array, int len)
{
	printf("Array: ");
	for(int i = 0; i < len; i++)
		printf("%d ", array[i]);
	printf("\n");
}
#endif


