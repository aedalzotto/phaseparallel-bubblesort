# phaseparallel-bubblesort
Distributed sort using MPI following the Phase Parallel model.

## Building

Compile with standard array size (1.000.000). It is highly recommended to 
optimize the code with -O3:
```
mpicc bubblesort.c -o bubblesort -std=c99 -O3
```

To choose different array size (N), define in compiler's command line
```
mpicc bubblesort.c -o bubblesort -std=c99 -O3 -DN=500000
```

To compile in DEBUG mode, also define in compiler's command line. The default 
array size for this mode is 40.
```
mpicc bubblesort.c -o bubblesort_debug -std=c99 -O3 -DDEBUG=1
```

To change the length of the converge phase, change its size divider by using:
```
mpicc bubblesort.c -o bubblesort -std=c99 -DCONV_DIV=4
```

To disable PREFER_MERGE optimization, where the bubblesort will run only in the 
first iteration and then the merge algorithm will run, use:
```
mpicc bubblesort.c -o bubblesort -std=c99 -DPREFER_MERGE=0
```

To disable the SKIP_CONVERGE optimization, where when a single not sorted neighbor
array is found the broadcast will be skipped, but all processes will try to 
converge even if already sorted, use: 
```
mpicc bubblesort.c -o bubblesort -std=c99 -DSKIP_CONVERGE=0
```

## Running

If you run with 1 processor, the result will be a sequential sorting. Please use
an odd number of processors.
```
$ mpirun -np NP bubblesort
```
