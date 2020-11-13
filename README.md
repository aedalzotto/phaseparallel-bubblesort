# phaseparallel-bubblesort
Distributed sort using MPI following the Phase Parallel model.

## Building

Compile with standard array size (1.000.000). It is highly recommended to 
optimize the code with -O3:
```
$ mpicc bubblesort.c -o bubblesort -std=c99 -O3
```

To choose different array size (N), define in compiler's command line
```
$ mpicc bubblesort.c -o bubblesort -std=c99 -O3 -DN=500000
```

To compile in DEBUG mode, also define in compiler's command line. The default 
array size for this mode is 40.
```
$ mpicc bubblesort.c -o bubblesort_debug -std=c99 -O3 -DDEBUG=1
```

To disable USE_COMBINE optimization, where the bubblesort will run only in the 
first iteration and then the combination algorithm will run, use:
```
mpicc bubblesort.c -o bubblesort -std=c99 -DUSE_COMBINE=0
```

To enable the FULL_CONVERGE optimization (not recommended), 
where when a single not sorted neighbor
array is found the broadcast will be skipped, but all processes will try to 
converge even if already sorted, use: 
```
mpicc bubblesort.c -o bubblesort -std=c99 -DFULL_CONVERGE=1
```

## Running

If you run with 1 processor, the result will be a sequential sorting. Please use
an odd number of processors.
```
$ mpirun -np NP bubblesort
```
