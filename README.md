This is an MPI parallelized version of a *very simple* Fortran
implementation of the [emcee](https://github.com/dfm/emcee)
algorithm.  See [fmc](https://github.com/dfm/fmc) for the
unparallelized version from which this is adapted.

You must have already installed some version of MPI (e.g. mvapich2,
openMPI).

See `src/test_mpi.f95` for a demo and run it as follows:

```
make
mpirun -np <np> bin/test_mpi
```

where `<np>` is the number (>1) of processors that you want to use.

The results will be printed to `stdout`.

This code is licensed under the terms of the MIT license (see `LICENSE`).
