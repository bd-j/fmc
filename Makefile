test: src/*.f95
	mpifort src/emcee_mpi.f95 src/test_mpi.f95 -o bin/test_mpi
