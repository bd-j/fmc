test: src/*.f95
	gfortran src/emcee.f95 src/test.f95 -o bin/test
