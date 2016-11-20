sequential:
	mpicxx -O2 InputFunctions.cpp SequentialCGM.cpp ConjugateGradientsMethod.cpp GridClass.cpp -o sequential
parallel:
	mpicxx -O2 InputFunctions.cpp MPICGM.cpp ConjugateGradientsMethod.cpp GridClass.cpp -o parallel
parallel_omp:
	mpixlcxx_r -O2  InputFunctions.cpp MPIOpenMPCGM.cpp ConjugateGradientsMethod.cpp GridClass.cpp -qsmp=omp -o parallel_omp
all: 
	mpicxx -O2 InputFunctions.cpp SequentialCGM.cpp ConjugateGradientsMethod.cpp GridClass.cpp -o sequential
	mpicxx -O2 InputFunctions.cpp MPICGM.cpp ConjugateGradientsMethod.cpp GridClass.cpp -o parallel
	mpixlcxx_r -O2  InputFunctions.cpp MPIOpenMPCGM.cpp ConjugateGradientsMethod.cpp GridClass.cpp -qsmp=omp -o parallel_omp