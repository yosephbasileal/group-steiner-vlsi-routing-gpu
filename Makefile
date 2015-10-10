CUDALDFLAGS=-L/share/apps/cuda/lib64 -lcudart 
all: twostar

floydWarshall.o: fw-cuda.cu 
	nvcc -c fw-cuda.cu -o fw-cuda.o

twostar-mpi.o: twostar-mpi.c 
	mpicc -std=c99 -c twostar-mpi.c -o twostar-mpi.o

utils.o: utils.c
	gcc -std=c99 -c utils.c -o utils.o

twostar: floydWarshall.o twostar-mpi.o utils.o 
	mpicc -o twostar fw-cuda.o utils.o twostar-mpi.o $(CUDALDFLAGS) -lstdc++
	
clean:
	rm -rf *.o *~ twostar
