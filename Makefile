CUDALDFLAGS=-L/share/apps/cuda/lib64 -lcudart 
all: onestar

floydWarshall.o: floydWarshall.cu 
	nvcc -c floydWarshall.cu -o floydWarshall.o

onestar-mpi.o: onestar-mpi.c 
	mpicc -std=c99 -c onestar-mpi.c -o onestar-mpi.o

onestar: floydWarshall.o onestar-mpi.o 
	mpicc -o onestar floydWarshall.o onestar-mpi.o $(CUDALDFLAGS) -lstdc++
	
clean:
	rm -rf *.o *~ onestar
