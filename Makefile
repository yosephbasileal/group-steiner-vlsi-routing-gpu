CUDALDFLAGS=-L/share/apps/cuda/lib64 -lcudart 
all: onestar

floydWarshall.o: floydWarshall.cu 
	nvcc -c floydWarshall.cu -o floydWarshall.o

onestar-mpi2.o: onestar-mpi2.c 
	mpicc -std=c99 -c onestar-mpi2.c -o onestar-mpi2.o

onestar: floydWarshall.o onestar-mpi2.o 
	mpicc -o onestar floydWarshall.o onestar-mpi2.o $(CUDALDFLAGS) -lstdc++
	
clean:
	rm -rf *.o *~ onestar
