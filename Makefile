CUDALDFLAGS=-L/share/apps/cuda/lib64 -lcudart 
all: twostar

floydWarshall.o: floydWarshall.cu 
	nvcc -c floydWarshall.cu -o floydWarshall.o

twostar-mpi.o: twostar-mpi.c 
	mpicc -std=c99 -c twostar-mpi.c -o twostar-mpi.o

twostar: floydWarshall.o twostar-mpi.o 
	mpicc -o twostar floydWarshall.o twostar-mpi.o $(CUDALDFLAGS) -lstdc++
	
clean:
	rm -rf *.o *~ twostar
