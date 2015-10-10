CUDALDFLAGS=-L/share/apps/cuda/lib64 -lcudart 
all: twostar

floydWarshall.o: floydWarshall.cu 
	nvcc -c floydWarshall.cu -o floydWarshall.o

twostar-mpi.o: twostar-mpi.c 
	mpicc -std=c99 -c twostar-mpi.c -o twostar-mpi.o

utils.o: utils.c
	gcc -std=c99 -c utils.c -o utils.o

twostar: floydWarshall.o twostar-mpi.o utils.o 
	mpicc -o twostar floydWarshall.o utils.o twostar-mpi.o $(CUDALDFLAGS) -lstdc++
	
clean:
	rm -rf *.o *~ twostar
