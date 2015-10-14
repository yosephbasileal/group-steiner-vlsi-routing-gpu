CUDALDFLAGS=-L/share/apps/cuda/lib64 -lcudart 
all: main

floydWarshall.o: fw-cuda.cu 
	nvcc -c fw-cuda.cu -o fw-cuda.o

main.o: main.c 
	mpicc -std=c99 -c main.c -o main.o

utils.o: utils.c
	gcc -std=c99 -c utils.c -o utils.o

main: floydWarshall.o main.o utils.o 
	mpicc -o main fw-cuda.o utils.o main.o $(CUDALDFLAGS) -lstdc++
	
clean:
	rm -rf *.o *~ main
