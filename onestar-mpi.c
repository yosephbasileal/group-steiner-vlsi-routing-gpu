//headers
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

//constant
#define INF 1061109567

//Function prototypes
void print(int *,int);                                //prints graph
void oneStarCost(int,int,int,int,int *,int *,int *);  //calculates rooted one star
void floydWarshall(int,int *,int *);                  //solves APSP problem

//main
int main(int argc, char *argv[])
{	
	//initialize MPI variables
	int numberOfProcessors;
	int procId;
	MPI_Status status;
	MPI_Request request;
	
	MPI_Init(&argc,&argv);  
	MPI_Comm_size(MPI_COMM_WORLD,&numberOfProcessors);  
	MPI_Comm_rank(MPI_COMM_WORLD,&procId);

	//graph variables
	unsigned int V;
	unsigned int E;
	unsigned int numTer;
	unsigned int numGroups;
	int * D, * G, * term, * groups;
	
	double starttime, endtime;

	//using proc 0
	if(procId == 0)  {
		unsigned int v1, v2, w;
		
		//get size of graph	
		scanf("%d %d", &V, &E);
		
		//allocate memory
		D = (int *) malloc(sizeof(int) * V * V);
		G = (int *) malloc(sizeof(int) * V * V);
		term = (int *) malloc(sizeof(int) * numTer);
		groups = (int *) malloc (sizeof(int) * numTer);
		
		//initialize graph to INF
		for(int i = 0; i < V; i++)
      	for(int j = 0; j < V; j++)
         	G[i * V + j] = INF;
	
		//read from file graph and number of terminals
		for(int e = 0; e < E; e++) {
			scanf("%d %d %d", &v1, &v2, &w);
			G[(v1 - 1) * V + (v2 - 1)] = w;
			G[(v2 - 1) * V + (v1 - 1)] = w;
		}
		
		//read from file terminals
		scanf("%d",&numTer);
		for(int i = 0; i < numTer; i++) {
			int v;
			scanf("%d", &v);
			term[i] = v - 1;
		}
		
		//read from file and groups
		scanf("%d", &numGroups);
		for(int i = 0; i < numGroups; i++) {
      for(int j = 0; j < numTer/numGroups; j++) {
        int v;
        scanf("%d", &v);
        groups[(i * (numTer/numGroups)) + j] = v - 1;
      }
    }
	}

	if(procId == 0) {
		starttime = MPI_Wtime();
		fw_gpu(V, G, D);
		endtime = MPI_Wtime();
		printf("Time spent: %lf \n", endtime - starttime);
		print(D,V);
	}

	MPI_Finalize();
	return 0;
}


//Prints 2D array
void print(int * G, int V) {
   for(int i = 0; i < V; i++, printf("\n")) {
      for(int j = 0; j < V; j++) {
         int out = G[i * V + j];
         if(out  == INF)
            printf("%3s " , "INF");
         else
            printf("%3d " , out );
      }
   }
   printf("\n");
}

/*
//One star algorithm, calculates one-star cost for each root and stores in 'oneStar' array
void oneStarCost(int V, int numTer, int src,  int numGroups, int * oneStar, int * groups, int * metClosure) {
	int currMin; // Current min in the group
	int curr; // Current cost being compared
	int vert;

	for(int i = 0; i < numGroups; i++) { //for each group
		currMin = INT_MAX;
		for(int j = 0; j < numTer/numGroups; j++) { //for each vertex in current group
			vert = groups[(i * (numTer/numGroups)) + j];
			if(vert == src) { //if current vertex is same as the root
				currMin = 0;
				break;
			}
			curr = metClosure[src * V + vert]; //get value form the metric closure

			if(curr < currMin) //find the minimum
				currMin = curr;
		}
		oneStar[src * numGroups + i] = currMin; //save minimum
	}
}
*/






//https://www.ccv.brown.edu/doc/mixing-mpi-and-cuda.html

//https://www.pdc.kth.se/resources/software/installed-software/mpi-libraries/cuda-and-mpi

//http://devblogs.nvidia.com/parallelforall/introduction-cuda-aware-mpi/
