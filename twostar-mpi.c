//headers
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <getopt.h>

//File headers
#include "twostar.h"
#include "floydWarshall.h"

//constant
#define INF 1061109567

//Global variables
bool gprint = false; // print graph and metric closure
bool debug = false;	// print more deatails for debugging

//Function prototypes
void print(int *,int);                                    //prints graph
void oneStarCost(int,int,int,int,int *,int *,int *,int);  //calculates rooted one star
void floydWarshall(int,int *,int *);                      //solves APSP problem
int *calcLaunchPar(int,int);                              //calculates how many roots for each compute node
bool validNumProc(int,int);                               //validates number of processes



//main
int main(int argc, char *argv[])
{	
	//initialize MPI variables
	int numProc;
	int procId;

	MPI_Status status;
	MPI_Request request;
	
	MPI_Init(&argc,&argv);  
	MPI_Comm_size(MPI_COMM_WORLD,&numProc);  
	MPI_Comm_rank(MPI_COMM_WORLD,&procId);
	
	//graph variables
	unsigned int V;
	unsigned int E;
	unsigned int numTer;
	unsigned int numGroups;
	int *D, *G, *term, *groups, *D_sub, *onestar, *onestar_sub;
	int overall_min;
	
	int * pars;
	int perParent;
	int perChild;
	
	double starttime, endtime;	
	
	//parent process
	if(!procId)  {
		int r;
		while ((r = getopt(argc, argv, "pd")) != -1) { //command line args
			switch(r)
			{
				case 'p':
					gprint = true;
					break;
				case 'd':
					debug = true;
					break;
				default:
					//printUsage();
					exit(1);
			}
		}

		unsigned int v1, v2, w;

		//get size of graph	
		scanf("%d %d", &V, &E);

		//allocate memory
		D = (int *) malloc(sizeof(int) * V * V);
		G = (int *) malloc(sizeof(int) * V * V);

		//initialize graph to INF
		for(int i = 0; i < V; i++)
				for(int j = 0; j < V; j++)
			   	G[i * V + j] = INF;

		//read graph and number of terminals
		for(int e = 0; e < E; e++) {
			scanf("%d %d %d", &v1, &v2, &w);
			G[(v1 - 1) * V + (v2 - 1)] = w;
			G[(v2 - 1) * V + (v1 - 1)] = w;
		}

		//read terminals
		scanf("%d",&numTer);
		term = (int *) malloc(sizeof(int) * numTer);
		groups = (int *) malloc (sizeof(int) * numTer);

		for(int i = 0; i < numTer; i++) {
			int v;
			scanf("%d", &v);
			term[i] = v - 1;
		}

		//read groups
		scanf("%d", &numGroups);
		for(int i = 0; i < numGroups; i++) {
			for(int j = 0; j < numTer/numGroups; j++) {
			  int v;
			  scanf("%d", &v);
			  groups[(i * (numTer/numGroups)) + j] = v - 1;
			}
		}
	
		//broadcast size variables
		MPI_Bcast(&V, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&numTer, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&numGroups, 1, MPI_INT, 0, MPI_COMM_WORLD);

		//buffer for combined onestar cost matrix
		onestar = (int *) malloc(sizeof(int) * V * numGroups);

		//broadbast groups
		MPI_Bcast(groups, numTer, MPI_INT, 0, MPI_COMM_WORLD);
	
		if (!validNumProc(V, numProc)) {
			printf("Error: More number of compute nodes than needed.\n");
			MPI_Finalize();
			return 0;
		}

		pars = calcLaunchPar(numProc, V);
		perParent = pars[1];
		perChild = pars[0];

		if(debug) {	
			printf("Number of vertices: %d\n", V);
			printf("Parent process gets: %d\n", pars[1]);
			printf("%d child processes get: %d\n", numProc - 1, pars[0]);
		}

		//construct metric closure
		fw_gpu(V, G, D);

		//broadcast metric closure
		MPI_Ibcast(D,V*V, MPI_INT, 0, MPI_COMM_WORLD,&request);

		if(gprint) {
			printf("Metric Closure:\n");
			print(D,V);
		}

		//reciveing buffer for distributing the metric closure
		D_sub = (int *) malloc(sizeof(int) * V * perChild);

		//buffer for sub onestar matrix in each process
		onestar_sub = (int *) malloc(sizeof(int) * perChild * numGroups);

		//Distribute one row of the metric closure to each proc at a time and construct one star
		int src;
		for(int i = 0; i < perChild; i++) {
			src = (perParent - perChild) + (i * numProc) + procId;
			MPI_Scatter(D + ((perParent - perChild) * V) + (i * (V * numProc)),V, MPI_INT, D_sub + (i * V),V, MPI_INT, 0, MPI_COMM_WORLD);			
			oneStarCost(V, numTer,src,numGroups,onestar_sub + (i * numGroups),groups,D_sub + (i * V),procId);
			MPI_Gather(onestar_sub + (i * numGroups),numGroups, MPI_INT, onestar + (src * numGroups),numGroups, MPI_INT, 0, MPI_COMM_WORLD);

			//Construct one star for any remaining roots that did not get divided among the processes
			for(int j = 0; j < (perParent - perChild); j++) {
				int srcc = j;
				oneStarCost(V, numTer,srcc,numGroups,onestar + (srcc * numGroups),groups,D + (j * V),procId);
			}
		}

		//broadbast onestar
		MPI_Bcast(onestar, V * numGroups, MPI_INT, 0, MPI_COMM_WORLD);
		
		if(debug) {
			for(int i = 0; i < V; i++) {
				 printf("Root%d:  ", i);
				 for(int j = 0; j < numGroups; j++)
					  printf("%3d  ", onestar[i * numGroups + j]);
				 printf("\n");
			}
		}
		
		//construct two star
		int MINIMUM = twostarwrapper(V,numGroups,perChild,perParent,numProc,procId,D,onestar);					
		MPI_Reduce(&MINIMUM,&overall_min,1,MPI_INT,MPI_MIN,0,MPI_COMM_WORLD);
		
		printf("\nOVERALL MINIMUM STEINER COST: %d\n", overall_min);
	}//end parent process
	
	
	//child processes
	if(procId) {
		//broadcast size variables
		MPI_Bcast(&V, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&numTer, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&numGroups, 1, MPI_INT, 0, MPI_COMM_WORLD);

		//allocate memory
		D = (int *) malloc(sizeof(int) * V * V);
		groups = (int *) malloc (sizeof(int) * numTer);

		//buffer for combined onestar cost matrix
		onestar = (int *) malloc(sizeof(int) * V * numGroups);

		//broadbast groups
		MPI_Bcast(groups, numTer, MPI_INT, 0, MPI_COMM_WORLD);

		//terminate if any excess processors
		if (!validNumProc(V, numProc)) {
			printf("Error: More number of compute nodes than needed.\n");
			MPI_Finalize();
			return 0;
		}

		//calculate number of roots per proc
		pars = calcLaunchPar(numProc, V);
		perParent = pars[1];
		perChild = pars[0];

		//broadcast metric closure
		MPI_Ibcast(D,V*V, MPI_INT, 0, MPI_COMM_WORLD,&request); //MPI-wait() should be added somewhere for safety

		//reciveing buffer for distributing the metric closure
		D_sub = (int *) malloc(sizeof(int) * V * perChild);

		//buffer for sub onestar matrix in each process
		onestar_sub = (int *) malloc(sizeof(int) * perChild * numGroups);

		int src;
		for(int i = 0; i < perChild; i++) {
			src = (perParent - perChild) + (i * numProc) + procId;
			MPI_Scatter(D + ((perParent - perChild) * V) + (i * (V * numProc)),V, MPI_INT, D_sub + (i * V),V, MPI_INT, 0, MPI_COMM_WORLD);		
			oneStarCost(V, numTer,src,numGroups,onestar_sub + (i * numGroups),groups,D_sub + (i * V),procId);
			MPI_Gather(onestar_sub + (i * numGroups),numGroups, MPI_INT, onestar + (src * numGroups),numGroups, MPI_INT, 0, MPI_COMM_WORLD);
		}
		
		//broadbast onestar
		MPI_Bcast(onestar, V * numGroups, MPI_INT, 0, MPI_COMM_WORLD);

		//construct two star
		int MINIMUM = twostarwrapper(V,numGroups,perChild,perParent,numProc,procId,D,onestar);					
		MPI_Reduce(&MINIMUM,&overall_min,1,MPI_INT,MPI_MIN,0,MPI_COMM_WORLD);
	}//end child processes		
	
	MPI_Finalize();
	return 0;
}

//calculates how many roots each proc gets
int * calcLaunchPar(int numProc, int V) {
	int * par = (int *) malloc(sizeof(int) * 2);
	
	if (numProc > V) {
		numProc = V;
	}
	
	int d = V/numProc;
	int r = V%numProc;
	
	if (r == 0) {
		par[0] = d;
		par[1] = d;
	} else {
		par[0] = d;
		par[1] = d + r;
	}

	return par;
}

//validates proc ID
bool validNumProc(int V, int numProc) {
	return (numProc > 0 && numProc <= V);
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

//One star algorithm, calculates one-star cost for each root and stores in 'oneStar' array
void oneStarCost(int V, int numTer, int src,  int numGroups, int * oneStar, int * groups, int * metClosure, int procId) {
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
			curr = metClosure[vert]; //get value form the metric closure
			if(curr < currMin) //find the minimum
				currMin = curr;
		}
		oneStar[i] = currMin; //save minimum
	}
}
