//headers
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>

//File headers
#include "partialStar.h"
#include "remSpanned.h"

//constant
#define INF 1061109567

//Global variables
bool gprint = false; // print graph and metric closure
bool debug = false;	// print more deatails for debugging

//Function prototypes
void print(int *,int);                                //prints graph
void oneStarCost(int,int,int,int,int *,int *,int *,int);  //calculates rooted one star
void floydWarshall(int,int *,int *);                  //solves APSP problem
int * calcLaunchPar(int,int);                         //calculates how many roots for each compute node
bool validNumProc(int,int);                           //validates number of processes

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
	
	int * pars;
	int perParent;
	int perChild;
	
	double starttime, endtime;	
	
	//proc 0
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
	
		//read from file graph and number of terminals
		for(int e = 0; e < E; e++) {
			scanf("%d %d %d", &v1, &v2, &w);
			G[(v1 - 1) * V + (v2 - 1)] = w;
			G[(v2 - 1) * V + (v1 - 1)] = w;
		}
		
		//read from file terminals
		scanf("%d",&numTer);
		term = (int *) malloc(sizeof(int) * numTer);
		groups = (int *) malloc (sizeof(int) * numTer);
		
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
	
		starttime = MPI_Wtime();
		//broadcast metric closure
		MPI_Bcast(D,V*V, MPI_INT, 0, MPI_COMM_WORLD);
	
		if(gprint) {
			printf("Metric Closure:\n");
			print(D,V);
		}
		
		//reciveing buffer for distributing the metric closure
		D_sub = (int *) malloc(sizeof(int) * V * perChild);
		
		//buffer for sub onestar matrix in each process
		onestar_sub = (int *) malloc(sizeof(int) * perChild * numGroups);
	
		double gatherTotal = 0;
		double e,s;
		
		//Distribute one row of the metric closure to each proc at a time
		for(int i = 0; i < perChild; i++) {
			
			s = MPI_Wtime();
			MPI_Scatter(D + ((perParent - perChild) * V) + (i * (V * numProc)),V, MPI_INT, D_sub + (i * V),V, MPI_INT, 0, MPI_COMM_WORLD);			
			int src = (perParent - perChild) + (i * numProc) + procId;
			//oneStarCost(V, numTer,src,numGroups,onestar + (src * numGroups),groups,D_sub + (i * V),procId);
			oneStarCost(V, numTer,src,numGroups,onestar_sub + (i * numGroups),groups,D_sub + (i * V),procId);
			int * p = onestar + (src * numGroups);
			//printf("Proc: %d, src: %d, onestar: %2d %2d %2d %2d %2d %2d\n\n", procId, src, p[0],p[1],p[2],p[3],p[4],p[5]);
			MPI_Gather(onestar_sub + (i * numGroups),numGroups, MPI_INT, onestar + (src * numGroups),numGroups, MPI_INT, 0, MPI_COMM_WORLD);
			e = MPI_Wtime();
			gatherTotal += e - s;
			
			
			

			for(int j = 0; j < (perParent - perChild); j++) {
				int srcc = j;
				oneStarCost(V, numTer,srcc,numGroups,onestar + (srcc * numGroups),groups,D + (j * V),procId);
				int * pp = onestar + (srcc * numGroups);
				//printf("Proc: %d, src: %d, onestar: %2d %2d %2d %2d %2d %2d\n\n", procId, srcc, pp[0],pp[1],pp[2],pp[3],pp[4],pp[5]);
			}
		}
		
		//broadbast onestar
		MPI_Bcast(onestar, V * numGroups, MPI_INT, 0, MPI_COMM_WORLD);
		
		printf("Data transfter time: %lf\n",gatherTotal);
		
		endtime = MPI_Wtime();
		printf("Total time %lfs \n", endtime - starttime);
		
		if(debug) {
			
			for(int i = 0; i < V; i++) {
		     printf("Root%d:  ", i);
		     for(int j = 0; j < numGroups; j++)
		        printf("%3d  ", onestar[i * numGroups + j]);
		     printf("\n");
		  }
		}
	}
	
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
	
		pars = calcLaunchPar(numProc, V);
		perParent = pars[1];
		perChild = pars[0];
		
		//broadcast metric closure
		MPI_Bcast(D,V*V, MPI_INT, 0, MPI_COMM_WORLD);
		
		//reciveing buffer for distributing the metric closure
		D_sub = (int *) malloc(sizeof(int) * V * perChild);
	
		//buffer for sub onestar matrix in each process
		onestar_sub = (int *) malloc(sizeof(int) * perChild * numGroups);
		
		for(int i = 0; i < perChild; i++) {
	
			MPI_Scatter(D + ((perParent - perChild) * V) + (i * (V * numProc)),V, MPI_INT, D_sub + (i * V),V, MPI_INT, 0, MPI_COMM_WORLD);		
			int src = (perParent - perChild) + (i * numProc) + procId;
			if(procId) {
				oneStarCost(V, numTer,src,numGroups,onestar_sub + (i * numGroups),groups,D_sub + (i * V),procId);
				int * p = onestar_sub + (i * numGroups);
				//printf("Proc: %d, src: %d, onestar: %2d %2d %2d %2d %2d %2d\n\n", procId, src, p[0],p[1],p[2],p[3],p[4],p[5]);
			}
			MPI_Gather(onestar_sub + (i * numGroups),numGroups, MPI_INT, onestar + (src * numGroups),numGroups, MPI_INT, 0, MPI_COMM_WORLD);
		}
		
		//broadbast onestar
		MPI_Bcast(onestar, V * numGroups, MPI_INT, 0, MPI_COMM_WORLD);
	}//child processes		
		
	
		
	printf("\nConstructing two star...\n");
	int TOTAL_COST;
	int MINIMUM = INT_MAX;
	int minRoot = INT_MAX;

	//Allocate memory for 2-star
	int * intermSet = (int *) malloc(sizeof(int) * (V + 1)); //index 0 shows how many intermediates are in the set
	int * groupIds = (int *) malloc(sizeof(int) * numGroups); //keeps track of groups already spanned
	int * newGroupIds = (int *) malloc(sizeof(int) * numGroups); //helper buffer when modifying groupIds
	int * partialStar1 = (int *) malloc(sizeof(int) * (2 + numGroups));//[0] -> interm v  [1] -> numGroups it spansthe rest are the groups V spans

	int remGroups; //number of groups not spanned yer
	int root;

	for(int i = 0; i < perChild; i++) {
		root = (perParent - perChild) + (i * numProc) + procId;
		
		TOTAL_COST = 0;
		remGroups = numGroups;
		for(int i = 0; i < numGroups;i++)
			groupIds[i] = i;
		intermSet[0] = 0;

		while(remGroups > 0) { //Find two star
			partialStar(root,groupIds,onestar,numGroups,remGroups, V, D + (root * V),partialStar1,intermSet); //find current partial star
		
			//if(root == 0) {
			//	printf("intem: %d NumGr: %d IDs ",partialStar1[0],partialStar1[1]);
			//	for(int i = 0; i < partialStar1[1]; i++)
			//		printf(" %d ", partialStar1[i+2]);
			//	printf("\n");
			//}
		
			intermSet[0] += 1; //increase counter for number of intermediates
			intermSet[intermSet[0]] = partialStar1[0]; //save current intermediate

			//printf("Inetrm: %d  NumGro: %d GroupIds: ", partialStar1[0], partialStar1[1]);	

			for(int i = 0; i < partialStar1[1]; i++) {
				TOTAL_COST += onestar[(partialStar1[0] * numGroups) + partialStar1[i+2]]; //intermediate to groups cost
				//printf(" %d ", partialStar1[i + 2]);
			}

			//printf("\n");

			TOTAL_COST += D[root * V + partialStar1[0]]; //root to intermediate cost

			remGroups -= partialStar1[1]; //update number of remaining groups
			if(remGroups == 0) //if all groups spanned
				break;

			remSpanned(groupIds,partialStar1,newGroupIds,remGroups,numGroups); //remove those already spanned
			groupIds = newGroupIds;
		}//rooted two star loop
		
		//printf("V: %d cost: %d\n",root,TOTAL_COST);
		if(TOTAL_COST < MINIMUM) { //update minimum
			MINIMUM = TOTAL_COST;
			minRoot = root;
		}
	}//all two star loop
	
	
	printf("\nMINIMUM STEINER COST: %d,   root: %d, proc ID: %d\n", MINIMUM, minRoot,procId);
	MPI_Finalize();
	return 0;
}

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
	
	//printf("Proc: %d, D: %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d\n", procId, metClosure[0],metClosure[1],metClosure[2],metClosure[3],metClosure[4],metClosure[5],metClosure[6],metClosure[7],metClosure[8],metClosure[9]);
	

	for(int i = 0; i < numGroups; i++) { //for each group
		currMin = INT_MAX;
		for(int j = 0; j < numTer/numGroups; j++) { //for each vertex in current group
			vert = groups[(i * (numTer/numGroups)) + j];
			if(vert == src) { //if current vertex is same as the root
				currMin = 0;
				break;
			}
			//curr = metClosure[src * V + vert]; //get value form the metric closure
			curr = metClosure[vert];
			if(curr < currMin) //find the minimum
				currMin = curr;
		}
		//oneStar[src * numGroups + i] = currMin; //save minimum
		oneStar[i] = currMin;
	}
}





//https://www.ccv.brown.edu/doc/mixing-mpi-and-cuda.html

//https://www.pdc.kth.se/resources/software/installed-software/mpi-libraries/cuda-and-mpi

//http://devblogs.nvidia.com/parallelforall/introduction-cuda-aware-mpi/
