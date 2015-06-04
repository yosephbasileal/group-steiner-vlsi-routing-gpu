//headers
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <getopt.h>

//constant
#define INF 1061109567

//File headers
#include "utils.h"
#include "readFile.h"
#include "onestar.h"
#include "twostar.h"
#include "floydWarshall.h"

//Global variables
bool gprint = false; // print graph and metric closure
bool debug = false;	// print more deatails for debugging

//main
int main(int argc, char *argv[])
{	
	//initialize MPI variables
	int numProc, procId;

	MPI_Status status;
	MPI_Request request;
	
	MPI_Init(&argc,&argv);  
	MPI_Comm_size(MPI_COMM_WORLD,&numProc);  
	MPI_Comm_rank(MPI_COMM_WORLD,&procId);
	
	//graph variables
	unsigned int V, E, numTer, numGroups;
	int *D, *G, *term, *groups, *D_sub, *onestar, *onestar_sub;
	
	int MINIMUM, overall_min;
	
	int *pars;
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
		
		//read from file and allocate memory
		readFile(&D, &G, &term, &groups, &V, &E, &numTer, &numGroups);
	
		//broadcast size variables
		MPI_Bcast(&V, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&numTer, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&numGroups, 1, MPI_INT, 0, MPI_COMM_WORLD);

		//buffer for combined onestar cost matrix
		onestar = (int *) malloc(sizeof(int) * V * numGroups);

		//broadbast groups
		MPI_Bcast(groups, numTer, MPI_INT, 0, MPI_COMM_WORLD);
	
		//validate number of processes
		if(!validNumProc(V, numProc)) {
			printf("Error: More number of compute nodes than needed.\n");
			MPI_Finalize();
			return 0;
		}

		//calculate roots per process
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

		//output metric closure
		if(gprint) print(D,V);

		//reciveing buffer for distributing the metric closure
		D_sub = (int *) malloc(sizeof(int) * V * perChild);

		//buffer for sub onestar matrix in each process
		onestar_sub = (int *) malloc(sizeof(int) * perChild * numGroups);

		//construct one star
		onestarWrapper(V,numTer,perChild,perParent,numProc,procId,numGroups,D,D_sub,onestar,onestar_sub,groups);

		//broadbast onestar
		MPI_Bcast(onestar, V * numGroups, MPI_INT, 0, MPI_COMM_WORLD);
		
		//output onestar
		if(debug) printOnestar(onestar,numGroups,V);
		
		//construct two star
		MINIMUM = twostarwrapper(V,numGroups,perChild,perParent,numProc,procId,D,onestar);					
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

		//validate number of processes
		if (!validNumProc(V, numProc)) {
			MPI_Finalize();
			return 0;
		}

		//calculate number of roots per process
		pars = calcLaunchPar(numProc, V);
		perParent = pars[1];
		perChild = pars[0];

		//broadcast metric closure
		MPI_Ibcast(D,V*V, MPI_INT, 0, MPI_COMM_WORLD,&request);

		//buffer for reciving rows of metric closure
		D_sub = (int *) malloc(sizeof(int) * V * perChild);

		//buffer for sub onestar matrix in each process
		onestar_sub = (int *) malloc(sizeof(int) * perChild * numGroups);

		//construct one star
		onestarWrapper(V,numTer,perChild,perParent,numProc,procId,numGroups,D,D_sub,onestar,onestar_sub,groups);
		
		//recieve onestar
		MPI_Bcast(onestar, V * numGroups, MPI_INT, 0, MPI_COMM_WORLD);

		//construct two star
		MINIMUM = twostarwrapper(V,numGroups,perChild,perParent,numProc,procId,D,onestar);		
		
		//get minimum of all			
		MPI_Reduce(&MINIMUM,&overall_min,1,MPI_INT,MPI_MIN,0,MPI_COMM_WORLD);
	}//end child processes		
	
	MPI_Finalize();
	return 0;
}
