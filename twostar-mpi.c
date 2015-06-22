//Library headers
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <getopt.h>

//Library for sched_getcpu()
#define _GNU_SOURCE 
#include <utmpx.h>

//Constant
#define INF 1061109567

//Global variables
bool gprint = false; // print graph and metric closure -o
bool debug = false;	// print more deatails for debugging -d
bool parallel = false; //construct metric closure in serial or parallel -p

//File headers
#include "lib/utils.h"
#include "lib/readFile.h"
#include "lib/onestar.h"
#include "lib/twostar.h"
#include "lib/floydSerial.h"

//Function declaration
int sched_getcpu(void);
void fw_gpu(const unsigned int n, const int * const G, int * const d);

//Main
int main(int argc, char *argv[])
{	
	//initialize MPI variables
	int numProc, procId;
	MPI_Status status;
	MPI_Request request;

	//initialize MPI environment
	MPI_Init(&argc,&argv);  
	MPI_Comm_size(MPI_COMM_WORLD,&numProc);  
	MPI_Comm_rank(MPI_COMM_WORLD,&procId);

	//graph variables
	unsigned int V, E, numTer, numGroups;
	int *D, *G, *term, *groups, *D_sub, *onestar, *onestar_sub;

	//solution variables
	int MINIMUM, overall_min;
	struct Solution solution;
	struct Solution minSolution;

	//variables for mapping roots to processes
	int *pars;
	int perParent;
	int perChild;

	//variables for calculating time
	double starttime, endtime;	

	//mapping of processes to CPU cores
	if(debug) {
		printf("ProcID = %d CpuID = %d\n", procId, sched_getcpu());
	}

	/*--------------------------------------------Parent process------------------------------------------------*/
	if(!procId)  {
		int r;
		while ((r = getopt(argc, argv, "odp")) != -1) { //command line args
			switch(r)
			{
				case 'o':
					gprint = true;
					break;
				case 'd':
					debug = true;
					break;
				case 'p':
					parallel = true;
					break;
				default:
					//printUsage();
					exit(1);
			}
		}

		//read graph from file and allocate memory
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

		//output for debug
		if(debug) {	
			printf("Number of vertices: %d\n", V);
			printf("Parent process gets: %d\n", pars[1]);
			printf("%d child processes get: %d\n", numProc - 1, pars[0]);
		}

		//construct metric closure
		if(parallel) {
			if(debug) {
				printf("Construction metric closure on the GPU...\n");			
			}
			fw_gpu(V, G, D);
		} else{
			if(debug) {
				printf("Construction metric closure on the CPU...\n");			
			}
			floydWarshall(V, G, D);
		}

		//broadcast metric closure - non-blocking
		MPI_Ibcast(D,V*V, MPI_INT, 0, MPI_COMM_WORLD,&request);

		//output for debug
		if(gprint) {
			print(D,V);
		}
		if(debug) {
			printTermGroups(numTer,numGroups,groups,term);
		}

		//reciveing buffer for distributing some rows of metric closure
		D_sub = (int *) malloc(sizeof(int) * V * perChild);

		//buffer for sub onestar matrix in each process
		onestar_sub = (int *) malloc(sizeof(int) * perChild * numGroups);

		//construct one star
		onestarWrapper(V,numTer,perChild,perParent,numProc,procId,numGroups,D,D_sub,onestar,onestar_sub,groups);

		//broadbast onestar
		MPI_Bcast(onestar, V * numGroups, MPI_INT, 0, MPI_COMM_WORLD);

		//output onestar
		if(debug) printOnestar(onestar,numGroups,V);

		//check if metric closure broadcast is done
		MPI_Wait(&request, &status);

		//construct two star
		twostarwrapper(V,numGroups,perChild,perParent,numProc,procId,D,onestar,&solution);

		//get minimum from all using reduction
		MPI_Reduce(&solution,&minSolution,1,MPI_2INT,MPI_MINLOC,0,MPI_COMM_WORLD);

		//ouput overall minimum cost
		printf("\nOVERALL MINIMUM STEINER COST: %d Root: %d\n\n", minSolution.cost, minSolution.root);

		//buildPrintSolution(minSolution.root,V,numGroups,D,onestar); //build complete solution with complete path - TODO
	
	}//end parent process
	

	/*--------------------------------------------Child processes------------------------------------------------*/
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

		//construct one star for assigned roots
		onestarWrapper(V,numTer,perChild,perParent,numProc,procId,numGroups,D,D_sub,onestar,onestar_sub,groups);

		//recieve onestar
		MPI_Bcast(onestar, V * numGroups, MPI_INT, 0, MPI_COMM_WORLD);

		//check if metric closure broadcast is done
		MPI_Wait(&request, &status);

		//construct two star
		twostarwrapper(V,numGroups,perChild,perParent,numProc,procId,D,onestar,&solution);

		//get minimum of all
		MPI_Reduce(&solution,&minSolution,1,MPI_2INT,MPI_MINLOC,0,MPI_COMM_WORLD);		

	}//end child processes

	MPI_Finalize();
	return 0;
}
