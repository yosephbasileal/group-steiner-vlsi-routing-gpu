//Library headers
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <getopt.h>
#include <string.h>
#include <time.h>

//Global variables
bool gprint = false; // print graph and metric closure -o
bool debug = false;	// print more deatails for debugging -d
bool serial = false; //construct metric closure in serial or parallel -n
bool stpFile = false; //an option to read standard stp file -t

//File headers
#include "lib/macros.h"
#include "lib/utils.h"
#include "lib/readFile.h"
#include "lib/readFile2.h"
#include "lib/floydSerial.h"
#include "lib/onestar.h"
#include "lib/twostar.h"
#include "lib/buildsolution.h"

void fw_gpu(const unsigned int n, const int * const G, int * const d, int * const p);

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
	int *D, *G, *P, *C, *term, *groups, *D_sub, *onestar, *onestar_sub, *onestar_V, *onestar_sub_V;

	//solution variables
	int MINIMUM, overall_min;
	struct Solution solution;
	struct Solution minSolution;
	struct TwoStar twostar;
	//struct Solution solutionTree

	//variables for mapping roots to processes
	int *pars;
	int perParent;
	int perChild;

	FILE * fout = fopen("twostar.out", "w+");

	//time variables
	clock_t start, end;
	double cpu_time_used;

	/*--------------------------------------------Parent process------------------------------------------------*/
	if(!procId)  {
		int r;
		while ((r = getopt(argc, argv, "odnt")) != -1) { //command line args
			switch(r)
			{
				case 'o':
					gprint = true;
					break;
				case 'd':
					debug = true;
					break;
				case 'n':
					serial = true;
					break;
				case 't':
					stpFile = true;
					break;
				default:
					//printUsage();
					exit(1);
			}
		}

		start = clock();
		printf("Reading graph from file...\n");
		//read graph from file and allocate memory
		if(stpFile) {
			readFile2(&D, &G, &P, &C, &term, &groups, &V, &E, &numTer, &numGroups);
		}
		else {
			readFile(&D, &G, &P, &term, &groups, &V, &E, &numTer, &numGroups);
		}
		
		fprintf(fout,"\nV: %d  E: %d\n",V,E);

		//broadcast size variables
		MPI_Bcast(&V, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&numTer, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&numGroups, 1, MPI_INT, 0, MPI_COMM_WORLD);

		//buffer for combined onestar cost matrix
		onestar = (int *) malloc(sizeof(int) * V * numGroups);
		onestar_V = (int *) malloc(sizeof(int) * V * numGroups);

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
			fprintf(fout, "\nParent process gets: %d process(es)\n", pars[1]);
			fprintf(fout, "%d child process(es) get: %d process(es)\n", numProc - 1, pars[0]);
		}

		//construct metric closure
		if(!serial) {
			printf("Constructing metric closure on the GPU...\n");
			fw_gpu(V, G, D, P);
		} else{
			if(debug) {
				printf("Constructing metric closure on the CPU...\n");			
			}
			//floydWarshall(V, G, D);
			floydWarshallWithPath(V,G,D,P);
		}

		//broadcast metric closure - non-blocking
		MPI_Ibcast(D,V*V, MPI_INT, 0, MPI_COMM_WORLD,&request);

		//output for debug
		if(gprint) {
			print(fout,G,V,"Graph");
			print(fout,D,V,"Metric Closure");
			print(fout,P,V,"Predecessors");
			printTerm(fout,numTer, term);
			printGroups(fout,numGroups, numTer, groups);
		}
		
		printf("Constructing onestar tress...\n");	
		//reciveing buffer for distributing some rows of metric closure
		D_sub = (int *) malloc(sizeof(int) * V * perChild);

		//buffer for sub onestar matrix in each process
		onestar_sub = (int *) malloc(sizeof(int) * perChild * numGroups);
		onestar_sub_V= (int *) malloc(sizeof(int) * perChild * numGroups);

		//construct one star
		onestarWrapper(V,numTer,perChild,perParent,numProc,procId,numGroups,D,D_sub,onestar,onestar_sub,onestar_V, onestar_sub_V,groups);

		//broadbast onestar
		MPI_Bcast(onestar, V * numGroups, MPI_INT, 0, MPI_COMM_WORLD);

		//output onestar
		if(gprint) {
			printOnestar(fout,onestar,numGroups,V,"One star cost");
			printOnestar(fout,onestar_V,numGroups,V,"One star vertices");
		}

		//check if metric closure broadcast is done
		MPI_Wait(&request, &status);
		
		printf("Constructing twostar trees...\n");	
		//construct two star
		twostarwrapper(V,numGroups,perChild,perParent,numProc,procId,D,onestar,&solution,&twostar);

		//get minimum from all using reduction
		MPI_Reduce(&solution,&minSolution,1,MPI_LONG_INT,MPI_MINLOC,0,MPI_COMM_WORLD);

		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		
		fprintf(fout, "\nTotal CPU time: %lf\n", cpu_time_used);

		//ouput overall minimum cost
		//if(!build) {		
			//printf("\nOVERALL MINIMUM STEINER COST: %d Root: %d\n\n", minSolution.cost, minSolution.root);
		//}
		

		//if(build) {
		printf("Post-processing tree to construct solution graph...\n");	
		buildWrapper(fout, minSolution,V,E,numGroups,P,G,D,C,onestar,onestar_V,term,numTer,perParent,perChild,numProc,procId,&twostar);
		//}
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
		onestar_V = (int *) malloc(sizeof(int) * V * numGroups);

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
		onestar_sub_V = (int *) malloc(sizeof(int) * perChild * numGroups);

		//construct one star for assigned roots
		onestarWrapper(V,numTer,perChild,perParent,numProc,procId,numGroups,D,D_sub,onestar,onestar_sub,onestar_V, onestar_sub_V,groups);

		//recieve onestar
		MPI_Bcast(onestar, V * numGroups, MPI_INT, 0, MPI_COMM_WORLD);

		//check if metric closure broadcast is done
		MPI_Wait(&request, &status);

		//construct two star
		twostarwrapper(V,numGroups,perChild,perParent,numProc,procId,D,onestar,&solution,&twostar);

		//get minimum of all
		MPI_Reduce(&solution,&minSolution,1,MPI_LONG_INT,MPI_MINLOC,0,MPI_COMM_WORLD);		

		buildWrapper(fout, minSolution,V,E,numGroups,P,G,D,C,onestar,onestar_V,term,numTer,perParent,perChild,numProc,procId,&twostar);
	}//end child processes

	MPI_Finalize();
	return 0;
}
