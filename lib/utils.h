#ifndef UTILS_H
#define URILS_H

//Prototypes
int sched_getcpu(void);
bool validNumProc(int V, int numProc);
int* calcLaunchPar(int numProc, int V);
int getProcId(int root, int perChild, int perParent, int numProc, int V);

void print(FILE * fout, int* G, int V, char* name);
void printTerm(FILE * fout, int numTer, int* terminals);
void printGroups(FILE * fout, int numGroups, int numTer, int* groups);
void printPartialStars(FILE * fout, int* partialStar1, int numGroups, int count);
void printOnestar(FILE * fout, int* onestar, int numGroups, int V, char* name);

void printTwoStarCost(FILE * fout, int root, int cost);
void printCpuID(FILE * fout, int procId);

void writetoFile(int* S, int* C, int V, char* filename);

int caclGraphCost(int* G, int V);
int isTerminal(int v, int numTer, int* terminals);
int countNonTerminals(int* G, int V, int numTer, int* terminals);
void copypartialStar(int* A, int * B, int numGroups, int count);

#endif
