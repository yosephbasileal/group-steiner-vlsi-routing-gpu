#ifndef UTILS_H
#define URILS_H

//Prototypes
int sched_getcpu(void);
bool validNumProc(int V, int numProc);
int* calcLaunchPar(int numProc, int V);
int getProcId(int root, int perChild, int perParent, int numProc, int V);

void print(int* G, int V, char* name);
void printTerm(int numTer, int* terminals);
void printGroups(int numGroups, int numTer, int* groups);
void printPartialStars(int* partialStar1, int numGroups, int count);
void printOnestar(int* onestar, int numGroups, int V, char* name);

void printTwoStarCost(int root, int cost);
void printCpuID(int procId);

int caclGraphCost(int* G, int V);
int isTerminal(int v, int numTer, int* terminals);
int countNonTerminals(int* G, int V, int numTer, int* terminals);
void copypartialStar(int* A, int * B, int numGroups, int count);

#endif
