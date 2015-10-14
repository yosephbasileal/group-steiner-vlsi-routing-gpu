#include "partial-star.h"
#include "remSpanned.h"

//Stores cost and root of solution twostar tree
struct Solution {
  long cost;
  int root;
};

//Stores all information about solution twostar tree
struct TwoStar {
  long cost;
  int root;
  int numPar;
  int* partialstars;
};

/* 
 * Constructs a twostar tree given a root
 */
long twoStar(int root,int * groupIds,int *onestar,int numGroups,int V, int *D, int *partialStar1,int *intermSet,int *newGroupIds, int * countPar) {    
  //initialize variables
  long TOTAL_COST = 0;
  int remGroups = numGroups;

  //initialize group IDs
  for(int i = 0; i < numGroups;i++)
    groupIds[i] = i;
  intermSet[0] = 0;

  //count partial stars
  int c = 0;

  //while there is a remaining group find a partial star
  while(remGroups > 0) {
    //find current partial star
    int * curParStar = partialStar1 + (c * (2 + numGroups));
    partialStar(root,groupIds,onestar,numGroups,remGroups, V, D + (root * V),curParStar,intermSet);

    //increate number of intermediates and save current intermediate
    intermSet[0] += 1;
    intermSet[intermSet[0]] = curParStar[0];

    //add intermediate to groups and root to intermediate costs
    for(int i = 0; i < curParStar[1]; i++) {
      TOTAL_COST += onestar[(curParStar[0] * numGroups) + curParStar[i+2]];
    }
    TOTAL_COST += D[root * V + curParStar[0]];

    //increment counter of partial stars    
    c++;

    //update number of remaining groups, stop there are none left
    remGroups -= curParStar[1];
    if(remGroups == 0) {
      break;
    }

    //remove groups that are already spanned
    remSpanned(groupIds,curParStar,newGroupIds,remGroups,numGroups);
    groupIds = newGroupIds;
  }//rooted two star loop
  
  *countPar = c;

  return TOTAL_COST;
}

/*
 * Contructs twostar trees with all possible roots and finds the one with minimum cost
 */
void twostarwrapper(int V, int numGroups, int perChild, int perParent, int numProc, int procId, int * D, int * onestar, struct Solution *solution, struct TwoStar * twostar) {
  //declare variables
  long TOTAL_COST;
  long MINIMUM = LONG_MAX;
  int minRoot = -1;
  int minNumPar;
  int countPar;

  //allocate memory for 2-star tree construction
  int * intermSet = (int *) malloc(sizeof(int) * (V + 1)); //index 0 shows how many intermediates are in the set
  int * groupIds = (int *) malloc(sizeof(int) * numGroups); //keeps track of groups already spanned
  int * newGroupIds = (int *) malloc(sizeof(int) * numGroups); //helper buffer when modifying groupIds
  int * partialStar1 = (int *) malloc(sizeof(int) * (2 + numGroups) * V);//[0] -> interm v,  [1] -> numGroups it spans, the rest are the groups V spans
  int * minParitalStar = (int *) malloc(sizeof(int) * (2 + numGroups) * V); //store partial stars of the solution

  //check if enough memory has been allocated
  if(!intermSet || !groupIds || !newGroupIds || !partialStar1) {
    printf("Memory allocation error inside twostar wrapper!\n");
    exit(1);
  }

  //current root
  int root;

  for(int i = 0; i < perChild; i++) {
    //get total cost of twostar tree for current root
    root = (perParent - perChild) + (i * numProc) + procId;
    TOTAL_COST = twoStar(root,groupIds,onestar,numGroups,V,D, partialStar1,intermSet,newGroupIds,&countPar);    
    //update minimum
    if(TOTAL_COST < MINIMUM) {
      MINIMUM = TOTAL_COST;
      minRoot = root;
      minNumPar = countPar;
      copypartialStar(partialStar1,minParitalStar,numGroups,countPar);
    }
  }

  //the first (perParent - perChild) roots are calculated by the first (perParent - perChild) processes
  if (procId < (perParent - perChild)) {
    //get total cost of twostar tree for current root
    root = procId;
    TOTAL_COST = twoStar(root,groupIds,onestar,numGroups,V,D, partialStar1,intermSet,newGroupIds,&countPar);
    //update minimum
    if(TOTAL_COST < MINIMUM) {
      MINIMUM = TOTAL_COST;
      minRoot = root;
      minNumPar = countPar;
      copypartialStar(partialStar1,minParitalStar,numGroups,countPar);
    }
  }

  //store solution tree
  twostar->root = minRoot;
  twostar->cost = MINIMUM;
  twostar->partialstars = minParitalStar;
  twostar->numPar = minNumPar;

  //Store solution cost and root
  solution->root = minRoot;
  solution->cost = MINIMUM;

  //free allocated space
  free(intermSet);
  free(groupIds);
  free(newGroupIds);
  free(partialStar1);
}
