#include "partial-star.h"
#include "remSpanned.h"

struct Solution {
  long cost;
  int root;
};

struct TwoStar {
  long cost;
  int root;
  int numPar;
  int* partialstars;
};

long twoStar(int root,int *groupIds,int *onestar,int numGroups,int remGroups,int V, int *D, int *partialStar1,int *intermSet,int *newGroupIds, int * countPar) {    
  long TOTAL_COST = 0;
  remGroups = numGroups;
  for(int i = 0; i < numGroups;i++)
    groupIds[i] = i;
  intermSet[0] = 0;

  int c = 0; //count partial stars

  //Find two star
  while(remGroups > 0) {
    int * curParStar = partialStar1 + (c * (2 + numGroups));
    partialStar(root,groupIds,onestar,numGroups,remGroups, V, D + (root * V),curParStar,intermSet); //find current partial star
    intermSet[0] += 1; //increase counter for number of intermediates
    intermSet[intermSet[0]] = curParStar[0]; //save current intermediate

    for(int i = 0; i < curParStar[1]; i++) {
      TOTAL_COST += onestar[(curParStar[0] * numGroups) + curParStar[i+2]]; //intermediate to groups cost
    }
    TOTAL_COST += D[root * V + curParStar[0]]; //root to intermediate cost

    //increment counter of partial stars    
    c++;

    remGroups -= curParStar[1]; //update number of remaining groups
    if(remGroups == 0) //if all groups spanned
      break;

    remSpanned(groupIds,curParStar,newGroupIds,remGroups,numGroups); //remove those already spanned
    groupIds = newGroupIds;
  }//rooted two star loop
  
  *countPar = c;

  return TOTAL_COST;
}

void twostarwrapper(int V, int numGroups, int perChild, int perParent, int numProc, int procId, int * D, int * onestar, struct Solution *solution, struct TwoStar * twostar) {
  long TOTAL_COST;
  long MINIMUM = LONG_MAX;
  int minRoot = -1;
  int minNumPar;
  int countPar;

  //Allocate memory for 2-star
  int * intermSet = (int *) malloc(sizeof(int) * (V + 1)); //index 0 shows how many intermediates are in the set
  int * groupIds = (int *) malloc(sizeof(int) * numGroups); //keeps track of groups already spanned
  int * newGroupIds = (int *) malloc(sizeof(int) * numGroups); //helper buffer when modifying groupIds
  int * partialStar1 = (int *) malloc(sizeof(int) * (2 + numGroups) * V);//[0] -> interm v,  [1] -> numGroups it spans, the rest are the groups V spans

  int * minParitalStar = (int *) malloc(sizeof(int) * (2 + numGroups) * V); //store partial stars of the solution

  if(!intermSet || !groupIds || !newGroupIds || !partialStar1) {
    printf("Memory allocation error inside twostar wrapper!\n");
    exit(1);
  }

  int remGroups; //number of groups not spanned yet
  int root;
  
  /*MPI_Win_lock(MPI_LOCK_EXCLUSIVE,0, 0, win);
  MPI_Rget(rootAvail,V,MPI_INT,0,0,V, MPI_INT ,win,&request);
  MPI_Wait(&request,&status);
  root = getNextAvailableRoot(rootAvail,V);
  MPI_Rput(rootAvail,V,MPI_INT,0,0,V, MPI_INT ,win,&request);
  MPI_Wait(&request,&status);
  MPI_Win_unlock(0, win);
  printf("Root: %d procId: %d\n",root,procId);*/

  for(int i = 0; i < perChild; i++) {
    root = (perParent - perChild) + (i * numProc) + procId;
    TOTAL_COST = twoStar(root,groupIds,onestar,numGroups,remGroups,V,D, partialStar1,intermSet,newGroupIds,&countPar);    
    if(TOTAL_COST < MINIMUM) { //update minimum
      MINIMUM = TOTAL_COST;
      minRoot = root;
      minNumPar = countPar;
      copypartialStar(partialStar1,minParitalStar,numGroups,countPar);
    }
  }

  //The first (perParent - perChild) roots are calculated by the parent process  
  if (procId < (perParent - perChild)) {
    //for(int j = 0; j < (perParent - perChild); j++) {
      root = procId;
      TOTAL_COST = twoStar(root,groupIds,onestar,numGroups,remGroups,V,D, partialStar1,intermSet,newGroupIds,&countPar);
      if(TOTAL_COST < MINIMUM) { //update minimum
        MINIMUM = TOTAL_COST;
        minRoot = root;
        minNumPar = countPar;
        copypartialStar(partialStar1,minParitalStar,numGroups,countPar);
      }
    //}
  }


  //store solution tree
  twostar->root = minRoot;
  twostar->cost = MINIMUM;
  twostar->partialstars = minParitalStar;
  twostar->numPar = minNumPar;

  //Store solution
  solution->root = minRoot;
  solution->cost = MINIMUM;

  //if(debug) {
  //printf("MINIMUM STEINER COST: %ld,   root: %d, proc ID: %d\n", MINIMUM, minRoot,procId);
  //}
  free(intermSet);
  free(groupIds);
  free(newGroupIds);
  free(partialStar1);
}


