#define MAX_NORM 147483640.0

//Stores an intermediate node, its norm and groups it is connected to
struct InterInfo
{  
  int * groupss; //groups spanned
  int interm; //intermediate vetex ID
  float norm; //norm of partial star
  int numGroups; //num of groups spanned
};

//Function prototypes
void sort(int, int, int*, int*, int, int);
float norm(int, int, int, int*, int, int*, int*, int, int);
void partialStarLoop(int, int*, int*, int, int, int, int*,struct InterInfo *);
int isAvailable(int, int *);
void bubble_sort2(int * groupIds, int n);

/*
 * Construct a partial star by inserting one intermediate node between root and any subset of remaining groups
 */
void partialStar(int src, int* groupIds, int* vertToGroup, int numGroups, int remGroups, int V, int* rTOv,int * groupsSpanned, int* intermSet ) {
  //rTOv is r^th row of the metric closure

  //array of stucture storing information about intermediates and groups they are connected to
  struct InterInfo * normV = (struct InterInfo *)  malloc(sizeof(struct InterInfo) * V); //array of structure
  for(int i = 0; i < V; i++) { //allocate memory
    normV[i].groupss =  (int *) malloc(sizeof(int) * remGroups);
  }

  //get all intermediates, their norms and groups spanned
  partialStarLoop(src,groupIds,vertToGroup,numGroups,remGroups,V,rTOv,normV); //return all v with thier norms and groups spanned

  //find an intermedate with minimum norm
  int firstTime = 1;
  int minV;
  float minNorm;
  for(int v = 0; v < V; v++) {
    if(v == src)
      continue;
    if(!isAvailable(v,intermSet))
      continue;
    if(firstTime) {
      minV = v;
      minNorm = normV[minV].norm;
      firstTime = 0;
      continue;
    }
    float currNorm = normV[v].norm;
    if(currNorm < minNorm) {
      minNorm = currNorm;
      minV = v;
    }
  }

  //store groups spanned by the intermediate with minimum norm
  int numSpanned = normV[minV].numGroups;
  groupsSpanned[0] = minV;
  groupsSpanned[1] = numSpanned;
  for(int i = 0; i < numSpanned; i++) 
    groupsSpanned[2 + i] = *(normV[minV].groupss + i);

  //free allocated memory
  for(int i = 0; i < V; i++) {
    free(normV[i].groupss);
  }
  free(normV);  
}

void partialStarLoop(int src, int *groupIds, int * vertToGroup, int numGroups, int remGroups, int V, int * rTOv, struct InterInfo* normV) {
  //for each candiate intermediate
  for(int v = 0; v < V; v++) {
    //current root cannot be an intermediate
    if(v == src) {
      continue;
    }

    //initialize varialbes
    float normm = MAX_NORM + 1.0;
    int minJ = remGroups;

    //if more than one group remaning, sort them
    if(remGroups > 1) {
      bubble_sort2(groupIds, remGroups);
      sort(src, v, groupIds, vertToGroup, numGroups, remGroups);
    }
    
    //find j (subset of remaining groups) that minimizes the norm of current intermediate
    for(int j = 1; j <= remGroups; j++) {
      //calculate current norm and update minimum
      float curr = norm(src,v,j,rTOv,numGroups,vertToGroup,groupIds,V,remGroups);
      if(curr < normm) {
        normm = curr;
        minJ = j;
      }
    }

    //save intermediten node, its norm and groups it is connected to, to a struct
    normV[v].interm = v;
    normV[v].norm = normm;
    normV[v].numGroups = minJ;
    for(int i = 0; i < minJ; i++) {
      *(normV[v].groupss + i) = groupIds[i];
    }
  }
}

/*
 * Calculates norm when a given intermediate is inserted between root and groups
 */
float norm(int src, int v, int j, int* rTOv, int numGroups, int* vertToGroup, int* groupIds, int V, int remGroups) {
  //initialize varables
  int sum1 = 0, sum2 = 0;
  
  for(int i = 0; i < j; i++) {
    //sum distance from intermediate to all groups
    sum1 += vertToGroup[v * numGroups + groupIds[i]];
    //sum distance from source to all groups
    sum2 += vertToGroup[src * numGroups + groupIds[i]];
  }

  //calculate norm
  float normm;
  if(sum2 == 0) {
    normm = MAX_NORM;
  } else {
    normm = ((float) rTOv[v] + (float) sum1) / (float) sum2;
  }

  return normm;
}

/*
 * Sorts groupIDs array based on ratio values
 */
void bubble_sort(float* v, int* groupIds, int n)
{
  int c, d;
  int t1;  
  float t2;
  for (c = 0 ; c < ( n - 1 ); c++) {
    for (d = 0 ; d < n - c - 1; d++) {
      if (v[d] > v[d+1]) {
        /* Swapping */
        t2 = v[d];
        v[d] = v[d+1];
        v[d+1] = t2;
        
        t1= groupIds[d];
        groupIds[d]   = groupIds[d+1];
        groupIds[d+1] = t1;  
      }
    }
  }
}

/*
 * Sorts groupIDs array
 */
void bubble_sort2(int * groupIds, int n)
{
  int c, d;
  int t1;  
  for (c = 0 ; c < ( n - 1 ); c++) {
    for (d = 0 ; d < n - c - 1; d++) {
      if (groupIds[d] > groupIds[d+1]) {
        /* Swapping */
        
        t1= groupIds[d];
        groupIds[d]   = groupIds[d+1];
        groupIds[d+1] = t1;  
      }
    }
  }
}

/*
 * Wrapper for sorting groups based on the ratio of how far they are from the intermediate and the root
 */
void sort(int src, int interm, int* groupIds, int* vertToGroup, int numGroups, int remGroups) {
  //stores ratio cost(r,Ni)/cost(v,Ni) for each group
  float * arr = (float *) malloc(sizeof(float) * remGroups);

  //calculate ratio for each group
  for(int i = 0; i < remGroups; i++) {
    if(vertToGroup[src * numGroups + i] == 0) {
      arr[i] = (float) INT_MAX;
      continue;
    }  
    arr[i] = ((float) vertToGroup[interm * numGroups + groupIds[i]]) / ((float) vertToGroup[src * numGroups + groupIds[i]]);    
  }
  
  //sort groups based on ratio
  bubble_sort(arr, groupIds, remGroups); //bubble sort

  free(arr);
}

/*
 * Checks if current vertex is not already part of the solution
 */
int isAvailable(int v, int * intermSet) { //checks if current vertex is not already part of the solution
  for(int i = 0; i < intermSet[0]; i++)
    if(v == intermSet[i+1])
      return 0;
  return 1;
}
