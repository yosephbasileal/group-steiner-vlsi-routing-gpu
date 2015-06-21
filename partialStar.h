#define MAX_NORM 147483640.0

struct InterInfo
{	
	int * groupss; //groups spanned
	int interm; //intermediate vetex ID
   float norm; //norm of partial star
   int numGroups; //num of groups spanned
};


void sort(int, int, int*, int*, int, int);
float norm(int, int, int, int*, int, int*, int*, int, int);
void partialStarLoop(int, int*, int*, int, int, int, int*,struct InterInfo *);
int isAvailable(int, int *);
void bubble_sort2(int * groupIds, int n);

void partialStar(int src, int* groupIds, int* vertToGroup, int numGroups, int remGroups, int V, int* rTOv,int * groupsSpanned, int* intermSet ) {
	//rTOv is r^th row of the metric closure
	struct InterInfo * normV = (struct InterInfo *)  malloc(sizeof(struct InterInfo) * V); //array of structure

	for(int i = 0; i < V; i++) { //allocate memory
		normV[i].groupss =  (int *) malloc(sizeof(int) * remGroups);
	}

	partialStarLoop(src,groupIds,vertToGroup,numGroups,remGroups,V,rTOv,normV); //return all v with thier norms and groups spanned

	//for(int i = 0; i < V; i++) 
	//	printf("V: %d norm: %f \n",i,normV[i].norm);
	int firstTime = 1;
	int minV;
	float minNorm;
	for(int v = 0; v < V; v++) { //find v with minimum norm
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

	int numSpanned = normV[minV].numGroups; //num of groups spanned by intermediate v

	groupsSpanned[0] = minV;
	groupsSpanned[1] = numSpanned;
	for(int i = 0; i < numSpanned; i++) 
		groupsSpanned[2 + i] = *(normV[minV].groupss + i); //copy groups spaneed to an array

	for(int i = 0; i < V; i++) { //deallocate memory
		free(normV[i].groupss);
	}
	free(normV);	
}

void partialStarLoop(int src, int *groupIds, int * vertToGroup, int numGroups, int remGroups, int V, int * rTOv, struct InterInfo* normV) {
	for(int v = 0; v < V; v++) { //for each V
		float normm;

		if(src  == 0 && (v == 8 || v == 7)) {
			//printf("rem: %d v: %d\n", remGroups,v);
			//printf("root: %d\n",src);
			//for(int i = 0; i < remGroups; i++) 
			//	printf("%d ", groupIds[i]);
			//printf("\n\n");
		}

		if(v == src) //if root
			continue;
		normm = MAX_NORM + 1.0;
		int minJ = remGroups;
	
		//if (src == 0) printf("Interm: %d\n",v);

		if(remGroups > 1) {
			bubble_sort2(groupIds, remGroups);
			sort(src, v, groupIds, vertToGroup, numGroups, remGroups); //sort groups
		}
		
		for(int j = 1; j <= remGroups; j++) { //find j that minimizes the norm of current v
			float curr = norm(src,v,j,rTOv,numGroups,vertToGroup,groupIds,V,remGroups);
			//if(src == 0) printf("Current norm: %f \n",curr);
			if(curr < normm) {
				normm = curr;
				minJ = j;
			}
		}

		if(src  == 0 && (v == 8 || v == 7)) {
			//printf("minNomr: %f\n",normm);
			//
		}
		//printf("minJ: %d\n", minJ);

		normV[v].interm = v; //save into struct
		normV[v].norm = normm;
		normV[v].numGroups = minJ;

		for(int i = 0; i < minJ; i++)  //copy the group spanned by v into an array
			*(normV[v].groupss + i) = groupIds[i];
	}
}

float norm(int src, int v, int j, int* rTOv, int numGroups, int* vertToGroup, int* groupIds, int V, int remGroups) {
	int sum1 = 0, sum2 = 0;
	
	/*if(v == 127) {
		printf("root: %d interm: %d j: %d\n",src,v,j);
	
		for(int i = 0; i < remGroups; i++) 
			printf("%d ", groupIds[i]);
		printf("\n");
	}*/

	for(int i = 0; i < j; i++) {
		//if(v == 117) printf("%dth loop117 sum1: %d sum2: %d\n",i,sum1,sum2);
		sum1 += vertToGroup[v * numGroups + groupIds[i]];
		sum2 += vertToGroup[src * numGroups + groupIds[i]];
	}
	//if(src == 0) printf("interm:%d     sum1: %d sum2: %d\n",v,sum1,sum2);
	float normm;
	if(sum2 == 0) {
		//printf("yeahhhhhh src: %d interm: %d \n",src,v);
		normm = MAX_NORM;
	} else {
		normm = ((float) rTOv[v] + (float) sum1) / (float) sum2;
	}
	return normm;
}

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

//this just sorts one array
void bubble_sort2(int * groupIds, int n)
{
  int c, d;
  int t1;  
  //float t2;
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

void sort(int src, int interm, int* groupIds, int* vertToGroup, int numGroups, int remGroups) {
	
	float * arr = (float *) malloc(sizeof(float) * remGroups); //stores cost(r,Ni)/cost(v,Ni) for each group
	if(arr == NULL) {
		return;
	}
	for(int i = 0; i < remGroups; i++) {
		if(vertToGroup[src * numGroups + i] == 0) {
			arr[i] = (float) INT_MAX;
			continue;
		}	
		arr[i] = ((float) vertToGroup[interm * numGroups + i]) / ((float) vertToGroup[src * numGroups + i]);
		//arr[i] = ((float) vertToGroup[src * numGroups + groupIds[i]]) / ((float) vertToGroup[interm * numGroups + groupIds[i]]);		
	}
	
	/*if(src == 0) {
		printf("\nbefore for interm %d: ", interm);
		for(int i = 0; i < remGroups; i++) {
			printf(" %d(%.12f) \n", groupIds[i],arr[i]);
		}
		printf("\n");
	}*/
	bubble_sort(arr, groupIds, remGroups); //bubble sort
	
	/*if(src == 0) {
		printf("sorted for interm %d: ", interm);
 		for(int i = 0; i < remGroups; i++) {
			printf(" %d(%.12f) \n", groupIds[i],arr[i]);
		}
		printf("\n");
	}*/
	free(arr);
}


int isAvailable(int v, int * intermSet) { //checks if current vertex is not already part of the solution
	for(int i = 0; i < intermSet[0]; i++)
		if(v == intermSet[i+1])
			return 0;
	return 1;
}
