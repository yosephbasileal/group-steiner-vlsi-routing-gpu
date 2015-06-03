#include "partialStar.h"
#include "remSpanned.h"

int twostar(int root,int *groupIds,int *onestar,int numGroups,int remGroups,int V, int *D, int *partialStar1,int *intermSet,int *newGroupIds) {		
	int TOTAL_COST = 0;
	remGroups = numGroups;
	for(int i = 0; i < numGroups;i++)
		groupIds[i] = i;
	intermSet[0] = 0;

	while(remGroups > 0) { //Find two star
		partialStar(root,groupIds,onestar,numGroups,remGroups, V, D + (root * V),partialStar1,intermSet); //find current partial star

		intermSet[0] += 1; //increase counter for number of intermediates
		intermSet[intermSet[0]] = partialStar1[0]; //save current intermediate

		for(int i = 0; i < partialStar1[1]; i++) {
			TOTAL_COST += onestar[(partialStar1[0] * numGroups) + partialStar1[i+2]]; //intermediate to groups cost
		}

		TOTAL_COST += D[root * V + partialStar1[0]]; //root to intermediate cost

		remGroups -= partialStar1[1]; //update number of remaining groups
		if(remGroups == 0) //if all groups spanned
			break;

		remSpanned(groupIds,partialStar1,newGroupIds,remGroups,numGroups); //remove those already spanned
		groupIds = newGroupIds;
	}//rooted two star loop
	
	return TOTAL_COST;
}


int twostarwrapper(int V, int numGroups, int perChild, int perParent, int numProc, int procId, int * D, int * onestar) { 
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
		TOTAL_COST = twostar(root,groupIds,onestar,numGroups,remGroups,V,D, partialStar1,intermSet,newGroupIds);		
		if(TOTAL_COST < MINIMUM) { //update minimum
			MINIMUM = TOTAL_COST;
			minRoot = root;
		}
	}//all two star loop
	
	if (!procId) {
		for(int j = 0; j < (perParent - perChild); j++) {
			root = j;
			TOTAL_COST = twostar(root,groupIds,onestar,numGroups,remGroups,V,D, partialStar1,intermSet,newGroupIds);
			if(TOTAL_COST < MINIMUM) { //update minimum
				MINIMUM = TOTAL_COST;
				//minRoot = root;
			}
		}//for any remaining roots - by proc 0 only	
	}
	
	printf("\nMINIMUM STEINER COST: %d,   root: %d, proc ID: %d\n", MINIMUM, minRoot,procId);
	return MINIMUM;
}
