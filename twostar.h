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
