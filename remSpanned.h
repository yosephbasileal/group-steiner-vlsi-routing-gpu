//Removes the alraedy spanned groups after each partial star is calculated
void remSpanned(int* groupIds, int* parStar, int* newGroupsIds, int remGroups, int numGroups) {
	int flag = 0;
	for(int i = 0, j = 0; i < numGroups; i++) { //for each group
		for(int k = 0; k < parStar[1]; k++) { //go through current partial star
			if(groupIds[i] == parStar[2+k]) { //if already spanned dont include in the newGroup
				flag = 1;
				break;
			}
		}
		if(flag) {
			flag = 0; 
			continue;
		}
		newGroupsIds[j++] = groupIds[i];
	}
}
