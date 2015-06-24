//Calculates one-star cost for a root and stores the result
void oneStarCost(int V, int numTer, int src,  int numGroups, int * oneStar, int * groups, int * metClosure, int procId) {
	int currMin; //current min in the group
	int curr; //current cost being compared
	int vert; //current vertex

	for(int i = 0; i < numGroups; i++) { //for each group
		currMin = INT_MAX;
		for(int j = 0; j < numTer/numGroups; j++) { //for each vertex in current group
			vert = groups[(i * (numTer/numGroups)) + j];
			if(vert == src) { //if current vertex is same as the root
				currMin = 0;
				break;
			}
			curr = metClosure[vert]; //get value form the metric closure
			if(curr < currMin) //find the minimum
				currMin = curr;
		}
		oneStar[i] = currMin; //save minimum
	}
}

void onestarWrapper(int V, int numTer, int perChild,int perParent,int numProc,int procId, int numGroups, int *D, int *D_sub, int *onestar, int *onestar_sub, int * groups) {
	//Distribute one row of the metric closure to each proc at a time and construct one star
	//The first (perParent - perChild) roots are calculated by the parent process
	//Starting from (perParent - perChild)+1'th root we do round robin between all processes including the parent
	int src;
	for(int i = 0; i < perChild; i++) {
		src = (perParent - perChild) + (i * numProc) + procId; //calculate root
		MPI_Scatter(D + ((perParent - perChild) * V) + (i * (V * numProc)),V, MPI_INT, D_sub + (i * V),V, MPI_INT, 0, MPI_COMM_WORLD);//recive a row
		oneStarCost(V, numTer,src,numGroups,onestar_sub + (i * numGroups),groups,D_sub + (i * V),procId);//construct a rooted one star
		MPI_Gather(onestar_sub + (i * numGroups),numGroups, MPI_INT, onestar + (src * numGroups),numGroups, MPI_INT, 0, MPI_COMM_WORLD);//send result to parent process

		//only for parent process
		if(!procId) {
			//Construct one star for any remaining roots that did not get divided among the processes
			for(int j = 0; j < (perParent - perChild); j++) {
				int srcc = j;
				oneStarCost(V, numTer,srcc,numGroups,onestar + (srcc * numGroups),groups,D + (j * V),procId);
			}
		}
	}
}
