void drawEdge(int i, int j, int V, int * S, int * G) {
	S[i * V + j] = G[i * V + j];
	S[j * V + i] = G[j * V + i];
}

int getPath(int i, int j, int V, int * path, int * P, int * G) {
	int count = 0;
	int final = reconstruct_path(V, i,j, P, G,path,&count);
	return count;
}

int buildSolution(int minRoot, int V, int numGroups, int * D, int * onestar,struct SolutionTree * solutionTree) {
	int TOTAL_COST;

	int * intermSet = (int *) malloc(sizeof(int) * (V + 1)); 
	int * groupIds = (int *) malloc(sizeof(int) * numGroups);
	int * newGroupIds = (int *) malloc(sizeof(int) * numGroups); 
	int * partialStar1 = (int *) malloc(sizeof(int) * (2 + numGroups) * V);

	if(!intermSet || !groupIds || !newGroupIds || !partialStar1) {
		printf("Error: Memory allocation error inside twostar wrapper!\n");
		exit(1);
	}

	int remGroups;
	int root = minRoot;

	int count = twostar(root,groupIds,onestar,numGroups,remGroups,V,D, partialStar1,intermSet,newGroupIds,true);		
	
	//Store solution
	solutionTree->root = minRoot;
	solutionTree->numPar = count;
	solutionTree->partialstars = partialStar1;
	

	free(intermSet);
	free(groupIds);
	free(newGroupIds);
	
	return count;
}

void buildWrapper(struct Solution minSolution, int V, int numGroups, int * P,int * G, int * D, int * onestar, int * onestar_V, int * terminals, int numTer ) {
	// variables
	struct SolutionTree solutionTree;
	int * partialStar1;
	int * S = (int *) malloc(sizeof(int) * V * V);

	//initialize graph to INF
	for(int i = 0; i < V; i++) {
		for(int j = 0; j < V; j++) {
			S[i * V + j] = INF;
		}
	}

	//construct two star
	int count = buildSolution(minSolution.root,V,numGroups,D,onestar,&solutionTree);
	int root = minSolution.root;	

	//output	
	partialStar1 = solutionTree.partialstars;
	
	printSolutionCost(minSolution.root,minSolution.cost);
	printPartialStars(partialStar1,numGroups,solutionTree.numPar);
	
	int * path = (int *) malloc(sizeof(int) * 2 * V);
	for(int i = 0; i < count; i++) { //for each partial star
		int * curStar = partialStar1 + (i * (2 + numGroups));
		int interm = curStar[0];
		int c = getPath(root,interm,V,path,P,G);
		//printf("INTERMEDIATE: %d \n", interm);
		//printf("\troot to inter count: %d\n", c);
		for(int j = 0; j < c; j++) {
			int s = path[j * 2 + 0];
			int d = path[j * 2 + 1];
			//printf("\t\t\tS: %d  D: %d\n",s,d);
			drawEdge(s, d , V, S, G);
		}
		//printf("\n");
		for(int j = 0; j < curStar[1]; j++) { //for each groupID
			int term = onestar_V[interm * numGroups + curStar[2+j]];
			int cc = getPath(interm, term,V,path,P,G);
			//printf("\tGROUP: %d TERM: %d \n", curStar[2+j], term);
			//printf("\t\tinterm to term count: %d\n", cc);
			for(int k = 0; k < cc; k++) {
				int s = path[k * 2 + 0];
				int d = path[k * 2 + 1];
				//printf("\t\t\tS: %d  D: %d\n",s,d);
				drawEdge(s, d , V, S, G);
			}
		}
		//printf("\n\n");
	}
	//print(S,V);
	printf("Final Graph Cost: %d\n", caclGraphCost(S,V));
	//printf("NonTer part of solution: %d\n",countNonTerminals(S,V,numTer, terminals));
}
