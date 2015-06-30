size_t getline(char **lineptr, size_t *n, FILE *stream);

void readFile2(int ** D, int ** G, int ** P, int ** C, int ** term, int ** groups, int * V, int * E, int * numTer, int * numGroups) {
	unsigned int v1, v2, w;
	int x, y;
	char *buffer = NULL;
	size_t size, s;

	//skip Comment section
	while((s = getline(&buffer, &size, stdin)) != -1) {
		if(strcmp(buffer,"Section Graph\n") == 0) {
			break;
		}
	}
	fscanf(stdin,"%s %d",buffer,V); //read # of Nodes
	fscanf(stdin,"%s %d",buffer,E); //read # of edges
	
	printf("V: %d  E: %d", *V, *E);
	
	//allocate memory
	*D = (int *) malloc(sizeof(int) * (*V) * (*V));//metric closure
	*G = (int *) malloc(sizeof(int) * (*V) * (*V));//graph
	*P = (int *) malloc(sizeof(int) * (*V) * (*V));//predecessors
	*C = (int *) malloc(sizeof(int) * (*V) * 2);//coordinates

	// Init Data for the graph G and p
	//memset(P, NONE, sizeof(int) * (*V) * (*V));

	//initialize graph to INF
	for(int i = 0; i < (*V); i++) {
		for(int j = 0; j < (*V); j++) {
			(*G)[i * (*V) + j] = INF;
			(*P)[i * (*V) + j] = NONE;
		}
	}
	
	//read graph and number of terminals
	for(int e = 0; e < (*E); e++) {
		fscanf(stdin,"%s %d %d %d", buffer, &v1, &v2, &w);
		(*G)[(v1 - 1) * (*V) + (v2 - 1)] = w;
		(*G)[(v2 - 1) * (*V) + (v1 - 1)] = w;
		if ((v1 - 1) != (v2 - 1)) {
			(*P)[(v1 - 1) * (*V) + (v2 - 1)] = (v1 - 1);
			(*P)[(v2 - 1) * (*V) + (v1 - 1)] = (v2 - 1);
		}
	}
	
	//skip Comment section
	while(getline(&buffer, &size, stdin) != -1) {
		if(strcmp(buffer,"Section Terminals\n") == 0) {
			break;
		}
	}
	
	fscanf(stdin,"%s %d",buffer,numTer); //read # of terminals
	
	//read terminals
	scanf("%d",numTer);
	*term = (int *) malloc(sizeof(int) * (*numTer));
	*groups = (int *) malloc (sizeof(int) * (*numTer));
	
	for(int i = 0; i < (*numTer); i++) {
		int v;
		fscanf(stdin, "%s %d", buffer, &v);
		(*term)[i] = v - 1;
	}
	
	int groupsSec = 0;
	//skip Comment section
	while(getline(&buffer, &size, stdin) != -1) {
		if(buffer[0] == 'S') {
			if(strcmp(buffer,"Section Groups\n") == 0) {
				groupsSec = 1;
			}
			if(strcmp(buffer,"Section Coordinates\n") == 0) {
				groupsSec = 0;
			}
			break;
		}
	}
	
	if(groupsSec) { //groups data in file
		//read groups
		fscanf(stdin,"%s %d",buffer,numGroups); //read # of groups
		for(int i = 0; i < (*numGroups); i++) {
			for(int j = 0; j < (*numTer)/(*numGroups); j++) {
				int v;
				fscanf(stdin,"%s %d",buffer, &v);
				(*groups)[(i * ((*numTer)/(*numGroups))) + j] = v - 1;
			}
		}
	}
	else { //no data, take each terminals as a group, one terminal per groups
		//copy groups from terminals
		*numGroups = *numTer;
		for(int i = 0; i < (*numGroups); i++) {
			for(int j = 0; j < (*numTer)/(*numGroups); j++) {
				int v = (*term)[i * ((*numTer)/(*numGroups)) + j];
				(*groups)[(i * ((*numTer)/(*numGroups))) + j] = v;
			}
		}
	}
	
	if(groupsSec) {
		while(getline(&buffer, &size, stdin) != -1) {
			if(strcmp(buffer,"Section Coordinates\n") == 0) {
				break;
			}
		}
	}
	
	//read coordinates
	
	for(int v = 0; v < (*V); v++) {
		fscanf(stdin,"%s %d %d %d", buffer, &v1, &x, &y);
		(*C)[v * 2 + 0] = x;
		(*C)[v * 2 + 1] = y;
	}
}
