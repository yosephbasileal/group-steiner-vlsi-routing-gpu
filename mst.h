//Header
#include "macros.h"

// A C / C++ program for Prim's Minimum Spanning Tree (MST) algorithm. 
// The program is for adjacency matrix representation of the graph
 
// A utility function to find the vertex with minimum key value, from
// the set of vertices not yet included in MST
int minKey(int key[], bool mstSet[], int V)
{
   // Initialize min value
   int min = INT_MAX, min_index;
 
   for (int v = 0; v < V; v++)
     if (mstSet[v] == false && key[v] < min)
         min = key[v], min_index = v;
 
   return min_index;
}

// Function to construct and print MST for a graph represented using adjacency
// matrix representation
void primMST(int * graph, int V, int * parent)
{
     int key[V];   // Key values used to pick minimum weight edge in cut
     bool mstSet[V];  // To represent set of vertices not yet included in MST
 
     // Initialize all keys as INFINITE
     for (int i = 0; i < V; i++)
        key[i] = INT_MAX, mstSet[i] = false;
 
     // Always include first 1st vertex in MST.
     key[0] = 0;     // Make key 0 so that this vertex is picked as first vertex
     parent[0] = -1; // First node is always root of MST 
 
     // The MST will have V vertices
     for (int count = 0; count < V-1; count++)
     {
        // Pick thd minimum key vertex from the set of vertices
        // not yet included in MST
        int u = minKey(key, mstSet,V);
 
        // Add the picked vertex to the MST Set
        mstSet[u] = true;
 
        // Update key value and parent index of the adjacent vertices of
        // the picked vertex. Consider only those vertices which are not yet
        // included in MST
        for (int v = 0; v < V; v++)
 
           // graph[u][v] is non zero only for adjacent vertices of m
           // mstSet[v] is false for vertices not yet included in MST
           // Update the key only if graph[u][v] is smaller than key[v]
          if (graph[u * V + v] && mstSet[v] == false && graph[u * V + v] <  key[v])
             parent[v]  = u, key[v] = graph[u * V + v];
     }
}

// A utility function to print the constructed MST stored in parent[]
int primMSTwrapper(int * graph, int V)
{
	int * parent = (int *) malloc(sizeof(int) * V); // Array to store constructed MST

	primMST(graph,V,parent);
	
	int * mst = (int *) malloc(sizeof(int) * V * V);
	for(int i = 0; i < V; i++) {
		for(int j = 0; j < V; j++) {
			mst[i * V + j] = INF;
		}
	}
   for (int i = 1; i < V; i++) {
		//printf("%d - %d    %d \n", parent[i], i, graph[i * V + parent[i]]);
		mst[parent[i] * V + i] = graph[i * V + parent[i]];
		mst[i * V + parent[i]] = graph[parent[i] * V + i];
	}
   
	for(int i = 0; i < V; i++) {
		for(int j =0; j < V; j++) {
			graph[i * V + j] = mst[i * V + j];
		}
	}
}


int getConnectedVert(int * S, int * connected, int V) {
	int count = 0;
	for(int i = 0; i < V; i++) {
		for(int j =0; j < V; j++) {
			if(i <= j) {
				continue;
			}
			int curr = S[i * V + j];
			if(curr == INF) {
				continue;
			}
			if(connected[i] == 0) {
				connected[i] = 1; //is connected, update
				count++;
			}
			if(connected[j] == 0) {	
				connected[j] = 1; //is connected, update
				count++;
			}	
		}
	}
	return count;
}

int getNewVertID(int * connected, int v) { //getID for mst graph from solution graph
	int newV = 0;	
	for(int i = 0; i < v; i++) { //count number of conected nodes before v
		if(connected[i] == 1) {
			newV++;
		}
	}
	return newV;
}

int getOldVertID(int * connected, int V,  int v) { //getID for solution graph from mst graph
	int oldV;
	int count;
	for(int i = 0; i < V; i++) { //count number of conected nodes before v
		if(connected[i] == 1) {
			count++;
		}
		if(count - 1 == v) {
			oldV = i;
			break;
		} 
	}
	return oldV;
}

void removeUnconnected(int * S, int * N, int * connected, int V, int newV) {
	for(int i = 0; i < V; i++) {
		for(int j =0; j < V; j++) {
			if(i <= j) {
				continue;
			}
			int curr = S[i * V + j];
			if(curr == INF) {
				continue;
			}
			int newI = getNewVertID(connected,i);
			int newJ = getNewVertID(connected,j);

			N[newI * newV + newJ] = curr;
			N[newJ * newV + newI] = curr;
		}
	}
}

void removeUnconnected2(int *N, int *G, int *connected, int V, int newV) {
	for(int i = 0; i < V; i++) {
		for(int j = 0; j < V; j++) {
			if(i <= j) {
				continue;
			}
			int curr = G[i * V + j];
			if(curr == INF) {
				continue;
			}
			if(connected[i] == 1 && connected[j] == 1) {
				int newI = getNewVertID(connected,i);
				int newJ = getNewVertID(connected,j);
				N[newI * newV + newJ] = curr;
				N[newJ * newV + newI] = curr;
			}
		}		
	}
}

