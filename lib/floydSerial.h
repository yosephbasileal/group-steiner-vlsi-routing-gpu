//Serial floyd warshall algorithm for APSP problem
void floydWarshall(int V, int *graph, int *dist ) {
 	int i, j, k;
	for (i = 0; i < V; i++)
  		for (j = 0; j < V; j++)
     		dist[i * V + j] = graph[i * V + j]; //copy graph to solution matrix
 
	for (k = 0; k < V; k++)
   	for (i = 0; i < V; i++) //pick source
      	for (j = 0; j < V; j++) //pick destinations one by one
         	if (dist[i * V + k] + dist[k * V + j] < dist[i * V + j])
         		dist[i * V + j] = dist[i * V + k] + dist[k * V + j];
}
