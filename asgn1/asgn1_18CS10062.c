/*
NAME : RADHIKA PATWARI
ROLL NO. : 18CS10062
Lab Assignment no. 1 of Algorithms 2 for Semester 5
*/

//including the C header files 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>

//structure for Edge
typedef struct EDGE{
	int y;
	int c;
	int f;
	struct EDGE* next;
}EDGE;

//structure for Vertex
typedef struct{
	int x;
	int n;
	EDGE* p;
}VERTEX;

//structure for Graph
typedef struct{
	int V;
	int E;
	VERTEX* H;
}GRAPH;

//function for adding an edge as a node to the adjancency list
void addEdge(GRAPH *graph, int x, int y, int c){

	if(graph == NULL)				//sanity check
		return;

	EDGE *edgeXY = (EDGE *)malloc(sizeof(EDGE));	//directed edge from x to y
	edgeXY->y = y;									//destination : y
	edgeXY->c = c;									//capacity : c
	edgeXY->f = 0;									//flow : 0

	EDGE *tempXY = (graph->H[x]).p;					//inserting edge node into list
	edgeXY->next = tempXY;
	(graph->H[x]).p = edgeXY;

	//inserting directed edge in opposite direction
	//however,it does not play any role in computing max flow
	/*
	EDGE *edgeYX = (EDGE *)malloc(sizeof(EDGE));	//directed edge from y to x
	edgeYX->y = x;									//destination : x
	edgeYX->c = 0;									//capacity : 0
	edgeYX->f = 0;									//flow : 0

	EDGE *tempYX = (graph->H[y]).p;					//inserting edge node into list
	edgeYX->next = tempYX;
	(graph->H[y]).p = edgeYX;*/
}
//Reading graph values from the given input file
GRAPH *ReadGraph(char *fname){           
	FILE *fileptr;
	fileptr = fopen(fname,"r");
	if(fileptr == NULL){
		printf("Error while opening file");
		return NULL;
	}

	int V, E, i, x, y, c, need;

	GRAPH *graph = (GRAPH *)malloc(sizeof(GRAPH));
	fscanf(fileptr, "%d %d", &V, &E);					//storing no. of vertices and edges

	graph->V = V;
	graph->E = E;
	graph->H = (VERTEX *)malloc((V+2)*sizeof(VERTEX));  //2 extra vertices for source and sink

	for(i=0;i<=V+1;i++){								//initialising adjancency list
		graph->H[i].x = i;								//vertex id
		graph->H[i].n = 0;								//vertex need
		graph->H[i].p = NULL;							//pointer to outgoiong edges of vertex
	}

	for(i=1;i<=V;i++){					
		fscanf(fileptr, "%d", &need);					//need of each vertex
		graph->H[i].n = need;
	}

	for(i=0;i<E;i++){
		fscanf(fileptr, "%d %d %d", &x, &y, &c);		//source,destination and capacity of edges
		addEdge(graph, x, y, c);
	}
	fclose(fileptr);									//closing file pointer
	return graph;
}
//Printing the graph with all edges having capacity not equal to 0
void PrintGraph(GRAPH graph){
	int V = graph.V;
	for(int i=1;i<=V;i++){
		printf("%d",graph.H[i].x);
		EDGE *curr = graph.H[i].p;
		while(curr != NULL){
			if(curr->y > 0 && curr->y <= V && curr->c !=0)
				printf(" -> (%d,%d,%d)",curr->y, curr->c, curr->f);
			curr = curr->next;
		}
		printf("\n");
	}
}
//finding minimum of 2 values
int min(int x, int y){
	return (x<y)?x:y;
}
//finding maximum of 2 values
int max(int x, int y){
	return (x>y)?x:y;
}

//Performing dfs for finding all possible augmenting paths 
//Finally storing smallest length augmenting path with maximum residual capacity
void dfs(GRAPH *graph, int src, int dest, int currLen, int **resCapacity,int currCapacity, int parent[], int currParent[], int **len, int visited[]){
	//base case
	if(currLen > (**len))		//ignoring paths of greater length than current path lenght
		return;

	if(src == dest){			//destination is reached
		if(currLen < (**len)){	//shorter augmenting path is obtained
			**len = currLen;		
			**resCapacity = currCapacity;
			for(int i=0;i<=(graph->V+1);i++){		//updating parents for final path 
				parent[i] = currParent[i];
			}
		}
		else if(currLen == (**len) && currCapacity>(**resCapacity)){ 
			**resCapacity = currCapacity;  //same length higher residual capacity path is obtained
			for(int i=0;i<=(graph->V+1);i++){
				parent[i] = currParent[i];			//updating parents for final path
			}	
		}
		return;
	}

	//recursive case
	visited[src] = 1;					 //setting visited as 1 
	EDGE *curr = graph->H[src].p;
	int tempCapacity,residual;
	while(curr != NULL){
		residual = curr->c - curr->f;				//residual capacity
		if(!visited[curr->y] && residual > 0){		//visit edges that are not visited and have residual capacity>0
			tempCapacity = min(currCapacity, curr->c);
			currParent[curr->y] = src;
			dfs(graph, curr->y, dest, currLen+1, resCapacity, tempCapacity, parent, currParent, len, visited);
		}
		curr = curr->next;							//moving to other edge in the linked list
	}
	visited[src] = 0;					//setting visited as 0
}

void dfsUtil(GRAPH *graph, int src, int dest, int *resCapacity, int parent[], int *len){

	int V=graph->V;
	int visited[V+2];
	int currParent[V+2];

	for(int i=0;i<=V+1;i++){
		visited[i] = 0;
		currParent[i] = -1;
	}
	//calling dfs for checking all possible paths between src and dest
	dfs(graph, src, dest, 0, &resCapacity, INT_MAX, parent, currParent, &len, visited);
}

//--------------------------------------------------------------------

struct QNode { 				//linked list node to store a queue entry 
    int key;
    struct QNode* next; 
}; 
  
struct Queue { 
    struct QNode *front, *rear; //front and rear end of queue
}; 
  
struct QNode* newNode(int key) 	//creating new linked list node
{ 
    struct QNode* temp = (struct QNode*)malloc(sizeof(struct QNode)); 
    temp->key = key; 
    temp->next = NULL; 
    return temp; 
} 
  
struct Queue* createQueue() 	//creating empty queue
{ 
    struct Queue* q = (struct Queue*)malloc(sizeof(struct Queue)); 
    q->front = q->rear = NULL; 
    return q; 
} 
  
void enQueue(struct Queue* q, int key) 		//pushing value into queue
{ 
    struct QNode* temp = newNode(key); 		//creating a new node
  
    if (q->rear == NULL) { 		//If queue is empty, then new node is front and rear both 
        q->front = q->rear = temp; 
        return; 
    } 
    q->rear->next = temp; 
    q->rear = temp; 
} 

bool isEmpty(struct Queue* q){     //check if queue is empty
	if(q->front == NULL)
		return true;
	return false;
} 

void deQueue(struct Queue* q) 				// remove a key from queue 
{ 
    if (isEmpty(q)) 				// If queue is empty, return NULL 
        return; 
 
    struct QNode* temp = q->front; 	//Store previous front and move front one node ahead 
    q->front = q->front->next; 
    if (q->front == NULL) 		// If front becomes NULL, then change rear also as NULL 
        q->rear = NULL; 
    free(temp); 
} 

int front(struct Queue* q){		   //return front element of queue
	return q->front->key;
}

//Performing bfs for finding all possible augmenting paths 
//Finally storing smallest length augmenting path with maximum residual capacity
void bfs(GRAPH *graph, int src, int dest,int *resCapacity,int parent[],int *len){
	int V = graph->V;
	int *visited;
	int *pathLen;
	int *resCap;

	visited = (int *)malloc((V+2)*sizeof(int));
	pathLen = (int *)malloc((V+2)*sizeof(int));
	resCap = (int *)malloc((V+2)*sizeof(int));

	EDGE *curr;
	int residual,vertex;
	for(int i=0;i<=V+1;i++){
		visited[i] = 0;
		pathLen[i] = 0;
		resCap[i] = 0;
	}

	struct Queue* q = createQueue();	//creating queue for performing bfs

	resCap[src] = INT_MAX;				//intialising arrays for source
	pathLen[src] = 0;
	enQueue(q, src);
	visited[src] = 1;					 //setting visited as 1 for source

	while(!isEmpty(q)){			
		vertex = front(q);
		deQueue(q);
		if(vertex == dest){				//checking if destination is reached
			*len = pathLen[dest];
			*resCapacity = resCap[dest];
			free(visited);
			free(pathLen);
			free(resCap);
			return;
		}
		curr = graph->H[vertex].p;
		while(curr != NULL){
			residual = curr->c - curr->f;				//residual capacity
			if(residual > 0){				//visit edges that have residual capacity>0
				if(!visited[curr->y]){		
					parent[curr->y] = vertex;
					visited[curr->y] = 1;
					resCap[curr->y] = min(residual,resCap[vertex]);
					pathLen[curr->y] = pathLen[vertex] + 1;
					enQueue(q, curr->y);
				}
				else{
					if(pathLen[vertex]+1 < pathLen[curr->y]){
						parent[curr->y] = vertex;
						resCap[curr->y] = min(residual,resCap[vertex]);
						pathLen[curr->y] = pathLen[vertex] + 1;
						enQueue(q, curr->y);
					}
					else if(pathLen[vertex]+1 == pathLen[curr->y]){
						parent[curr->y] = vertex;
						resCap[curr->y] = max(resCap[curr->y],min(residual,resCap[vertex]));
						enQueue(q, curr->y);
					}
				}
			}
			
			curr = curr->next;							//moving to other edge in the linked list
		}
	}
	free(visited);
	free(pathLen);
	free(resCap);
}
//Finding shortest path augmenting path of maximum residual capacity between src and dest
bool findShortestPath(GRAPH *graph, int src, int dest, int parent[], int *path_flow){

	int len = INT_MAX,resCapacity = INT_MIN;

	//calling dfs for checking all possible paths between src and dest
	dfsUtil(graph, src, dest, &resCapacity, parent, &len);

	//calling bfs for checking all possible paths between src and dest
	//bfs(graph, src, dest, &resCapacity, parent, &len);
	
	if(len == INT_MAX)		//no path is present between src and dest
		return false;

	*path_flow = resCapacity;	//storing residual capacity required for augmenting the obtained path
	return true;
}

//updating flow values of edge between u and v with weight
void updateResidualGraph(GRAPH *graph, int u, int v, int weight){
	EDGE *curr = graph->H[u].p;

	while(curr != NULL){
		if(curr->y == v){
			curr->f = curr->f + weight;
			return;
		}
		curr = curr->next;
	}
}
//Ford Fulkerson method for computing max flow between source and destination
void ComputeMaxFlow(GRAPH *graph, int s, int d){

	if(graph == NULL)		//sanity check
		return;

	int V = graph->V;
	int parent[V+2];
	int max_flow = 0,u,v,path_flow;

	//finding shortest path with maximum residual capacity
	while(findShortestPath(graph, s, d, parent, &path_flow)){  
        // update residual capacities of the edges and reverse edges along the path 
        for (v=d; v != s; v=parent[v]){        //augmenting path
            u = parent[v];
            updateResidualGraph(graph, u, v, path_flow);
            updateResidualGraph(graph, v, u, -1*path_flow); 
        } 
        // Add path flow to overall flow 
        max_flow += path_flow;
	}
	printf("Maximum flow obtained : %d\n",max_flow);
}
//finding need based flow using Max-flow algorithm
void NeedBasedFlow(GRAPH *graph){
	int V = graph->V;

	int totProduction = 0,totConsumption = 0,maxFlow = 0;

	for(int i=1;i<=V;i++){			//creating modified graph with single source and destination
		if(graph->H[i].n < 0){							//checking procedure
			addEdge(graph,0,i,-1*(graph->H[i].n));		//directed edge from source:0 to i
			totProduction = totProduction + (-1*(graph->H[i].n));  //updating total production
		}
		else if(graph->H[i].n > 0){						//checking consumers
			addEdge(graph,i,V+1,graph->H[i].n);			//directed edge from i to destination:V+1
			totConsumption = totConsumption + graph->H[i].n;	   //updating total totConsumption
		}
	}

	ComputeMaxFlow(graph,0,V+1);	//computing max flow between source and sink

	EDGE *curr = graph->H[0].p;	
	while(curr!=NULL){				//storing max flow
		maxFlow+=curr->f;
		curr = curr->next;
	}
	//checking if need-based flow constraints are satisfied
	if(maxFlow == totConsumption && maxFlow == totProduction)
		printf("It is possible to assign a need based flow to the graph\n");
	else
		printf("It is not possible to assign a need based flow to the graph\n");
}
//Main() function
int main(){
	char fname[20];
	printf("Enter the file name : ");
	scanf("%s",fname);
	
	GRAPH *graph = ReadGraph(fname);
	PrintGraph(*graph);
	
	printf("\n\n----------------------------------------\n\n");
	int src,dest;
	printf("Enter the index of source : ");
	scanf("%d",&src);
	printf("Enter the index of sink : ");
	scanf("%d",&dest);
	ComputeMaxFlow(graph,src,dest);
	printf("Graph after computing max flow with src:%d and sink:%d is \n",src,dest);
	PrintGraph(*graph);

	printf("\n\n----------------------------------------\n\n");
	GRAPH *new_graph = ReadGraph(fname);
	NeedBasedFlow(new_graph);
	printf("Graph after the need based flow computation: \n");
	PrintGraph(*new_graph);

	return 0;
}