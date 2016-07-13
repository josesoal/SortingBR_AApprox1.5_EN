/*
 ============================================================================
 Project     : 1.5 Approximation Algorithm for sorting by reversals
 Author      : Jose Luis Soncco Alvarez and Mauricio Ayala Rincon
                Universidade de Brasilia (UnB) - Brazil
 Version     : 1.0
 
 Last Modified On: January 23, 2012.
 ============================================================================
 */

#ifndef GRAPHS_H_
#define GRAPHS_H_

#define MAX 202 /** Maximum number of vertices + 2 **/

#define BLUE 0
#define RED 1

#define UNORIENTED 0
#define ORIENTED 1

//#define BLACK 0
//#define GREY 1

typedef struct {  
	int adj[MAX][MAX];             /* adjacency matrix */
	int nadj[MAX];                 /* numbero of adjacencies of each node */
	int color[MAX];                /* color of the vertices of the graph */
	int component[MAX];            /* indicates the leader vertex of each connected component */
	int componentOrientation[MAX]; /* orientation of component */
	int size;                      /* number of vertices of the graph (indexed from 0 to size -1) */
	int isolatedVertices;          /* number of isolated vertices (blue) */
}reversalGraph;

typedef struct {                    /* indices of arrays are the indices of the permutation */
    int blackAdj[MAX][2];           /* adjacency matrix of black-edges */ 
    int greyAdj[MAX][2];            /* adjacency matrix of grey-edges */ 
    int nblackAdj[MAX];             /* number of black adjacencies */ 
    int ngreyAdj[MAX];              /* number of grey adjacencies */      
}breakpointGraph;

typedef struct {
    breakpointGraph cycle[MAX];   
    int nCycles; 
}unorientedDecomposition; 

typedef struct {
    /* matrix linking grey and black edges (save values and not positions) */ 
    int vertexEdge[MAX][3]; 
    /* the index represents the vertex (the grey edge) */ 
    /* column 1 represents the vertex + 1 */ 
    /* column 2 represents the black edge related with the vertex */ 
    /* column 3 represents the black edge related with the vertex + 1  */ 
}orientedCycle; 

typedef struct {
    orientedCycle cycle[MAX];
    int nCycles;
}cycleDecomposition; 

typedef struct {
    int size;                   /* size of permutation */
    int pi[MAX];                /* permutation */
    int piInverse[MAX];         /* inverse permutation */
}permutation;

typedef struct{
	int left;
	int right;
}interval; 

typedef struct {
    int p1,p2,p3,p4;    /* represent black edges (p1<p2 and p3<p4) */
    int c1,c2,c3,c4;    /* represent grey edges (c1<c2 and c3<c4) */
    /* p1 is not necessarily < than p3 */
    /* c1 is not necessarily < than c3 */
}twoCycle;

typedef struct{
	int* vertexReversal; 
	interval* intervals;
	int size;
}reversals;

#endif
