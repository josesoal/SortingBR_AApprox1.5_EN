/*
 ============================================================================
 Project     : 1.5 Approximation Algorithm for sorting by reversals
 Author      : Jose Luis Soncco Alvarez and Mauricio Ayala Rincon
               Universidade de Brasilia (UnB) - Brazil
 Version     : 1.0
 
 Last Modified On: January 23, 2012.
 ============================================================================
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "eliminate.h"
#include "reversalGraph.h"
#include "cycleDecomposition.h"
#include "breakpointGraph.h"
#include "main.h"

/**
* Return the sequence of reversals
**/
void eliminateVertices(permutation *p, reversals* revs){
	int i,j,k;
    int r=0;
    int verticesComponent[MAX];
    char marcado2[MAX];
    int answer;    
    
    reversalGraph copy;    
    revs->vertexReversal = (int*)malloc(sizeof(int)*MAX);//(int*) para g++
    revs->intervals = (interval*)malloc(sizeof(interval)*MAX);//(interval*) para g++
    
    int isolatedVertices = -1;
    int size = 0;
    
    /* while vertices are not isolated  */
    while(isolatedVertices < size){
                    
        breakpointGraph bpg;
        reversalGraph g;	
        unorientedDecomposition D;
        cycleDecomposition C;
            
        p->size -= 2;  
        buildBreakpointGraph(p, &bpg); //--from breakpointGraph.c       
        generateCycleDecomposition(&bpg, p,&D,&C); //--from cycleDecomposition.c
        p->size += 2;   
        buildReversalGraph(&C, p, &g); //--from reversalGraph.c
        
        preProcess(&g);
    
        for(i=0;i<g.size;i++){
            if(g.componentOrientation[i]== UNORIENTED){         
                applyVertexReversal(&g,i,&C);
                revs->vertexReversal[r]=i;        
                updateCycleDecomposition(p,&C,r,i,revs,g.color[i]);
                r++;   
            }
        }
        
        memset(verticesComponent,0,sizeof(verticesComponent));
        
        //take a component with adjacencies
        for(i=0;i<g.size;i++){
            if(g.nadj[i]!=0){
                break;
            }
        }
        
        //recover all vertices that are in the same component of vertex i
        for(j=0;j<g.size;j++){
            if(g.component[j] == g.component[i]){
                verticesComponent[j] = 1;
            }
        }
        
        //for all vertex in the same component of vertex i     
        for(j=0;j<g.size;j++){
            if(g.component[j]==g.component[i] && g.color[j]==RED){
                copy = g;
                applyVertexReversal(&g,j,&C);
                answer = ORIENTED;
                for(k=0;k<g.size;k++){
                    if(verticesComponent[k]==1){
                        memset(marcado2,0,sizeof(marcado2));
                        if(verifyOrientation(&g,k,marcado2)==UNORIENTED){
                            answer = UNORIENTED;
                            break;
                        }
                    }
                }
                if(answer==ORIENTED){							
                    revs->vertexReversal[r]=j;
                    updateCycleDecomposition(p,&C,r,j,revs,RED);						
                    r++;
                    break;
                }
                else{
                    g = copy;
                }
            }
        }//end for
            
        isolatedVertices = g.isolatedVertices; 
        size = g.size; 
                   
        
    }//fin while
    revs->size = r;
}

void updateCycleDecomposition(permutation *p, cycleDecomposition *C,int reversalPosition, int vertice, 
                                reversals* revs,int color){
	
	int j,current,grey1, grey2, temp, black1, black2, right,twoCycle, a, b, v;
	
    /* update the cycle decomposition */
	if (color == RED){
        /* find to which cycle the vertex belongs to */
		grey1 = p->piInverse[vertice];
		grey2 = p->piInverse[vertice+1];
		
		current = -1;//initializing
        
        for(j=0; j<C->nCycles; j++){
            if (C->cycle[j].vertexEdge[vertice][0] == vertice+1){
                current = j;
                break;
            }
        }
        
        /* verifiy if black edges are at the right or left of the grey edge */	
		if (grey2 < grey1){
			temp = grey1;
			grey1 = grey2;
			grey2 = temp;
		}		
		
		black1 = chooseBlackEdge(grey1, grey2, current, C,p);
		black2 = chooseBlackEdge(grey2, grey1, current, C,p);		
		
		right = 0;//is the left
        
		if ((black1 > grey1)&&(black2 > grey2))
			right = 1;//is the right

        /* verify if the current cycle is a 2-cycle */
		twoCycle=0;		
        if (p->pi[black1] < p->pi[black2]){
            if (p->pi[black1]+1 == p->pi[black2])
                twoCycle =1;            
        }
        else{//p->pi[black1] > p->pi[black2]
            if (p->pi[black2]+1 == p->pi[black1])
                twoCycle =1;
        }
        
        /* generate the reversal interval (between two black edgees) */
		if (right){
			a = black1;
			b = grey2;
		}
		else{//left
			a = grey1;
			b = black2;
		}	
        
        /* MODIFY CURRENT CYCLE DECOMPOSITION */
        
        //eliminate grey edge (grey1-grey2)
        if (p->pi[grey1] < p->pi[grey2]){
            C->cycle[current].vertexEdge[p->pi[grey1]][0] = -1;
            C->cycle[current].vertexEdge[p->pi[grey1]][1] = -1;
            C->cycle[current].vertexEdge[p->pi[grey1]][2] = -1;
        }            
        else{//p->pi[grey1] > p->pi[grey2]
            C->cycle[current].vertexEdge[p->pi[grey2]][0] = -1;
            C->cycle[current].vertexEdge[p->pi[grey2]][1] = -1;
            C->cycle[current].vertexEdge[p->pi[grey2]][2] = -1;
        }

        //update the vertex v that contains the black edge black2-grey2 with the black edge black2-black1
        v = p->pi[black2];
        if (v != p->size-1 && C->cycle[current].vertexEdge[v][0] == v+1){//verify if black2 represents a v vertex
            if (C->cycle[current].vertexEdge[v][1] == p->pi[grey2]){//verify if vertex v (black2) is related with black-edge grey2                 
                C->cycle[current].vertexEdge[v][1] = p->pi[black1];//update vertex v with the black-edge black2-black1
            }
        }
        
        if (v != 0 && C->cycle[current].vertexEdge[v-1][0] == v){//verify if black2 is part of vertex v-1
            if (C->cycle[current].vertexEdge[v-1][2] == p->pi[grey2]){//verify if vertex v (black2) is related with black-edge grey2
                
                C->cycle[current].vertexEdge[v-1][2] = p->pi[black1];//update vertex v with a black-edge black2-black1
            }
        }
        
        //update the vertex v that contains the black-egde black1-grey1 with the black-edge black1-black2 
        v = p->pi[black1];
        if (v != p->size-1 && C->cycle[current].vertexEdge[v][0] == v+1){//verify if black1 represents a v vertex
            if (C->cycle[current].vertexEdge[v][1] == p->pi[grey1]){//verify if vertex v (black1) is related with black edge grey1
                
                C->cycle[current].vertexEdge[v][1] = p->pi[black2];//update the vertex v with a black edge black1-black2  
            }
        }
        
        if (v != 0 && C->cycle[current].vertexEdge[v-1][0] == v){//verify if black1 is part of a vertex v-1 
            if (C->cycle[current].vertexEdge[v-1][2] == p->pi[grey1]){//verify if vertex v (black1) is related with black-edge grey1
                
                C->cycle[current].vertexEdge[v-1][2] = p->pi[black2];//update vertex v with a black-edge black1-black2
            }
        }
        
        //if it is a 2-cycle, eliminate the remaining edge (black1-black2)
        if (twoCycle){ 
            if (p->pi[black1] < p->pi[black2]){
                C->cycle[current].vertexEdge[p->pi[black1]][0] = -1;
                C->cycle[current].vertexEdge[p->pi[black1]][1] = -1;
                C->cycle[current].vertexEdge[p->pi[black1]][2] = -1;
            }            
            else{//p->pi[black1] > p->pi[black2]
                C->cycle[current].vertexEdge[p->pi[black2]][0] = -1;
                C->cycle[current].vertexEdge[p->pi[black2]][1] = -1;
                C->cycle[current].vertexEdge[p->pi[black2]][2] = -1;
            }
        }
        
    }
    else{//color == BLUE
        
        /* FIND TO WHICH CYCLE THE VERTEX BELONGS TO */
        
        grey1 = p->piInverse[vertice];
		grey2 = p->piInverse[vertice+1];
		
		current = -1;//initializing 
        
        for(j=0; j<C->nCycles; j++){
            if (C->cycle[j].vertexEdge[vertice][0] == vertice+1){
                current = j;
                break;
            }
        }
        
        /* GENERATE REVERSAL INTERVAL (BETWEEN TWO BLACK-EDGES) */
        
        //put grey1 less than grey2 
        if (grey2 < grey1){
			temp = grey1;
			grey1 = grey2;
			grey2 = temp;
		}
        
		black1 = chooseBlackEdge(grey1, grey2, current, C,p);//returns the position of the black-edge related with grey1
		black2 = chooseBlackEdge(grey2, grey1, current, C,p);//returns the position of the black-edge related with grey2
        
        //reversal interval
        if (black1 < grey1 && black2 > grey2){
            a = grey1;
            b = grey2; 
        }
        else if (black1 > grey1 && black2 < grey2){
            a = black1;
            b = black2;
        }
        
        /* MODIFY THE CURRENT CYCLE DECOMPOSITION */

        //update the two black-edges (g1-b1 and g2-b2) by (g1-b2 and g2-b1), related with vertex v (g1-g2) and v+1
        
        //choose if g1 or g2 is the label of the vertex v
        if (grey1 == p->piInverse[vertice]){
            v = p->pi[grey1];
            
            //update the edges
            C->cycle[current].vertexEdge[v][1] = p->pi[black2];// (g1-b1) by (g1-b2)
            C->cycle[current].vertexEdge[v][2] = p->pi[black1];// (g2-b2) by (g2-b1)
            
        }else{
            v = p->pi[grey2];
            
            //update the edges
            C->cycle[current].vertexEdge[v][1] = p->pi[black1];// (g2-b2) by (g2-b1) 
            C->cycle[current].vertexEdge[v][2] = p->pi[black2];// (g1-b1) by (g1-b2)
            
        }
        
        //update the black-edge (b1-g1) por (b1-g2), related with vertex v (that have as fixed point b1)
        v = p->pi[black1];
        if (v != p->size-1 && C->cycle[current].vertexEdge[v][0] == v+1){//verify if black1 represents a vertex v  
            if (C->cycle[current].vertexEdge[v][1] == p->pi[grey1]){//verify if vertex v (black1) is related with black-edge grey1 

                C->cycle[current].vertexEdge[v][1] = p->pi[grey2];//update vertex v with a black edge black1-grey2   
            }
        }
        
        if (v != 0 && C->cycle[current].vertexEdge[v-1][0] == v){//verify if black1 is part of a vertex v-1 
            if (C->cycle[current].vertexEdge[v-1][2] == p->pi[grey1]){//verify if vertex v (black1) is related with black-edge grey1 
                
                C->cycle[current].vertexEdge[v-1][2] = p->pi[grey2];//update the vertex v with the black-egde black1-grey2 
            }
        }
        
        //update the black-edge (b2-g2) by (b2-g1), related with vertex v (that have as fixed point b2)
        v = p->pi[black2];
        if (v != p->size-1 && C->cycle[current].vertexEdge[v][0] == v+1){//verify if black2 represents a vertex v
            if (C->cycle[current].vertexEdge[v][1] == p->pi[grey2]){//verify if vertex v (black2) is related with black-edge grey2
                
                C->cycle[current].vertexEdge[v][1] = p->pi[grey1];//update vertex v with a black-edge black2-grey1  
            }
        }
        
        if (v != 0 && C->cycle[current].vertexEdge[v-1][0] == v){//verify if black2 is part of a vertex v-1   
            if (C->cycle[current].vertexEdge[v-1][2] == p->pi[grey2]){//verify if vertex v (black2) is related with black-edge grey2
                
                C->cycle[current].vertexEdge[v-1][2] = p->pi[grey1];//update vertex v with a black-edge black2-grey1 
            }
        }
        
        
    } 
	
    revs->intervals[reversalPosition].left = a;
    revs->intervals[reversalPosition].right = b;		
    
	applyReversal(p,a,b);	
    
}

/* apply the reversal corresponding to vertex v and make the modifications
in the reversal graph, return the number of oriented components */
void applyVertexReversal(reversalGraph* g, int v ,cycleDecomposition *C){
    int i,j;
    int component;
    int marcado[MAX];
    /**O(n^2)**/
    /* if we have (v,i) and (v,j), (i,j) is flipped */
    for(i=0;i<g->size;i++){
        for(j=i+1;j<g->size;j++){

            if(g->adj[v][i] == 1 && g->adj[v][j]==1){
                flipEdge(g,i,j,C);
            }
        }
    }

    /* flip the color of all the adjacencies to the reversal represented by v */
    for(i=0;i<g->size;i++){
        if(g->adj[v][i]==1){
            flipColor(g,i);
        }
    }

    /* if red, becomes an isolated vertex */
    if(g->color[v] == RED){
        flipColor(g,v);
        /* isolate the vertex */
        for(i=0;i<g->size;i++){
            g->adj[v][i]=0;
        }
        g->nadj[v]=0;
        for(i=0;i<g->size;i++){
			if(g->adj[i][v]==1){
				g->adj[i][v]=0;
				g->nadj[i]--;
			}
		}
        memset(marcado,0,sizeof(marcado));
        component = g->component[v];
        /* for each vertex of the same component of v, update */
        for(i=0;i<g->size;i++){
            if(g->component[i] == component){
                marcado[i] = 1;
                g->component[i] = i;
            }
        }
        for(i=0;i<g->size;i++){
            if(marcado[i]){
                dfsReversal(g,i,g->component[i]);
            }
        }
    }
    else{//g->color[v] == BLUE
        /* for each vertex of the same component of v, update the orientation (to ORIENTED) */
        component = g->component[v];
        for(i=0;i<g->size;i++){
            if(g->component[i] == component){
                g->componentOrientation[i] = ORIENTED;
            }
        }
        
    }
    
    g->isolatedVertices=0;
    for(i=0;i<g->size;i++){
		if(g->nadj[i]==0){
			g->isolatedVertices++;
		}
	}	
}

void flipEdge(reversalGraph* g,int u,int v,cycleDecomposition *C){
    int j, current, specialCase;
    
    //find to which cycle the vertex belongs to
    current = -1;//initializing
    for(j=0; j<C->nCycles; j++){
        if (C->cycle[j].vertexEdge[u][0] == u+1){
            current = j;
            break;
        }
    }
    
    //verify the special case of lemma 4.1
    specialCase = 0;
    if (C->cycle[current].vertexEdge[v][0] == v+1){//if the other vertex is in the same cycle 
        if (u < v){
            if (u+1 == v && // verify if they are consecutives and 
                C->cycle[current].vertexEdge[u][1] == v+1 && // exists a black edge between u and v+1, and  
                C->cycle[current].vertexEdge[v][2] == u){ // exists a black edge between v+1 and u                
                specialCase = 1;
            }
        }
        else{ //u > v
            if (v+1 == u && // verify if they are consecutives and
                C->cycle[current].vertexEdge[v][1] == u+1 && // exists a black edge between v and u+1, and   
                C->cycle[current].vertexEdge[u][2] == v){ // exists a black edge between u+1 and v 
                specialCase = 1;
            }
        }
    }
     
    //flip the edges if it is not the special case
    if (!specialCase){    
        if(g->adj[u][v]==1){
            g->adj[u][v]=0;
            g->adj[v][u]=0;
            g->nadj[u]--;
            g->nadj[v]--;
        }
        else{
            g->adj[u][v]=1;
            g->adj[v][u]=1;
            g->nadj[u]++;
            g->nadj[v]++;
        }
    }
}

void flipColor(reversalGraph* g,int v){
    g->color[v] = g->color[v]^1; //invert from 0 to 1 or 1 to 0. 
}

/* apply a depth-first-search and label the vertices that 
    belongs to the same connected component of the leader */
void dfsReversal(reversalGraph* g,int v,int lider){
    int i;
    g->component[v] = lider;
    for(i=0;i<g->size;i++){
        if(g->adj[v][i]==1 && g->component[i]!=lider){
            dfsReversal(g,i,lider);
        }
    }
}

int verifyOrientation(reversalGraph *g,int v,char marcado[]){
	int i;

	marcado[v]=1;
    if(g->nadj[v]==0){
        return(ORIENTED);
    }
    else{
        if(g->color[v]==RED){
            return(ORIENTED);
        }
        for(i=0;i<g->size;i++){
            if(g->adj[v][i] && marcado[i]==0){
                if(verifyOrientation(g,i,marcado)==ORIENTED){
                    return(ORIENTED);
                }
            }
        }
        return(UNORIENTED);
    }
}

/* pre-process the connected and not-connected components  */
void preProcess(reversalGraph* g){
    int i,j;
    int orientation;
    char marcado[MAX];
    char marcado2[MAX];
    for(i=0;i<g->size;i++){
        g->component[i]=i;
    }

    /* indicates the connected components */
    for(i=0;i<g->size;i++){
        dfsReversal(g,i,g->component[i]);
    }

	memset(marcado,0,sizeof(marcado));
    for(i=0;i<g->size;i++){
		if(!marcado[i]){
			marcado[i]=1;
			memset(marcado2,0,sizeof(marcado2));
			orientation = verifyOrientation(g,i,marcado2);
			g->componentOrientation[i] = orientation;
			for(j=0;j<g->size;j++){
				if(g->component[j] == g->component[i]){
					g->componentOrientation[j] = g->componentOrientation[i];
					marcado[j]=1;
				}
			}
		}
    }

    /* count the number of isolated blue vertices */
    g->isolatedVertices=0;
    for(i=0;i<g->size;i++){
        if(g->nadj[i]==0){
            g->isolatedVertices++;
        }
    }
}

void applyReversal(permutation* p,int a,int b){
	int aux,i;
	for(i=0;i<(b-a+1)/2;i++){ 
        aux = p->pi[a+i]; 
        p->pi[a+i] = p->pi[b-i]; 
        p->pi[b-i] = aux; 

		p->piInverse[p->pi[a+i]] = a+i;
        p->piInverse[p->pi[b-i]] = b-i;      		
    }
}

//calculate the difference of breakpoint if the reversal of vertex v BLUE would be applied
int difBreakpoints(permutation *p, int v){
	int before = 0, after = 0, pos_a, pos_b, temp;
	
    pos_a = p->piInverse[v];
    pos_b = p->piInverse[v+1];
    
	if (pos_a > pos_b){
        temp = pos_b;
        pos_b = pos_a;
        pos_a = temp;
    }	
	
	//numero de breakpoints before
	if (abs(p->pi[pos_a - 1] - p->pi[pos_a]) > 1)
		before++;
	if (abs(p->pi[pos_b] - p->pi[pos_b + 1]) > 1)
		before++;	
	//numero de breakpoints after
	if (abs(p->pi[pos_a - 1] - p->pi[pos_b]) > 1)
		after++;
	if (abs(p->pi[pos_a] - p->pi[pos_b + 1]) > 1)
		after++;	
	//diferencia de breakpoint
	return after - before;	
}





