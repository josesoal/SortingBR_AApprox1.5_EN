/*
 ============================================================================
 Project     : 1.5 Approximation Algorithm for sorting by reversals
 Author      : Jose Luis Soncco Alvarez and Mauricio Ayala Rincon
               Universidade de Brasilia (UnB) - Brazil
 Version     : 1.0
 
 Last Modified On: January 23, 2012.
 ============================================================================
 */

#include <stdlib.h>
#include <string.h>
#include "reversalGraph.h"

void buildReversalGraph(cycleDecomposition *C, permutation *p, reversalGraph *revGraph){
    int i, j, k, f, g;

	/* initialize variables */
	revGraph->isolatedVertices=0;
	revGraph->size=0;
	memset(revGraph->adj,0,sizeof(revGraph->adj));
	memset(revGraph->nadj,0,sizeof(revGraph->nadj));

	/* generate isolated vertices */
	for(i=0; i<p->size-1; i++){
		if (abs(p->pi[i] - p->pi[i+1]) == 1){//it is not a black edge
			revGraph->color[p->pi[i]] = BLUE;
			revGraph->isolatedVertices++;
		}
	}	
	
	/* generate red and blue vertices */	
    for(i=0; i<C->nCycles; i++){
        for(j=0; j<p->size-1; j++){
            if (C->cycle[i].vertexEdge[j][0] == j+1)
                calculateVertexColor(j, i, C, p, revGraph);
        }    
    }

	revGraph->size += revGraph->isolatedVertices;
   
   	/* generate edges */ 	
   	//compare inside the cycles	
    for(i=0; i<C->nCycles; i++){				
        
		for(j=0; j<p->size-2; j++){
			if (C->cycle[i].vertexEdge[j][0] == j+1){//if vertex j forms an edge                 
				
                for(k=j+1; k<p->size-1; k++){
                    if (C->cycle[i].vertexEdge[k][0] == k+1){//if vertex j forms an edge 
                        createEdge(j, k, i, i, C, p, revGraph);
                        
                    }
                }
                
            }
        }
    }
    
	//compare among cycles
    for(f=0; f<C->nCycles-1; f++){	
		for(g=f+1; g<C->nCycles; g++){
            
            for(j=0; j<p->size-1; j++){
                if (C->cycle[f].vertexEdge[j][0] == j+1){//if vertex j forms an edge               
                    
                    for(k=0; k<p->size-1; k++){
                        if (C->cycle[g].vertexEdge[k][0] == k+1){//if vertex k forms an edge 
                            
                            createEdge(j, k, f, g, C, p, revGraph);
                            
                        }
                    }
                    
                }
            }            
            
        }
    }
       
	
}


void createEdge(int vertice_i,int vertice_j, 
        int i, int j,cycleDecomposition *C,permutation *p, reversalGraph *revGraph){ 
     
    int a,b,lgu,rgu,lgv,rgv,lbu,rbu,lbv,rbv,blackAdj;
	//compare if lg(u)<lg(v)<rg(u)<rg(v) o lg(v)<lg(u)<rg(v)<rg(u)
    if (p->piInverse[vertice_i] < p->piInverse[vertice_i+1]){
        lgu = p->piInverse[vertice_i]; 
        rgu = p->piInverse[vertice_i+1]; 
    }
    else{
        lgu = p->piInverse[vertice_i+1]; 
        rgu = p->piInverse[vertice_i]; 
    }
    
    if (p->piInverse[vertice_j] < p->piInverse[vertice_j+1]){
        lgv = p->piInverse[vertice_j]; 
        rgv = p->piInverse[vertice_j+1]; 
    }
    else{
        lgv = p->piInverse[vertice_j+1]; 
        rgv = p->piInverse[vertice_j]; 
    }
						
		
	
	if (((lgu<lgv)&&(lgv<rgu)&&(rgu<rgv))||
		((lgv<lgu)&&(lgu<rgv)&&(rgv<rgu)))
	{
		//create edge
		a = rgu;
		if (p->pi[lgu] < p->pi[rgu]) a = lgu;							
		b = rgv;
		if (p->pi[lgv] < p->pi[rgv]) b = lgv;
								
		revGraph->adj[p->pi[a]][p->pi[b]] = 1;
		revGraph->adj[p->pi[b]][p->pi[a]] = 1;
		
		revGraph->nadj[p->pi[a]]++;
		revGraph->nadj[p->pi[b]]++;

		//printf("lgu=%d,rgu=%d,lgv=%d,rgv=%d",lgu,rgu,lgv,rgv);
	}				
	else{ 
		//comparar se lb(u)<lb(v)<rb(u)<rb(v) o lb(v)<lb(u)<rb(v)<rb(u)
		blackAdj = chooseBlackEdge(lgu, rgu, i, C,p);				
		if (lgu < blackAdj) lbu = lgu;
		else lbu = blackAdj;
		
		blackAdj = chooseBlackEdge(rgu, lgu, i, C,p);					
		if (rgu < blackAdj) rbu = rgu;
		else rbu = blackAdj;
		
		blackAdj = chooseBlackEdge(lgv, rgv, j, C,p);
		if (lgv < blackAdj) lbv = lgv;
		else lbv = blackAdj;
		
		blackAdj = chooseBlackEdge(rgv, lgv, j, C,p); 
		if (rgv < blackAdj) rbv = rgv;
		else rbv = blackAdj;
		
		if (((lbu<lbv)&&(lbv<rbu)&&(rbu<rbv))||
			((lbv<lbu)&&(lbu<rbv)&&(rbv<rbu)))
		{
			//create edge								
			a = rgu;
			if (p->pi[lgu] < p->pi[rgu]) a = lgu;							
			b = rgv;
			if (p->pi[lgv] < p->pi[rgv]) b = lgv;
									
			revGraph->adj[p->pi[a]][p->pi[b]] = 1;
			revGraph->adj[p->pi[b]][p->pi[a]] = 1;
			
			revGraph->nadj[p->pi[a]]++;
			revGraph->nadj[p->pi[b]]++;
		}															
	}	
}

void calculateVertexColor(int vertice, int i, cycleDecomposition *C, permutation *p, reversalGraph *revGraph){
	int grey1,grey2,seta1,seta2,pos,blackAdj1, blackAdj2;
    grey1 = p->piInverse[vertice];
    grey2 = p->piInverse[vertice+1];
	
    /* Note: head = 1, and tail = 0 */

	seta1 = 0;//tail
	seta2 = 0;//tail
	
	blackAdj1 = chooseBlackEdge(grey1,grey2,i, C,p);//para grey1
	blackAdj2 = chooseBlackEdge(grey2,grey1,i, C,p);//para grey2
	
	if (grey1 > blackAdj1)
		seta1 = 1;//head
	if (grey2 > blackAdj2)
		seta2 = 1;//head
	
	if (seta1 != seta2){//head with tail
		pos = grey2;
		if (p->pi[grey1] < p->pi[grey2]) pos = grey1;
		
		revGraph->color[p->pi[pos]] = BLUE;
	}
	else{//head with head, or tail with tail
		pos = grey2;
		if (p->pi[grey1] < p->pi[grey2]) pos = grey1;
		
		revGraph->color[p->pi[pos]] = RED;
	}

	revGraph->size++; 
}


int chooseBlackEdge(int grey1, int grey2, int i,cycleDecomposition *C,permutation *p){
	
    int blackEdge;
    blackEdge = -1;
    //choose black edge according to orientation of grey edge
    if (p->pi[grey1] < p->pi[grey2]){
        blackEdge = p->piInverse[C->cycle[i].vertexEdge[p->pi[grey1]][1]];//black edge related with grey1
    }
    else{//p->pi[grey1] > p->pi[grey2]
        blackEdge = p->piInverse[C->cycle[i].vertexEdge[p->pi[grey2]][2]];//black edge related with grey2
    }  
    
    return blackEdge;    
}


