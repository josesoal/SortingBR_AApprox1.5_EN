/*
 ============================================================================
 Project     : 1.5 Approximation Algorithm for sorting by reversals
 Author      : Jose Luis Soncco Alvarez and Mauricio Ayala Rincon
               Universidade de Brasilia (UnB) - Brazil
 Version     : 1.0
 
 Last Modified On: January 23, 2012.
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "breakpointGraph.h"

void buildBreakpointGraph(const permutation *p, \
	breakpointGraph *b)
{
	int i;
	
	memset(b->nblackAdj, 0, sizeof(b->nblackAdj));
    memset(b->ngreyAdj, 0, sizeof(b->ngreyAdj));
    memset(b->blackAdj, 0, sizeof(b->blackAdj));
    memset(b->greyAdj, 0, sizeof(b->greyAdj));
	
	//Create black edges
	for(i=0;i<=p->size;i++){ 
		if (abs(p->pi[i] - p->pi[i+1]) > 1){
			b->blackAdj[i][b->nblackAdj[i]] = i+1; 
			b->blackAdj[i+1][b->nblackAdj[i+1]] = i; 
			b->nblackAdj[i]++;
			b->nblackAdj[i+1]++;
		}			
	}
	
	//Create grey edges 
	for(i=0;i<=p->size;i++){ 
		if (abs(p->piInverse[i] - p->piInverse[i+1]) > 1){
			b->greyAdj[p->piInverse[i]][b->ngreyAdj[p->piInverse[i]]] = p->piInverse[i+1]; 
			b->greyAdj[p->piInverse[i+1]][b->ngreyAdj[p->piInverse[i+1]]] = p->piInverse[i]; 
			b->ngreyAdj[p->piInverse[i]]++;
			b->ngreyAdj[p->piInverse[i+1]]++;			
		}
	}
}

