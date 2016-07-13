/*
 ============================================================================
 Project     : 1.5 Approximation Algorithm for sorting by reversals
 Author      : Jose Luis Soncco Alvarez and Mauricio Ayala Rincon
               Universidade de Brasilia (UnB) - Brazil
 Version     : 1.0
 
 Last Modified On: January 23, 2012.
 ============================================================================
 */

#ifndef _REVERSAL_GRAPH_H_
#define _REVERSAL_GRAPH_H_

#include "graphs.h"

void buildReversalGraph(cycleDecomposition *C, permutation *p, reversalGraph *revGraph);
void createEdge(int vertice_i,int vertice_j, int i, int j,cycleDecomposition *C,permutation *p, reversalGraph *revGraph);
void calculateVertexColor(int vertice, int i, cycleDecomposition *C, permutation *p, reversalGraph *revGraph);
int chooseBlackEdge(int grey1, int grey2, int i,cycleDecomposition *C,permutation *p);

#endif
