/*
 ============================================================================
 Project     : 1.5 Approximation Algorithm for sorting by reversals
 Author      : Jose Luis Soncco Alvarez and Mauricio Ayala Rincon
               Universidade de Brasilia (UnB) - Brazil
 Version     : 1.0
 
 Last Modified On: January 23, 2012.
 ============================================================================
 */

#ifndef _ELIMINACAO_H_
#define _ELIMINACAO_H_

#include "graphs.h"

void eliminateVertices(permutation *p,reversals* revs);
void updateCycleDecomposition(permutation* p,cycleDecomposition *C,int reversalPosition, int vertice,reversals* revs, int color);
void applyVertexReversal(reversalGraph* g, int v,cycleDecomposition *C);

void flipEdge(reversalGraph* g,int u,int v,cycleDecomposition *C);
void flipColor(reversalGraph* g,int v);
void dfsReversal(reversalGraph* g,int v,int lider);
int verifyOrientation(reversalGraph *g,int v,char marcado[]);
void preProcess(reversalGraph* g);

void applyReversal(permutation* p,int a,int b);
int difBreakpoints(permutation *p, int v);
#endif
