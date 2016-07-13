/*
 ============================================================================
 Project     : 1.5 Approximation Algorithm for sorting by reversals
 Author      : Jose Luis Soncco Alvarez and Mauricio Ayala Rincon
               Universidade de Brasilia (UnB) - Brazil
 Version     : 1.0
 
 Last Modified On: January 23, 2012.
 ============================================================================
 */

#ifndef _CYCLE_DECOMPOSITION_H_
#define _CYCLE_DECOMPOSITION_H_

#include "graphs.h"

void generateCycleDecomposition(breakpointGraph *b, permutation *p, unorientedDecomposition *D, cycleDecomposition *C);
int checkExistence2Cycle(twoCycle c[],int n, twoCycle aux);
void makeCyclesOriented(unorientedDecomposition *D, permutation *p, cycleDecomposition *C);
void eliminateBlackAdj(breakpointGraph *cycle,int black1, int black2);
void eliminateGreyAdj(breakpointGraph *cycle,int grey1, int grey2);

void showBreakpointGraph(orientedCycle *b,permutation *p);
#endif
