/*
 ============================================================================
 Project     : 1.5 Approximation Algorithm for sorting by reversals
 Author      : Jose Luis Soncco Alvarez and Mauricio Ayala Rincon
               Universidade de Brasilia (UnB) - Brazil
 Version     : 1.0
 
 Last Modified On: January 23, 2012.
 ============================================================================
 */
#ifndef _MAIN_H_
#define _MAIN_H_

#include "graphs.h"

void showResults(permutation *p, reversals *rev, double secs);
double timeval_diff(struct timeval *a, struct timeval *b);
int numberBreakpoints(permutation *p);
void showPermutation(const int *p, int size);
void generateInversePermutation(permutation *p);
void showReversalGraph(const reversalGraph *revGraph);
#endif
