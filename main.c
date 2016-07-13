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
#include <time.h>
#include <sys/time.h>//for measuring time
#include "main.h"
#include "breakpointGraph.h"
#include "cycleDecomposition.h"
#include "reversalGraph.h"
#include "eliminate.h"

int main(int argc, char** argv)
{	
    struct timeval t_ini, t_fin;
    double msecs;
    
    gettimeofday(&t_ini, NULL);//-------------Takes initial TIME
    
	permutation p;
	reversals rev;	
	int i;

	/* Read length of permutation */
    scanf("%d",&p.size);

	/* Read elements of the permutation */
    p.pi[0]=0;
    p.pi[p.size+1]=p.size+1;
    for(i=1;i<=p.size;i++){
        scanf("%d",&p.pi[i]);
    }    

    /* Generate a inverse permutation */
    generateInversePermutation(&p);
    
    /* Run main algorithm: elimination of vertices in a reversal graph */
    p.size += 2;
    eliminateVertices(&p,&rev);//--from eliminate.c
    p.size -= 2;
    
    gettimeofday(&t_fin, NULL);//-------------Takes final TIME    
    msecs = timeval_diff(&t_fin, &t_ini);//Calculate total time
    
    /* Show the number of reversals found */        
    showResults(&p,&rev,msecs);
    
    return 0;
}

void showResults(permutation *p, reversals *rev, double msecs){
    printf("\n1.5 Approximation Algorithm for Sorting Unsigned Permutations\n");
    printf("\tNumber of Reversals: %d\n",rev->size);
    printf("\tTotal Time: %.16g miliseconds\n", msecs * 1000.0);
    printf("\n");
}


/*void showAllResults(permutation *p, reversals *rev, double secs){
	int i;    
	printf("Final permutation:"); showPermutation(p->pi, p->size);
	printf("Number of breakpoints:%d\n",numberBreakpoints(p));
	printf("Number of reversals: %d\n",rev->size);
	printf("Reversals(positions):"); 
	for(i=0; i<rev->size-1; i++)
		printf("(%d,%d),",rev->intervals[i].left,rev->intervals[i].right);
	printf("(%d,%d)\n",rev->intervals[i].left,rev->intervals[i].right);
	printf("Total time: %.16g milliseconds\n", secs * 1000.0);
    printf("\n");
}*/

double timeval_diff(struct timeval *a, struct timeval *b)
{
    return
    (double)(a->tv_sec + (double)a->tv_usec/1000000) -
    (double)(b->tv_sec + (double)b->tv_usec/1000000);
}

int numberBreakpoints(permutation *p){
	int i,nroBp=0;
	for(i=0; i<=p->size; i++){
		if (abs(p->pi[i] - p->pi[i+1]) > 1)
			nroBp++;		
	}
	return nroBp;
}

void showPermutation(const int *p, int size){
int i;
	for (i=0 ; i<= size+1; i++)
		printf("%d ", p[i]);
	printf("\n");
}

void generateInversePermutation(permutation *p){
int i;
	for (i=0 ; i<= p->size+1; i++)
		p->piInverse[p->pi[i]]=i;
}

void showReversalGraph(const reversalGraph *revGraph){
    int i, j;
    
    printf("\nREVERSAL GRAPH R(C)\n");
    
    for(i=0;i<=revGraph->size;i++){
		printf("Color of v[%d] = ",i);
		if(revGraph->color[i]==RED){
			printf("RED\n");
		}
		else{
			printf("BLUE\n");
		}
	}	
    
    printf("\t Adjacency Matrix:\n");
    printf("\t");
    for(j=0; j<= revGraph->size; j++){
		printf("%d\t",j);
	}
	printf("\n");
    for(i=0; i <= revGraph->size; i++){
		printf("%d\t",i);
		for(j=0; j<= revGraph->size; j++){
			printf("%d\t",revGraph->adj[i][j]);
		}
		printf("\n");
	} 
    
    printf("Number of vertices:%d\n", revGraph->size);
    printf("Number of isolated vertices:%d\n", revGraph->isolatedVertices);
    
    printf("Degree of incidence of each vertex:\n");
    for(i=0;i<=revGraph->size;i++){
		printf("Vertex %d with degree %d\n",i, revGraph->nadj[i]);
	}
}



