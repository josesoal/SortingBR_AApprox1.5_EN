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
#include "cycleDecomposition.h"

#include <string>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <cassert>
#include <boost/graph/max_cardinality_matching.hpp>
using namespace boost;
typedef adjacency_list<vecS, vecS, undirectedS> my_graph; 


void generateCycleDecomposition(breakpointGraph *b, 
		permutation *p, unorientedDecomposition *D, cycleDecomposition *C)
{
	int i, j, k, l, blackAdj, greyAdj1, greyAdj2,temp;
	
	int num2cycles=0;
	twoCycle cycles2[MAX];
	int num2cyclesNoBlackEdges=0;
	int cycle2NoBlackEdges[MAX];//save the indices of 2cycles not sharing black edges
	int num2cyclesNoGreyEdges=0;
	int cycle2NoGreyEdges[MAX];//save the indices of 2cycles not sharing grey edges        
	twoCycle auxCiclo2;    
    
    /* Perform the decomposition compound of just 2-cycles */
	//for each vertex, verify if it is part of a 2-cycle
    for(i=0;i<=p->size+1;i++){
    	//for each black edge starting at node i
        for(j=0;j<b->nblackAdj[i];j++){//for each black adjacency of node i
           blackAdj = b->blackAdj[i][j];//black adjacency j of node i

           for(k=0;k<b->ngreyAdj[i];k++){//for each grey adjacency of node i  
               for(l=0;l<b->ngreyAdj[blackAdj];l++){//for each grey adjacency of blackAdj
                   greyAdj1 = b->greyAdj[i][k];// grey adjacency k of node i
                   greyAdj2 = b->greyAdj[blackAdj][l];// grey adjacency l of node blackAdj
                   
                   if (abs(greyAdj1-greyAdj2)==1){ /* if are adjacents */
					   if (((b->nblackAdj[greyAdj1] == 1)&&(b->blackAdj[greyAdj1][0] == greyAdj2))||
						   ((b->nblackAdj[greyAdj1] == 2)&&(b->blackAdj[greyAdj1][1] == greyAdj2))){					   
						  
					   	   //check if this 2-cycle already exists in the set of 2-cycles
						   if(i<blackAdj){
							   auxCiclo2.p1=i;  auxCiclo2.p2=blackAdj;
						   }
						   else{
							   auxCiclo2.p1=blackAdj;  auxCiclo2.p2=i;
						   }

						   if(greyAdj1 < greyAdj2){
							   auxCiclo2.p3=greyAdj1;  auxCiclo2.p4=greyAdj2;
						   }
						   else{
								auxCiclo2.p3=greyAdj2;  auxCiclo2.p4=greyAdj1;
						   }

						   if(i<greyAdj1){
							   auxCiclo2.c1=i;  auxCiclo2.c2=greyAdj1;
						   }
						   else{
							   auxCiclo2.c1=greyAdj1;  auxCiclo2.c2=i;
						   }

						   if(blackAdj<greyAdj2){
							   auxCiclo2.c3=blackAdj;  auxCiclo2.c4=greyAdj2;
						   }
						   else{
							   auxCiclo2.c3=greyAdj2;  auxCiclo2.c4=blackAdj;
						   }

						   if(checkExistence2Cycle(cycles2,num2cycles,auxCiclo2)==0){
							   cycles2[num2cycles++]=auxCiclo2;
						   }
					   }
				   }
               }//End for
           }//End for
        }
    }
    
	/* create the matching graph F(pi) */	
    const int n_vertices = MAX;
	my_graph g(n_vertices);
	
    for(i=0;i<num2cycles;i++){
		add_edge(cycles2[i].p1, cycles2[i].p3, g);
	}

	/* Find the maximum cardinality matching M, by eliminating
		2-cycles that share black edges */	    
	std::vector<graph_traits<my_graph>::vertex_descriptor> mate(n_vertices);
	checked_edmonds_maximum_cardinality_matching(g, &mate[0]);
	
	graph_traits<my_graph>::vertex_iterator vi, vi_end;
	for(tie(vi,vi_end) = vertices(g); vi != vi_end; ++vi){
		if (mate[*vi] != graph_traits<my_graph>::null_vertex() && *vi < mate[*vi]){
			//std::cout << "{" << *vi << ", " << mate[*vi] << "}" << std::endl;//--------
			for(i=0;i<num2cycles;i++){
				if ((cycles2[i].p1 == (int)*vi)&&(cycles2[i].p3 == (int)mate[*vi])){
					cycle2NoBlackEdges[num2cyclesNoBlackEdges] = i;
					num2cyclesNoBlackEdges++;
				}
			}
		}
	}
    
	//std::cout << std::endl << "Found a matching of size " << matching_size(g, &mate[0]) << std::endl;
	//std::cout << "The matching is:" << std::endl;
	
	/* create the ladder graph L(M) */
	my_graph g2(n_vertices);	
	
	int c,d;
	for(i=0;i<num2cyclesNoBlackEdges;i++){
		if (p->pi[cycles2[cycle2NoBlackEdges[i]].c1] < p->pi[cycles2[cycle2NoBlackEdges[i]].c2])
			c = cycles2[cycle2NoBlackEdges[i]].c1;
		else
			c = cycles2[cycle2NoBlackEdges[i]].c2;
			
		if (p->pi[cycles2[cycle2NoBlackEdges[i]].c3] < p->pi[cycles2[cycle2NoBlackEdges[i]].c4])
			d = cycles2[cycle2NoBlackEdges[i]].c3;
		else
			d = cycles2[cycle2NoBlackEdges[i]].c4;
		
		add_edge(c, d, g2);
		
	}	
		
	/* eliminate 2-cycles that share grey edges */
	std::vector<graph_traits<my_graph>::vertex_descriptor> mate2(n_vertices);
	checked_edmonds_maximum_cardinality_matching(g2, &mate2[0]);
	
	//graph_traits<my_graph>::vertex_iterator vi, vi_end;
	for(tie(vi,vi_end) = vertices(g2); vi != vi_end; ++vi){
		if (mate2[*vi] != graph_traits<my_graph>::null_vertex() && *vi < mate2[*vi]){			
			//std::cout << "*{" << *vi << ", " << mate2[*vi] << "}" << std::endl;//----------
			for(i=0;i<num2cycles;i++){
				if (p->pi[cycles2[i].c1] < p->pi[cycles2[i].c2])
					c = cycles2[i].c1;
				else
					c = cycles2[i].c2;
				
				if (p->pi[cycles2[i].c3] < p->pi[cycles2[i].c4])
					d = cycles2[i].c3;
				else
					d = cycles2[i].c4;
				
				if (c > d){
					temp = c;
					c = d;
					d = temp;
				}				
					
				if (c == (int)*vi && d == (int)mate2[*vi]){					
					cycle2NoGreyEdges[num2cyclesNoGreyEdges] = i;
					num2cyclesNoGreyEdges++;								
				}
			}
				
		}
	}
	
	/* eliminate the 2-cycles of the breakpoint graph "b" */
	for(i=0;i<num2cyclesNoGreyEdges;i++){
		j = cycle2NoGreyEdges[i];
		//eliminate black edges p1 and p2
		if (b->nblackAdj[cycles2[j].p1] == 1){
			b->nblackAdj[cycles2[j].p1]--;
			b->blackAdj[cycles2[j].p1][0] = 0;			
		}
		else if (b->nblackAdj[cycles2[j].p1] == 2){
			b->nblackAdj[cycles2[j].p1]--;
			if (b->blackAdj[cycles2[j].p1][0] == cycles2[j].p2){
				b->blackAdj[cycles2[j].p1][0] = b->blackAdj[cycles2[j].p1][1]; 				
			}
			b->blackAdj[cycles2[j].p1][1] = 0;			
		}
		
		if (b->nblackAdj[cycles2[j].p2] == 1){
			b->nblackAdj[cycles2[j].p2]--;
			b->blackAdj[cycles2[j].p2][0] = 0;			
		}
		else if (b->nblackAdj[cycles2[j].p2] == 2){
			b->nblackAdj[cycles2[j].p2]--;
			if (b->blackAdj[cycles2[j].p2][0] == cycles2[j].p1){
				b->blackAdj[cycles2[j].p2][0] = b->blackAdj[cycles2[j].p2][1]; 				
			}
			b->blackAdj[cycles2[j].p2][1] = 0;			
		}
		//eliminate black edges p3 and p4
		if (b->nblackAdj[cycles2[j].p3] == 1){
			b->nblackAdj[cycles2[j].p3]--;
			b->blackAdj[cycles2[j].p3][0] = 0;			
		}
		else if (b->nblackAdj[cycles2[j].p3] == 2){
			b->nblackAdj[cycles2[j].p3]--;
			if (b->blackAdj[cycles2[j].p3][0] == cycles2[j].p4){
				b->blackAdj[cycles2[j].p3][0] = b->blackAdj[cycles2[j].p3][1]; 				
			}
			b->blackAdj[cycles2[j].p3][1] = 0;			
		}
		
		if (b->nblackAdj[cycles2[j].p4] == 1){
			b->nblackAdj[cycles2[j].p4]--;
			b->blackAdj[cycles2[j].p4][0] = 0;			
		}
		else if (b->nblackAdj[cycles2[j].p4] == 2){
			b->nblackAdj[cycles2[j].p4]--;
			if (b->blackAdj[cycles2[j].p4][0] == cycles2[j].p3){
				b->blackAdj[cycles2[j].p4][0] = b->blackAdj[cycles2[j].p4][1]; 				
			}
			b->blackAdj[cycles2[j].p4][1] = 0;			
		}			
		//eliminate grey edges c1 and c2
		if (b->ngreyAdj[cycles2[j].c1] == 1){
			b->ngreyAdj[cycles2[j].c1]--;
			b->greyAdj[cycles2[j].c1][0] = 0;			
		}
		else if (b->ngreyAdj[cycles2[j].c1] == 2){
			b->ngreyAdj[cycles2[j].c1]--;
			if (b->greyAdj[cycles2[j].c1][0] == cycles2[j].c2){
				b->greyAdj[cycles2[j].c1][0] = b->greyAdj[cycles2[j].c1][1]; 				
			}
			b->greyAdj[cycles2[j].c1][1] = 0;			
		}
		
		if (b->ngreyAdj[cycles2[j].c2] == 1){
			b->ngreyAdj[cycles2[j].c2]--;
			b->greyAdj[cycles2[j].c2][0] = 0;			
		}
		else if (b->ngreyAdj[cycles2[j].c2] == 2){
			b->ngreyAdj[cycles2[j].c2]--;
			if (b->greyAdj[cycles2[j].c2][0] == cycles2[j].c1){
				b->greyAdj[cycles2[j].c2][0] = b->greyAdj[cycles2[j].c2][1]; 				
			}
			b->greyAdj[cycles2[j].c2][1] = 0;			
		}
		//eliminate grey edges c3 and c4
		if (b->ngreyAdj[cycles2[j].c3] == 1){
			b->ngreyAdj[cycles2[j].c3]--;
			b->greyAdj[cycles2[j].c3][0] = 0;			
		}
		else if (b->ngreyAdj[cycles2[j].c3] == 2){
			b->ngreyAdj[cycles2[j].c3]--;
			if (b->greyAdj[cycles2[j].c3][0] == cycles2[j].c4){
				b->greyAdj[cycles2[j].c3][0] = b->greyAdj[cycles2[j].c3][1]; 				
			}
			b->greyAdj[cycles2[j].c3][1] = 0;			
		}
		
		if (b->ngreyAdj[cycles2[j].c4] == 1){
			b->ngreyAdj[cycles2[j].c4]--;
			b->greyAdj[cycles2[j].c4][0] = 0;			
		}
		else if (b->ngreyAdj[cycles2[j].c4] == 2){
			b->ngreyAdj[cycles2[j].c4]--;
			if (b->greyAdj[cycles2[j].c4][0] == cycles2[j].c3){
				b->greyAdj[cycles2[j].c4][0] = b->greyAdj[cycles2[j].c4][1]; 				
			}
			b->greyAdj[cycles2[j].c4][1] = 0;			
		}	
	}
	
	//showBreakpointGraph(b,p);//--		

	/* perform the decomposition with remaining edges 
		and add it into the cycle decomposition */
	int startingPoint, found, point1, point2;
	
	D->nCycles = 0;
	found = 1;
	while(found){
		memset(D->cycle[D->nCycles].nblackAdj, 0, sizeof(D->cycle[D->nCycles].nblackAdj));
		memset(D->cycle[D->nCycles].ngreyAdj, 0, sizeof(D->cycle[D->nCycles].ngreyAdj));
		memset(D->cycle[D->nCycles].blackAdj, 0, sizeof(D->cycle[D->nCycles].blackAdj));
		memset(D->cycle[D->nCycles].greyAdj, 0, sizeof(D->cycle[D->nCycles].ngreyAdj));

		//find starting point
		found = 0;
		for(i=0;i<=p->size+1;i++){
			if (b->ngreyAdj[i] > 0){
				startingPoint = i;
				found = 1;
				break;
			}
		}	
		
		if (found){
			point1 = startingPoint;
			point2 = -1; 
			while (point2 != startingPoint){
				//add grey edge
				point2 = b->greyAdj[point1][b->ngreyAdj[point1]-1];
				
				D->cycle[D->nCycles].greyAdj[point1][D->cycle[D->nCycles].ngreyAdj[point1]] = point2;			
				D->cycle[D->nCycles].ngreyAdj[point1]++;			
				
				D->cycle[D->nCycles].greyAdj[point2][D->cycle[D->nCycles].ngreyAdj[point2]] = point1;			
				D->cycle[D->nCycles].ngreyAdj[point2]++;						
				
				//eliminate grey edge
				if (b->ngreyAdj[point1] == 1){
					b->ngreyAdj[point1]--;
					b->greyAdj[point1][0] = 0;				
				}
				else if (b->ngreyAdj[point1] == 2){
					b->ngreyAdj[point1]--;
					if (b->greyAdj[point1][0] == point2){
						b->greyAdj[point1][0] = b->greyAdj[point1][1];
					}
					b->greyAdj[point1][1] = 0;	
				}
				
				if (b->ngreyAdj[point2] == 1){
					b->ngreyAdj[point2]--;
					b->greyAdj[point2][0] = 0;				
				}
				else if (b->ngreyAdj[point2] == 2){
					b->ngreyAdj[point2]--;
					if (b->greyAdj[point2][0] == point1){
						b->greyAdj[point2][0] = b->greyAdj[point2][1];
					}
					b->greyAdj[point2][1] = 0;	
				}

				//add black edge			
				point1 = point2;
				point2 = b->blackAdj[point1][b->nblackAdj[point1]-1];
				
				D->cycle[D->nCycles].blackAdj[point1][D->cycle[D->nCycles].nblackAdj[point1]] = point2;
				D->cycle[D->nCycles].nblackAdj[point1]++;
				
				D->cycle[D->nCycles].blackAdj[point2][D->cycle[D->nCycles].nblackAdj[point2]] = point1;
				D->cycle[D->nCycles].nblackAdj[point2]++;			
				
				//eliminate black edge
				if (b->nblackAdj[point1] == 1){
					b->nblackAdj[point1]--;
					b->blackAdj[point1][0] = 0;				
				}
				else if (b->nblackAdj[point1] == 2){
					b->nblackAdj[point1]--;
					if (b->blackAdj[point1][0] == point2){
						b->blackAdj[point1][0] = b->blackAdj[point1][1];
					}
					b->blackAdj[point1][1] = 0;	
				}
				
				if (b->nblackAdj[point2] == 1){
					b->nblackAdj[point2]--;
					b->blackAdj[point2][0] = 0;				
				}
				else if (b->nblackAdj[point2] == 2){
					b->nblackAdj[point2]--;
					if (b->blackAdj[point2][0] == point1){
						b->blackAdj[point2][0] = b->blackAdj[point2][1];
					}
					b->blackAdj[point2][1] = 0;	
				}			
				point1 = point2;			
			}
			
			D->nCycles++;
		}
	}	

	/* add the 2-cycles into the cycle decomposition */
	for(i=0;i<num2cyclesNoGreyEdges;i++){		
		memset(D->cycle[D->nCycles].nblackAdj, 0, sizeof(D->cycle[D->nCycles].nblackAdj));
		memset(D->cycle[D->nCycles].ngreyAdj, 0, sizeof(D->cycle[D->nCycles].ngreyAdj));
		memset(D->cycle[D->nCycles].blackAdj, 0, sizeof(D->cycle[D->nCycles].blackAdj));
		memset(D->cycle[D->nCycles].greyAdj, 0, sizeof(D->cycle[D->nCycles].ngreyAdj));
		
		j = cycle2NoGreyEdges[i];
				
		//add black edges 
		D->cycle[D->nCycles].nblackAdj[cycles2[j].p1]++;
		D->cycle[D->nCycles].nblackAdj[cycles2[j].p2]++;
		D->cycle[D->nCycles].blackAdj[cycles2[j].p1][0] = cycles2[j].p2;
		D->cycle[D->nCycles].blackAdj[cycles2[j].p2][0] = cycles2[j].p1;
		
		D->cycle[D->nCycles].nblackAdj[cycles2[j].p3]++;
		D->cycle[D->nCycles].nblackAdj[cycles2[j].p4]++;
		D->cycle[D->nCycles].blackAdj[cycles2[j].p3][0] = cycles2[j].p4;
		D->cycle[D->nCycles].blackAdj[cycles2[j].p4][0] = cycles2[j].p3;		
		
		//add grey edges	
		D->cycle[D->nCycles].ngreyAdj[cycles2[j].c1]++;
		D->cycle[D->nCycles].ngreyAdj[cycles2[j].c2]++;
		D->cycle[D->nCycles].greyAdj[cycles2[j].c1][0] = cycles2[j].c2;
		D->cycle[D->nCycles].greyAdj[cycles2[j].c2][0] = cycles2[j].c1;
		
		D->cycle[D->nCycles].ngreyAdj[cycles2[j].c3]++;
		D->cycle[D->nCycles].ngreyAdj[cycles2[j].c4]++;
		D->cycle[D->nCycles].greyAdj[cycles2[j].c3][0] = cycles2[j].c4;
		D->cycle[D->nCycles].greyAdj[cycles2[j].c4][0] = cycles2[j].c3;		
		 
		D->nCycles++; 
	}
	
	//generate the cycle orientation
    makeCyclesOriented(D, p, C);
    
}

int checkExistence2Cycle(twoCycle c[],int n, twoCycle aux){
    int i;

    for(i=0;i<n;i++){
        if(aux.p1 == c[i].p1 && aux.p2 == c[i].p2 &&
            aux.p3==c[i].p3 && aux.p4==c[i].p4){
            return(1);
        }
        else if(aux.p3 == c[i].p1 && aux.p4 == c[i].p2 &&
                aux.p1 == c[i].p3 && aux.p2 == c[i].p4){
            return(1);
        }
    }
    return(0);
}

void makeCyclesOriented(unorientedDecomposition *D, permutation *p, cycleDecomposition *C){
    
    int i,j,startingPoint, black1,black2,grey1,grey2,tempGrey0,tempGrey1,diffBlacks;
    int aGrey,pBlack,blackEnd, found;
    breakpointGraph cycle;   
    
    C->nCycles = 0;
    for(i=0; i<D->nCycles; i++){
        
        cycle = D->cycle[i];
        found = 1;
        
        while (found){           
            //find starting point
            startingPoint = -1;
            found = 0;
            for(j=0; j<=p->size+1; j++){
                if(cycle.nblackAdj[j] > 0){
                    startingPoint = j; 
                    found = 1;
                    break;
                }
            }
            
            if (found){                                
                memset(C->cycle[C->nCycles].vertexEdge, -1, sizeof(C->cycle[C->nCycles].vertexEdge));
                //find orientation of current cycle, and initial black and grey
                black1 = startingPoint;
                black2 = cycle.blackAdj[startingPoint][0];
                blackEnd = cycle.blackAdj[startingPoint][0];
                
                grey1 = black2;
                if (cycle.ngreyAdj[grey1]==1){
                    grey2 = cycle.greyAdj[grey1][0];            
                }
                else{//==2
                    diffBlacks = black1 - black2;//direction of the cycle
                    
                    tempGrey0 = cycle.greyAdj[grey1][0];
                    tempGrey1 = cycle.greyAdj[grey1][1];        
                    grey2 = tempGrey0; 
                    
                    if ((cycle.nblackAdj[tempGrey1] == 1)&&
                                (tempGrey1 - cycle.blackAdj[tempGrey1][0] != diffBlacks)){
                         grey2 = tempGrey1; 
                     }
                     else if (cycle.nblackAdj[tempGrey0] == 2){
                         grey2 = tempGrey0; 
                     }
                     else if (cycle.nblackAdj[tempGrey1] == 2){
                         grey2 = tempGrey1;  
                     }
                }
                                
             
                if (p->pi[grey1] < p->pi[grey2]){
                    aGrey = p->pi[grey1];
                    pBlack = p->pi[black1];
                    C->cycle[C->nCycles].vertexEdge[aGrey][0] = aGrey + 1;//add grey edge, aGrey-aGrey+1
                    C->cycle[C->nCycles].vertexEdge[aGrey][1] = pBlack;//add black edge related with aGrey
                }
                else{
                    aGrey = p->pi[grey2];
                    pBlack = p->pi[black1];
                    C->cycle[C->nCycles].vertexEdge[aGrey][0] = aGrey + 1;//add grey edge, aGrey-aGrey+1
                    C->cycle[C->nCycles].vertexEdge[aGrey][2] = pBlack;//add black edge related with aGrey+1
                }
                     
                eliminateBlackAdj(&cycle,black1,black2);
                eliminateGreyAdj(&cycle,grey1,grey2);
                
                //find the other black and grey edges
                while(grey2 != startingPoint){
                    diffBlacks = black1 - black2;//direction of the cycle            
                    //find black edge 
                    black1 = grey2;
                    if (cycle.nblackAdj[black1] == 1){
                        black2 = cycle.blackAdj[black1][0];
                    }
                    else{ // == 2
                        if (black1-cycle.blackAdj[black1][0] != diffBlacks){
                            black2 = cycle.blackAdj[black1][0];                           
                        }
                        else{
                            black2 = cycle.blackAdj[black1][1]; 
                        }
                    }
                    
                    //find grey edge
                    grey1 = black2;
                    if (cycle.ngreyAdj[grey1]==1){
                        grey2 = cycle.greyAdj[grey1][0];            
                    }
                    else{//==2
                        diffBlacks = black1 - black2;//direction of the cycle
                        
                        tempGrey0 = cycle.greyAdj[grey1][0];
                        tempGrey1 = cycle.greyAdj[grey1][1];        
                        grey2 = tempGrey0; 
                        
                        if ((cycle.nblackAdj[tempGrey1] == 1)&&
                                    (tempGrey1 - cycle.blackAdj[tempGrey1][0] != diffBlacks)){
                             grey2 = tempGrey1; 
                         }
                         else if (cycle.nblackAdj[tempGrey0] == 2){
                             grey2 = tempGrey0; 
                         }
                         else if (cycle.nblackAdj[tempGrey1] == 2){
                             grey2 = tempGrey1;  
                         }
                    }
                    
                    if (p->piInverse[aGrey] == black1){
                        C->cycle[C->nCycles].vertexEdge[aGrey][1] = p->pi[black2];//add black edge related with aGrey             
                    }
                    else{
                        C->cycle[C->nCycles].vertexEdge[aGrey][2] = p->pi[black2];//add black edge related with aGrey+1
                    }
                                        
                    
                    if (p->pi[grey1] < p->pi[grey2]){
                        aGrey = p->pi[grey1];
                        pBlack = p->pi[black1];
                        C->cycle[C->nCycles].vertexEdge[aGrey][0] = aGrey + 1;//add grey edge, aGrey-aGrey+1
                        C->cycle[C->nCycles].vertexEdge[aGrey][1] = pBlack;//add black edge related with aGrey
                    }
                    else{
                        aGrey = p->pi[grey2];
                        pBlack = p->pi[black1];
                        C->cycle[C->nCycles].vertexEdge[aGrey][0] = aGrey + 1;//add grey edge, aGrey-aGrey+1
                        C->cycle[C->nCycles].vertexEdge[aGrey][2] = pBlack;//add black edge related with aGrey+1
                    }
                    
                    eliminateBlackAdj(&cycle,black1,black2);
                    eliminateGreyAdj(&cycle,grey1,grey2);          
                    
                }//END-WHILE
                
                if (p->piInverse[aGrey] == grey2){
                    C->cycle[C->nCycles].vertexEdge[aGrey][1] = p->pi[blackEnd];//add black edge related with aGrey             
                }
                else{
                    C->cycle[C->nCycles].vertexEdge[aGrey][2] = p->pi[blackEnd];//add black edge related with aGrey+1
                }
                                
                C->nCycles++;
            }//end if
        }//end while
		
	}//end for
}

void eliminateBlackAdj(breakpointGraph *cycle,int black1, int black2){
    if (cycle->nblackAdj[black1] == 1){
        cycle->nblackAdj[black1]--;
        cycle->blackAdj[black1][0] = 0;				
    }
    else if (cycle->nblackAdj[black1] == 2){
        cycle->nblackAdj[black1]--;
        if (cycle->blackAdj[black1][0] == black2){
            cycle->blackAdj[black1][0] = cycle->blackAdj[black1][1];
        }
        cycle->blackAdj[black1][1] = 0;	
    }
    
    if (cycle->nblackAdj[black2] == 1){
        cycle->nblackAdj[black2]--;
        cycle->blackAdj[black2][0] = 0;				
    }
    else if (cycle->nblackAdj[black2] == 2){
        cycle->nblackAdj[black2]--;
        if (cycle->blackAdj[black2][0] == black1){
            cycle->blackAdj[black2][0] = cycle->blackAdj[black2][1];
        }
        cycle->blackAdj[black2][1] = 0;	
    }
    
}

void eliminateGreyAdj(breakpointGraph *cycle,int grey1, int grey2){
    if (cycle->ngreyAdj[grey1] == 1){
        cycle->ngreyAdj[grey1]--;
        cycle->greyAdj[grey1][0] = 0;				
    }
    else if (cycle->ngreyAdj[grey1] == 2){
        cycle->ngreyAdj[grey1]--;
        if (cycle->greyAdj[grey1][0] == grey2){
            cycle->greyAdj[grey1][0] = cycle->greyAdj[grey1][1];
        }
        cycle->greyAdj[grey1][1] = 0;	
    }
    
    if (cycle->ngreyAdj[grey2] == 1){
        cycle->ngreyAdj[grey2]--;
        cycle->greyAdj[grey2][0] = 0;				
    }
    else if (cycle->ngreyAdj[grey2] == 2){
        cycle->ngreyAdj[grey2]--;
        if (cycle->greyAdj[grey2][0] == grey1){
            cycle->greyAdj[grey2][0] = cycle->greyAdj[grey2][1];
        }
        cycle->greyAdj[grey2][1] = 0;	
    }
    
}

void showBreakpointGraph(orientedCycle *b,permutation *p){
	int i;
    
    printf("grey-black edges:\n");	
	for(i=0;i<=p->size;i++)
		printf("%*d",3,i);
	printf("\n");
    for(i=0;i<=p->size;i++){
		printf("%*d",3,b->vertexEdge[i][0]);
	}
    printf("\n");
    for(i=0;i<=p->size;i++){
		printf("%*d",3,b->vertexEdge[i][1]);
	}
    printf("\n");
    for(i=0;i<=p->size;i++){
		printf("%*d",3,b->vertexEdge[i][2]);
	}
    printf("\n");
}





