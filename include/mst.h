/*
 Petar 'PetarV' Velickovic
 Algorithm: MST-TSP approximation
*/

#ifndef MST
#define MST

#include <vector>

#define MAX_N 5001

/*
 This is a polynomial time 2-approximation algorithm to the Traveling
 Salesman Problem, making advantage of a minimal spanning tree (MST)
 subroutine. It proceeds as follows:
    1. Compute the MST of the input graph;
    2. Perform a preorder traversal of the MST;
    3. Traverse vertices in that exact order, closing the circuit at the end.
*/

std::pair<double, std::vector<std::pair<int, int> > > tsp_mst(int n, double adj[MAX_N][MAX_N]);

#endif
