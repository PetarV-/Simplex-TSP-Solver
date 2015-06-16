/*
 Petar 'PetarV' Velickovic
 Linear Programming TSP Solver (Dantzig-Fulkerson-Johnson)
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <algorithm>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <complex>

#include <simplex.h>

#define MAX_N 5001

typedef long long lld;
typedef unsigned long long llu;
using namespace std;

char path_adj[150], path_graph[150];
FILE *f_adj, *f_graph;

char cmd[150];

int n;
double adj[MAX_N][MAX_N];

int simp_n;
double c[MAX_N * MAX_N];

// Constraints
vector<vector<double> > Ap;
vector<double> bp;

inline void print_delimiter()
{
    printf("---------------------------------------------------------------\n");
}

inline int encode_edge(int i, int j)
{
    // Precondition: i > j
    assert(i > j);
    return (((i - 1) * (i - 2)) >> 1) + (j - 1);
}

inline pair<int, int> decode_edge(int x)
{
    int i = 2, j = 1;
    
    while (((i * (i - 1)) >> 1) <= x) i++;
    j = x - (((i - 1) * (i - 2)) >> 1) + 1;
    
    return make_pair(i, j);
}

inline void dump_edge(int x, double val)
{
    fprintf(f_graph, "\\draw[edge] (%d) to node[lab]{%g} (%d);\n", decode_edge(x).first, val, decode_edge(x).second);
}

int main()
{
    printf("Linear Programming TSP Solver, implemented by Petar Veličković.\n");
    printf("Special thanks to Thomas Sauerwald!\n");
    print_delimiter();
    
    printf("Enter the path to the file containing the adjacency matrix:\n");
    scanf("%s", path_adj);
    
    if ((f_adj = fopen(path_adj, "r")) == NULL)
    {
        printf("Error: Adjacency matrix file could not be opened!\n");
        return 1;
    }
    
    printf("Processing adjacency matrix...\n");
    if (fscanf(f_adj, "%d", &n) == 0)
    {
        printf("Error: Improper adjacency matrix format!\n");
        return 2;
    }
    
    printf("The graph has %d nodes.\n", n);
    
    for (int i=2;i<=n;i++)
    {
        for (int j=1;j<i;j++)
        {
            if (fscanf(f_adj, "%lf", &adj[i][j]) == 0)
            {
                printf("Error: Improper adjacency matrix format!\n");
                return 3;
            }
        }
    }
    
    fclose(f_adj);
    
    printf("Adjacency matrix successfully read!\n");
    print_delimiter();
    
    printf("Generating basic constraints...\n");
    
    simp_n = (n * (n-1)) >> 1; // the number of edges
    
    // Generating objective function
    for (int i=0;i<simp_n;i++)
    {
        c[i] = -adj[decode_edge(i).first][decode_edge(i).second];
    }
    
    // Generating constraints
    for (int i=1;i<=n;i++)
    {
        vector<double> constr_1(simp_n, 0.0), constr_2(simp_n, 0.0);
        
        for (int j=1;j<=n;j++)
        {
            if (j < i)
            {
                constr_1[encode_edge(i, j)] = -1.0;
                constr_2[encode_edge(i, j)] = 1.0;
            }
            else if (j > i)
            {
                constr_1[encode_edge(j, i)] = -1.0;
                constr_2[encode_edge(j, i)] = 1.0;
            }
        }
        
        Ap.push_back(constr_1);
        bp.push_back(2.0);
        
        Ap.push_back(constr_2);
        bp.push_back(-2.0);
    }
    
    for (int i=0;i<simp_n;i++)
    {
        vector<double> constr(simp_n, 0.0);
        constr[i] = -1.0;
        
        Ap.push_back(constr);
        bp.push_back(1.0);
    }
    
    printf("Constraints generated!\n");
    print_delimiter();
    
    printf("Enter the path to the file containing the graph:\n");
    scanf("%s", path_graph);
    print_delimiter();
    
    while (true)
    {
        printf("Enter one of the following:\n");
        printf("- SOLVE : to launch the Simplex Algorithm on the constraints given thus far;\n");
        printf("- REM_LOOP N x_1 x_2 ... x_N : to add a loop-removing constraint for the subset of N nodes x_1, x_2, ..., x_N;\n");
        printf("- SET x_1 x_2 v : to add constraints that set the edge between x_1 and x_2 to v (where v is either 0 or 1);\n");
        printf("- EXIT : to stop the program.\n");
        scanf("%s", cmd);
        
        if (strcmp(cmd, "SOLVE") == 0)
        {
            int constr_n = Ap.size();
            
            double **A = new double*[constr_n];
            for (int i=0;i<constr_n;i++)
            {
                A[i] = new double[simp_n];
                for (int j=0;j<simp_n;j++)
                {
                    A[i][j] = Ap[i][j];
                }
            }
            
            double *b = new double[constr_n];
            for (int i=0;i<constr_n;i++)
            {
                b[i] = bp[i];
            }
            
            Simplex *s = new Simplex(simp_n, constr_n, A, b, c, 0);
            pair<vector<double>, double> ret = s -> simplex();
            
            delete s;
            for (int i=0;i<constr_n;i++) delete[] A[i];
            delete[] A;
            delete[] b;
            
            printf("Simplex subroutine finished!\n");
            
            if (isinf(ret.second))
            {
                if (ret.first[0] == -1) printf("Objective function unbounded!\n");
                else if (ret.first[0] == -2) printf("Linear program infeasible!\n");
            }
            
            else
            {
                printf("Objective value: %lf\n", ret.second);
                
                //for (int i=0;i<simp_n;i++) if (ret.first[i] != 0) cout << ret.first[i] << endl;
        
                if ((f_graph = fopen(path_graph, "w")) == NULL)
                {
                    printf("Error: Graph file could not be opened!\n");
                    return 4;
                }
                
                for (int i=0;i<simp_n;i++)
                {
                    if (ret.first[i] > 0) dump_edge(i, ret.first[i]);
                }
            
                fclose(f_graph);
            
                printf("Results written to the graph file!\n");
            }
            
            print_delimiter();
        }
        
        else if (strcmp(cmd, "REM_LOOP") == 0)
        {
            int num;
            scanf("%d", &num);
            
            vector<int> vals(num);
            for (int i=0;i<num;i++) scanf("%d", &vals[i]);
            
            vector<double> constr(simp_n, 0.0);
            
            for (int i=0;i<num;i++)
            {
                for (int j=i+1;j<num;j++)
                {
                    if (vals[i] > vals[j]) constr[encode_edge(vals[i], vals[j])] = -1.0;
                    else if (vals[i] < vals[j]) constr[encode_edge(vals[j], vals[i])] = -1.0;
                }
            }
            
            double val = num - 1.0;
            
            Ap.push_back(constr);
            bp.push_back(val);
            
            printf("Loop-removing constraint added!\n");
            print_delimiter();
        }
        
        else if (strcmp(cmd, "SET") == 0)
        {
            // x(i, j) = v
            int x, y, v;
            scanf("%d%d%d", &x, &y, &v);
            
            vector<double> constr1(simp_n, 0.0), constr2(simp_n, 0.0);
            
            if (x > y)
            {
                constr1[encode_edge(x, y)] = -1.0;
                constr2[encode_edge(x, y)] = 1.0;
            }
            else if (x < y)
            {
                constr1[encode_edge(y, x)] = -1.0;
                constr2[encode_edge(y, x)] = 1.0;
            }
            
            Ap.push_back(constr1);
            bp.push_back(v);
            
            if (v == 1)
            {
                Ap.push_back(constr2);
                bp.push_back(-1.0);
            }
            
            printf("Value setting constraints added!\n");
            print_delimiter();
        }
        
        else if (strcmp(cmd, "EXIT") == 0)
        {
            break;
        }
        
        else
        {
            printf("Incorrectly entered command! Try again.\n");
            print_delimiter();
        }
    }
    
    return 0;
}