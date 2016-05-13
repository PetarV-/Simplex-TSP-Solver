/*
 Petar 'PetarV' Velickovic
 Linear Programming TSP Solver (Dantzig-Fulkerson-Johnson)
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
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

#include <mst.h>
#include <simplex.h>

#define MAX_N 5001
#define EPS 1e-5

typedef long long lld;
typedef unsigned long long llu;
typedef unsigned int uint;
using namespace std;

char path_adj[150], path_demo[150], file_graph[150], file_map[150];
char path_graph[150], path_map[150], path_pdf[150];
FILE *f_adj, *f_graph, *f_lp;

char cmd[150];
char cmd_map[150];
char cmd_prev[150];

int n;
double adj[MAX_N][MAX_N];

int simp_n;
double c[MAX_N * MAX_N];

// Constraints
vector<vector<double> > Ap;
vector<double> bp;

int curr_pos = 0;
stack<int> positions;

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

inline void dump_lp(int n, int m, double **A, double *b, double *c, double v)
{
    fprintf(f_lp, "minimise %.2lf ", v);
    for (int i=0;i<n;i++)
    {
        if (fabs(c[i]) > EPS)
        {
            fprintf(f_lp, "%c %.2lf * x_{%02d, %02d} ", (c[i] >= 0.0 ? '+' : '-'), fabs(c[i]), decode_edge(i).first, decode_edge(i).second);
        }
    }
    fprintf(f_lp, "\nsubject to\n");
    for (int i=0;i<m;i++)
    {
        fprintf(f_lp, "s_{%03d} = %.2lf ", i, b[i]);
        for (int j=0;j<n;j++)
        {
            if (fabs(A[i][j]) > EPS)
            {
                fprintf(f_lp, "%c %.2lf * x_{%02d, %02d} ", (A[i][j] >= 0.0 ? '+' : '-'), fabs(A[i][j]), decode_edge(j).first, decode_edge(j).second);
            }
        }
        fprintf(f_lp, "\n");
    }
}

inline void dump_title(double obj, int n, int m, int iter)
{
    fprintf(f_graph, "\\node[] (title) at (10, 25) {\\Huge Objective value: $%lf$, $%d$ variables, $%d$ constraints, $%d$ iterations};\n", obj, n, m, iter);
}

inline void dump_edge(int x, double val)
{
    fprintf(f_graph, "\\draw[edge,%s] (%d) to node[lab]{%g} (%d);\n", ((fabs(1.0 - val) < EPS) ? "black" : "red"), decode_edge(x).first, val, decode_edge(x).second);
}

inline void write_full_edge(int x, int y)
{
    fprintf(f_graph, "\\draw[edge] (%d) to node[lab]{%d} (%d);\n", x, 1, y);
}

char *rep_ext(char *name, const char *ext)
{
    char *ret = new char[150];
    strcpy(ret, name);
    int ii = strlen(ret) - 1;
    while (ii >= 0 && ret[ii] != '.') ii--;
    assert(ii >= 0);
    int l = strlen(ext);
    for (int i=0;i<=l;i++) // takes care of '\0' at the end as well
    {
        ret[ii + i] = ext[i];
    }
    
    return ret;
}

int main()
{
    srand(time(NULL));
    
    printf("Linear Programming TSP Solver, implemented by Petar Veličković.\n");
    printf("Thanks to Thomas Sauerwald for the visualisation data!\n");
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
    
    curr_pos = Ap.size();
    
    printf("Constraints generated!\n");
    print_delimiter();
    
    printf("Enter the path to the folder containing the demo .tex files:\n");
    scanf("%s", path_demo);
    print_delimiter();
    
    printf("Enter the name of the file containing the graph:\n");
    scanf("%s", file_graph);
    print_delimiter();
    
    strcpy(path_graph, path_demo);
    strcat(path_graph, file_graph);
    
    printf("Enter the name of the file containing the map:\n");
    scanf("%s", file_map);
    print_delimiter();
    
    strcpy(path_map, path_demo);
    strcat(path_map, file_map);
    
    strcpy(path_pdf, path_demo);
    strcat(path_pdf, rep_ext(file_map, ".pdf"));
    
    
    while (true)
    {
        printf("Enter one of the following:\n");
        printf("- SOLVE : to launch the Simplex Algorithm on the constraints given thus far;\n");
        printf("- REM_LOOP N x_1 x_2 ... x_N : to add a loop-removing constraint for the subset of N nodes x_1, x_2, ..., x_N;\n");
        printf("- REM_LOOP_RNG lo hi : to add a loop-removing constraint for the contiguous interval of nodes [lo, hi];\n");
        printf("- SET x_1 x_2 v : to add constraints that set the edge between x_1 and x_2 to v (where v is either 0 or 1);\n");
        printf("- UNDO N : to remove the N previously generated constraint sets;\n");
        printf("- APPROX_MST : to run the MST-based 2-approximation algorithm;\n");
        printf("- EXIT : to stop the program.\n");
        scanf("%s", cmd);
        
        if (strcmp(cmd, "SOLVE") == 0)
        {
            int constr_n = curr_pos;
            
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
            
            if ((f_lp = fopen("lp.txt", "w")) == NULL)
            {
                printf("Error: LP dump file could not be opened\n");
                return 4;
            }
            
            dump_lp(simp_n, constr_n, A, b, c, 0);
            
            fclose(f_lp);
            
            Simplex *s = new Simplex(simp_n, constr_n, A, b, c, 0);
            auto ret = s -> simplex();
            
            while (std::isnan(get<1>(ret)))
            {
                delete s;
                s = new Simplex(simp_n, constr_n, A, b, c, 0);
                ret = s -> simplex();
            }
            
            delete s;
            for (int i=0;i<constr_n;i++) delete[] A[i];
            delete[] A;
            delete[] b;
            
            printf("Simplex subroutine finished!\n");
            
            vector<double> xs = get<0>(ret);
            double val = get<1>(ret);
            int iters = get<2>(ret);
            
            if (std::isinf(val))
            {
                if (xs[0] == -1) printf("Objective function unbounded!\n");
                else if (xs[0] == -2) printf("Linear program infeasible!\n");
            }
            
            else
            {
                printf("Objective value: %lf\n", val);
                
                //for (int i=0;i<simp_n;i++) if (ret.first[i] != 0) cout << ret.first[i] << endl;
        
                if ((f_graph = fopen(path_graph, "w")) == NULL)
                {
                    printf("Error: Graph file could not be opened!\n");
                    return 4;
                }
                
                dump_title(val, simp_n, constr_n, iters);
                
                for (int i=0;i<simp_n;i++)
                {
                    if (xs[i] > 0) dump_edge(i, xs[i]);
                }
            
                fclose(f_graph);
            
                printf("Results written to the graph file! Compiling the map...\n");
                
                sprintf(cmd_map, "(cd %s && exec pdflatex %s &> /dev/null)", path_demo, file_map);
                
                system(cmd_map);
                
                printf("Map compiled!\n");
#ifdef __APPLE__
                // These shell instructions will work only on OS X
                sprintf(cmd_prev, "killall qlmanage &> /dev/null; qlmanage -p %s &> /dev/null &", path_pdf);
                
                printf("%s\n", cmd_prev);
                
                system(cmd_prev);
#endif
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
            
            if (curr_pos == (int)Ap.size())
            {
                Ap.push_back(constr);
                bp.push_back(val);
            }
            else
            {
                Ap[curr_pos] = constr;
                bp[curr_pos] = val;
            }
            
            positions.push(curr_pos);
            curr_pos++;
            
            printf("Loop-removing constraint added!\n");
            print_delimiter();
        }
        
        else if (strcmp(cmd, "REM_LOOP_RNG") == 0)
        {
            int lo, hi;
            scanf("%d%d", &lo, &hi);
            
            assert(hi >= lo);
            
            int num = hi - lo + 1;
            vector<int> vals(num);
            for (int i=0;i<num;i++) vals[i] = lo + i;
            
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
            
            if (curr_pos == (int)Ap.size())
            {
                Ap.push_back(constr);
                bp.push_back(val);
            }
            else
            {
                Ap[curr_pos] = constr;
                bp[curr_pos] = val;
            }
            
            positions.push(curr_pos);
            curr_pos++;
            
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
            
            positions.push(curr_pos);
            
            if (curr_pos == (int)Ap.size())
            {
                Ap.push_back(constr1);
                bp.push_back(v);
            }
            else
            {
                Ap[curr_pos] = constr1;
                bp[curr_pos] = v;
            }
            
            curr_pos++;
            
            if (v == 1)
            {
                if (curr_pos == (int)Ap.size())
                {
                    Ap.push_back(constr2);
                    bp.push_back(-1.0);
                }
                else
                {
                    Ap[curr_pos] = constr2;
                    bp[curr_pos] = -1.0;
                }
                
                curr_pos++;
            }
            
            printf("Value setting constraints added!\n");
            print_delimiter();
        }
        
        else if (strcmp(cmd, "UNDO") == 0)
        {
            int steps;
            scanf("%d", &steps);
            
            bool done = true;
            
            while (steps--)
            {
                if (positions.empty())
                {
                    printf("There's nothing to undo!\n");
                    print_delimiter();
                    done = false;
                    break;
                }
            
                curr_pos = positions.top();
                positions.pop();
            }
            
            if (!done) continue;
            
            printf("Previous generated constraint set(s) removed!\n");
            print_delimiter();
        }
        
        else if (strcmp(cmd, "APPROX_MST") == 0)
        {
            auto mst_sol = tsp_mst(n, adj);
            
            printf("Objective value: %lf\n", mst_sol.first);
            
            //for (int i=0;i<simp_n;i++) if (ret.first[i] != 0) cout << ret.first[i] << endl;
            
            if ((f_graph = fopen(path_graph, "w")) == NULL)
            {
                printf("Error: Graph file could not be opened!\n");
                return 4;
            }
            
            for (uint i=0;i<mst_sol.second.size();i++)
            {
                write_full_edge(mst_sol.second[i].first, mst_sol.second[i].second);
            }
            
            fclose(f_graph);
            
            printf("Results written to the graph file! Compiling the map...\n");
            
            sprintf(cmd_map, "(cd %s && exec pdflatex %s &> /dev/null)", path_demo, file_map);
            
            system(cmd_map);
            
            printf("Map compiled!\n");
#ifdef __APPLE__
            // These shell instructions will work only on OS X
            sprintf(cmd_prev, "killall qlmanage &> /dev/null; qlmanage -p %s &> /dev/null &", path_pdf);
            
            printf("%s\n", cmd_prev);
            
            system(cmd_prev);
#endif
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
