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
#include <tuple>

#include <simplex.h>

#define EPS 1e-3

typedef long long lld;
typedef unsigned long long llu;
using namespace std;

int total_steps = 0;

Simplex::Simplex(int n, int m, double **A, double *b, double *c, double v) : n(n), m(m), v(v)
{
    this -> A = new double*[m];
    for (int i=0;i<m;i++)
    {
        this -> A[i] = new double[n+1]; // allocate +1 because of aux. LP
        for (int j=0;j<n;j++)
        {
            this -> A[i][j] = A[i][j];
        }
    }
    
    this -> b = new double[m];
    for (int i=0;i<m;i++)
    {
        this -> b[i] = b[i];
    }
    
    this -> c = new double[n+1]; // allocate +1 because of aux. LP
    for (int i=0;i<n;i++)
    {
        this -> c[i] = c[i];
    }
    
    this -> N = new int[n+1]; // allocate +1 because of aux. LP
    this -> B = new int[m];
}

Simplex::~Simplex()
{
    for (int i=0;i<m;i++) delete[] A[i];
    delete[] A;
    delete[] b;
    delete[] c;
    delete[] N;
    delete[] B;
}

// pivot yth variable around xth constraint
void Simplex::pivot(int x, int y)
{
    //printf("Pivoting variable %d around constraint %d.\n", y, x);
    total_steps++;
    
    // first rearrange the x-th row
    for (int j=0;j<n;j++)
    {
        if (j != y)
        {
            A[x][j] /= -A[x][y];
        }
    }
    b[x] /= -A[x][y];
    A[x][y] = 1.0 / A[x][y];
    
    // now rearrange the other rows
    for (int i=0;i<m;i++)
    {
        if (i != x)
        {
            for (int j=0;j<n;j++)
            {
                if (j != y)
                {
                    A[i][j] += A[i][y] * A[x][j];
                }
            }
            b[i] += A[i][y] * b[x];
            A[i][y] *= A[x][y];
        }
    }
    
    // now rearrange the objective function
    for (int j=0;j<n;j++)
    {
        if (j != y)
        {
            c[j] += c[y] * A[x][j];
        }
    }
    v += c[y] * b[x];
    c[y] *= A[x][y];
    
    // finally, swap the basic & nonbasic variable
    swap(B[x], N[y]);
}

// Run a single iteration of the simplex algorithm.
// Returns: 0 if OK, 1 if STOP, -1 if UNBOUNDED
int Simplex::iterate_simplex()
{
    /*
     printf("--------------------\n");
     printf("State:\n");
     printf("Maximise: ");
     for (int j=0;j<n;j++) printf("%lfx_%d + ", c[j], N[j]);
     printf("%lf\n", v);
     printf("Subject to:\n");
     for (int i=0;i<m;i++)
     {
         for (int j=0;j<n;j++) printf("%lfx_%d + ", A[i][j], N[j]);
         printf("%lf = x_%d\n", b[i], B[i]);
     }
    */
    // getchar(); // uncomment this for debugging purposes!
    
    double v_prev = v;
    
    vector<int> vars;
    int best_var = -1;
    for (int j=0;j<n;j++)
    {
        if (c[j] > 0)
        {
            vars.push_back(j);
            /* Bland's Rule
            if (best_var == -1 || N[j] < ind)
            {
                ind = N[j];
                best_var = j;
            }
            */
        }
    }
    if (vars.empty()) return 1;
    
    best_var = vars[rand() % vars.size()];
    
    double max_constr = INFINITY;
    int best_constr = -1;
    for (int i=0;i<m;i++)
    {
        if (A[i][best_var] < 0)
        {
            double curr_constr = -b[i] / A[i][best_var];
            if (curr_constr < max_constr)
            {
                max_constr = curr_constr;
                best_constr = i;
            }
        }
    }
    if (std::isinf(max_constr)) return -1;
    else pivot(best_constr, best_var);
    
    // underflow avoiding: round all entries that are at most 1e-3 away from an integer
    for (int i=0;i<m;i++)
    {
        for (int j=0;j<n;j++)
        {
            if (fabs(ceil(A[i][j]) - A[i][j]) < EPS) A[i][j] = ceil(A[i][j]);
            else if (fabs(floor(A[i][j]) - A[i][j]) < EPS) A[i][j] = floor(A[i][j]);
        }
    }
    
    for (int i=0;i<m;i++)
    {
        if (fabs(ceil(b[i]) - b[i]) < EPS) b[i] = ceil(b[i]);
        if (fabs(floor(b[i]) - b[i]) < EPS) b[i] = floor(b[i]);
    }
    for (int j=0;j<n;j++)
    {
        if (fabs(ceil(c[j]) - c[j]) < EPS) c[j] = ceil(c[j]);
        if (fabs(floor(c[j]) - c[j]) < EPS) c[j] = floor(c[j]);
    }
    
    if (fabs(ceil(v) - v) < EPS) v = ceil(v);
    if (fabs(floor(v) - v) < EPS) v = floor(v);
    
    if (fabs(v - v_prev) > EPS) cout << "Iteration " << total_steps << ": " << v << endl;
    
    return 0;
}

// (Possibly) converts the LP into a slack form with a feasible basic solution.
// Returns 0 if OK, -1 if INFEASIBLE
int Simplex::initialise_simplex()
{
    total_steps = 0;
    
    int k = -1;
    double min_b = -1;
    for (int i=0;i<m;i++)
    {
        if (k == -1 || b[i] < min_b)
        {
            k = i;
            min_b = b[i];
        }
    }
    
    if (b[k] >= 0) // basic solution feasible!
    {
        for (int j=0;j<n;j++) N[j] = j;
        for (int i=0;i<m;i++) B[i] = n + i;
        return 0;
    }
    
    // generate auxiliary LP
    n++;
    for (int j=0;j<n;j++) N[j] = j;
    for (int i=0;i<m;i++) B[i] = n + i;
    
    // store the objective function
    double *c_old = new double[n];
    for (int j=0;j<n-1;j++) c_old[j] = c[j];
    double v_old = v;
    
    // aux. objective function
    c[n-1] = -1;
    for (int j=0;j<n-1;j++) c[j] = 0;
    v = 0;
    // aux. coefficients
    for (int i=0;i<m;i++) A[i][n-1] = 1;
    
    // perform initial pivot
    pivot(k, n - 1);
    
    // now solve aux. LP
    int code;
    while (!(code = iterate_simplex()));
    
    assert(code == 1); // aux. LP cannot be unbounded!!!
    
    if (v != 0)
    {
        delete[] c_old;
        return -1; // infeasible!
    }
    
    int z_basic = -1;
    for (int i=0;i<m;i++)
    {
        if (B[i] == n - 1)
        {
            z_basic = i;
            break;
        }
    }
    
    // if x_n basic, perform one degenerate pivot to make it nonbasic
    if (z_basic != -1) pivot(z_basic, n - 1);
    
    int z_nonbasic = -1;
    for (int j=0;j<n;j++)
    {
        if (N[j] == n - 1)
        {
            z_nonbasic = j;
            break;
        }
    }
    assert(z_nonbasic != -1);
    
    for (int i=0;i<m;i++)
    {
        A[i][z_nonbasic] = A[i][n-1];
    }
    swap(N[z_nonbasic], N[n - 1]);
    
    n--;
    for (int j=0;j<n;j++) if (N[j] > n) N[j]--;
    for (int i=0;i<m;i++) if (B[i] > n) B[i]--;
    
    for (int j=0;j<n;j++) c[j] = 0;
    v = v_old;
    
    for (int j=0;j<n;j++)
    {
        bool ok = false;
        for (int jj=0;jj<n;jj++)
        {
            if (j == N[jj])
            {
                c[jj] += c_old[j];
                ok = true;
                break;
            }
        }
        if (ok) continue;
        for (int i=0;i<m;i++)
        {
            if (j == B[i])
            {
                for (int jj=0;jj<n;jj++)
                {
                    c[jj] += c_old[j] * A[i][jj];
                }
                v += c_old[j] * b[i];
                break;
            }
        }
    }
    
    delete[] c_old;
    return 0;
}

// Runs the simplex algorithm to optimise the LP.
// Returns a vector of -1s if unbounded, -2s if infeasible.
tuple<vector<double>, double, int> Simplex::simplex()
{
    printf("Running the Simplex algorithm with %d variables and %d constraints.\n", n, m);
    
    if (initialise_simplex() == -1)
    {
        return {vector<double>(n + m, -2), INFINITY, total_steps};
    }
    
    int code;
    while (!(code = iterate_simplex()));
    
    if (code == -1) return {vector<double>(n + m, -1), INFINITY, total_steps};
    
    vector<double> ret;
    ret.resize(n + m);
    for (int j=0;j<n;j++)
    {
        ret[N[j]] = 0;
    }
    for (int i=0;i<m;i++)
    {
        ret[B[i]] = b[i];
    }
    
    printf("Finished in %d iterations.\n", total_steps);
    
    return {ret, v, total_steps};
}
