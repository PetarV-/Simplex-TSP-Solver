/*
 Petar 'PetarV' Velickovic
 Algorithm: Simplex Algorithm
*/

#ifndef SIMPLEX
#define SIMPLEX

#include <tuple>
#include <vector>

/*
 The Simplex algorithm aims to solve a linear program - optimising a linear function subject
 to linear constraints. As such it is useful for a very wide range of applications.
 
 N.B. The linear program has to be given in *slack form*, which is as follows:
 maximise
 c_1 * x_1 + c_2 * x_2 + ... + c_n * x_n + v
 subj. to
 a_11 * x_1 + a_12 * x_2 + ... + a_1n * x_n + b_1 = s_1
 a_21 * x_1 + a_22 * x_2 + ... + a_2n * x_n + b_2 = s_2
 ...
 a_m1 * x_1 + a_m2 * x_2 + ... + a_mn * x_n + b_m = s_m
 and
 x_1, x_2, ..., x_n, s_1, s_2, ..., s_m >= 0
 
 Every linear program can be translated into slack form; the parameters to specify are:
 - the number of variables, n, and the number of constraints, m;
 - the matrix A = [[A_11, A_12, ..., A_1n], ..., [A_m1, A_m2, ..., A_mn]];
 - the vector b = [b_1, b_2, ..., b_m];
 - the vector c = [c_1, c_2, ..., c_n] and the constant v.
 
 Complexity:    O(m^(n/2)) worst case
                O(n + m) average case (common)
*/

class Simplex
{
private:
    int n, m;
    double **A, *b, *c, v;
    int *N, *B; // nonbasic & basic
    
public:
    Simplex(int n, int m, double **A, double *b, double *c, double v);
    
    ~Simplex();
    
    // pivot yth variable around xth constraint
    void pivot(int x, int y);
    
    // Run a single iteration of the simplex algorithm.
    // Returns: 0 if OK, 1 if STOP, -1 if UNBOUNDED
    int iterate_simplex();
    
    // (Possibly) converts the LP into a slack form with a feasible basic solution.
    // Returns 0 if OK, -1 if INFEASIBLE
    int initialise_simplex();

    // Runs the simplex algorithm to optimise the LP.
    // Returns a vector of -1s if unbounded, -2s if infeasible.
    std::tuple<std::vector<double>, double, int> simplex();
};

#endif
