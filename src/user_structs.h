#ifndef USER_STRUCTS_H
#define USER_STRUCTS_H

#include <stdbool.h>
#include <mkl.h>

typedef struct {
    int i;
    int j;
    int k;
    int l;
    int m;
    int type;
    double value;
    double y;
} Cut;

typedef struct {
    int *i; // i >= j
    int *j; // i >= j
    double *val;
    int nnz;
} Sparse_Matrix;

typedef struct {
    int number_of_facilities;
    int *lengths_of_facilities;
	int **weights;
	long **C; // ganzzahlig, d.h. mit Faktor 4 // TODO
	double K; // Konstante
	long K_times_two;
    long evaluations;
} Instance_Data;

typedef struct {
    double *Q; // Objective matrix in DENSE format
    Sparse_Matrix Qs; // TODO ist das noetig?
    int dimension_of_Q;
    Sparse_Matrix *Bs; // list of sparse matrices for the equality constraints
    double *b; // right-hand-side vector of equality constraints
    int number_of_equations; // number of equality constraints
    double *X; // Stores current X (primal solution)
    double fixed_objective_value;
    int *index_shift;
    long evaluations;
    VSLStreamStatePtr rng_stream;
    int number_of_cuts;
    Cut *current_cuts;
    Cut *new_cuts;
    
    
    // fuer Heuristiken
    bool *temp_solution;
    bool *temp_solution_fixed_variables;
    int *permutation;
    int *to_do_list;
    long *distances_times_two;
    int *random_permutation;
    int *most_fractional_permutation;
    double *sca;
    double *vector_orthogonal_to_random_hyperplane;
    int *number_of_facilities_on_the_right;

    /* L-BFGS-B variables */
	double gradInorm; // norm of the gradient of f for inequality constraints
	double gradEnorm; // norm of the gradient of f for equality constraints
	double *g; // gradient
	double f; // function value
	double *y; // dual variables
	double *binf; // lower bounds on the variables y
	double *bsup; // upper bounds on the variables y
	int *nbd; // indicates which variables are bounded
	double *wa;  // Double workspace for L-BFGS-B
	int *iwa; // Integer workspace for L-BFGS-B

	/* dsyevr function variables (projection) */
	int M; // number of eigenvalues
	double *W; // contains the eigenvalues
	double *Z; // contains the eigenvectors
	int *ISUPPZ; // dim = 2*max(1,M), support of the eigvecs. in Z
	double *WORK; // dim = LWORK
	int *IWORK; // dim = LIWORK
	int sizeWORK; // size of WORK
	int sizeIWORK; // size of IWORK
} Subproblem_Data;

#endif
