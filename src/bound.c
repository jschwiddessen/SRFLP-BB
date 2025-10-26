#include "constants_and_macros.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mkl.h>

#include "bound.h"
#include "user_structs.h"

// computes gradient and function values for x
/*************** Simulator (evaluates function and gradient) ******************/
void sim(Subproblem *subproblem, long double alpha) {

    // ***************************************************************************
    // ""hole"" alle wichtigen Sachen
    Subproblem_Data *subproblem_data = subproblem->subproblem_data;
    int N = subproblem_data->dimension_of_Q;
    double *Q = subproblem_data->Q;
    double *X = subproblem_data->X;
    double *ff = &(subproblem_data->f);
    // ***************************************************************************

    int N2 = N * N;

    int incx = 1;
    int incy = 1;
    dcopy(&N2, Q, &incx, X, &incy);

    B(true, subproblem_data, alpha);
    A(true, subproblem_data, alpha);

    // Overwrites X with X_+
    projSDP(subproblem_data); // returns *ff = 0.5*|X_+|_F^2

    //double alphainv = 1.0 / alpha;
    *ff = (*ff / alpha); // now *ff = 0.5*|X_+|_F^2/alpha

    B(false, subproblem_data, alpha);
    A(false, subproblem_data, alpha);

    *ff += alpha * (0.5L * N2);
	
    try_to_update_local_dualbound(subproblem, (subproblem_data->fixed_objective_value) + (*ff));
    (subproblem_data->evaluations)++;

}


/***********************  projSDP  ************************/
void projSDP(Subproblem_Data *subproblem_data) {

    // ***************************************************************************
    // ""hole"" alle wichtigen Sachen
    int N = subproblem_data->dimension_of_Q;
    double *X = subproblem_data->X;
    double *WORK = subproblem_data->WORK;
    int sizeWORK = subproblem_data->sizeWORK;
    int sizeIWORK = subproblem_data->sizeIWORK;
    int *M = &(subproblem_data->M);
    double *W = subproblem_data->W;
    double *Z = subproblem_data->Z;
    int *ISUPPZ = subproblem_data->ISUPPZ;
    int *IWORK = subproblem_data->IWORK;
    double *ff = &(subproblem_data->f);
    // ***************************************************************************

    double borneinf = 0.0;

    char UPLO = 'L';
    int LDX = N;

    // Operator norm for X to get a upper bound of the spectral radius of X
    // |X|_2 \leq \sqrt{ |X|_1 |X|_inf }  (Holder's inequality)
    //          = |X|_1 = |X|_inf  (since X is symmetric)
    //
    // Frobenius norm is also an upper bound on the spectral radius of X:
    //      |X|_2 <= |X|_F
    char NORM = 'I';
    double norminf = dlansy(&NORM, &UPLO, &N, X, &LDX, WORK);
    NORM = 'F';
    double normfro = dlansy(&NORM, &UPLO, &N, X, &LDX, WORK);

    // bornesup = min(norminf, normfro)
    double bornesup = (norminf < normfro) ? norminf : normfro;

    // Ensure that borneinf <= bornesup.
    if (bornesup < borneinf) {
        bornesup = 2.0 * borneinf;
    }

    bornesup += 1e-8;

    /* Compute the positive eigenvalues and associated eigenvectors of X.
     *
     * The M columns of Z will contain the orthonormal eigenvectors of the
     * matrix X corresponding to the positive eigenvalues, the i-th column
     * of Z holding the eigenvector associated with W[i].
     */

    char JOBZ = 'V';
    char RANGE = 'V';
    double VL = borneinf;
    double VU = bornesup;
    int IL = 0;
    int IU = 0;
    double ABSTOL = 1e-42;
    int LDZ = N;
    int LWORK = sizeWORK;
    int LIWORK = sizeIWORK;
    int INFO;
    //printf("Eigenwert Start\n");
    dsyevr(&JOBZ, &RANGE, &UPLO, &N, X, &LDX, &VL, &VU, &IL, &IU, &ABSTOL,
            M, W, Z, &LDZ, ISUPPZ, WORK, &LWORK, IWORK, &LIWORK, &INFO);
    //printf("Eigenwert Ende\n");


    // Check if the eigensolver failed (i.e., if INFO != 0)

    if (INFO) {
        fprintf(stderr,
                "Error: eigenvalue computation failed (INFO = %d)\n",
                INFO);
        exit(1);
    }

    // Compute ff = 0.5*||X_+||^2 = 0.5*||W||^2
    int INCW = 1;
    double normW = dnrm2(M, W, &INCW);
    *ff = 0.5 * (normW * normW);

    // Compute Z = Z*Diag(W)^{1/2}
    int INCZ = 1;
    for (int j = 0; j < (*M); j++) {
        // Scale jth column of Z by sqrt(W[j])
        double temp = sqrt(W[j]);
        dscal(&N, &temp, Z + j * N, &INCZ);
    }

    char TRANS = 'N';
    double ALPHA = 1.0;
    double BETA = 0.0;

    /* X = ALPHA*Z*Z^T + BETA*X = Z*Z^T
     * Only writes to lower-triangular part of X.
     * When *M = 0, we will obtain X = 0.
     */
    dsyrk(&UPLO, &TRANS, &N, M, &ALPHA, Z, &LDZ, &BETA, X, &LDX);

}



void A(bool transposed, Subproblem_Data *subproblem_data, long double alpha) {

	//double alphainv = 1.0 / alpha;

    int N = subproblem_data->dimension_of_Q;
    int number_of_equations = subproblem_data->number_of_equations;
    int number_of_cuts = subproblem_data->number_of_cuts;
    double *ff = &(subproblem_data->f);
    double *yy = subproblem_data->y;
    double *gg = subproblem_data->g;
    double *X = subproblem_data->X;
    Cut *current_cuts = subproblem_data->current_cuts;

    if (true == transposed) {

        for (int ineq = 0; ineq < number_of_cuts; ineq++) {

            int type = current_cuts[ineq].type;
            int ii = current_cuts[ineq].i;
            int jj = current_cuts[ineq].j;
            int kk = current_cuts[ineq].k;
            int ll = current_cuts[ineq].l;
            int mm = current_cuts[ineq].m;
            
            long double temp = 0.0;
                     
            if (type <= 4) {
				temp = 0.5L * SCALE_FACTOR_OF_TRIANGLE_INEQUALITIES * yy[number_of_equations + ineq];
			} else if (type <= 20) {
				temp = 0.5L * SCALE_FACTOR_OF_PENTAGONAL_INEQUALITIES * yy[number_of_equations + ineq];
			}
			
            switch (type) {
                case 1:
                    X[ii + jj * N] += temp;
                    X[ii + kk * N] += temp;
                    X[jj + kk * N] += temp;
                    break;
                case 2:
                    X[ii + jj * N] += temp;
                    X[ii + kk * N] -= temp;
                    X[jj + kk * N] -= temp;
                    break;
                case 3:
                    X[ii + jj * N] -= temp;
                    X[ii + kk * N] += temp;
                    X[jj + kk * N] -= temp;
                    break;
                case 4:
                    X[ii + jj * N] -= temp;
                    X[ii + kk * N] -= temp;
                    X[jj + kk * N] += temp;
                    break;
				case 5:
					X[ii + jj * N] += temp;
					X[ii + kk * N] += temp;
					X[ii + ll * N] += temp;
					X[ii + mm * N] += temp;
					X[jj + kk * N] += temp;
					X[jj + ll * N] += temp;
					X[jj + mm * N] += temp;
					X[kk + ll * N] += temp;
					X[kk + mm * N] += temp;
					X[ll + mm * N] += temp;
					break;
				case 6:
					X[ii + jj * N] += temp;
					X[ii + kk * N] += temp;
					X[ii + ll * N] += temp;
					X[ii + mm * N] -= temp;
					X[jj + kk * N] += temp;
					X[jj + ll * N] += temp;
					X[jj + mm * N] -= temp;
					X[kk + ll * N] += temp;
					X[kk + mm * N] -= temp;
					X[ll + mm * N] -= temp;
					break;
				case 7:
					X[ii + jj * N] += temp;
					X[ii + kk * N] += temp;
					X[ii + ll * N] -= temp;
					X[ii + mm * N] += temp;
					X[jj + kk * N] += temp;
					X[jj + ll * N] -= temp;
					X[jj + mm * N] += temp;
					X[kk + ll * N] -= temp;
					X[kk + mm * N] += temp;
					X[ll + mm * N] -= temp;
					break;
				case 8:
					X[ii + jj * N] += temp;
					X[ii + kk * N] += temp;
					X[ii + ll * N] -= temp;
					X[ii + mm * N] -= temp;
					X[jj + kk * N] += temp;
					X[jj + ll * N] -= temp;
					X[jj + mm * N] -= temp;
					X[kk + ll * N] -= temp;
					X[kk + mm * N] -= temp;
					X[ll + mm * N] += temp;
					break;
				case 9:
					X[ii + jj * N] += temp;
					X[ii + kk * N] -= temp;
					X[ii + ll * N] += temp;
					X[ii + mm * N] += temp;
					X[jj + kk * N] -= temp;
					X[jj + ll * N] += temp;
					X[jj + mm * N] += temp;
					X[kk + ll * N] -= temp;
					X[kk + mm * N] -= temp;
					X[ll + mm * N] += temp;
					break;
				case 10:
					X[ii + jj * N] += temp;
					X[ii + kk * N] -= temp;
					X[ii + ll * N] += temp;
					X[ii + mm * N] -= temp;
					X[jj + kk * N] -= temp;
					X[jj + ll * N] += temp;
					X[jj + mm * N] -= temp;
					X[kk + ll * N] -= temp;
					X[kk + mm * N] += temp;
					X[ll + mm * N] -= temp;
					break;
				case 11:
					X[ii + jj * N] += temp;
					X[ii + kk * N] -= temp;
					X[ii + ll * N] -= temp;
					X[ii + mm * N] += temp;
					X[jj + kk * N] -= temp;
					X[jj + ll * N] -= temp;
					X[jj + mm * N] += temp;
					X[kk + ll * N] += temp;
					X[kk + mm * N] -= temp;
					X[ll + mm * N] -= temp;
					break;
				case 12:
					X[ii + jj * N] += temp;
					X[ii + kk * N] -= temp;
					X[ii + ll * N] -= temp;
					X[ii + mm * N] -= temp;
					X[jj + kk * N] -= temp;
					X[jj + ll * N] -= temp;
					X[jj + mm * N] -= temp;
					X[kk + ll * N] += temp;
					X[kk + mm * N] += temp;
					X[ll + mm * N] += temp;
					break;
				case 13:
					X[ii + jj * N] -= temp;
					X[ii + kk * N] += temp;
					X[ii + ll * N] += temp;
					X[ii + mm * N] += temp;
					X[jj + kk * N] -= temp;
					X[jj + ll * N] -= temp;
					X[jj + mm * N] -= temp;
					X[kk + ll * N] += temp;
					X[kk + mm * N] += temp;
					X[ll + mm * N] += temp;
					break;
				case 14:
					X[ii + jj * N] -= temp;
					X[ii + kk * N] += temp;
					X[ii + ll * N] += temp;
					X[ii + mm * N] -= temp;
					X[jj + kk * N] -= temp;
					X[jj + ll * N] -= temp;
					X[jj + mm * N] += temp;
					X[kk + ll * N] += temp;
					X[kk + mm * N] -= temp;
					X[ll + mm * N] -= temp;
					break;
				case 15:
					X[ii + jj * N] -= temp;
					X[ii + kk * N] += temp;
					X[ii + ll * N] -= temp;
					X[ii + mm * N] += temp;
					X[jj + kk * N] -= temp;
					X[jj + ll * N] += temp;
					X[jj + mm * N] -= temp;
					X[kk + ll * N] -= temp;
					X[kk + mm * N] += temp;
					X[ll + mm * N] -= temp;
					break;
				case 16:
					X[ii + jj * N] -= temp;
					X[ii + kk * N] += temp;
					X[ii + ll * N] -= temp;
					X[ii + mm * N] -= temp;
					X[jj + kk * N] -= temp;
					X[jj + ll * N] += temp;
					X[jj + mm * N] += temp;
					X[kk + ll * N] -= temp;
					X[kk + mm * N] -= temp;
					X[ll + mm * N] += temp;
					break;
				case 17:
					X[ii + jj * N] -= temp;
					X[ii + kk * N] -= temp;
					X[ii + ll * N] += temp;
					X[ii + mm * N] += temp;
					X[jj + kk * N] += temp;
					X[jj + ll * N] -= temp;
					X[jj + mm * N] -= temp;
					X[kk + ll * N] -= temp;
					X[kk + mm * N] -= temp;
					X[ll + mm * N] += temp;
					break;
				case 18:
					X[ii + jj * N] -= temp;
					X[ii + kk * N] -= temp;
					X[ii + ll * N] += temp;
					X[ii + mm * N] -= temp;
					X[jj + kk * N] += temp;
					X[jj + ll * N] -= temp;
					X[jj + mm * N] += temp;
					X[kk + ll * N] -= temp;
					X[kk + mm * N] += temp;
					X[ll + mm * N] -= temp;
					break;
				case 19:
					X[ii + jj * N] -= temp;
					X[ii + kk * N] -= temp;
					X[ii + ll * N] -= temp;
					X[ii + mm * N] += temp;
					X[jj + kk * N] += temp;
					X[jj + ll * N] += temp;
					X[jj + mm * N] -= temp;
					X[kk + ll * N] += temp;
					X[kk + mm * N] -= temp;
					X[ll + mm * N] -= temp;
					break;
				case 20:
					X[ii + jj * N] -= temp;
					X[ii + kk * N] -= temp;
					X[ii + ll * N] -= temp;
					X[ii + mm * N] -= temp;
					X[jj + kk * N] += temp;
					X[jj + ll * N] += temp;
					X[jj + mm * N] += temp;
					X[kk + ll * N] += temp;
					X[kk + mm * N] += temp;
					X[ll + mm * N] += temp;
					break;
            }
        }

    } else {

        for (int ineq = 0; ineq < number_of_cuts; ineq++) {

            int type = current_cuts[ineq].type;
            int ii = current_cuts[ineq].i;
            int jj = current_cuts[ineq].j;
            int kk = current_cuts[ineq].k;
            int ll = current_cuts[ineq].l;
            int mm = current_cuts[ineq].m;
          
            if (type <= 4) {
				*ff += SCALE_FACTOR_OF_TRIANGLE_INEQUALITIES * yy[number_of_equations + ineq];
			} else if (type <= 20) {
				*ff += 2.0L * SCALE_FACTOR_OF_PENTAGONAL_INEQUALITIES * yy[number_of_equations + ineq];
			}
            
			long double temp = 0.0;
						
            switch (type) {
                case 1:
                    temp = + X[ii + jj * N] + X[ii + kk * N] + X[jj + kk * N];
                    break;
                case 2:
                    temp = + X[ii + jj * N] - X[ii + kk * N] - X[jj + kk * N];
                    break;
                case 3:
                    temp = - X[ii + jj * N] + X[ii + kk * N] - X[jj + kk * N];
                    break;
                case 4:
                    temp = - X[ii + jj * N] - X[ii + kk * N] + X[jj + kk * N];
                    break;
				case 5:
					temp = + X[ii + jj * N] + X[ii + kk * N] + X[ii + ll * N] + X[ii + mm * N] + X[jj + kk * N] + X[jj + ll * N] + X[jj + mm * N] + X[kk + ll * N] + X[kk + mm * N] + X[ll + mm * N];
					break;
				case 6:
					temp = + X[ii + jj * N] + X[ii + kk * N] + X[ii + ll * N] - X[ii + mm * N] + X[jj + kk * N] + X[jj + ll * N] - X[jj + mm * N] + X[kk + ll * N] - X[kk + mm * N] - X[ll + mm * N];
					break;
				case 7:
					temp = + X[ii + jj * N] + X[ii + kk * N] - X[ii + ll * N] + X[ii + mm * N] + X[jj + kk * N] - X[jj + ll * N] + X[jj + mm * N] - X[kk + ll * N] + X[kk + mm * N] - X[ll + mm * N];
					break;
				case 8:
					temp = + X[ii + jj * N] + X[ii + kk * N] - X[ii + ll * N] - X[ii + mm * N] + X[jj + kk * N] - X[jj + ll * N] - X[jj + mm * N] - X[kk + ll * N] - X[kk + mm * N] + X[ll + mm * N];
					break;
				case 9:
					temp = + X[ii + jj * N] - X[ii + kk * N] + X[ii + ll * N] + X[ii + mm * N] - X[jj + kk * N] + X[jj + ll * N] + X[jj + mm * N] - X[kk + ll * N] - X[kk + mm * N] + X[ll + mm * N];
					break;
				case 10:
					temp = + X[ii + jj * N] - X[ii + kk * N] + X[ii + ll * N] - X[ii + mm * N] - X[jj + kk * N] + X[jj + ll * N] - X[jj + mm * N] - X[kk + ll * N] + X[kk + mm * N] - X[ll + mm * N];
					break;
				case 11:
					temp = + X[ii + jj * N] - X[ii + kk * N] - X[ii + ll * N] + X[ii + mm * N] - X[jj + kk * N] - X[jj + ll * N] + X[jj + mm * N] + X[kk + ll * N] - X[kk + mm * N] - X[ll + mm * N];
					break;
				case 12:
					temp = + X[ii + jj * N] - X[ii + kk * N] - X[ii + ll * N] - X[ii + mm * N] - X[jj + kk * N] - X[jj + ll * N] - X[jj + mm * N] + X[kk + ll * N] + X[kk + mm * N] + X[ll + mm * N];
					break;
				case 13:
					temp = - X[ii + jj * N] + X[ii + kk * N] + X[ii + ll * N] + X[ii + mm * N] - X[jj + kk * N] - X[jj + ll * N] - X[jj + mm * N] + X[kk + ll * N] + X[kk + mm * N] + X[ll + mm * N];
					break;
				case 14:
					temp = - X[ii + jj * N] + X[ii + kk * N] + X[ii + ll * N] - X[ii + mm * N] - X[jj + kk * N] - X[jj + ll * N] + X[jj + mm * N] + X[kk + ll * N] - X[kk + mm * N] - X[ll + mm * N];
					break;
				case 15:
					temp = - X[ii + jj * N] + X[ii + kk * N] - X[ii + ll * N] + X[ii + mm * N] - X[jj + kk * N] + X[jj + ll * N] - X[jj + mm * N] - X[kk + ll * N] + X[kk + mm * N] - X[ll + mm * N];
					break;
				case 16:
					temp = - X[ii + jj * N] + X[ii + kk * N] - X[ii + ll * N] - X[ii + mm * N] - X[jj + kk * N] + X[jj + ll * N] + X[jj + mm * N] - X[kk + ll * N] - X[kk + mm * N] + X[ll + mm * N];
					break;
				case 17:
					temp = - X[ii + jj * N] - X[ii + kk * N] + X[ii + ll * N] + X[ii + mm * N] + X[jj + kk * N] - X[jj + ll * N] - X[jj + mm * N] - X[kk + ll * N] - X[kk + mm * N] + X[ll + mm * N];
					break;
				case 18:
					temp = - X[ii + jj * N] - X[ii + kk * N] + X[ii + ll * N] - X[ii + mm * N] + X[jj + kk * N] - X[jj + ll * N] + X[jj + mm * N] - X[kk + ll * N] + X[kk + mm * N] - X[ll + mm * N];
					break;
				case 19:
					temp = - X[ii + jj * N] - X[ii + kk * N] - X[ii + ll * N] + X[ii + mm * N] + X[jj + kk * N] + X[jj + ll * N] - X[jj + mm * N] + X[kk + ll * N] - X[kk + mm * N] - X[ll + mm * N];
					break;
				case 20:
					temp = - X[ii + jj * N] - X[ii + kk * N] - X[ii + ll * N] - X[ii + mm * N] + X[jj + kk * N] + X[jj + ll * N] + X[jj + mm * N] + X[kk + ll * N] + X[kk + mm * N] + X[ll + mm * N];
					break;
            }
            
			if (type <= 4) {
				gg[number_of_equations + ineq] = (temp / alpha + 1.0L) * SCALE_FACTOR_OF_TRIANGLE_INEQUALITIES;
			} else if (type <= 20) {
				gg[number_of_equations + ineq] = (temp / alpha + 2.0L) * SCALE_FACTOR_OF_PENTAGONAL_INEQUALITIES;
			}
        }

    } // end if

}


/***************** B *********************/
void B(bool transposed, Subproblem_Data *subproblem_data, long double alpha) {

    //double alphainv = 1.0 / alpha;

    int N = subproblem_data->dimension_of_Q;
    int number_of_equations = subproblem_data->number_of_equations;
    Sparse_Matrix *Bs = subproblem_data->Bs;
    double *b = subproblem_data->b;
    double *ff = &(subproblem_data->f);
    double *yy = subproblem_data->y;
    double *gg = subproblem_data->g;
    double *X = subproblem_data->X;

    if (true == transposed) {
        for (int eq = 0; eq < number_of_equations; eq++) { // for each equality constraint
            long double temp = yy[eq];

            // Using sparse representation
            for (int entry = 0; entry < Bs[eq].nnz; entry++) { // for each nonzero entry
                int ii = Bs[eq].i[entry];
                int jj = Bs[eq].j[entry];
                double dd = Bs[eq].val[entry];
                if (ii == jj) {
                    X[ii + jj * N] -= temp * dd;
                } else {
                    X[ii + jj * N] -= temp * dd;
                    X[jj + ii * N] -= temp * dd;
                }
            }
        }

    } else {

        for (int eq = 0; eq < number_of_equations; eq++) { // for each equality constraint

            *ff += b[eq] * yy[eq];

            gg[eq] = 0.0;
            for (int entry = 0; entry < Bs[eq].nnz; entry++) { // for each nonzero entry
                int ii = Bs[eq].i[entry];
                int jj = Bs[eq].j[entry];
                double dd = Bs[eq].val[entry];
                if (ii == jj) {
                    gg[eq] -= dd * X[ii + jj * N];
                } else {
                    gg[eq] -= 2.0 * dd * X[ii + jj * N];
                }
            }
            gg[eq] /= alpha;
            gg[eq] += b[eq];
        }

    } // end if
}
