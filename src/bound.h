#ifndef BOUND_H
#define BOUND_H

#include <stdbool.h>

#include "user_structs.h"
#include "b_and_b_framework.h"

/*
extern double dnrm2_(int *n, double *x, int *incx);
extern void dscal_(int *n, double *da, double *x, int *incx);
extern void dcopy_(int *n, double *dx, int *incx, double *dy, int *incy);
extern void dsyevr_(char *JOBZ, char *RANGE, char *UPLO, int *N, double *A, int *LDA, double *VL, double *VU, int *IL, int *IU, double *ABSTOL, int *M, double *W, double *Z, int *LDZ, int *ISUPPZ, double *WORK, int *LWORK, int *IWORK, int *LIWORK, int *INFO);
extern void dsyrk_(char *UPLO, char *TRANS, int *N, int *K, double *ALPHA, double *A, int *LDA, double *BETA, double *C, int *LDC);
extern double dlansy_(char *NORM, char *UPLO, int *N, double *A, int *LDA, double *WORK);
*/

void sim(Subproblem *subproblem, long double alpha);
void projSDP(Subproblem_Data *subproblem_data);
void A(bool transposed, Subproblem_Data *subproblem_data, long double alpha);
void B(bool transposed, Subproblem_Data *subproblem_data, long double alpha);

#endif
