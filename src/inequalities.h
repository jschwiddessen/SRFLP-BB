#ifndef INEQUALITIES_H
#define INEQUALITIES_H

#include "b_and_b_framework.h"

int getViolatedCuts(Subproblem *subproblem);
void updateInequalities(Subproblem *subproblem, int *number_of_subtracted_cuts, int *number_of_added_cuts);
double eval_ineq(double *XX, int N, int type, int ii, int jj, int kk, int ll, int mm);


#endif

