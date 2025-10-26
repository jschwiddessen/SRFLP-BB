#ifndef BOUNDING_PROCEDURE_H
#define BOUNDING_PROCEDURE_H

#include <stdlib.h>

#include "srflp.h"

void SDPbound(Subproblem *subproblem);
int calllbfgsb(Subproblem *subproblem, long double alpha, double tol, int *nbit);

void setulb_(
        int *n, // number of variables
        int *m, // size of the limited memory
        double *x, // current iterate (length n)
        double *l, // lower bounds (length n)
        double *u, // upper bounds (length n)
        int *nbd, // indicates which vars are bounded
        // nbd[i] = 0 : x[i] unbounded
        // nbd[i] = 1 : x[i] has only a lower bound
        // nbd[i] = 2 : x[i] has both lower and upper bounds
        // nbd[i] = 3 : x[i] has only an upper bound
        double *f, // value of function at the point x
        double *g, // value of the gradient at the point x (length n)
        double *factr, // termination tolerance
        // factr = 1.d+12 for low accuracy;
        //         1.d+7  for moderate accuracy;
        //         1.d+1  for extremely high accuracy.
        double *pgtol, // projected gradient tolerance
        // (suppress this termination test by pgtol = 0)
        double *wa, // workspace (length (2MAXMIUM_NUMBER_OF_MEMORY_CORRECTIONS + 4)nmax + 12MAXMIUM_NUMBER_OF_MEMORY_CORRECTIONS^2 + 12MAXMIUM_NUMBER_OF_MEMORY_CORRECTIONS)
        int *iwa, // workspace (length 3nmax)
        char *task, // character array (length 60)
        // 'START'    : when starting
        // 'FG'        : user must evaluate function f and gradient g
        // 'NEW_X'    : user can decide whether to continue or stop
        // 'CONV'    : termination test has been satisfied
        // 'ABNO'    : abnormal termination
        // 'ERROR'    : error in input parameters
        // 'STOP'    : set by user to stop L-BFGS-B
        int *iprint, // set level of output
        // iprint<0    no output is generated;
        // iprint=0    print only one line at the last iteration;
        // 0<iprint<99 print also f and |proj g| every iprint iterations;
        // iprint=99   print details of every iteration except n-vectors;
        // iprint=100  print also the changes of active set and final x;
        // iprint>100  print details of every iteration including x and g;
        char *csave, // character array (length 60)
        int *lsave, // logicial array (length 4)
        int *isave, // integer array (length 44)
        double *dsave // double array (length 29)
        );

#endif
