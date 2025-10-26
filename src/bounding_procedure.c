#include "constants_and_macros.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <mkl.h>
#include <assert.h>

#include "b_and_b_framework.h"
#include "user_structs.h"
#include "bound.h"
#include "bounding_procedure.h"
#include "build_subproblem_data.h"
#include "heuristics.h"
#include "inequalities.h"

// es wird angenommen, dass alle impliziten Fixierungen im Teilproblem vorgenommen worden sind
// soll Funktion double oder void zurueckgeben?
// jedes mal, wenn "sim" aufgerufen wird, versuchen duale Schranke zu aktualisieren
void SDPbound(Subproblem *subproblem) {
		
	assert(true == check_instance_consistency(subproblem->instance));

	Subproblem_Data *subproblem_data = subproblem->subproblem_data;

	// pruefe, ob ein Blatt des BB-Baums vorliegt
    if (1 == subproblem_data->dimension_of_Q) {

        Instance *instance = subproblem->instance;
        double value = calculate_objective_value(instance, subproblem->fixed_values);

        try_to_update_local_and_global_primalbound_and_solution(subproblem, subproblem->fixed_values, value);

        try_to_update_local_dualbound(subproblem, subproblem_data->fixed_objective_value);

        return;
    }

	// Parameter
    long double alpha = ALPHA_MAX;
    double tol = TOLERANCE_MAX;
    
    subproblem_data->number_of_cuts = 0;

    // initial dual vector y = 0
    int ylength = subproblem_data->number_of_equations + subproblem_data->number_of_cuts;
    for (int i = 0; i < ylength; i++) {
        (subproblem_data->y)[i] = 0.0;
    }

    sim(subproblem, alpha); // bestimmt Funktionswert und Gradient

    int count = 0; // major iteration counter
    int nbitalpha = 0; // current number of iterations for a given alpha
    int nbit; // number of iterations in BFGS

	bool done = false;

    // Main loop
    while (!done) {

        // Update iteration counter
        count++;
        nbitalpha++;

        // Call BFGS solver
	    bool BFGSsuccess = calllbfgsb(subproblem, alpha, tol, &nbit);
	    
	    done = !BFGSsuccess || (alpha == ALPHA_MIN) || (true ==  can_subproblem_be_pruned(subproblem));

		if (!done && ((3 == HEURISTIC_FREQUENCY) || (2 == HEURISTIC_FREQUENCY && 1 == nbitalpha))) {
			run_heuristics(subproblem);
		}

	    done = !BFGSsuccess || (alpha == ALPHA_MIN) || (true ==  can_subproblem_be_pruned(subproblem));
	    
	    int number_of_subtracted_cuts = 0;
	    int number_of_added_cuts = 0;
	    if (!done && 1 == WITH_CUTS) {
	    	updateInequalities(subproblem, &number_of_subtracted_cuts, &number_of_added_cuts);
    	}

        if (!done && MINIMUM_NUMBER_OF_CUTS_PER_ITERATION > number_of_added_cuts) {
        	if (1 == HEURISTIC_FREQUENCY) {
        		run_heuristics(subproblem);
        	}
        	//printf("alpha = %Lf und dualbound = %f\n", alpha, subproblem->instance->instance_data->K - subproblem->local_dualbound);
            nbitalpha = 0;
            alpha *= SCALE_ALPHA;
            alpha = (alpha < ALPHA_MIN) ? ALPHA_MIN : alpha;
            tol *= SCALE_TOLERANCE;
            tol = (tol < TOLERANCE_MIN) ? TOLERANCE_MIN : tol;          
        }

    } // end main while loop

    try_to_update_local_dualbound(subproblem, (subproblem_data->fixed_objective_value) + (subproblem_data->f));
    
    if (0 < HEURISTIC_FREQUENCY) {
		run_heuristics(subproblem);
	}

}



/************************ calllbfgs *********************************/
int calllbfgsb(Subproblem *subproblem, long double alpha, double tol, int *nbit) {

	//**************************************************
	// "hole" alle wichtigen Sachen
	Subproblem_Data *subproblem_data = subproblem->subproblem_data;
	int N = subproblem_data->dimension_of_Q;
	double *g = subproblem_data->g;
	double *f = &(subproblem_data->f);
	double *Z = subproblem_data->Z;
	double *X = subproblem_data->X;
	int *M = &(subproblem_data->M);
	double *y = subproblem_data->y;
	double *binf = subproblem_data->binf;
	double *bsup = subproblem_data->bsup;
	int *nbd = subproblem_data->nbd;
	double *wa = subproblem_data->wa;
	int *iwa = subproblem_data->iwa;
	
	Cut *current_cuts = subproblem_data->current_cuts;


    char task[60], csave[60];
    int lsave[4], isave[44];
    double dsave[29];
    int status = 1;

    int iprint = -1;
    int mem = MAXMIUM_NUMBER_OF_MEMORY_CORRECTIONS; // memory corrections

    // Compute the number of variables in BFGS
    int number_of_equations = subproblem_data->number_of_equations;
    int number_of_cuts = subproblem_data->number_of_cuts;
    int nvars = number_of_equations + number_of_cuts;
    
    // Initialize y
    for (int ineq = 0; ineq < number_of_cuts; ineq++) {
        y[number_of_equations + ineq] = current_cuts[ineq].y;
    }

    double factr = 0.0; // CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH
    double pgtol = 0.0;

    // evaluate f,g with updated y, alpha, and Inequalities (stored in vector Cuts)
    sim(subproblem, alpha); // bestimmt Funktionswert und Gradient

    /*
     * Compute the bounds for y
     */


    // equalities
    for (int i = 0; i < number_of_equations; i++) {
        nbd[i] = 0; // y[i] unbounded
    }
    
    //  inequalities
    for (int i = number_of_equations; i < nvars; i++) {
        nbd[i] = 1; // y[i] bounded from below
        binf[i] = 0.0; // lower bound is zero
    }


    strcpy(task, "START");
    for (int i = 5; i < 60; i++) {
    	task[i] = ' ';
    }
    int StopBFGS = 0;
    *nbit = 0; // counts number of function/gradient evaluates
    while (!StopBFGS) {
        factr = 0.0; // CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH
        pgtol = 0.0; // NORM OF PROJECTED GRADIENT <= PGTOL

        /* This calls the main L-BFGS-B function */
        setulb_(&nvars, &mem, y, binf, bsup, nbd, f, g, &factr, &pgtol, wa, iwa, task, &iprint, csave, lsave, isave, dsave);

        if (0 == strncmp(task, "FG", 2)) { // L-BFGS-B requesting new f and g
            sim(subproblem, alpha); // bestimmt Funktionswert und Gradient
            (*nbit)++;
        } else if (0 == strncmp(task, "NEW_X", 5)) { // L-BFGS-B found new x
        
            // test if we should stop

            // aktualisiere "gradEnorm"
            double gradEnorm = 0.0;
            for (int i = 0; i < number_of_equations; i++) {
                gradEnorm = (gradEnorm < fabs(g[i])) ? fabs(g[i]) : gradEnorm;
            }
			subproblem_data->gradEnorm = gradEnorm;
			
			// aktualisiere "gradInorm"
			double gradInorm = 0.0;
            for (int i = number_of_equations; i < nvars; i++) {
                double mintemp = (g[i] > 0.0) ? 0.0 : g[i];
                gradInorm = (gradInorm < fabs(mintemp)) ? fabs(mintemp) : gradInorm;
            }
           	gradInorm /= SCALE_FACTOR_OF_TRIANGLE_INEQUALITIES;
			subproblem_data->gradInorm = gradInorm;

            if ((gradEnorm < tol && gradInorm < tol && *nbit >= 100) || (true == can_subproblem_be_pruned(subproblem)) || ((*nbit) >= MAXMIUM_NUMBER_OF_EVALUATIONS_PER_ITERATION)) {
                strcpy(task, "STOP");
                for (int i = 4; i < 60; i++) {
                	task[i] = ' ';
                }
            }
            
        }
        else { // L-BFGS-B has terminated
            StopBFGS = 1;
            if (0 != strncmp(task, "STOP", 4)) {
                status = 0;
            }
        }
    }


    // Scale: X = X/alpha
    double alphainv = (double) 1.0L / alpha;
    int N2 = N * N;
    int incx = 1;
    dscal(&N2, &alphainv, X, &incx);  // X = X/alpha


    // Update Z so that X = Z*Z'
    double sca = sqrt(alphainv);
    int NM = N * (*M);
    int incz = 1;
    dscal(&NM, &sca, Z, &incz);  // Z = Z/sqrt(alpha)

    return status;

}
