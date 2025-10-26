#include "constants_and_macros.h"

#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <float.h>
#include <stdio.h>

#include "inequalities.h"
#include "user_structs.h"
#include "miscellaneous_functions.h"
#include "b_and_b_framework.h"


int getViolatedCuts(Subproblem *subproblem) {

	Subproblem_Data *subproblem_data = subproblem->subproblem_data;
	double *X = subproblem_data->X;
	int N = subproblem_data->dimension_of_Q;
	Cut *new_cuts = subproblem_data->new_cuts;

    int LeastViolatedIneq = 0;
    double LeastViolatedIneqValue = - DBL_MAX;

    int size = 0;
    for (int type = 1; type <= 4; type++) {
        for (int ii = 0; ii < N;  ii++) {
            for (int jj = 0; jj < ii; jj++) {
                for (int kk = 0; kk < jj; kk++) {
                    double testineqvalue = eval_ineq(X, N, type, ii, jj, kk, -1, -1);
                    if (testineqvalue < GAP_CUTS) {
                        if (size < MAXIMUM_NUMBER_OF_CUTS_PER_ITERATION) {
                            new_cuts[size].type = type;
                            new_cuts[size].i = ii;
                            new_cuts[size].j = jj;
                            new_cuts[size].k = kk;
                            new_cuts[size].l = -1;
                            new_cuts[size].m = -1;
                            new_cuts[size].value = testineqvalue;
                            if (0 == size || testineqvalue > LeastViolatedIneqValue) {
                                LeastViolatedIneqValue = testineqvalue;
                                LeastViolatedIneq = size;
                            }
                            size++;
                        } else if (testineqvalue < LeastViolatedIneqValue) {
                            new_cuts[LeastViolatedIneq].type = type;
                            new_cuts[LeastViolatedIneq].i = ii;
                            new_cuts[LeastViolatedIneq].j = jj;
                            new_cuts[LeastViolatedIneq].k = kk;
                            new_cuts[LeastViolatedIneq].l = -1;
                            new_cuts[LeastViolatedIneq].m = -1;
                            new_cuts[LeastViolatedIneq].value = testineqvalue;
                            LeastViolatedIneqValue = - DBL_MAX;
                            for (int ListCount = 0; ListCount < size; ListCount++) {
                                if (new_cuts[ListCount].value > LeastViolatedIneqValue) {
                                    LeastViolatedIneqValue = new_cuts[ListCount].value;
                                    LeastViolatedIneq = ListCount;
                                }
                            }
                        }
                    }
                } // kk loop
            } // jj loop
        } // ii loop
    } // type loop
    
    
    Instance *instance = subproblem->instance;
    Instance_Data *instance_data = instance->instance_data;
    int number_of_facilities = instance_data->number_of_facilities;
    int *index_shift = subproblem_data->index_shift;
    bool *fixed_variables = subproblem->fixed_variables;
    
    for (int f1 = 0; f1 < number_of_facilities; f1++) {
    	for (int f2 = 0; f2 < f1; f2++) {
    		for (int f3 = 0; f3 < f2; f3++) {
    			for (int f4 = 0; f4 < f3; f4++) {
    				for (int f5 = 0; f5 < f4; f5++) {
    					for (int f6 = 0; f6 < number_of_facilities; f6++) {
    						if (f6 == f1 || f6 == f2 || f6 == f3 || f6 == f4 || f6 == f5) continue;
    						int pair1 = inMa(number_of_facilities, f1, f6);
    						int pair2 = inMa(number_of_facilities, f2, f6);
    						int pair3 = inMa(number_of_facilities, f3, f6);
    						int pair4 = inMa(number_of_facilities, f4, f6);
    						int pair5 = inMa(number_of_facilities, f5, f6);
    						if (fixed_variables[pair1] || fixed_variables[pair2] || fixed_variables[pair3] || fixed_variables[pair4] || fixed_variables[pair5]) continue;
    						int ii = index_shift[pair1];
    						int jj = index_shift[pair2];
    						int kk = index_shift[pair3];
    						int ll = index_shift[pair4];
    						int mm = index_shift[pair5];
						    for (int type = 5; type <= 20; type++) {

								double testineqvalue = eval_ineq(X, N, type, ii, jj, kk, ll, mm);
								if (testineqvalue < GAP_CUTS) {
									if (size < MAXIMUM_NUMBER_OF_CUTS_PER_ITERATION) {
									    new_cuts[size].type = type;
									    new_cuts[size].i = ii;
									    new_cuts[size].j = jj;
									    new_cuts[size].k = kk;
									    new_cuts[size].l = ll;
									    new_cuts[size].m = mm;
									    new_cuts[size].value = testineqvalue;
									    if (0 == size || testineqvalue > LeastViolatedIneqValue) {
									        LeastViolatedIneqValue = testineqvalue;
									        LeastViolatedIneq = size;
									    }
									    size++;
									} else if (testineqvalue < LeastViolatedIneqValue) {
									    new_cuts[LeastViolatedIneq].type = type;
									    new_cuts[LeastViolatedIneq].i = ii;
									    new_cuts[LeastViolatedIneq].j = jj;
									    new_cuts[LeastViolatedIneq].k = kk;
									    new_cuts[LeastViolatedIneq].l = ll;
									    new_cuts[LeastViolatedIneq].m = mm;
									    new_cuts[LeastViolatedIneq].value = testineqvalue;
									    LeastViolatedIneqValue = - DBL_MAX;
									    for (int ListCount = 0; ListCount < size; ListCount++) {
									        if (new_cuts[ListCount].value > LeastViolatedIneqValue) {
									            LeastViolatedIneqValue = new_cuts[ListCount].value;
									            LeastViolatedIneq = ListCount;
									        }
									    }
									}
								}

							} // type loop
    					}
    				}
    			}
    		}
    	}
    }
    
    
    
    
    

    return size;

}


/******************** Update Inequalities ***************************/
void updateInequalities(Subproblem *subproblem, int *number_of_subtracted_cuts, int *number_of_added_cuts) {

	Subproblem_Data *subproblem_data = subproblem->subproblem_data;
    int N = subproblem_data->dimension_of_Q;
    int number_of_equations = subproblem_data->number_of_equations;
    //int number_of_cuts = subproblem_data->number_of_cuts;
   
    Cut *current_cuts = subproblem_data->current_cuts;
    Cut *new_cuts = subproblem_data->new_cuts;
    double *X = subproblem_data->X;
    double *y = subproblem_data->y;
    double gradInorm = subproblem_data->gradInorm;
    
    int number_of_new_cuts = getViolatedCuts(subproblem);
    
    int subtracted = 0;
    int yindex = number_of_equations;
    int next_ineq = 0;
    for (int ineq = 0; ineq < subproblem_data->number_of_cuts; ineq++) {
    	int type = current_cuts[ineq].type;
        int ii = current_cuts[ineq].i;
        int jj = current_cuts[ineq].j;
        int kk = current_cuts[ineq].k;
        int ll = current_cuts[ineq].l;
        int mm = current_cuts[ineq].m;

        double ineqvalue = eval_ineq(X, N, type, ii, jj, kk, ll, mm);

        current_cuts[ineq].y = y[yindex];
        yindex++;

        if (current_cuts[ineq].y < 1e-8 && ineqvalue > gradInorm / 10.0) {
            subtracted++;
        } else {
            current_cuts[next_ineq].type = current_cuts[ineq].type;
            current_cuts[next_ineq].i    = current_cuts[ineq].i;
            current_cuts[next_ineq].j    = current_cuts[ineq].j;
            current_cuts[next_ineq].k    = current_cuts[ineq].k;
            current_cuts[next_ineq].l    = current_cuts[ineq].l;
            current_cuts[next_ineq].m    = current_cuts[ineq].m;
            current_cuts[next_ineq].y    = current_cuts[ineq].y;
            next_ineq++;
        }
    }
    subproblem_data->number_of_cuts -= subtracted;

    // Add List to Cuts
    int added = 0;
    for (int ListCount = 0; ListCount < number_of_new_cuts; ListCount++) { 
        if (MAXIMUM_NUMBER_OF_CUTS == next_ineq) break;

        // Check if inequality is already included in Cuts
        bool found_ineq = false;
        for (int ineq = 0; ineq < subproblem_data->number_of_cuts; ineq++) {
            if (current_cuts[ineq].type == new_cuts[ListCount].type &&
                	current_cuts[ineq].i == new_cuts[ListCount].i &&
                	current_cuts[ineq].j == new_cuts[ListCount].j &&
                	current_cuts[ineq].k == new_cuts[ListCount].k &&
                	current_cuts[ineq].l == new_cuts[ListCount].l &&
                	current_cuts[ineq].m == new_cuts[ListCount].m) {
                found_ineq = true;
            }
        }
        
        if (false == found_ineq) {
            current_cuts[next_ineq].type = new_cuts[ListCount].type;
            current_cuts[next_ineq].i = new_cuts[ListCount].i;
            current_cuts[next_ineq].j = new_cuts[ListCount].j;
            current_cuts[next_ineq].k = new_cuts[ListCount].k;
            current_cuts[next_ineq].l = new_cuts[ListCount].l;
            current_cuts[next_ineq].m = new_cuts[ListCount].m;
            current_cuts[next_ineq].y = 0.0;
            next_ineq++;
            added++;
        }
    }
    
    subproblem_data->number_of_cuts += added;
    
	*number_of_subtracted_cuts = subtracted;
	*number_of_added_cuts = added;
}


/************************ eval_ineq *********************************/
double eval_ineq(double *XX, int N, int type, int ii, int jj, int kk, int ll, int mm) {
    double ineqvalue = 0.0;

    // compute the inequality
    switch (type) {
        case 1:
            ineqvalue = + XX[ii + jj * N] + XX[ii + kk * N] + XX[jj + kk * N] + 1.0;
            break;
        case 2:
            ineqvalue = + XX[ii + jj * N] - XX[ii + kk * N] - XX[jj + kk * N] + 1.0;
            break;
        case 3:
            ineqvalue = - XX[ii + jj * N] + XX[ii + kk * N] - XX[jj + kk * N] + 1.0;
            break;
        case 4:
            ineqvalue = - XX[ii + jj * N] - XX[ii + kk * N] + XX[jj + kk * N] + 1.0;
            break;
		case 5:
			ineqvalue = + XX[ii + jj * N] + XX[ii + kk * N] + XX[ii + ll * N] + XX[ii + mm * N] + XX[jj + kk * N] + XX[jj + ll * N] + XX[jj + mm * N] + XX[kk + ll * N] + XX[kk + mm * N] + XX[ll + mm * N] + 2.0;
			break;
		case 6:
			ineqvalue = + XX[ii + jj * N] + XX[ii + kk * N] + XX[ii + ll * N] - XX[ii + mm * N] + XX[jj + kk * N] + XX[jj + ll * N] - XX[jj + mm * N] + XX[kk + ll * N] - XX[kk + mm * N] - XX[ll + mm * N] + 2.0;
			break;
		case 7:
			ineqvalue = + XX[ii + jj * N] + XX[ii + kk * N] - XX[ii + ll * N] + XX[ii + mm * N] + XX[jj + kk * N] - XX[jj + ll * N] + XX[jj + mm * N] - XX[kk + ll * N] + XX[kk + mm * N] - XX[ll + mm * N] + 2.0;
			break;
		case 8:
			ineqvalue = + XX[ii + jj * N] + XX[ii + kk * N] - XX[ii + ll * N] - XX[ii + mm * N] + XX[jj + kk * N] - XX[jj + ll * N] - XX[jj + mm * N] - XX[kk + ll * N] - XX[kk + mm * N] + XX[ll + mm * N] + 2.0;
			break;
		case 9:
			ineqvalue = + XX[ii + jj * N] - XX[ii + kk * N] + XX[ii + ll * N] + XX[ii + mm * N] - XX[jj + kk * N] + XX[jj + ll * N] + XX[jj + mm * N] - XX[kk + ll * N] - XX[kk + mm * N] + XX[ll + mm * N] + 2.0;
			break;
		case 10:
			ineqvalue = + XX[ii + jj * N] - XX[ii + kk * N] + XX[ii + ll * N] - XX[ii + mm * N] - XX[jj + kk * N] + XX[jj + ll * N] - XX[jj + mm * N] - XX[kk + ll * N] + XX[kk + mm * N] - XX[ll + mm * N] + 2.0;
			break;
		case 11:
			ineqvalue = + XX[ii + jj * N] - XX[ii + kk * N] - XX[ii + ll * N] + XX[ii + mm * N] - XX[jj + kk * N] - XX[jj + ll * N] + XX[jj + mm * N] + XX[kk + ll * N] - XX[kk + mm * N] - XX[ll + mm * N] + 2.0;
			break;
		case 12:
			ineqvalue = + XX[ii + jj * N] - XX[ii + kk * N] - XX[ii + ll * N] - XX[ii + mm * N] - XX[jj + kk * N] - XX[jj + ll * N] - XX[jj + mm * N] + XX[kk + ll * N] + XX[kk + mm * N] + XX[ll + mm * N] + 2.0;
			break;
		case 13:
			ineqvalue = - XX[ii + jj * N] + XX[ii + kk * N] + XX[ii + ll * N] + XX[ii + mm * N] - XX[jj + kk * N] - XX[jj + ll * N] - XX[jj + mm * N] + XX[kk + ll * N] + XX[kk + mm * N] + XX[ll + mm * N] + 2.0;
			break;
		case 14:
			ineqvalue = - XX[ii + jj * N] + XX[ii + kk * N] + XX[ii + ll * N] - XX[ii + mm * N] - XX[jj + kk * N] - XX[jj + ll * N] + XX[jj + mm * N] + XX[kk + ll * N] - XX[kk + mm * N] - XX[ll + mm * N] + 2.0;
			break;
		case 15:
			ineqvalue = - XX[ii + jj * N] + XX[ii + kk * N] - XX[ii + ll * N] + XX[ii + mm * N] - XX[jj + kk * N] + XX[jj + ll * N] - XX[jj + mm * N] - XX[kk + ll * N] + XX[kk + mm * N] - XX[ll + mm * N] + 2.0;
			break;
		case 16:
			ineqvalue = - XX[ii + jj * N] + XX[ii + kk * N] - XX[ii + ll * N] - XX[ii + mm * N] - XX[jj + kk * N] + XX[jj + ll * N] + XX[jj + mm * N] - XX[kk + ll * N] - XX[kk + mm * N] + XX[ll + mm * N] + 2.0;
			break;
		case 17:
			ineqvalue = - XX[ii + jj * N] - XX[ii + kk * N] + XX[ii + ll * N] + XX[ii + mm * N] + XX[jj + kk * N] - XX[jj + ll * N] - XX[jj + mm * N] - XX[kk + ll * N] - XX[kk + mm * N] + XX[ll + mm * N] + 2.0;
			break;
		case 18:
			ineqvalue = - XX[ii + jj * N] - XX[ii + kk * N] + XX[ii + ll * N] - XX[ii + mm * N] + XX[jj + kk * N] - XX[jj + ll * N] + XX[jj + mm * N] - XX[kk + ll * N] + XX[kk + mm * N] - XX[ll + mm * N] + 2.0;
			break;
		case 19:
			ineqvalue = - XX[ii + jj * N] - XX[ii + kk * N] - XX[ii + ll * N] + XX[ii + mm * N] + XX[jj + kk * N] + XX[jj + ll * N] - XX[jj + mm * N] + XX[kk + ll * N] - XX[kk + mm * N] - XX[ll + mm * N] + 2.0;
			break;
		case 20:
			ineqvalue = - XX[ii + jj * N] - XX[ii + kk * N] - XX[ii + ll * N] - XX[ii + mm * N] + XX[jj + kk * N] + XX[jj + ll * N] + XX[jj + mm * N] + XX[kk + ll * N] + XX[kk + mm * N] + XX[ll + mm * N] + 2.0;
			break;
    }

    return ineqvalue;
}

