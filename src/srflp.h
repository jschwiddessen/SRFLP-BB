#ifndef SRFLP_H
#define SRFLP_H

#include "user_structs.h"
#include "b_and_b_framework.h"


// Funktionen, die der "Benutzer" selbst implementieren muss
bool is_solution_feasible(Instance *instance, bool *sol);
double calculate_objective_value(Instance *instance, bool *sol);
void compute_and_update_dualbound_and_run_heuristics(Subproblem *subproblem);
void print_optimal_solution(Instance *instance);

#endif
