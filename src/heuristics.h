#ifndef HEURISTICS_H
#define HEURISTICS_H

#include <mkl.h>

#include "b_and_b_framework.h"

void generate_random_vector_on_unit_sphere(MKL_INT dimension, double *vector, VSLStreamStatePtr stream);
void run_heuristics(Subproblem *subproblem);
void make_boolean_solution_feasible(int number_of_facilities, bool *solution, int *number_of_facilities_on_the_right, int *permutation, VSLStreamStatePtr stream);
void run_goemans_williamson_heuristic_with_implicit_fixations(Subproblem *subproblem, double *sca, int *permutation_for_heuristic);
void shuffle_permutation(int *permutation, int length, VSLStreamStatePtr stream);
void run_goemans_williamson_heuristic_without_implicit_fixations(Subproblem *subproblem, double *sca);


#endif

