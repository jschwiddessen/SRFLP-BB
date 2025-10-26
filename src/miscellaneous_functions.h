#ifndef DIVERSE_FUNCTIONS_H
#define DIVERSE_FUNCTIONS_H

#include <stdbool.h>

#include "b_and_b_framework.h"


double calculate_objective_value_from_permutation(Instance *instance, int *permutation, long *distances_times_two);
void compute_permutation_from_boolean_solution(int number_of_facilities, bool *solution, int *permutation, int *number_of_facilities_on_the_right);
int inMa(int n, int one, int two);
bool check_implicit_fixations(int number_of_facilities, bool *fixed_variables, bool *fixed_values);
void fix_implicitly(int number_of_facilities, bool *fixed_variables, bool *fixed_values, int *to_do_list, int just_fixed);

#endif
