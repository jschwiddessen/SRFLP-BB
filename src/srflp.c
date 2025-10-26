#include "constants_and_macros.h"

#include <limits.h>
#include <assert.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

#include "user_structs.h"
#include "b_and_b_framework.h"
#include "srflp.h"
#include "miscellaneous_functions.h"
#include "build_subproblem_data.h"
#include "instance_data.h"
#include "bounding_procedure.h"


// teste alle 3-Kreis-Gleichungen
bool is_solution_feasible(Instance *instance, bool *solution) {
	int number_of_facilities = (instance->instance_data)->number_of_facilities;
	for (int i = 0; i < number_of_facilities; i++) {
		for (int j = i + 1; j < number_of_facilities; j++) {
			for (int k = j + 1; k < number_of_facilities; k++) {
			
                int pair_i_j = inMa(number_of_facilities, i, j);
                int pair_i_k = inMa(number_of_facilities, i, k);
                int pair_j_k = inMa(number_of_facilities, j, k);
			
            	int right_hand_side = 0;
            	right_hand_side += (2 * solution[pair_i_j] - 1) * (2 * solution[pair_j_k] - 1);
            	right_hand_side -= (2 * solution[pair_i_j] - 1) * (2 * solution[pair_i_k] - 1);
            	right_hand_side -= (2 * solution[pair_i_k] - 1) * (2 * solution[pair_j_k] - 1);
			
            	if (-1 != right_hand_side) {
            		return false;
            	}
							
			}
		}
	}
	return true;
}



// es wird angenommen, dass die Loesung zulaessig ist
double calculate_objective_value(Instance *instance, bool *sol) {
	long **C = (instance->instance_data)->C;
	int number_of_variables = instance->number_of_variables;
	long return_value = 0;
	for (int i = 0; i < number_of_variables; i++) {
		for (int j = i + 1; j < number_of_variables; j++) {
			if (sol[i] == sol[j]) {
				return_value += C[i][j];
			} else {
				return_value -= C[i][j];
			}
		}
	}
	return (return_value / 2.0);
}

void compute_and_update_dualbound_and_run_heuristics(Subproblem *subproblem) {
	SDPbound(subproblem);
}

void print_optimal_solution(Instance *instance) {
	Instance_Data *instance_data = instance->instance_data;
	int number_of_facilities = instance_data->number_of_facilities;
	int permutation[number_of_facilities];
	int number_of_facilities_on_the_right[number_of_facilities];
	compute_permutation_from_boolean_solution(number_of_facilities, instance->best_global_solution, permutation, number_of_facilities_on_the_right);

	printf("The instance has been solved to global optimality.\n");

	printf("Optimal value: %.1f\n\n", instance_data->K - instance->global_primalbound);
	printf("Optimal ordering:\n\n");
	printf("%d", permutation[0] + 1);
	for (int i = 1; i < number_of_facilities; i++) {
		printf(" %d", permutation[i] + 1);
	}
	printf("\n\n");
	printf("Solution file has been written.\n\n");
}
