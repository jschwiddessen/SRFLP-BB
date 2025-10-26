#include "constants_and_macros.h"

#include <limits.h>
#include <stdbool.h>
#include <stdio.h>

#include "b_and_b_framework.h"
#include "initial_heuristics.h"
#include "srflp.h"


void run_initial_primal_heuristics(Instance *instance) {

	// verwende Identitaet: setze alles auf "true"
	int number_of_variables = instance->number_of_variables;
	bool solution[number_of_variables];
	
	for (int i = 0; i < number_of_variables; i++) {
		solution[i] = true;
	}
	
	double value = calculate_objective_value(instance, solution);

	update_initial_primalbound_and_solution(instance, solution, value);
}


void run_initial_dual_heuristics(Instance *instance) {

	// verwende betragsmaessig groessten Eintrag von C
	Instance_Data *instance_data = instance->instance_data;

	long **C = instance_data->C;
	int number_of_variables = instance->number_of_variables;
	long C_abs_max = LONG_MIN;
		
	for (int i = 0; i < number_of_variables; i++) {
		for (int j = i + 1; j < number_of_variables; j++) {
			if (+ C[i][j] > C_abs_max) {
				C_abs_max = C[i][j];
			}
			if (- C[i][j] > C_abs_max) {
				C_abs_max = - C[i][j];
			}
		}
	}
	
	int number_of_facilities = instance_data->number_of_facilities;
	
	// jede Zeile von C hat maximal 2 * n - 1 Eintraege, die ungleich 0 sind
	long upper_bound = C_abs_max * number_of_variables * (2 * number_of_facilities - 1);
	
	update_initial_dualbound(instance, upper_bound / 4.0);
}
