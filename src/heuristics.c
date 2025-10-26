#include "constants_and_macros.h"

#include <math.h>
#include <mkl.h>
#include <assert.h>
#include <stdio.h>
#include <float.h>

#include "heuristics.h"
#include "b_and_b_framework.h"
#include "user_structs.h"
#include "srflp.h"
#include "miscellaneous_functions.h"


void generate_random_vector_on_unit_sphere(MKL_INT dimension, double *vector, VSLStreamStatePtr stream) {
	double norm = - DBL_MAX;
	MKL_INT inc = 1;
	do {
		vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, stream, dimension, vector, 0.0, 1.0);	
		norm = dnrm2(&dimension, vector, &inc);
	} while ((1e-42) >= norm);
	double inorm = 1.0 / norm;
	dscal(&dimension, &inorm, vector, &inc);
}


void run_heuristics(Subproblem *subproblem) {

    if (0 == HEURISTIC_LEVEL) return; // Loesungen nur an Blaettern


    Subproblem_Data *subproblem_data = subproblem->subproblem_data;
    
    Instance *instance = subproblem->instance;
    
    Instance_Data *instance_data = instance->instance_data;
    
    int number_of_facilities = instance_data->number_of_facilities;

    VSLStreamStatePtr rng_stream = subproblem_data->rng_stream;
    int number_of_variables = instance->number_of_variables;
    int dimension_of_Q = subproblem_data->dimension_of_Q;
    
    int M = subproblem_data->M;
    bool *fixed_variables = subproblem->fixed_variables;
    double *Z = subproblem_data->Z;
    double *X = subproblem_data->X;
    int *index_shift = subproblem_data->index_shift;
    
    
    double *vector_orthogonal_to_random_hyperplane = subproblem_data->vector_orthogonal_to_random_hyperplane;

	double *sca = subproblem_data->sca;
	
	int *random_permutation = subproblem_data->random_permutation;
	int *most_fractional_permutation = subproblem_data->most_fractional_permutation;
	
	for (int i = 0; i < number_of_variables; i++) {
		random_permutation[i] = i;
		most_fractional_permutation[i] = i;
		sca[i] = fabs(X[dimension_of_Q - 1 + index_shift[i] * dimension_of_Q]);
	}
	
	// bestimme "most_fractional_permutation"
	// verwende "sca" zur Hilfe
	for (int i = 0; i < number_of_variables; i++) {
		for (int j = i + 1; j < number_of_variables; j++) {
			if (sca[i] < sca[j]) {
				double sca_temp = sca[i];
				sca[i] = sca[j];
				sca[j] = sca_temp;
				int perm_temp = most_fractional_permutation[i];
				most_fractional_permutation[i] = most_fractional_permutation[j];
				most_fractional_permutation[j] = perm_temp;
			}
		}
	}
	
	long number_of_hyperplanes = 0;
	switch (HEURISTIC_LEVEL) {
		case 1:
			number_of_hyperplanes = 1;
			break;
		case 2:
			number_of_hyperplanes = number_of_facilities;
			break;
		case 3:
			number_of_hyperplanes = dimension_of_Q;
			break;
		case 4:
			number_of_hyperplanes = dimension_of_Q;
			break;
	}
	
    // Goemans/Williamson random hyperplane rounding
    for (long count = 0; count < number_of_hyperplanes; count++) {

        // Compute the "random hyperplane"
        generate_random_vector_on_unit_sphere((MKL_INT) M, vector_orthogonal_to_random_hyperplane, rng_stream);
        
        // berechne "sca", um es fuer verschiedene Heuristiken verwenden zu koennen
        int index = 0;
        for (int i = 0; i < number_of_variables; i++) {
        	if (false == fixed_variables[i]) {
        		sca[i] = 0.0;
        		for (int j = 0; j < M; j++) {
        			sca[i] += vector_orthogonal_to_random_hyperplane[j] * Z[j * dimension_of_Q + index];
        		}
        		index++;
        	}
        }
		// Heuristik, die zulaessige Loesungen liefert, welche aber nicht unbedingt im Teilproblem zulaessig sind
		run_goemans_williamson_heuristic_without_implicit_fixations(subproblem, sca);
       	// Heuristik, welche zulaessige Loesungen fuer das Teilproblem bestimmt
        shuffle_permutation(random_permutation, number_of_variables, rng_stream);      
        run_goemans_williamson_heuristic_with_implicit_fixations(subproblem, sca, random_permutation);
        run_goemans_williamson_heuristic_with_implicit_fixations(subproblem, sca, most_fractional_permutation);
              
    }
}


void run_goemans_williamson_heuristic_without_implicit_fixations(Subproblem *subproblem, double *sca) {
	Subproblem_Data *subproblem_data = subproblem->subproblem_data;
	Instance *instance = subproblem->instance;
	Instance_Data *instance_data = instance->instance_data;
	VSLStreamStatePtr rng_stream = subproblem_data->rng_stream;
	int number_of_facilities = instance_data->number_of_facilities;
	int number_of_variables = instance->number_of_variables;
	bool *fixed_variables = subproblem->fixed_variables;
	bool *fixed_values = subproblem->fixed_values;
	bool *temp_solution = subproblem_data->temp_solution;
	int *number_of_facilities_on_the_right = subproblem_data->number_of_facilities_on_the_right;
	int *permutation = subproblem_data->permutation;
	long *distances_times_two = subproblem_data->distances_times_two;
	
    // beide Moeglichkeiten, um die Menge der Variablen zu partitionieren
    for (int side_of_hyperplane = 0; side_of_hyperplane <= 1; side_of_hyperplane++) {

    	// ohne implizite Fixierungen (liefert zulaessige Loesungen, welche NICHT unbedingt zulaessig im Teilproblem sind)
	    for (int i = 0; i < number_of_variables; i++) {
	        if (true == fixed_variables[i]) {
	            temp_solution[i] = fixed_values[i];
	        } else {
	            if (sca[i] < 0) {
	                temp_solution[i] = side_of_hyperplane;
	            } else {
	                temp_solution[i] = 1 - side_of_hyperplane;
	            }	            
	        }
	    }
		make_boolean_solution_feasible(number_of_facilities, temp_solution, number_of_facilities_on_the_right, permutation, rng_stream);
		double value = calculate_objective_value_from_permutation(instance, permutation, distances_times_two);
		assert(true == is_solution_feasible(instance, temp_solution));
		assert(value == calculate_objective_value(instance, temp_solution));
		try_to_update_local_and_global_primalbound_and_solution(subproblem, temp_solution, value);
	}
}


void run_goemans_williamson_heuristic_with_implicit_fixations(Subproblem *subproblem, double *sca, int *permutation_for_heuristic) {
	Subproblem_Data *subproblem_data = subproblem->subproblem_data;
	Instance *instance = subproblem->instance;
	Instance_Data *instance_data = instance->instance_data;
	int number_of_facilities = instance_data->number_of_facilities;
	int number_of_variables = instance->number_of_variables;
	bool *fixed_variables = subproblem->fixed_variables;
	bool *fixed_values = subproblem->fixed_values;
	bool *temp_solution = subproblem_data->temp_solution;
	bool *temp_solution_fixed_variables = subproblem_data->temp_solution_fixed_variables;
	int *permutation = subproblem_data->permutation;
	int *to_do_list = subproblem_data->to_do_list;
	long *distances_times_two = subproblem_data->distances_times_two;
	int *number_of_facilities_on_the_right = subproblem_data->number_of_facilities_on_the_right;
	
	 // beide Moeglichkeiten, um die Menge der Variablen zu partitionieren
	for (int side_of_hyperplane = 0; side_of_hyperplane <= 1; side_of_hyperplane++) {
		// gehe die gegebene Permutation einmal von links und einmal von rechts durch
		for (int direction = 0; direction <= 1; direction++) {
			
			// "Initialisierung"
			for (int i = 0; i < number_of_variables; i++) {
				if (true == fixed_variables[i]) {
					temp_solution_fixed_variables[i] = true;
					temp_solution[i] = fixed_values[i];
				} else {
					temp_solution_fixed_variables[i] = false;
				}
			}
			
			// Durchfuehrung
			for (int i = 0; i < number_of_variables; i++) {
				int variable = (0 == direction) ? permutation_for_heuristic[i] : permutation_for_heuristic[number_of_variables - 1 - i];
				if (false == temp_solution_fixed_variables[variable]) {
					temp_solution_fixed_variables[variable] = true;
			        if (sca[variable] < 0) {
			            temp_solution[variable] = side_of_hyperplane;
			        } else {
			            temp_solution[variable] = 1 - side_of_hyperplane;
			        }
			        fix_implicitly(number_of_facilities, temp_solution_fixed_variables, temp_solution, to_do_list, variable);
				}
			}
			
			// Auswertung
			assert(true == is_solution_feasible(instance, temp_solution));
			assert(true == is_solution_feasible_in_subproblem(subproblem, temp_solution));
			compute_permutation_from_boolean_solution(number_of_facilities, temp_solution, permutation, number_of_facilities_on_the_right);
			double value = calculate_objective_value_from_permutation(instance, permutation, distances_times_two);
			assert(value == calculate_objective_value(instance, temp_solution));
			try_to_update_local_and_global_primalbound_and_solution(subproblem, temp_solution, value);
			
		}
	        
	}
	
}


void make_boolean_solution_feasible(int number_of_facilities, bool *solution, int *number_of_facilities_on_the_right, int *permutation, VSLStreamStatePtr stream) {
	// "Initialisierung"
	for (int i = 0; i < number_of_facilities; i++) {
		number_of_facilities_on_the_right[i] = 0; 
		permutation[i] = i;
	}
	int number_of_variables = (number_of_facilities * number_of_facilities - number_of_facilities) / 2;
	for (int i = 0; i < number_of_variables; i++) {
		int f1 = ceil(0.5 * (2 * number_of_facilities - 3 - sqrt(4 * number_of_facilities * number_of_facilities - 4 * number_of_facilities - 7 - 8 * i)));
		int f2 = f1 * (f1 - 2 * number_of_facilities + 3) / 2 + i + 1;
		if (true == solution[i]) { // f1 links von f2
			number_of_facilities_on_the_right[f1] += 1;
		} else {
			number_of_facilities_on_the_right[f2] += 1;
		}
	}
	// permutiere die Paare (permutation[i], number_of_facilities_on_the_right[i]) zufaellig
	for (int i = number_of_facilities - 1; i >= 1; i--) {
		int j;
		viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 1, &j, 0, i + 1);
		int temp = permutation[i];
		permutation[i] = permutation[j];
		permutation[j] = temp;
		temp = number_of_facilities_on_the_right[i];
		number_of_facilities_on_the_right[i] = number_of_facilities_on_the_right[j];
		number_of_facilities_on_the_right[j] = temp;
	}
	// sortiere ABSTEIGEND nach number_of_facilities_on_the_right
	for (int i = 0; i < number_of_facilities; i++) {
		for (int j = i + 1; j < number_of_facilities; j++) {
			if (number_of_facilities_on_the_right[i] < number_of_facilities_on_the_right[j]) {
				int temp = permutation[i];
				permutation[i] = permutation[j];
				permutation[j] = temp;
				temp = number_of_facilities_on_the_right[i];
				number_of_facilities_on_the_right[i] = number_of_facilities_on_the_right[j];
				number_of_facilities_on_the_right[j] = temp;
			}
		}
	}
	// mache aus der Permutation wieder eine binaere Loesung und speichere sie in "solution"
	for (int i = 0; i < number_of_facilities; i++) {
		for (int j = i + 1; j < number_of_facilities; j++) {
			int index = inMa(number_of_facilities, permutation[i], permutation[j]);
			if (permutation[i] < permutation[j]) {
				solution[index] = true;
			} else {
				solution[index] = false;
			}
		}
	}
}

void shuffle_permutation(int *permutation, int length, VSLStreamStatePtr stream) {
	for (int i = length - 1; i >= 1; i--) {
		int j;
		viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 1, &j, 0, i + 1);
		int temp = permutation[i];
		permutation[i] = permutation[j];
		permutation[j] = temp;
	}
}
