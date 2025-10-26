#include "constants_and_macros.h"

#include "assert.h"
#include "math.h"
#include "float.h"

#include "branching.h"
#include "b_and_b_framework.h"
#include "miscellaneous_functions.h"


void add_initial_subproblems_to_bb_tree(Instance *instance) {
	
	Subproblem *subproblem = create_subproblem(instance);
	
	if (0 != SYMMETRY_BREAKING) {
		int number_of_variables = instance->number_of_variables;
		int symmetry_breaking_variable = 0;
		bool symmetry_breaking_value = true;
		if (2 == SYMMETRY_BREAKING) {
			symmetry_breaking_variable = rand() % number_of_variables;
			symmetry_breaking_value = rand() % 2;
		}
		subproblem->fixed_variables[symmetry_breaking_variable] = true;
		subproblem->fixed_values[symmetry_breaking_variable] = symmetry_breaking_value;
	}
	
	if (0 == INITIAL_SUBPROBLEMS) {
		update_local_primalbound_and_solution_and_insert_subproblem_into_bb_tree(subproblem);
	} else {
	
	    (instance->number_of_created_subproblems)--;
    	(instance->number_of_left_subproblems)--;
		
		create_and_add_subproblems_to_bb_tree_edge(subproblem);

		clean_up_subproblem(subproblem);
		free(subproblem);
	}
}

void create_and_add_subproblems_to_bb_tree_edge(Subproblem *subproblem_father) {
    Instance *instance = subproblem_father->instance;
    Instance_Data *instance_data = instance->instance_data;
    bool *fixed_variables = subproblem_father->fixed_variables;
    bool *fixed_values = subproblem_father->fixed_values;
    int number_of_facilities = instance_data->number_of_facilities;
    int number_of_facilities_on_the_left[number_of_facilities];
    int number_of_facilities_on_the_right[number_of_facilities];
    for (int i = 0; i < number_of_facilities; i++) {
    	number_of_facilities_on_the_left[i] = 0;
    	number_of_facilities_on_the_right[i] = 0;
    }
    
  	for (int f1 = 0; f1 < number_of_facilities; f1++) {
  		for (int f2 = f1 + 1; f2 < number_of_facilities; f2++) {
  			int pair = inMa(number_of_facilities, f1, f2);
  			if (true == fixed_variables[pair]) {
  				if (true == fixed_values[pair]) {
  					number_of_facilities_on_the_right[f1]++;
  					number_of_facilities_on_the_left[f2]++;
  				} else {
  					number_of_facilities_on_the_right[f2]++;
  					number_of_facilities_on_the_left[f1]++;
  				}
  			}
  		}
  	}
  	if (rand() % 2) {
	  	// wir wollen links am Rand fixieren
	  	// bestimme die Position, an der man fixieren muss (0, 1, ..., n - 1)
		int position = -1;
		for (int candidate_position = 0; candidate_position < number_of_facilities; candidate_position++) {
			bool found = true;
			for (int f = 0; f < number_of_facilities; f++) {
				if (((number_of_facilities - 1 - candidate_position) == number_of_facilities_on_the_right[f]) && (candidate_position == number_of_facilities_on_the_left[f])) {
				 	found = false;
					break;
				}
			}
			if (true == found) {
				position = candidate_position;
				break;
			}
		}
		
		// position gibt nun die "Position" an
		// an dieser Position kommen nun nur Maschinen infrage, die
		// GENAU "position viele" Maschinen links von sich liegen haben
		
		int number_of_variables = (number_of_facilities * number_of_facilities - number_of_facilities) / 2;
		int to_do_list[number_of_variables];
		
		for (int f1 = 0; f1 < number_of_facilities; f1++) {
			if ((position == number_of_facilities_on_the_left[f1])) {
				Subproblem *child = create_child_from_father(subproblem_father);
				for (int f2 = 0; f2 < number_of_facilities; f2++) {
					if (f2 == f1) continue;
					int pair = inMa(number_of_facilities, f1, f2);
					if (false == (child->fixed_variables)[pair]) {
						(child->fixed_variables)[pair] = true;
						if (f1 < f2) {
							(child->fixed_values)[pair] = true;
						} else {
							(child->fixed_values)[pair] = false;
						}
						fix_implicitly(number_of_facilities, child->fixed_variables, child->fixed_values, to_do_list, pair);
						update_local_primalbound_and_solution_and_insert_subproblem_into_bb_tree(child);
					}
				}
			}
		}
	} else {
	  	// wir wollen rechts am Rand fixieren
	  	// bestimme die Position, an der man fixieren muss (n - 1, n - 2, ... 1, 0)
		int position = -1;
		for (int candidate_position = number_of_facilities - 1; candidate_position >= 0; candidate_position--) {
			bool found = true;
			for (int f = 0; f < number_of_facilities; f++) {
				if (((candidate_position) == number_of_facilities_on_the_left[f]) && (number_of_facilities - 1 - candidate_position == number_of_facilities_on_the_right[f])) {
				 	found = false;
					break;
				}
			}
			if (true == found) {
				position = candidate_position;
				break;
			}
		}
		
		// position gibt nun die "Position" an
		// an dieser Position kommen nun nur Maschinen infrage, die
		// GENAU "position viele" Maschinen rechts von sich liegen haben
		
		int number_of_variables = (number_of_facilities * number_of_facilities - number_of_facilities) / 2;
		int to_do_list[number_of_variables];
		
		for (int f1 = 0; f1 < number_of_facilities; f1++) {
			if ((number_of_facilities - 1 - position == number_of_facilities_on_the_right[f1])) {
				Subproblem *child = create_child_from_father(subproblem_father);
				for (int f2 = 0; f2 < number_of_facilities; f2++) {
					if (f2 == f1) continue;
					int pair = inMa(number_of_facilities, f1, f2);
					if (false == (child->fixed_variables)[pair]) {
						(child->fixed_variables)[pair] = true;
						if (f1 < f2) {
							(child->fixed_values)[pair] = false;
						} else {
							(child->fixed_values)[pair] = true;
						}
						fix_implicitly(number_of_facilities, child->fixed_variables, child->fixed_values, to_do_list, pair);
						update_local_primalbound_and_solution_and_insert_subproblem_into_bb_tree(child);
					}
				}
			}
		}
	}
	
}



void create_and_add_subproblems_to_bb_tree(Subproblem *subproblem_father) {
    
    if (2 == BRANCHING) {
    	create_and_add_subproblems_to_bb_tree_edge(subproblem_father);
    	return;
    }
    
    Instance *instance = subproblem_father->instance;
    int branching_variable = -1;
    
    if (1 == BRANCHING) {
    	branching_variable = choose_branching_variable(subproblem_father);
    } else {
    	branching_variable = determine_most_fractional_branching_variable(subproblem_father);
    }
    
    // binaeres Branching
	assert(branching_variable >= 0 && branching_variable < instance->number_of_variables);
    assert(false == (subproblem_father->fixed_variables)[branching_variable]); // die Branching-Variable darf nicht bereits fixiert sein
    for (int fixierung = 0; fixierung <= 1; fixierung++) {

        // erzeuge Kind, das alle noetigen problemUNSPEZIFISCHEN Daten erhaelt
        Subproblem *child = create_child_from_father(subproblem_father);
        // setze das Branching
        (child->fixed_variables)[branching_variable] = true;
        (child->fixed_values)[branching_variable] = fixierung;

		Instance_Data *instance_data = instance->instance_data;

		int number_of_facilities = instance_data->number_of_facilities;
		int number_of_variables = (number_of_facilities * number_of_facilities - number_of_facilities) / 2;
		int to_do_list[number_of_variables];
		fix_implicitly(number_of_facilities, child->fixed_variables, child->fixed_values, to_do_list, branching_variable);

        update_local_primalbound_and_solution_and_insert_subproblem_into_bb_tree(child);
    }
}

int choose_branching_variable(Subproblem *subproblem) {
	Instance *instance = subproblem->instance;
	Instance_Data *instance_data = instance->instance_data;
	Subproblem_Data *subproblem_data = subproblem->subproblem_data;
	int number_of_facilities = instance_data->number_of_facilities;
	int permutation[number_of_facilities];
	int number_of_facilities_on_the_right[number_of_facilities];
	bool *solution = subproblem->best_local_solution;
	double *X = subproblem_data->X;
	bool *fixed_variables = subproblem->fixed_variables;
	int *index_shift = subproblem_data->index_shift;
	int dimension_of_X = subproblem_data->dimension_of_Q;
	// da durch Heuristiken in jedem Teilproblem zulaessige Loesungen gefunden werden
	if (false ==  is_solution_feasible_in_subproblem(subproblem, solution)) {
		return determine_most_fractional_branching_variable(subproblem);
	}
	assert(true ==  is_solution_feasible_in_subproblem(subproblem, solution));
	assert(- DBL_MAX != subproblem->local_primalbound);
	compute_permutation_from_boolean_solution(number_of_facilities, solution, permutation, number_of_facilities_on_the_right);
	
	int candidate = -1;
	double distance = DBL_MAX;
	for (int i = 0; i < number_of_facilities - 1; i++) {
		// gucke Paar (permutation[i], permutation[i + 1]) an
		int pair = inMa(number_of_facilities, permutation[i], permutation[i + 1]);
		if (true == fixed_variables[pair]) continue;
		pair = index_shift[pair];
		if (fabs(X[dimension_of_X - 1 + pair * dimension_of_X]) < distance) {
			candidate = pair;
			distance = fabs(X[dimension_of_X - 1 + pair * dimension_of_X]);
		}
	}
	if (-1 == candidate) {
		return determine_most_fractional_branching_variable(subproblem);
	}
	
	int number_of_variables = instance->number_of_variables;
	for (int i = 0; i < number_of_variables; i++) {
		if (false == fixed_variables[i] && candidate == index_shift[i]) {
			return i;
		}
	}
	assert(false);
	return -1;

}

int determine_most_fractional_branching_variable(Subproblem *subproblem) {
	Subproblem_Data *subproblem_data = subproblem->subproblem_data;
	int dimension_of_X = subproblem_data->dimension_of_Q;
	double *X = subproblem_data->X;
	int candidate = -1;
	double distance = DBL_MAX;
	for (unsigned int i = 0; i < dimension_of_X - 1; i++) {
		if (fabs(X[dimension_of_X - 1 + i * dimension_of_X]) < distance) {
			candidate = i;
			distance = fabs(X[dimension_of_X - 1 + i * dimension_of_X]);
		}
	}
	// finde nun "richtigen" Index
	bool *fixed_variables = subproblem->fixed_variables;
	int *index_shift = subproblem_data->index_shift;

	int number_of_variables = (subproblem->instance)->number_of_variables;
	for (int i = 0; i < number_of_variables; i++) {
		if (false == fixed_variables[i] && candidate == index_shift[i]) {
			return i;
		}
	}
	assert(false);
	return -1;
}
