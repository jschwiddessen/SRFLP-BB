#include "constants_and_macros.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "miscellaneous_functions.h"
#include "b_and_b_framework.h"
#include "user_structs.h"



double calculate_objective_value_from_permutation(Instance *instance, int *permutation, long *distances_times_two) {
	Instance_Data *instance_data = instance->instance_data;
	int number_of_facilities = instance_data->number_of_facilities;
	long K_times_two = instance_data->K_times_two; // neu!
	int *lengths_of_facilities = instance_data->lengths_of_facilities;
	int **weights = instance_data->weights;
	// berechne Abstaende von Maschine permutation[0] zu allen anderen
	distances_times_two[0] = 0; // Maschine permutation[0] hat Abstand 0 zu sich selbst
	for (int i = 1; i < number_of_facilities; i++) {
		distances_times_two[i] = distances_times_two[i - 1] + lengths_of_facilities[permutation[i]] + lengths_of_facilities[permutation[i - 1]];
	}
	// fuer i < j gilt nun: 2 * d_permutation[i]_permutation[j] = distances_times_two[j] - distances_times_two[i]
	long value = 0;
	for (int i = 0; i < number_of_facilities; i++) {
		for (int j = i + 1; j < number_of_facilities; j++) {
			value += (distances_times_two[j] - distances_times_two[i]) * weights[permutation[i]][permutation[j]];
		}
	}
	// nun gilt: K - x = value / 2.0 bzw. K_times_two / 2.0 - x = value / 2.0
	// => x = (K_times_two - value) / 2.0
	K_times_two -= value;
	return (K_times_two / 2.0);
}


// Loesung muss zulaessig sein!
void compute_permutation_from_boolean_solution(int number_of_facilities, bool *solution, int *permutation, int *number_of_facilities_on_the_right) {
	int number_of_variables = (number_of_facilities * number_of_facilities - number_of_facilities) / 2;
	for (int i = 0; i < number_of_facilities; i++) {
		number_of_facilities_on_the_right[i] = 0; // "Initialisierung"
	}
	for (int i = 0; i < number_of_variables; i++) {
    	int f1 = ceil(0.5 * (2 * number_of_facilities - 3 - sqrt(4 * number_of_facilities * number_of_facilities - 4 * number_of_facilities - 7 - 8 * i)));
    	int f2 = f1 * (f1 - 2 * number_of_facilities + 3) / 2 + i + 1;
    	if (true == solution[i]) { // f1 links von f2
    		number_of_facilities_on_the_right[f1] += 1;
    	} else {
    		number_of_facilities_on_the_right[f2] += 1;
    	}
	}
	for (int i = 0; i < number_of_facilities; i++) {
		permutation[number_of_facilities - 1 - number_of_facilities_on_the_right[i]] = i;
	}
}

int inMa(int number_of_facilities, int f1, int f2) {
	assert(number_of_facilities > 0);
	assert(f1 >= 0 && f1 < number_of_facilities);
	assert(f2 >= 0 && f2 < number_of_facilities);
	assert(f1 != f2);
	if (f1 < f2) {
		return (f1 * (2 * number_of_facilities - f1 - 1) / 2 + f2 - f1 - 1);
	} else {
		return (f2 * (2 * number_of_facilities - f2 - 1) / 2 + f1 - f2 - 1);
	}
	assert(false);
	return -1;
}

// liefert "true" zurueck, wenn alle impliziten Fixierungen durchgefuehrt wurden bzw. die 3-Kreis-Gleichungen gelten
bool check_implicit_fixations(int number_of_facilities, bool *fixed_variables, bool *fixed_values) {
	assert(number_of_facilities > 0);
	assert(NULL != fixed_variables);
	assert(NULL != fixed_values);
    for (int i = 0; i < number_of_facilities; i++) {
        for (int j = i + 1; j < number_of_facilities; j++) {
            for (int k = j + 1; k < number_of_facilities; k++) {
                int fixed = 0;
                int pair_i_j = inMa(number_of_facilities, i, j);
                int pair_i_k = inMa(number_of_facilities, i, k);
                int pair_j_k = inMa(number_of_facilities, j, k);
                if (true == fixed_variables[pair_i_j]) fixed++;
                if (true == fixed_variables[pair_i_k]) fixed++;
                if (true == fixed_variables[pair_j_k]) fixed++;
                if (2 == fixed) {
                	// pruefe "Transitivitaet"
                	if (true == fixed_variables[pair_i_j] && true == fixed_variables[pair_i_k]) {
                		if (fixed_values[pair_i_j] != fixed_values[pair_i_k]) {
                			return false;
                		}
                	}
                	if (true == fixed_variables[pair_i_j] && true == fixed_variables[pair_j_k]) {
                		if (fixed_values[pair_i_j] == fixed_values[pair_j_k]) {
                			return false;
                		}
                	}
                	if (true == fixed_variables[pair_i_k] && true == fixed_variables[pair_j_k]) {
                		if (fixed_values[pair_i_k] != fixed_values[pair_j_k]) {
                			return false;
                		}
                	}
                }
                if (3 == fixed) {
                	int right_hand_side = 0;
                	right_hand_side += (2 * fixed_values[pair_i_j] - 1) * (2 * fixed_values[pair_j_k] - 1);
                	right_hand_side -= (2 * fixed_values[pair_i_j] - 1) * (2 * fixed_values[pair_i_k] - 1);
                	right_hand_side -= (2 * fixed_values[pair_i_k] - 1) * (2 * fixed_values[pair_j_k] - 1);
                	if (-1 != right_hand_side) {
                		return false;
                	}
                }
            }
        }
    }
    return true;
}


void fix_implicitly(int number_of_facilities, bool *fixed_variables, bool *fixed_values, int *to_do_list, int just_fixed) {
	assert(number_of_facilities > 0);
	assert(just_fixed >= 0 && just_fixed < (number_of_facilities * number_of_facilities - number_of_facilities) / 2);
	assert(NULL != fixed_variables);
	assert(NULL != fixed_values);
	assert(NULL != to_do_list);

    to_do_list[0] = just_fixed; // just_fixed wurde gerade fixiert, d.h. es muss zur to_do_list-Liste hinzugefuegt werden
    int head = 0;
    int tail = 1;
    while (head != tail) {
    	assert(true == fixed_variables[to_do_list[head]]);
    	assert(head < tail);
    	// verwende "magische Formel"
    	int f1 = ceil(0.5 * (2 * number_of_facilities - 3 - sqrt(4 * number_of_facilities * number_of_facilities - 4 * number_of_facilities - 7 - 8 * to_do_list[head])));
    	int f2 = f1 * (f1 - 2 * number_of_facilities + 3) / 2 + to_do_list[head] + 1;
    	assert(f1 < f2);
    	int pair_f1_f2 = to_do_list[head];
    	assert(pair_f1_f2 == inMa(number_of_facilities, f1, f2));
		// gehe alle 3-Kreis-Gleichungen durch
    	for (int f3 = 0; f3 < f1; f3++) {
    		// f3 < f1 < f2
    		int pair_f3_f1 = inMa(number_of_facilities, f3, f1);
    		int pair_f3_f2 = inMa(number_of_facilities, f3, f2);
    	    if (true == fixed_variables[pair_f3_f1] && fixed_values[pair_f3_f1] == fixed_values[pair_f1_f2] && false == fixed_variables[pair_f3_f2]) {
    	    	fixed_variables[pair_f3_f2] = true;
    	    	fixed_values[pair_f3_f2] = fixed_values[pair_f1_f2];
    	    	to_do_list[tail] = pair_f3_f2;
    	    	tail += 1;
    	    } else if (true == fixed_variables[pair_f3_f2] && fixed_values[pair_f3_f2] != fixed_values[pair_f1_f2] && false == fixed_variables[pair_f3_f1]) {
    	    	fixed_variables[pair_f3_f1] = true;
    	    	fixed_values[pair_f3_f1] = fixed_values[pair_f3_f2];
    	    	to_do_list[tail] = pair_f3_f1;
    	    	tail += 1;
    	    }
    	}
    	for (int f3 = f1 + 1; f3 < f2; f3++) {
    		// f1 < f3 < f2
    		int pair_f1_f3 = inMa(number_of_facilities, f1, f3);
    		int pair_f3_f2 = inMa(number_of_facilities, f3, f2);
    	    if (true == fixed_variables[pair_f3_f2] && fixed_values[pair_f3_f2] != fixed_values[pair_f1_f2] && false == fixed_variables[pair_f1_f3]) {
    	    	fixed_variables[pair_f1_f3] = true;
    	    	fixed_values[pair_f1_f3] = fixed_values[pair_f1_f2];
    	    	to_do_list[tail] = pair_f1_f3;
    	    	tail += 1;
    	    } else if (true == fixed_variables[pair_f1_f3] && fixed_values[pair_f1_f3] != fixed_values[pair_f1_f2] && false == fixed_variables[pair_f3_f2]) {
    	    	fixed_variables[pair_f3_f2] = true;
    	    	fixed_values[pair_f3_f2] = fixed_values[pair_f1_f2];
    	    	to_do_list[tail] = pair_f3_f2;
    	    	tail += 1;
    	    }
    	}
    	for (int f3 = f2 + 1; f3 < number_of_facilities; f3++) {
    		// f1 < f2 < f3
    		int pair_f1_f3 = inMa(number_of_facilities, f1, f3);
    		int pair_f2_f3 = inMa(number_of_facilities, f2, f3);
    	    if (true == fixed_variables[pair_f2_f3] && fixed_values[pair_f2_f3] == fixed_values[pair_f1_f2] && false == fixed_variables[pair_f1_f3]) {
    	    	fixed_variables[pair_f1_f3] = true;
    	    	fixed_values[pair_f1_f3] = fixed_values[pair_f1_f2];
    	    	to_do_list[tail] = pair_f1_f3;
    	    	tail += 1;
    	    } else if (true == fixed_variables[pair_f1_f3] && fixed_values[pair_f1_f3] != fixed_values[pair_f1_f2] && false == fixed_variables[pair_f2_f3]) {
    	    	fixed_variables[pair_f2_f3] = true;
    	    	fixed_values[pair_f2_f3] = fixed_values[pair_f1_f3];
    	    	to_do_list[tail] = pair_f2_f3;
    	    	tail += 1;
    	    }
    	}
    	head += 1;
    }

    assert(true == check_implicit_fixations(number_of_facilities, fixed_variables, fixed_values));
}
