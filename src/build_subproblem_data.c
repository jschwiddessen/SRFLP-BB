#include "constants_and_macros.h"

#include <stdbool.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <mkl.h>

#include "user_structs.h"
#include "build_subproblem_data.h"
#include "miscellaneous_functions.h"
#include "b_and_b_framework.h"



int count_fixed_variables(int number_of_variables, bool *fixed_variables) {
	assert(number_of_variables > 0);
	assert(NULL != fixed_variables);
	int number_of_fixed_variables = 0;
	for (int i = 0; i < number_of_variables; i++) {
		if (true == fixed_variables[i]) {
			number_of_fixed_variables++;
		}
	}
	return number_of_fixed_variables;
}


// es wird angenommen, dass alle impliziten Fixierungen vorgenommen wurden;
// liefert die Anzahl der 3-Kreis-Gleichungen zurueck, welche NICHT trivialerweise erfuellt sind,
// also zum Problem hinzugefuegt werden muessen;
// Hinweis: sind zwei der drei Variablen in einer 3-Kreis-Gleichung fixiert, so ist sie
// trivialerweise erfuellt!
int count_three_cycle_equations(int number_of_facilities, bool *fixed_variables) {
	assert(number_of_facilities > 0);
	assert(NULL != fixed_variables);
    int count = 0;
    for (int i = 0; i < number_of_facilities; i++) {
        for (int j = i + 1; j < number_of_facilities; j++) {
            for (int k = j + 1; k < number_of_facilities; k++) {
                int fixed = 0;
                if (true == fixed_variables[inMa(number_of_facilities, i, j)]) fixed++;
                if (true == fixed_variables[inMa(number_of_facilities, i, k)]) fixed++;
                if (true == fixed_variables[inMa(number_of_facilities, j, k)]) fixed++;
                if (fixed <= 1) count++;
            }
        }
    }
    return count;
}


int *generate_index_shift(int number_of_variables, bool *fixed_variables) {
    int *index_shift;
    calloc_makro(index_shift, number_of_variables, int);
    index_shift[0] = 0;
    for (int i = 1; i < number_of_variables; i++) {
        if (fixed_variables[i - 1]) index_shift[i] = index_shift[i - 1];
        else index_shift[i] = index_shift[i - 1] + 1;
    }

    for (int i = 1; i < number_of_variables; i++) {
        assert(index_shift[i] == (index_shift[i - 1] + 1 - fixed_variables[i - 1])); // TODO Assertion entfernen
    }

    return index_shift;
}

// es wird angenommen, dass alle IMPLIZITEN FIXIERUNGEN bereits vorgenommen worden sind!
// Achtung: C ist ganzzahlig und entspricht demnach in etwa 4 * Q_0!
Subproblem_Data *construct_problem(int number_of_facilities, long **C, bool *fixed_variables, bool *fixed_values, int id) {

    Subproblem_Data *subproblem_data;
    calloc_makro(subproblem_data, 1, Subproblem_Data);

	int number_of_variables = (number_of_facilities * number_of_facilities - number_of_facilities) / 2;
	int number_of_fixed_variables = count_fixed_variables(number_of_variables, fixed_variables);
	assert(number_of_fixed_variables <= number_of_variables);

	int dimension_of_Q = number_of_variables - number_of_fixed_variables + 1;
	// setze dimension_of_Q
	subproblem_data->dimension_of_Q = dimension_of_Q;

    // stelle sicher, dass alle impliziten Fixierungen vorgenommen wurden
    assert(true == check_implicit_fixations(number_of_facilities, fixed_variables, fixed_values));

    int number_of_three_cycle_equations = count_three_cycle_equations(number_of_facilities, fixed_variables);

    int *index_shift = generate_index_shift(number_of_variables, fixed_variables);
    // setze index_shift
    subproblem_data->index_shift = index_shift;

    int number_of_equations = dimension_of_Q + number_of_three_cycle_equations;
    // setze number_of_equations
    subproblem_data->number_of_equations = number_of_equations;

    Sparse_Matrix *Bs;
    calloc_makro(Bs, number_of_equations, Sparse_Matrix);
    // setze Bs
	subproblem_data->Bs = Bs;
    double *b;
    calloc_makro(b, number_of_equations, double);
    // setze b
	subproblem_data->b = b;

    // diag(X) = e
    for (int eq = 0; eq < dimension_of_Q; eq++) {
        b[eq] = 1.0;
        Bs[eq].nnz = 1;
        calloc_makro(Bs[eq].i, 1, int);
        calloc_makro(Bs[eq].j, 1, int);
        calloc_makro(Bs[eq].val, 1, double);
        Bs[eq].i[0] = eq;
        Bs[eq].j[0] = eq;
        Bs[eq].val[0] = 1.0;
    }

    // 3-Kreis-Gleichungen werden in der Form +/-1/2, +/-1/2, +/-1/2 = -1.0 abgespeichert
    int eq = dimension_of_Q;
	// TODO andere "Repraesentierung" kann gewaehlt werden; was ist mit "Normierung" ?
    for (int i = 0; i < number_of_facilities; i++) {
        for (int j = i + 1; j < number_of_facilities; j++) {
            for (int k = j + 1; k < number_of_facilities; k++) {
                int fixed = 0;
                if (true == fixed_variables[inMa(number_of_facilities, i, j)]) fixed++;
                if (true == fixed_variables[inMa(number_of_facilities, i, k)]) fixed++;
                if (true == fixed_variables[inMa(number_of_facilities, j, k)]) fixed++;
                if (fixed <= 1) {
                    b[eq] = - 1.0;
                    Bs[eq].nnz = 3;
                    calloc_makro(Bs[eq].i, 3, int);
        			calloc_makro(Bs[eq].j, 3, int);
        			calloc_makro(Bs[eq].val, 3, double);
                    // X_ij_jk - X_ij_ik - X_ik_jk = -1
                    int pair_ij = inMa(number_of_facilities, i, j);
                    int pair_ik = inMa(number_of_facilities, i, k);
                    int pair_jk = inMa(number_of_facilities, j, k);
                    if (true == fixed_variables[pair_ij]) {
                        assert(false == fixed_variables[pair_jk]);
                        assert(false == fixed_variables[pair_ik]);
                        Bs[eq].i[2] = dimension_of_Q - 1;
                        Bs[eq].j[2] = index_shift[pair_jk];
                        Bs[eq].val[2] = ((true == fixed_values[pair_ij]) ? (+ 0.5) : (- 0.5));
                        Bs[eq].i[1] = dimension_of_Q - 1;
                        Bs[eq].j[1] = index_shift[pair_ik];
                        Bs[eq].val[1] = ((true == fixed_values[pair_ij]) ? (- 0.5) : (+ 0.5));
                        Bs[eq].i[0] = index_shift[pair_jk];
                        Bs[eq].j[0] = index_shift[pair_ik];
                        Bs[eq].val[0] = - 0.5;
                    } else if (true == fixed_variables[pair_ik]) {
                        assert(false == fixed_variables[pair_ij]);
                        assert(false == fixed_variables[pair_jk]);
                        Bs[eq].i[0] = index_shift[pair_jk];
                        Bs[eq].j[0] = index_shift[pair_ij];
                        Bs[eq].val[0] = + 0.5;
                        Bs[eq].i[1] = dimension_of_Q - 1;
                        Bs[eq].j[1] = index_shift[pair_ij];
                        Bs[eq].val[1] = ((true == fixed_values[pair_ik]) ? (- 0.5) : (+ 0.5));
                        Bs[eq].i[2] = dimension_of_Q - 1;
                        Bs[eq].j[2] = index_shift[pair_jk];
                        Bs[eq].val[2] = ((true == fixed_values[pair_ik]) ? (- 0.5) : (+ 0.5));
                    } else if (true == fixed_variables[pair_jk]) {
                        assert(false == fixed_variables[pair_ij]);
                        assert(false == fixed_variables[pair_ik]);
                        Bs[eq].i[1] = dimension_of_Q - 1;
                        Bs[eq].j[1] = index_shift[pair_ij];
                        Bs[eq].val[1] = ((true == fixed_values[pair_jk]) ? (+ 0.5) : (- 0.5));
                        Bs[eq].i[0] = index_shift[pair_ik];
                        Bs[eq].j[0] = index_shift[pair_ij];
                        Bs[eq].val[0] = - 0.5;
                        Bs[eq].i[2] = dimension_of_Q - 1;
                        Bs[eq].j[2] = index_shift[pair_ik];
                        Bs[eq].val[2] = ((true == fixed_values[pair_jk]) ? (- 0.5) : (+ 0.5));
                    } else {
                        // nichts fixiert
                        Bs[eq].i[1] = index_shift[pair_jk];
                        Bs[eq].j[1] = index_shift[pair_ij];
                        Bs[eq].val[1] = + 0.5;
                        Bs[eq].i[0] = index_shift[pair_ik];
                        Bs[eq].j[0] = index_shift[pair_ij];
                        Bs[eq].val[0] = - 0.5;
                        Bs[eq].i[2] = index_shift[pair_jk];
                        Bs[eq].j[2] = index_shift[pair_ik];
                        Bs[eq].val[2] = - 0.5;
                    }
                    eq++;
                }
            }
        }
    }

	long fixed_objective_value = 0;
	for (int i = 0; i < number_of_variables; i++) {
		for (int j = i + 1; j < number_of_variables; j++) {
			if (fixed_variables[i] && fixed_variables[j]) {
				if (fixed_values[i] == fixed_values[j]) {
					fixed_objective_value += C[i][j];
				} else {
					fixed_objective_value -= C[i][j];
				}
			}
		}
	}
	// teile nur durch 2 statt 4, da nur eine Dreiecksmatrix verwendet wurde
	subproblem_data->fixed_objective_value = fixed_objective_value / 2.0; // TODO Wert richtig?

	// Q_long ist "Hilfsmatrix"; Speicher wird fuer gesamte Matrix allokiert, aber nur als untere Dreiecksmatrix verwendet
    long *Q_long;
    calloc_makro(Q_long, dimension_of_Q * dimension_of_Q, long);

    for (int i = 0; i < number_of_variables; i++) {
        for (int j = 0; j < i; j++) {
            if (!fixed_variables[i] && !fixed_variables[j]) {
                Q_long[index_shift[i] + dimension_of_Q * index_shift[j]] = C[i][j];
            } else if (fixed_variables[i] && !fixed_variables[j]) {
                Q_long[dimension_of_Q - 1 + dimension_of_Q * index_shift[j]] += ((true == fixed_values[i]) ? (+ C[i][j]) : (- C[i][j]));
            } else if (!fixed_variables[i] && fixed_variables[j]) {
                Q_long[dimension_of_Q - 1 + dimension_of_Q * index_shift[i]] += ((true == fixed_values[j]) ? (+ C[i][j]) : (- C[i][j]));
            }
        }
    }


    int count_nonzeros = 0;
    double *Q;
    calloc_makro(Q, dimension_of_Q * dimension_of_Q, double);
	// setze Q
	subproblem_data->Q = Q;

    for (int i = 0; i < dimension_of_Q; i++) {
        for (int j = 0; j <= i; j++) {
            if (0 == Q_long[i + dimension_of_Q * j]) {
                Q[i + dimension_of_Q * j] = 0.0;
                Q[j + dimension_of_Q * i] = 0.0;
            } else {
                count_nonzeros++;
                Q[i + dimension_of_Q * j] = Q_long[i + dimension_of_Q * j] / 4.0; // TODO: Korrekter Wert?
                Q[j + dimension_of_Q * i] = Q_long[i + dimension_of_Q * j] / 4.0;
            }
        }
    }

	// Qs
    (subproblem_data->Qs).nnz = count_nonzeros;
    calloc_makro((subproblem_data->Qs).i, count_nonzeros, int);
    calloc_makro((subproblem_data->Qs).j, count_nonzeros, int);
    calloc_makro((subproblem_data->Qs).val, count_nonzeros, double);

    int counter = 0;
    for (int i = 0; i < dimension_of_Q; i++) {
        for (int j = 0; j < i; j++) {
            if (0 != Q_long[i + dimension_of_Q * j]) {
                (subproblem_data->Qs).i[counter] = i;
                (subproblem_data->Qs).j[counter] = j;
                (subproblem_data->Qs).val[counter] = Q[i + dimension_of_Q * j];
                counter++;
            }
        }
    }

    free(Q_long);


    // initialisiere alle anderen Sachen

    size_t maximum_problem_size = number_of_equations + MAXIMUM_NUMBER_OF_CUTS;

    size_t wa_length = (2 * MAXMIUM_NUMBER_OF_MEMORY_CORRECTIONS + 5) * maximum_problem_size + 11 * MAXMIUM_NUMBER_OF_MEMORY_CORRECTIONS * MAXMIUM_NUMBER_OF_MEMORY_CORRECTIONS + 8 * MAXMIUM_NUMBER_OF_MEMORY_CORRECTIONS;
	size_t iwa_length = 3 * maximum_problem_size;

	// number_of_cuts
	subproblem_data->number_of_cuts = 0;
	// evaluations
	subproblem_data->evaluations = 0;
	// f
	subproblem_data->f = 0.0;
	// M
	subproblem_data->M = 0;
	// rng_stream
	vslNewStream(&(subproblem_data->rng_stream), VSL_BRNG_MT19937, time(NULL) + id);
	
	// current_cuts
	calloc_makro(subproblem_data->current_cuts, MAXIMUM_NUMBER_OF_CUTS, Cut);
	// new_cuts
	calloc_makro(subproblem_data->new_cuts, MAXIMUM_NUMBER_OF_CUTS_PER_ITERATION, Cut);
	
	// temp_solution
	calloc_makro(subproblem_data->temp_solution, number_of_variables, bool);
	// temp_solution_fixed_variables
	calloc_makro(subproblem_data->temp_solution_fixed_variables, number_of_variables, bool);
	// permutation
	calloc_makro(subproblem_data->permutation, number_of_facilities, int);
	// to_do_list
	calloc_makro(subproblem_data->to_do_list, number_of_variables, int);
	// distances_times_two
	calloc_makro(subproblem_data->distances_times_two, number_of_facilities, long);
	// random_permutation
	calloc_makro(subproblem_data->random_permutation, number_of_variables, int);
	// sca
	calloc_makro(subproblem_data->sca, number_of_variables, double);
	// vector_orthogonal_to_random_hyperplane
	calloc_makro(subproblem_data->vector_orthogonal_to_random_hyperplane, dimension_of_Q, double);
	// number_of_facilities_on_the_right
	calloc_makro(subproblem_data->number_of_facilities_on_the_right, number_of_facilities, int);
	// most_fractional_permutation
	calloc_makro(subproblem_data->most_fractional_permutation, number_of_variables, int);
	
	// X
	calloc_makro(subproblem_data->X, dimension_of_Q * dimension_of_Q, double);
	// g
	calloc_makro(subproblem_data->g, maximum_problem_size, double);
	// y
	calloc_makro(subproblem_data->y, maximum_problem_size, double);
	// binf
	calloc_makro(subproblem_data->binf, maximum_problem_size, double);
	// bsup
	calloc_makro(subproblem_data->bsup, maximum_problem_size, double);
	// nbd
	calloc_makro(subproblem_data->nbd, maximum_problem_size, int);
	// wa
	calloc_makro(subproblem_data->wa, wa_length, double);
	// iwa
	calloc_makro(subproblem_data->iwa, iwa_length, int);
	// sizeWORK
	subproblem_data->sizeWORK = 26 * dimension_of_Q;
	// sizeIWORK
	subproblem_data->sizeIWORK = 10 * dimension_of_Q;
	// WORK
	calloc_makro(subproblem_data->WORK, subproblem_data->sizeWORK, double);
	// IWORK
	calloc_makro(subproblem_data->IWORK, subproblem_data->sizeIWORK, int);
	// W
	calloc_makro(subproblem_data->W, dimension_of_Q, double);
	// Z
	calloc_makro(subproblem_data->Z, dimension_of_Q * dimension_of_Q, double);
	// ISUPPZ
	calloc_makro(subproblem_data->ISUPPZ, 2 * dimension_of_Q, int);

	return subproblem_data;

}


void clean_up_subproblem_data(Subproblem_Data *subproblem_data) {
	//if (NULL == subproblem_data) return;
	free(subproblem_data->Q);
	free((subproblem_data->Qs).i);
	free((subproblem_data->Qs).j);
	free((subproblem_data->Qs).val);
	for (int i = 0; i < subproblem_data->number_of_equations; i++) {
		free(((subproblem_data->Bs)[i]).i);
		free(((subproblem_data->Bs)[i]).j);
		free(((subproblem_data->Bs)[i]).val);
	}
	free(subproblem_data->Bs);
	free(subproblem_data->b);
	free(subproblem_data->X);
	free(subproblem_data->index_shift);
	free(subproblem_data->g);
	free(subproblem_data->y);
	free(subproblem_data->binf);
	free(subproblem_data->bsup);
	free(subproblem_data->nbd);
	free(subproblem_data->wa);
	free(subproblem_data->iwa);
	free(subproblem_data->W);
	free(subproblem_data->Z);
	free(subproblem_data->ISUPPZ);
	free(subproblem_data->WORK);
	free(subproblem_data->IWORK);
	vslDeleteStream(&(subproblem_data->rng_stream));
	free(subproblem_data->temp_solution);
	free(subproblem_data->temp_solution_fixed_variables);
	free(subproblem_data->permutation);
	free(subproblem_data->to_do_list);
	free(subproblem_data->distances_times_two);
	free(subproblem_data->random_permutation);
	free(subproblem_data->sca);
	free(subproblem_data->vector_orthogonal_to_random_hyperplane);
	free(subproblem_data->number_of_facilities_on_the_right);
	free(subproblem_data->most_fractional_permutation);
	free(subproblem_data->current_cuts);
	free(subproblem_data->new_cuts);
}
