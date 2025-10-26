#ifndef BUILD_SUBPROBLEM_DATA_H
#define BUILD_SUBPROBLEM_DATA_H


int count_fixed_variables(int number_of_variables, bool *fixed_variables);
int count_three_cycle_equations(int number_of_facilities, bool *fixed_variables);
int *generate_index_shift(int number_of_variables, bool *fixed_variables);
Subproblem_Data *construct_problem(int number_of_facilities, long **C, bool *fixed_variables, bool *fixed_values, int id);
void clean_up_subproblem_data(Subproblem_Data *subproblem_data);

#endif
