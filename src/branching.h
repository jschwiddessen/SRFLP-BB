#ifndef BRANCHING_H
#define BRANCHING_H

#include "b_and_b_framework.h"

void add_initial_subproblems_to_bb_tree(Instance *instance);
void create_and_add_subproblems_to_bb_tree_edge(Subproblem *subproblem_father);
int determine_most_fractional_branching_variable(Subproblem *subproblem);
int choose_branching_variable(Subproblem *subproblem);
void create_and_add_subproblems_to_bb_tree(Subproblem *subproblem_father);





#endif

