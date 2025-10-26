#ifndef B_AND_B_FRAMEWORK_H
#define B_AND_B_FRAMEWORK_H

#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <sys/time.h>

#include "tree.h"
#include "user_structs.h"

typedef struct {
   long number_of_created_subproblems; // Anzahl der bisher erzeugten Teilprobleme
   long number_of_evaluated_subproblems; // Anzahl der bisher untersuchten Teilprobleme ( die Teilprobleme, fuer die ein pthread zur Bound-Berechnung gestartet wurde)
   long number_of_left_subproblems; // Anzahl der Teilprobleme, fuer die noch kein pthread gestartet wurde bzw. die gerade untersucht werden (Anzahl im BB-Baum + "currently_evaluated_subproblems")
   double global_dualbound; // globale dualbound (beste bekannte Schranke)
   double global_primalbound; // globale primalbound (Wert der besten bekannten zulaessigen Loesung)
   bool *best_global_solution; // Array fuer die beste bekannte zulaessige Loesung
   int number_of_variables; // Anzahl der binaeren Variablen
   struct timeval start_t; // Startzeit der Optimierung
   clock_t cpu_t;
   Instance_Data *instance_data; // konkrete Daten zur Instanz
   const char *filename;
} Instance;

typedef struct {
   Instance *instance; // jedes Teilproblem "kennt" die Instanz
   long id; // einzigartige ID des Teilproblems; Nummerierung start bei 0
   double local_dualbound; // dualbound fuer dieses Teilproblem; wird beim Erzeugen vom Vater "geerbt"; wird zum Sortieren der BB_Tree_Node's verwendet
   double local_primalbound; // Wert der besten bekannten zulaessigen Loesung in diesem Teilproblem
   bool *best_local_solution; // Array fuer die entsprechende beste bekannte zulaessige Loesung in diesem Teilproblem
   bool *fixed_variables; // Array, welches die fixierten Variablen markiert
   bool *fixed_values; // Array fuer den entsprechenden Wert einer fixierten Variable
   Subproblem_Data *subproblem_data; // konkrete Daten und benoetigte Variablen fuer das Teilproblem
   bool priority;
} Subproblem;

typedef struct BB_Tree_Node { // ("struct BB_Tree_Node" ist wichtig)
   RB_ENTRY(BB_Tree_Node) entry;
   Subproblem *subproblem;
} BB_Tree_Node;

// fertige Funktionen, die der "Benutzer" NICHT verwenden sollte:
int subproblem_compare(BB_Tree_Node *node1, BB_Tree_Node *node2);
void *evaluate_subproblem(void *arg);
void set_up_instance(Instance *instance);
void prune_bb_tree(Instance *instance, double new_global_primal_bound);
void update_global_dualbound(Instance *instance);
void print_terminal_output(Instance *instance);
void *print_terminal_output_in_loop(void *arg);
bool is_solution_feasible_in_subproblem(Subproblem *subproblem, bool *solution);
void update_local_primalbound_and_local_solution(Subproblem *subproblem, bool *new_sol, double new_val);
void update_global_primalbound_and_global_solution_and_prune_bb_tree(Instance *instance, bool *new_sol, double new_val);
void clean_up_subproblem(Subproblem *subproblem);
void clean_up_instance(Instance *instance);
void write_best_global_solution_to_file(Instance *instance);
bool check_instance_consistency(Instance *instance);
bool check_subproblem_consistency(Subproblem *subproblem);

// fertige Funktionen, die der "Benutzer" in seinen Funktionen aufrufen sollte:

// um zu ueberpruefen, ob die Bound-Berechnung abgebrochen werden kann
bool can_subproblem_be_pruned(Subproblem *subproblem);
// erzeugt aus dem Vater ein Kind und initialisiert alle "generischen" Variablen
Subproblem *create_child_from_father(Subproblem *subproblem_father);
// erzeugt ein Subproblem und initialisiert alle "generischen" Variablen
Subproblem *create_subproblem(Instance *instance);
// muss aufgerufen werden, wenn eine neue beste zulaessige Loesung gefunden wurde
// soll aufgerufen werden, wenn eine neue dualbound fuer das Teilproblem bekannt ist (ist nicht schlimm, wenn sie schlechter ist als die alte)
void try_to_update_local_dualbound(Subproblem *subproblem, double new_dualbound);
void update_local_primalbound_and_solution_and_insert_subproblem_into_bb_tree(Subproblem *subproblem);
void try_to_update_local_and_global_primalbound_and_solution(Subproblem *subproblem, bool *solution, double value);
void update_initial_dualbound(Instance *instance, double upper_bound);
void update_initial_primalbound_and_solution(Instance *instance, bool *sol, double value);
void set_number_of_variables(Instance *instance, int number_of_variables);

/* redirect_output.f */
void redirect_output_(); // Redirect the Fortran output to /dev/null

#endif

