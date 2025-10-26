#include "constants_and_macros.h"

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <limits.h>
#include <sys/time.h>

#include "tree.h"
#include "b_and_b_framework.h"
#include "srflp.h"
#include "instance_data.h"
#include "build_subproblem_data.h"
#include "miscellaneous_functions.h"
#include "initial_heuristics.h"
#include "parameters.h"
#include "branching.h"

//#define NDEBUG // Assertions werden nur geprueft, wenn auskommentiert


/*	lege maximale Anzahl an gleichzeitig laufenden Threads fest;
	entspricht der maximalen Anzahl an Teilproblemen, die gleichzeitig bearbeitet werden koennen;
	sollte nicht groesser sein als die tatsaechlich zur Verfuegung stehenden Threads
*/


//*********************************************************************************************************************************
//*********************************************************************************************************************************

/*	mutex, der verwendet wird, damit nur ein Thread (main-Thread und pthread) gleichzeitig
    den BB-Baum manipulieren und die Variable "currently_evaluated_subproblems" veraendern
    kann
*/
static pthread_mutex_t main_mutex = PTHREAD_MUTEX_INITIALIZER;

/*	Condition Variable fuer: der main-Thread legt sich schlafen, wenn er keinen weiteren pthread mehr
	starten kann/sollte; er wartet dann, bis ein pthread sich beendet und ihm ein
	Signal zum Aufwachen sendet
*/
static pthread_cond_t main_cond = PTHREAD_COND_INITIALIZER;

/*	mutex, um dafuer zu sorgen, dass nur ein pthread gleichzeitig die globale
	primalbound und die entsprechende Loesung (Belegung der Variablen) veraendern
	kann
*/
static pthread_mutex_t solution_mutex = PTHREAD_MUTEX_INITIALIZER;


// Array von Pointern auf diejenigen Teilprobleme, die gerade von Threads untersucht werden
static Subproblem *currently_evaluated_subproblems[MAXIMUM_NUMBER_OF_CONCURRENTLY_RUNNING_THREADS];

/*	globale Variable, die die Anzahl der aktuell laufenden Threads nachhaelt;
	jeder Thread hat Zugriff hierdrauf und MUSS diese Variable selbst veraendern;
	main-Thread inkrementiert ausschliesslich und jeder pthread dekrementiert
	bevor er sich beendet
*/
static int number_of_currently_running_threads = 0;

//*********************************************************************************************************************************
//*********************************************************************************************************************************

RB_HEAD(bb_tree, BB_Tree_Node) head_bb_tree = RB_INITIALIZER(&head_bb_tree);
RB_PROTOTYPE_STATIC(bb_tree, BB_Tree_Node, entry, subproblem_compare) // Achtung: "static"!
RB_GENERATE_STATIC(bb_tree, BB_Tree_Node, entry, subproblem_compare) // Achtung: "static"!

// ****************************************************************************************
// ****************************************************************************************

int main(int argc, char *argv[]) {

	redirect_output_();

	printf("\n");

//while (true) {
	if (argc != 2) {
        printf("\n\tusage:\t%s <instance_file>\n\n", argv[0]);
        return EXIT_FAILURE;
    }

    // set seed for random number generator
    srand((unsigned) time(NULL));
    // spaeter: erzeuge fuer jedes Teilproblem einen eigenen Zufallszahlengenerator

    Instance instance;

	instance.filename = argv[1];

    instance.number_of_variables = - INT_MAX; // um spaeter zu testen, ob der Benutzer die Anzahl der Variablen gesetzt hat
    read_and_initialize_given_instance_data(&instance, argv[1]);
	assert(instance.number_of_variables > 0);

    // setzt voraus, dass "instance.number_of_variables" den richtigen Wert hat
    set_up_instance(&instance);

	print_parameters();

	printf("\nRunning heuristics ...\n");
    run_initial_primal_heuristics(&instance);
	run_initial_dual_heuristics(&instance);

	printf("Starting with initial primalbound of %.1f\n", (instance.instance_data)->K - (instance.global_primalbound));
	printf("Starting with initial dualbound of %.2f\n\n", (instance.instance_data)->K - (instance.global_dualbound));

	printf("K=%.1f\n", instance.instance_data->K);

    add_initial_subproblems_to_bb_tree(&instance);

    for (int i = 0; i < MAXIMUM_NUMBER_OF_CONCURRENTLY_RUNNING_THREADS; i++) {
        currently_evaluated_subproblems[i] = NULL;
    }

    // da keine besonderen pthread-Attribute gesetzt werden, verwende immer denselben pthread
    pthread_t thread;

	// fuer staendigen Terminal-Output
	pthread_t terminal_output_thread;

    pthread_mutex_lock(&main_mutex);

	pthread_create(&terminal_output_thread, NULL, print_terminal_output_in_loop, &instance);

    while (true) {
        if (true == RB_EMPTY(&head_bb_tree) && 0 == number_of_currently_running_threads) {
            // fertig, optimale Loesung gefunden!
            break;
        }
        while (!RB_EMPTY(&head_bb_tree) && number_of_currently_running_threads < MAXIMUM_NUMBER_OF_CONCURRENTLY_RUNNING_THREADS) {
            BB_Tree_Node *node = RB_MAX(bb_tree, &head_bb_tree); // Achtung: MAX vorausgesetzt
            RB_REMOVE(bb_tree, &head_bb_tree, node);
            Subproblem *subproblem = node->subproblem;
            free(node); // der Knoten wird nicht mehr gebraucht, nur das Teilproblem ist noch wichtig
            // speichere einen Pointer auf das Teilproblem im globalen Array "currently_evaluated_subproblems"
            for (int i = 0; i < MAXIMUM_NUMBER_OF_CONCURRENTLY_RUNNING_THREADS; i++) {
                if (NULL == currently_evaluated_subproblems[i]) {
                    currently_evaluated_subproblems[i] = subproblem;
                    break;
                }
                assert(i < MAXIMUM_NUMBER_OF_CONCURRENTLY_RUNNING_THREADS - 1); // spaetestens im letzten Schleifendurchlauf muss ein "freier Platz" gefunden worden sein
            } // end i
            number_of_currently_running_threads++;
            assert(number_of_currently_running_threads > 0 && number_of_currently_running_threads <= MAXIMUM_NUMBER_OF_CONCURRENTLY_RUNNING_THREADS);
            pthread_create(&thread, NULL, evaluate_subproblem, subproblem); // starte den Thread
            pthread_detach(thread); // niemand muss auf irgendwen warten
        }
        // main-Thread legt sich schlafen, bis er wieder aufgeweckt wird
        pthread_cond_wait(&main_cond, &main_mutex);
    } // end main loop


	pthread_cancel(terminal_output_thread);
	pthread_mutex_unlock(&main_mutex);
	pthread_join(terminal_output_thread, NULL);

	// mache einen letzten Terminal-Output

	pthread_mutex_lock(&main_mutex);
	print_terminal_output(&instance);
	pthread_mutex_unlock(&main_mutex);
	write_best_global_solution_to_file(&instance);
	printf("=============================================================================================================================================\n\n");

	print_optimal_solution(&instance);

	//printf("Time spent in heuristics:\n");
	//printf("Time spent in separation:\n");

	// raeume Instanz auf und gib Speicher frei
	clean_up_instance_data(instance.instance_data);
	clean_up_instance(&instance);


//}
    return EXIT_SUCCESS;
}

void clean_up_instance(Instance *instance) {
	assert(NULL != instance);
	free(instance->instance_data);
	free(instance->best_global_solution);
}

// uebergebene Loesung muss zulaessig sein und den angegebenen Zielfunktionswert besitzen
void update_initial_primalbound_and_solution(Instance *instance, bool *solution, double value) {
	assert(NULL != instance);
	assert(NULL != solution);
	assert(true == is_solution_feasible(instance, solution));
	assert(value == calculate_objective_value(instance, solution));
	if (value > instance->global_primalbound) {
		instance->global_primalbound = value;
		memcpy(instance->best_global_solution, solution, (instance->number_of_variables) * sizeof(bool));
		write_best_global_solution_to_file(instance);
	}
}

void set_number_of_variables(Instance *instance, int number_of_variables) {
	assert(NULL != instance);
	assert(- INT_MAX == instance->number_of_variables);
	assert(number_of_variables > 0);
	instance->number_of_variables = number_of_variables;
}

void update_initial_dualbound(Instance *instance, double dualbound) {
	assert(NULL != instance);
	if (dualbound < instance->global_dualbound) {
		instance->global_dualbound = dualbound;
	}
}


void update_local_primalbound_and_solution_and_insert_subproblem_into_bb_tree(Subproblem *subproblem) {
	assert(NULL != subproblem);
	if (true == is_solution_feasible_in_subproblem(subproblem, subproblem->best_local_solution)) {
		subproblem->local_primalbound = calculate_objective_value(subproblem->instance, subproblem->best_local_solution);
	} else {
		subproblem->local_primalbound = - DBL_MAX;
	}
	BB_Tree_Node *node;
	calloc_makro(node, 1, BB_Tree_Node);
	node->subproblem = subproblem;
	assert(true == check_subproblem_consistency(subproblem));
	RB_INSERT(bb_tree, &head_bb_tree, node);	
}


int subproblem_compare(BB_Tree_Node *node1, BB_Tree_Node *node2) {
	assert(NULL != node1);
	assert(NULL != node2);
    // vergleiche zunaechst anhand der jeweiligen dualbounds
	if ((node1->subproblem)->local_dualbound < (node2->subproblem)->local_dualbound) return -1;
	if ((node1->subproblem)->local_dualbound > (node2->subproblem)->local_dualbound) return +1;
	// waehle dann das Teilproblem mit der hoechsten local_primalbound (beste Loesung des Vaters ist im Kind zulaessig)
	if ((node1->subproblem)->priority < (node2->subproblem)->priority) return -1;
	if ((node1->subproblem)->priority > (node2->subproblem)->priority) return +1;
	// loese ein "Unentschieden" mithilfe der ID auf
	if ((node1->subproblem)->id < (node2->subproblem)->id) return -1;
	if ((node1->subproblem)->id > (node2->subproblem)->id) return +1;
    assert(true);
	return 0; // sollte niemals eintreten
}


/*  diese Funktion wird aufgerufen, nachdem global_primalbound aktualisiert
    worden ist; es werden alle Teilprobleme aus dem BB-Baum entfernt, in denen keine
    bessere Loesung mehr liegen kann; der main_mutex muss beim Aufruf gelocked sein!
*/
void prune_bb_tree(Instance *instance, double new_global_primalbound) {
	assert(NULL != instance);
	assert(instance->global_dualbound >= new_global_primalbound);
	assert(instance->global_primalbound == new_global_primalbound);
	// kein anderer Thread darf gleichzeitig den BB-Baum manipulieren
    BB_Tree_Node *var;
    BB_Tree_Node *nxt;
    for (var = RB_MIN(bb_tree, &head_bb_tree); NULL != var; var = nxt) {
        nxt = RB_NEXT(bb_tree, &head_bb_tree, var);
		Subproblem *subproblem = var->subproblem;
		assert(true == check_subproblem_consistency(subproblem));
        // pruefe, ob subproblem gepruned werden kann
		if (((subproblem->local_dualbound) + 1e-6) < new_global_primalbound + 1.0) {
            // verringere die Anzahl der noch zu untersuchenden Teilprobleme
            ((subproblem->instance)->number_of_left_subproblems)--;
            assert((subproblem->instance)->number_of_left_subproblems >= 0);
            // entferne Teilproblem aus dem BB-Baum, "raeume" es auf und gib den Speicher frei
            RB_REMOVE(bb_tree, &head_bb_tree, var);
			// nicht notwendig, da Subproblem_Data noch nicht erzeugt wurde
			//clean_up_subproblem_data((var->subproblem)->subproblem_data);
            clean_up_subproblem(subproblem);
			free(subproblem);
            free(var);
        }
    }
}



void try_to_update_local_dualbound(Subproblem *subproblem, double dualbound) {
	assert(true == check_subproblem_consistency(subproblem));
	assert(NULL != subproblem);
	assert(dualbound >= subproblem->local_primalbound);
	if (dualbound < (subproblem->local_dualbound)) {
		subproblem->local_dualbound = dualbound;
	}
	assert(true == check_subproblem_consistency(subproblem));
}

void clean_up_subproblem(Subproblem *subproblem) {
	assert(true == check_subproblem_consistency(subproblem));
	assert(NULL != subproblem);
	// Instanz wird natuerlich nicht "gefreed"
    free(subproblem->fixed_variables);
    free(subproblem->fixed_values);
    free(subproblem->best_local_solution);
	free(subproblem->subproblem_data);
}

void try_to_update_local_and_global_primalbound_and_solution(Subproblem *subproblem, bool *solution, double value) {
	assert(true == check_subproblem_consistency(subproblem));
	assert(NULL != subproblem);
	assert(NULL != solution);
	Instance *instance = subproblem->instance;
	assert(value == calculate_objective_value(instance, solution));
	assert(value <= instance->global_dualbound);
	assert((false == is_solution_feasible_in_subproblem(subproblem, solution)) || (value <= (subproblem->local_dualbound)));
	// Loesung muss global zulaessig sein, aber NICHT unbedingt lokal!
	assert(true == is_solution_feasible(instance, solution));
	update_global_primalbound_and_global_solution_and_prune_bb_tree(instance, solution, value);

	if (true == is_solution_feasible_in_subproblem(subproblem, solution)) {
		update_local_primalbound_and_local_solution(subproblem, solution, value);
	}
	assert(true == check_subproblem_consistency(subproblem));
}

void update_local_primalbound_and_local_solution(Subproblem *subproblem, bool *solution, double value) {
	assert(true == check_subproblem_consistency(subproblem));
	assert(NULL != subproblem);
	assert(NULL != solution);
    Instance *instance = subproblem->instance;
    assert(value == calculate_objective_value(instance, solution));
    assert(value <= instance->global_primalbound); // da vorher bereits aktualisiert
    assert(value <= instance->global_dualbound);
    assert(value <= subproblem->local_dualbound);
	assert(true == is_solution_feasible(instance, solution));
	assert(true == is_solution_feasible_in_subproblem(subproblem, solution));
    if (value > subproblem->local_primalbound) { // Achtung: MAX vorausgesetzt
        subproblem->local_primalbound = value;
        memcpy(subproblem->best_local_solution, solution, (instance->number_of_variables) * sizeof(bool));
    }
    assert(true == check_subproblem_consistency(subproblem));
}

bool is_solution_feasible_in_subproblem(Subproblem *subproblem, bool *solution) {
	assert(NULL != subproblem);
	assert(NULL != solution);
	Instance *instance = subproblem->instance;
	assert(true == is_solution_feasible(instance, solution));
	int number_of_variables = instance->number_of_variables;
	bool *fixed_variables = subproblem->fixed_variables;
	bool *fixed_values = subproblem->fixed_values;
	for (int i = 0; i < number_of_variables; i++) {
		if ((true == fixed_variables[i]) && (solution[i] != fixed_values[i])) {
			return false;;
		}
	}
	return true;
}



void write_best_global_solution_to_file(Instance *instance) {
	assert(NULL != instance);
	Instance_Data *instance_data = instance->instance_data;
	int number_of_facilities = instance_data->number_of_facilities;
	int permutation[number_of_facilities];
	int number_of_facilities_on_the_right[number_of_facilities];
	compute_permutation_from_boolean_solution(number_of_facilities, instance->best_global_solution, permutation, number_of_facilities_on_the_right);

	FILE *fptr;
	char filename[5 + strlen(instance->filename)];
	sprintf(filename, "%s.sol", instance->filename);
	fptr = fopen(filename, "a");

	if (NULL == fptr) return;

	fprintf(fptr, "%s value: %.1f\n", (instance->global_primalbound == instance->global_dualbound) ? "optimal" : "objective", instance_data->K - instance->global_primalbound);
	fprintf(fptr, "%d", permutation[0] + 1);
	for (int i = 1; i < number_of_facilities; i++) {
		fprintf(fptr, " %d", permutation[i] + 1);
	}
	fprintf(fptr, "\n\n");

	if (instance->global_primalbound == instance->global_dualbound) {
		for (int i = 0; i < number_of_facilities; i++) {
			fprintf(fptr, "---");
		}
		fprintf(fptr, "\n\n");
	}

	fclose(fptr);
}



/*  diese Funktion darf NICHT aufgerufen werden, wenn der aufrufende Thread bereits den main_mutex gelocked hat,
	da dieser in der Funktion prune_bb_tree() bereits gelocked wird
*/
void update_global_primalbound_and_global_solution_and_prune_bb_tree(Instance *instance, bool *solution, double value) {
    /*  ACHTUNG: zwar wird die Loesung und deren Wert von subproblem selbst uebergeben, aber wir koennen nicht annehmen,
        dass diese Loesung auch zulaessig in diesem Teilproblem ist!!! (manche Heuristiken erzeugen Loesungen, die nicht
        unbedingt im jeweiligen Teilproblem zulaessig sind)
    */
    assert(NULL != instance);
    assert(NULL != solution);
    // pruefe, ob der Zielfunktionswert richtig uebergeben worden ist
    assert(value == calculate_objective_value(instance, solution));
	assert(true == is_solution_feasible(instance, solution));
	assert(value <= (instance->global_dualbound));
    // pruefe, ob die neue Loesung wirklich besser ist
    if (value > instance->global_primalbound) { // Achtung: MAX vorausgesetzt
        // sorge dafuer, dass nur ein Thread gleichzeitig den folgenden Codeabschnitt betreten kann
        pthread_mutex_lock(&solution_mutex); // lock solution_mutex
        // pruefe ein letztes Mal, ob die neue Loesung wirklich besser ist
        if (value > instance->global_primalbound) { // Achtung: MAX vorausgesetzt
            // aktualisiere global_primalbound
            instance->global_primalbound = value;
            // speichere die neue Loesung
            memcpy(instance->best_global_solution, solution, (instance->number_of_variables) * sizeof(bool));
            // versuche, Teilprobleme aus dem BB-Baum rauszuwerfen
			pthread_mutex_lock(&main_mutex);
            prune_bb_tree(instance, value);
			print_terminal_output(instance); // darf aufgerufen werden, weil vorher der main_mutex gelocked wurde

			write_best_global_solution_to_file(instance);

			pthread_mutex_unlock(&main_mutex); // unlock main_mutex
        }
		// Terminal-Output, da global_primalbound aktualisiert worden ist
        pthread_mutex_unlock(&solution_mutex); // unlock solution_mutex
    }
}

void *print_terminal_output_in_loop(void *arg) {
	Instance *instance = (Instance *) arg;
	while (true) {
		sleep(TERMINAL_OUTPUT_FREQUENCY_IN_SECONDS);
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, NULL);
		pthread_mutex_lock(&main_mutex);
		print_terminal_output(instance);
		pthread_mutex_unlock(&main_mutex);
		pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
	}
}

void *evaluate_subproblem(void *arg) {

    Subproblem *subproblem = (Subproblem *) arg;
	Instance *instance = subproblem->instance;
	Instance_Data *instance_data = instance->instance_data;

	Subproblem_Data *subproblem_data = construct_problem(instance_data->number_of_facilities, instance_data->C, subproblem->fixed_variables, subproblem->fixed_values, subproblem->id);

	subproblem->subproblem_data = subproblem_data;
	assert(true == check_subproblem_consistency(subproblem));

    // hier wird die Hauptarbeit verrichtet
    compute_and_update_dualbound_and_run_heuristics(subproblem);
    
    assert(true == check_subproblem_consistency(subproblem));

	bool can_be_pruned = can_subproblem_be_pruned(subproblem); // boolscher Wert muss zwischengespeichert werden, da can_subproblem_be_pruned() nicht im main_mutex-lock aufgerufen werden darf!!

    pthread_mutex_lock(&main_mutex);
    // pruefe, ob gebrancht werden muss
    if (!can_be_pruned) {
        // Branching
        create_and_add_subproblems_to_bb_tree(subproblem);
    }

    (instance->number_of_evaluated_subproblems)++;
    (instance->number_of_left_subproblems)--;

    assert(instance->number_of_evaluated_subproblems <= instance->number_of_created_subproblems);
    assert(instance->number_of_left_subproblems >= 0);
    assert((instance->number_of_evaluated_subproblems) + (instance->number_of_left_subproblems) <= instance->number_of_created_subproblems);
    assert((instance->global_primalbound) >= (subproblem->local_primalbound));
    assert((instance->global_dualbound) >= (subproblem->local_dualbound));

    // "loesche" Teilproblem aus "currently_evaluated_subproblems"
    for (int i = 0; i < MAXIMUM_NUMBER_OF_CONCURRENTLY_RUNNING_THREADS; i++) {
        if (currently_evaluated_subproblems[i] == subproblem) {
            currently_evaluated_subproblems[i] = NULL;
            break;
        }
        assert(i < MAXIMUM_NUMBER_OF_CONCURRENTLY_RUNNING_THREADS - 1); // spaetestens im letzten Schleifendurchlauf muss das Teilproblem gefunden worden sein
    }

	instance_data->evaluations += subproblem_data->evaluations;

	clean_up_subproblem_data(subproblem->subproblem_data);
    clean_up_subproblem(subproblem);
	free(subproblem);

    number_of_currently_running_threads--;

    assert(number_of_currently_running_threads >= 0 && number_of_currently_running_threads < MAXIMUM_NUMBER_OF_CONCURRENTLY_RUNNING_THREADS);
    // wecke main-Thread auf
    pthread_cond_signal(&main_cond);
    pthread_mutex_unlock(&main_mutex);
    return NULL;
}

// diese Funktion darf unter keinen Umstaenden aufgerufen werden, wenn der main_mutex gelocked ist!
bool can_subproblem_be_pruned(Subproblem *subproblem) {
	assert(true == check_subproblem_consistency(subproblem));
    Instance *instance = subproblem->instance;
    assert((instance->global_primalbound) >= (subproblem->local_primalbound)); // Maximierungsproblem
    assert((instance->global_dualbound) >= (subproblem->local_dualbound)); // Maximierungsproblem
    bool return_value = false;

    if (((subproblem->local_dualbound) + 1e-6) < instance->global_primalbound + 1.0) { // Maximierungsproblem
        /*  fuer den unwahrscheinlichen Fall, dass global_primalbound gerade von
            einem anderen Thread veraendert wird und der Wert inkonsistent ist,
            setze einen Lock und frage erneut ab
        */
        pthread_mutex_lock(&solution_mutex);
        if (((subproblem->local_dualbound) + 1e-6) < instance->global_primalbound + 1.0) { // Maximierungsproblem
            // Teilproblem kann "gepruned" werden
            return_value = true;
        }
        pthread_mutex_unlock(&solution_mutex);
    }
    
    assert(true == check_subproblem_consistency(subproblem));

    return return_value;
}

/*  diese Funktion erzeugt aus dem gegebenen "Vaterknoten" einen "Kindknoten";
    alle in jedem Fall notwendigen, generischen Variablen werden korrekt initialisiert und
    Informationen vom Vaterknoten in den Kindknoten uebernommen;
    subproblem_data wird jedoch NICHT initialisiert
*/
Subproblem *create_child_from_father(Subproblem *subproblem_father) {
	//assert(true == check_subproblem_consistency(subproblem_father));
    Instance *instance = subproblem_father->instance;
    Subproblem *subproblem_child;
    calloc_makro(subproblem_child, 1, Subproblem);
    // mache das Kind mit der Instanz betraut
    subproblem_child->instance = instance;
    // vergebe dem Kind eine einzigartige ID
    subproblem_child->id = instance->number_of_created_subproblems;
    // erhoehe die Anzahl der bisher erzeugten Teilprobleme
    (instance->number_of_created_subproblems)++;
    // erhoehe die Anzahl der verbleibenden Teilprobleme
    (instance->number_of_left_subproblems)++;
    // das Kind "erbt" die dualbound des Vaterknotens (schlechter kann es nicht werden)
    subproblem_child->local_dualbound = subproblem_father->local_dualbound;
    // "initialisiere" die lokale primalbound des Kindes
    // da noch keine zulaessige Loesung fuer dieses Teilproblem bekannt ist, ist der Wert der schlechtestmoegliche
    subproblem_child->local_primalbound = - DBL_MAX;
    subproblem_child->priority = false;
    // allokiere Speicher fuer die beste zulaessige Loesung in diesem Teilproblem (bool *best_local_solution)
    calloc_makro(subproblem_child->best_local_solution, instance->number_of_variables, bool);
    // allokiere Speicher fuer das Array, welches angibt, welche Variablen in diesem Teilproblem fixiert sind
    calloc_makro(subproblem_child->fixed_variables, instance->number_of_variables, bool);
    // kopiere die Informationen ueber die fixierten Variaben aus dem Vaterknoten in den Kindknoten
    memcpy(subproblem_child->fixed_variables, subproblem_father->fixed_variables, (instance->number_of_variables) * sizeof(bool));
    // allokiere Speicher fuer das Array, welches die fixierten Werte angibt
    calloc_makro(subproblem_child->fixed_values, instance->number_of_variables, bool);
    // kopiere die Informationen ueber die fixierten Werte aus dem Vaterknoten in den Kindknoten
    memcpy(subproblem_child->fixed_values, subproblem_father->fixed_values, (instance->number_of_variables) * sizeof(bool));
	// "versuche" die beste Loesung des Vaters ins Kind zu uebernehmen:
	memcpy(subproblem_child->best_local_solution, subproblem_father->best_local_solution, (instance->number_of_variables) * sizeof(bool));
	if (true == is_solution_feasible_in_subproblem(subproblem_child, subproblem_child->best_local_solution)) {
		subproblem_child->local_primalbound = calculate_objective_value(instance, subproblem_child->best_local_solution);
		subproblem_child->priority = true;
	}
	assert(true == check_subproblem_consistency(subproblem_child));
    return subproblem_child;
}


Subproblem *create_subproblem(Instance *instance) {
    Subproblem *subproblem;
    calloc_makro(subproblem, 1, Subproblem);
    // mache das Teilproblem mit der Instanz betraut
    subproblem->instance = instance;
    // vergebe eine einzigartige ID
    subproblem->id = instance->number_of_created_subproblems;
    // erhoehe die Anzahl der bisher erzeugten Teilprobleme
    (instance->number_of_created_subproblems)++;
    // erhoehe die Anzahl der verbleibenden Teilprobleme
    (instance->number_of_left_subproblems)++;
    // initialisiere die local_dualbound
    subproblem->local_dualbound = instance->global_dualbound;
    // "initialisiere" die lokale primalbound
    subproblem->local_primalbound = - DBL_MAX;
    subproblem->priority = false;
    // allokiere Speicher fuer das Array, welches angibt, welche Variablen in diesem Teilproblem fixiert sind
    // initialisiere alles mit 0 (=false)
    calloc_makro(subproblem->fixed_variables, instance->number_of_variables, bool);
    // allokiere Speicher fuer das Array, welches die fixierten Werte angibt;
    calloc_makro(subproblem->fixed_values, instance->number_of_variables, bool);
    // allokiere Speicher fuer die beste zulaessige Loesung in diesem Teilproblem (bool *best_local_solution)
    calloc_makro(subproblem->best_local_solution, instance->number_of_variables, bool);
    // kopiere die bisher beste Loesung der Instanz in das Teilproblem
    memcpy(subproblem->best_local_solution, instance->best_global_solution, (instance->number_of_variables) * sizeof(bool));
    if (true == is_solution_feasible_in_subproblem(subproblem, subproblem->best_local_solution)) {
		subproblem->local_primalbound = calculate_objective_value(instance, subproblem->best_local_solution);
		subproblem->priority = true;
	}
	assert(true == check_subproblem_consistency(subproblem));
    return subproblem;
}

// Maximierungsproblem
void set_up_instance(Instance *instance) {
	assert(instance->number_of_variables > 0);
	instance->number_of_created_subproblems = 0;
	instance->number_of_evaluated_subproblems = 0;
	instance->number_of_left_subproblems = 0;
	instance->global_dualbound = DBL_MAX; // Maximierungsproblem
	instance->global_primalbound = - DBL_MAX; // Maximierungsproblem
	calloc_makro(instance->best_global_solution, instance->number_of_variables, bool);
	gettimeofday(&(instance->start_t), NULL);
	instance->cpu_t = clock();
}

// diese Funktion darf nur aufgerufen werden, wenn der aufrufende Thread den main_mutex gelocked hat

// ACHTUNG: es koennte Probleme geben, wenn die dualbound eines Teilproblems gerade aktualisiert wird...
// Wie kann man das umgehen? braucht man dafuer einen weiteren Lock???

void update_global_dualbound(Instance *instance) {
    double bound = instance->global_primalbound;
	// leider koennte jetzt bound < instance->global_primalbound gelten
    BB_Tree_Node *var;
    BB_Tree_Node *nxt;
    // gehe alle Teilprobleme im aktuellen BB-Baum durch
    for (var = RB_MIN(bb_tree, &head_bb_tree); var != NULL; var = nxt) {
        nxt = RB_NEXT(bb_tree, &head_bb_tree, var);
        double temp = (var->subproblem)->local_dualbound;
        if (temp > bound) { // Achtung: MAX vorausgesetztclean_up_instance(&instance);
            bound = temp;
        }
    }
    // gehe alle Teilprobleme, die gerade untersucht werden, durch
    for (int i = 0; i < MAXIMUM_NUMBER_OF_CONCURRENTLY_RUNNING_THREADS; i++) {
        Subproblem *subproblem = currently_evaluated_subproblems[i];
        if (NULL == subproblem) continue; // tritt auf, wenn nicht alle Threads arbeiten
        double temp = subproblem->local_dualbound;
        if (temp > bound) {
            bound = temp;
        }
    }
    assert(bound <= instance->global_dualbound); // es kann sich nicht verschlechtert haben
    instance->global_dualbound = bound;
}

// diese Funktion darf nur aufgerufen werden, wenn der aufrufende Thread den main_mutex gelocked hat
void print_terminal_output(Instance *instance) {
	assert(NULL != instance);
	Instance_Data *instance_data = instance->instance_data;
	// berechne zunaechst die aktuelle globale dualbound
    update_global_dualbound(instance);
	// berechne Anzahl der evaluations der aktuell laufenden Threads
	long evaluations = 0;
	for (int i = 0; i < MAXIMUM_NUMBER_OF_CONCURRENTLY_RUNNING_THREADS; i++) {
		Subproblem *subproblem = currently_evaluated_subproblems[i];
		if (NULL == subproblem) continue; // tritt auf, wenn nicht alle Threads arbeiten
		if (NULL == subproblem->subproblem_data) continue; // Thread laeuft schon, aber hat subproblem_data noch nicht erzeugt
		assert((subproblem->subproblem_data)->evaluations >= 0);
		evaluations += (subproblem->subproblem_data)->evaluations;
	}
	// addiere die #evaluations der schon fertigen Teilprobleme
	evaluations += instance_data->evaluations;

    double K = instance_data->K;
    double dualbound = K - (instance->global_dualbound);
    double primalbound = K - (instance->global_primalbound);
    assert(primalbound >= dualbound); // Wir wollen eigentlich MINimieren!
	static int count = 0;
	if (0 == count) {
		printf("=============================================================================================================================================\n");
		printf("|%16s|%18s|%14s|%14s|%15s|%21s|%17s|%17s|\n",
			"time      ",
			"CPU time     ",
			"explored   ",
			"unexplored  ",
			"evaluations  ",
			"dualbound      ",
			"primalbound   ",
			"gap       ");
		printf("=============================================================================================================================================\n");
	}
	count++;
	count %= 20;

	struct timeval now_t;
	gettimeofday(&now_t, NULL);



	if (dualbound <= 0.0) {
		printf("|%14.2fs |%16.2fs |%13ld |%13ld |%14lu |%20.4f |%16.1f |%16s |\n",
			now_t.tv_sec - instance->start_t.tv_sec + ((now_t.tv_usec - instance->start_t.tv_usec) / (1e6)), // Laufzeit in Sekunden seit Start
			(clock() - (instance->cpu_t)) / ((double) CLOCKS_PER_SEC),
			instance->number_of_evaluated_subproblems,
			instance->number_of_left_subproblems,
			evaluations,
			dualbound,
			primalbound,
			"-----    ");
	} else if (1000.0 <= ((primalbound - dualbound) / dualbound) * 100.0){
		printf("|%14.2fs |%16.2fs |%13ld |%13ld |%14lu |%20.4f |%16.1f |%16s |\n",
			now_t.tv_sec - instance->start_t.tv_sec + ((now_t.tv_usec - instance->start_t.tv_usec) / (1e6)), // Laufzeit in Sekunden seit Start
			(clock() - (instance->cpu_t)) / ((double) CLOCKS_PER_SEC),
			instance->number_of_evaluated_subproblems,
			instance->number_of_left_subproblems,
			evaluations,
			dualbound,
			primalbound,
			"LARGE    ");
	} else {
		printf("|%14.2fs |%16.2fs |%13ld |%13ld |%14lu |%20.4f |%16.1f |%15.8f%% |\n",
	    	now_t.tv_sec - instance->start_t.tv_sec + ((now_t.tv_usec - instance->start_t.tv_usec) / (1e6)), // Laufzeit in Sekunden seit Start
			(clock() - (instance->cpu_t)) / ((double) CLOCKS_PER_SEC),
	    	instance->number_of_evaluated_subproblems,
	    	instance->number_of_left_subproblems,
			evaluations,
	    	dualbound,
	   	 	primalbound,
	    	((primalbound - dualbound) / dualbound) * 100.0);
	}
	
}


// es werden keine Tests mit Daten gemacht, die schrittweise veraendert werden
bool check_instance_consistency(Instance *instance) {
	if (NULL == instance) return false;
	// die "number_of"-Variablen werden nicht getestet
	if (instance->global_primalbound > instance->global_dualbound) return false;
	if (NULL == instance->best_global_solution) return false;
	if (1 > instance->number_of_variables) return false;
	if (NULL == instance->filename) return false;
	if (NULL == instance->instance_data) return false;
	return true;
}

bool check_subproblem_consistency(Subproblem *subproblem) {
	assert(NULL != subproblem);
	assert(NULL != (subproblem->instance));
	assert(0 <= subproblem->id);
	assert(subproblem->local_dualbound <= (subproblem->instance)->global_dualbound);
	assert(subproblem->local_primalbound <= (subproblem->instance)->global_primalbound);
	assert(subproblem->local_primalbound <= (subproblem->instance)->global_dualbound);
	assert(subproblem->local_primalbound <= subproblem->local_dualbound);
	assert(NULL != subproblem->best_local_solution);
	assert(NULL != subproblem->fixed_variables);
	assert(NULL != subproblem->fixed_values);
	assert(((false == is_solution_feasible_in_subproblem(subproblem, subproblem->best_local_solution)) && ((- DBL_MAX) == subproblem->local_primalbound)) || ((true == is_solution_feasible_in_subproblem(subproblem, subproblem->best_local_solution)) && ((- DBL_MAX) != subproblem->local_primalbound)));
	assert((true == is_solution_feasible_in_subproblem(subproblem, subproblem->best_local_solution) && true == is_solution_feasible(subproblem->instance, subproblem->best_local_solution)) || (false == is_solution_feasible_in_subproblem(subproblem, subproblem->best_local_solution)));
	assert(((- DBL_MAX) == subproblem->local_primalbound) || (true == is_solution_feasible_in_subproblem(subproblem, subproblem->best_local_solution)));
	assert((true == is_solution_feasible_in_subproblem(subproblem, subproblem->best_local_solution) && subproblem->local_primalbound == calculate_objective_value(subproblem->instance, subproblem->best_local_solution)) || (false == is_solution_feasible_in_subproblem(subproblem, subproblem->best_local_solution)));	
	return true;
}
