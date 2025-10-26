#include "constants_and_macros.h"

#include "parameters.h"


void print_parameters() {

	if (false == PRINT_PARAMETERS) return;

	printf("\nParameters:\n");
	printf("Max. Threads: %d\n", MAXIMUM_NUMBER_OF_CONCURRENTLY_RUNNING_THREADS);
	printf("Output frequency: %d\n", TERMINAL_OUTPUT_FREQUENCY_IN_SECONDS);
	printf("With cuts: %d\n", WITH_CUTS);
	printf("Gap cuts: %e\n", GAP_CUTS);
	printf("Scale factor triangle inequalities: %Lf\n", SCALE_FACTOR_OF_TRIANGLE_INEQUALITIES);
	printf("Max. cuts: %d\n", MAXIMUM_NUMBER_OF_CUTS);
	printf("Max. cuts per iteration: %d\n", MAXIMUM_NUMBER_OF_CUTS_PER_ITERATION);
	printf("Min. cuts per iteration: %d\n", MINIMUM_NUMBER_OF_CUTS_PER_ITERATION);
	printf("Memory corrections: %d\n", MAXMIUM_NUMBER_OF_MEMORY_CORRECTIONS);
	printf("Max. alpha: %e\n", ALPHA_MAX);
	printf("Min. alpha: %e\n", ALPHA_MIN);
	printf("Scale alpha: %f\n", SCALE_ALPHA);
	printf("Max. tol: %e\n", TOLERANCE_MAX);
	printf("Min. tol: %e\n", TOLERANCE_MIN);
	printf("Scale tol: %f\n", SCALE_TOLERANCE);
	printf("Heuristic level: %d\n", HEURISTIC_LEVEL);
	printf("Heuristic frequency: %d\n", HEURISTIC_FREQUENCY);
	printf("Max. evals: %d\n", MAXMIUM_NUMBER_OF_EVALUATIONS_PER_ITERATION);
	printf("Symmetry breaking: %d\n", SYMMETRY_BREAKING);
	printf("Initial subproblems: %d\n", INITIAL_SUBPROBLEMS);
	printf("Branching: %d\n", BRANCHING);
}
