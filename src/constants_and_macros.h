#ifndef CONSTANTS_AND_MACROS_H
#define CONSTANTS_AND_MACROS_H

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#define PRINT_PARAMETERS true

#define GAP_CUTS (- 5e-2)

#define SCALE_FACTOR_OF_TRIANGLE_INEQUALITIES (1.0L / (1.0L + sqrt(3.0L)))

//#define SCALE_FACTOR_OF_TRIANGLE_INEQUALITIES (0.5)

#define SCALE_FACTOR_OF_PENTAGONAL_INEQUALITIES (1.0L / (1.0L + sqrt(5.0L)))

#define SYMMETRY_BREAKING (2)
// 0: keine Symmetriebrechung
// 1: Maschine 1 muss links von Maschine 2 liegen
// 2: zufaellige Ordnungsvariable wird auf zufaelligen Wert fixiert

#define INITIAL_SUBPROBLEMS (0)
// 0: nur (ein) root node
// 1: sofort am Rand fixieren

#define BRANCHING (1)
// 0: most-fractional Branching
// 1: most-fractional Branching auf "benachbarten Maschinen"
// 2: am Rand fixieren

#define WITH_CUTS 0
// 0: keine Cuts
// 1: mit Cuts

#define MAXIMUM_NUMBER_OF_CUTS (1000000)

#define MAXIMUM_NUMBER_OF_CUTS_PER_ITERATION (5000)

#define MINIMUM_NUMBER_OF_CUTS_PER_ITERATION (500)

#define MAXMIUM_NUMBER_OF_MEMORY_CORRECTIONS (10)

#define MAXIMUM_NUMBER_OF_CONCURRENTLY_RUNNING_THREADS (8)

#define TERMINAL_OUTPUT_FREQUENCY_IN_SECONDS (10)

#define MAXMIUM_NUMBER_OF_EVALUATIONS_PER_ITERATION (10000)
// breche L-BFGS-B Solver ab, wenn er "ein neues x gefunden hat" und schon so viele evaluations gefordert hat

#define ALPHA_MAX (1.0)

#define SCALE_ALPHA (0.5)

#define ALPHA_MIN (1e-5)

#define TOLERANCE_MAX (1e-1)

#define SCALE_TOLERANCE (0.95)

#define TOLERANCE_MIN (1e-2)

#define HEURISTIC_LEVEL (3)
// 0: keine Heuristiken (bis auf initiale Heuristik nur primale Loesungen an Blaettern)
// 1: sehr wenig Aufwand (ausreichend, wenn Problem am root node geloest werden soll)
// 2: wenig Aufwand
// 3: hoher Aufwand
// 4: maximaler Aufwand

#define HEURISTIC_FREQUENCY (3)
// 0: keine Heuristiken
// 1: nur, wenn alpha GLEICH sofort verringert wird oder Bound-Berechnung abgeschlossen ist
// 2: nur, wenn alpha GERADE EBEN verringert worden ist oder Bound-Berechnung abgeschlossen ist
// 3: vor jedem Separieren oder wenn Bound-Berechnung abgeschlossen ist


#define calloc_makro(var, size, type) \
    var = (type *) calloc((size) , sizeof(type)); \
    if (NULL == var) { \
        fprintf(stderr, "\n\nError: Not enough memory allocating "#var" variable in %s line %d\n\n", __FILE__, __LINE__); \
        exit(EXIT_FAILURE); \
    }


#define free_makro(var) \
    if (NULL != var) { \
        free(var); \
    }



#endif

