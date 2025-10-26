#include "constants_and_macros.h"

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#include "user_structs.h"
#include "instance_data.h"
#include "miscellaneous_functions.h"


void read_and_initialize_given_instance_data(Instance *instance, const char *filename) {
	Instance_Data *instance_data;
	calloc_makro(instance_data, 1, Instance_Data);
	instance->instance_data = instance_data;
	read_instance_data(instance_data, filename);
	construct_C_and_K(instance_data);
	int number_of_facilities = instance_data->number_of_facilities;
	int number_of_variables = (number_of_facilities * number_of_facilities - number_of_facilities) / 2;
	set_number_of_variables(instance, number_of_variables);
	instance_data->evaluations = 0;
}

void construct_C_and_K(Instance_Data *instance_data) {

	int number_of_facilities = instance_data->number_of_facilities;
	int *lengths_of_facilities = instance_data->lengths_of_facilities;
	int **weights = instance_data->weights;

	long sum_of_lengths = 0;
	for (int i = 0; i < number_of_facilities; i++) {
		sum_of_lengths += lengths_of_facilities[i];
	}

	long sum_of_weights = 0;
	for (int i = 0; i < number_of_facilities; i++) {
		for (int j = i + 1; j < number_of_facilities; j++) {
			sum_of_weights += weights[i][j];
		}
	}

	long product_of_sum_of_lengths_and_sum_of_weights = sum_of_lengths * sum_of_weights;
	// setze K_times_two
	instance_data->K_times_two = product_of_sum_of_lengths_and_sum_of_weights;
	// setze K
	instance_data->K = product_of_sum_of_lengths_and_sum_of_weights / 2.0;

	int number_of_variables = (number_of_facilities * number_of_facilities - number_of_facilities) / 2;

	long **C;
	calloc_makro(C, number_of_variables, long *);
	// setze C
	instance_data->C = C;

	for (int i = 0; i < number_of_variables; i++) {
		calloc_makro(C[i], number_of_variables, long);
	}

	for (int i = 0; i < number_of_facilities; i++) {
		for (int j = i + 1; j < number_of_facilities; j++) {
			for (int k = 0; k < i; k++) {
				C[inMa(number_of_facilities, k, i)][inMa(number_of_facilities, k, j)] += weights[i][j] * lengths_of_facilities[k];
				C[inMa(number_of_facilities, k, j)][inMa(number_of_facilities, k, i)] += weights[i][j] * lengths_of_facilities[k];
			}
			for (int k = i + 1; k < j; k++) {
				C[inMa(number_of_facilities, i, k)][inMa(number_of_facilities, k, j)] -= weights[i][j] * lengths_of_facilities[k];
				C[inMa(number_of_facilities, k, j)][inMa(number_of_facilities, i, k)] -= weights[i][j] * lengths_of_facilities[k];
			}
			for (int k = j + 1; k < number_of_facilities; k++) {
				C[inMa(number_of_facilities, i, k)][inMa(number_of_facilities, j, k)] += weights[i][j] * lengths_of_facilities[k];
				C[inMa(number_of_facilities, j, k)][inMa(number_of_facilities, i, k)] += weights[i][j] * lengths_of_facilities[k];
			}
		}
	}

}


void check_instance_data(Instance_Data *instance_data) {

	int number_of_facilities = instance_data->number_of_facilities;
	int *lengths_of_facilities = instance_data->lengths_of_facilities;
	int **weights = instance_data->weights;
	
	if (1 > number_of_facilities) {
		fprintf(stderr, "\nERROR: The number of facilities must be at least 1.\n");
		exit(EXIT_FAILURE);
	}
	
	for (int i = 0; i < number_of_facilities; i++) {
		if (0 > lengths_of_facilities[i]) {
			fprintf(stderr, "\nERROR: Negative length of a facility detected.\n");
			exit(EXIT_FAILURE);
		}
	}

	// Diagonaleintraege des Transportaufkommens muessen gleich 0 sein
	for (int i = 0; i < number_of_facilities; i++) {
		if (0 != weights[i][i]) {
			fprintf(stderr, "\nERROR: Nonzero entry on the diagonal of the weight matrix detected.\n");
			exit(EXIT_FAILURE);
		}
	}
	
	bool only_upper_triangular_matrix = true;
	bool symmetric_matrix = true;
	for (int i = 0; i < number_of_facilities; i++) {
		for (int j = i + 1; j < number_of_facilities; j++) {
			if (weights[i][j] != weights[j][i]) {
				symmetric_matrix = false;
			}
			if (0 != weights[j][i]) {
				only_upper_triangular_matrix = false;
			}
		}
	}
	
	if (false == symmetric_matrix && false == only_upper_triangular_matrix) {
		fprintf(stderr, "\nERROR: The weights must be given either as a upper triangular matrix or as a symmetric matrix.\n");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < number_of_facilities; i++) {
		if (0 == lengths_of_facilities[i]) {
			printf("WARNING: Zero length of a facility detected.\n");
			break;
		}
	}

	if (true == only_upper_triangular_matrix && true == symmetric_matrix) {
		printf("WARNING: All weights are zero. Optimization is obsolete.\n");
	} else if (true == only_upper_triangular_matrix) {
		printf("Input weights are given as an upper triangular matrix.\n");
	} else {
		printf("Input weights are given as a symmetric matrix.\n");
	}
	
	
	printf("\nNumber of facilities: %d\n", number_of_facilities);
	int min_length = INT_MAX;
	int max_length = INT_MIN;
	for (int i = 0; i < number_of_facilities; i++) {
		min_length = (lengths_of_facilities[i] < min_length) ? lengths_of_facilities[i] : min_length;
		max_length = (lengths_of_facilities[i] > max_length) ? lengths_of_facilities[i] : max_length;
	}
	printf("Lengths of facilities range: [%d, %d]\n", min_length, max_length);
	
	int min_weight = INT_MAX;
	int max_weight = INT_MIN;
	for (int i = 0; i < number_of_facilities; i++) {
		for (int j = i + 1; j < number_of_facilities; j++) {
			min_weight = (weights[i][j] < min_weight ) ? weights[i][j] : min_weight;
			max_weight = (weights[i][j] > max_weight) ? weights[i][j] : max_weight;
		}
	}
	printf("Weights range: [%d, %d]\n", min_weight, max_weight);
}


// es wird automatisch dafuer gesorgt, dass das Transportaufkommen ("weights") als symmetrische Matrix gespeichert wird
// es muss in "filename" komplett oder als obere Dreiecksmatrix gegeben sein
void read_instance_data(Instance_Data *instance_data, const char *filename) {

	// versuche die Datei "filename" lesend zu oeffnen
	FILE *fp = fopen(filename, "r");
	if (NULL == fp) {
		fprintf(stderr, "\nERROR: Could not open input file %s.\n", filename);
		exit(EXIT_FAILURE);
	}

	printf("Instance name: %s\n", filename);

	// verwende einen hinreichend grossen Buffer (auch fuer 100 Maschinen)
	char buffer[4096];

	// lese erste Zeile, d.h. die Anzahl der Maschinen ein
	fgets(buffer, 4096, fp);
	int number_of_facilities = atoi(buffer);
	
	if (1 > number_of_facilities) {
		fprintf(stderr, "\nERROR: The number of facilities must be at least 1.\n");
		exit(EXIT_FAILURE);
	}
	
	// setze number_of_facilities
	instance_data->number_of_facilities = number_of_facilities;

	// allokiere Speicher fuer das Array mit den Maschinenlaengen
	int *lengths_of_facilities;
	calloc_makro(lengths_of_facilities, number_of_facilities, int);
	// setze lengths_of_facilities
	instance_data->lengths_of_facilities = lengths_of_facilities;

	// lese Maschinenlaengen ein
	fgets(buffer, 4096, fp);
	char delimiter[] = " ,;\t"; // Trennzeichen in der Instanz
	char *ptr; // zeigt immer auf den naechsten Token

	ptr = strtok(buffer, delimiter);

	for (int i = 0; i < number_of_facilities; i++) {
		lengths_of_facilities[i] = atoi(ptr);
		ptr = strtok(NULL, delimiter);
	}

	// allokiere Speicher fuer das 2D-Array mit dem Transportaufkommen
	int **weights;
	calloc_makro(weights, number_of_facilities, int *);
	// setze weights
	instance_data->weights = weights;
	for (int i = 0; i < number_of_facilities; i++) {
		calloc_makro(weights[i], number_of_facilities, int);
	}

	// lese das Transportaufkommen in das Array weights ein
	for (int row = 0; row < number_of_facilities; row++) {
		fgets(buffer, 4096, fp);
		ptr = strtok(buffer, delimiter);
		for (int column = 0; column < number_of_facilities; column++) {
			weights[row][column] = atoi(ptr);
			ptr = strtok(NULL, delimiter);
		}
	}

	// gib Fehler und Warnungen aus
	check_instance_data(instance_data);


	// stelle sicher, dass das Transportaufkommen eine symmetrische Matrix ist
	for (int i = 0; i < number_of_facilities; i++) {
		for (int j = i + 1; j < number_of_facilities; j++) {
			weights[j][i] = weights[i][j];
		}
	}
	
	fclose(fp); // Datei schliessen

}


void clean_up_instance_data(Instance_Data *instance_data) {
	int number_of_facilities = instance_data->number_of_facilities;
	free(instance_data->lengths_of_facilities);
	for (int i = 0; i < number_of_facilities; i++) {
		free((instance_data->weights)[i]);
	}
	free(instance_data->weights);
	int number_of_variables = (number_of_facilities * number_of_facilities - number_of_facilities) / 2;
	for (int i = 0; i < number_of_variables; i++) {
		free((instance_data->C)[i]);
	}
	free(instance_data->C);
}
