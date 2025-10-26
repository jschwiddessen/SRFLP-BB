#ifndef INSTANCE_DATA_H
#define INSTANCE_DATA_H

#include "b_and_b_framework.h"
#include "user_structs.h"

void construct_C_and_K(Instance_Data *instance_data);
void read_instance_data(Instance_Data *instance_data, const char *filename);
void clean_up_instance_data(Instance_Data *instance_data);
void read_and_initialize_given_instance_data(Instance *instance, const char *filename);
void check_instance_data(Instance_Data *instance_data);


#endif
