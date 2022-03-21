/* header file containing the declaration of the functions defined
in common_usage.c

Program Written by : Ranadip Pal
Last Modified      : 15th July, 2005
*/


#include "parameters.h"


int dis2dec(char *x, int n);  /* binary or ternary to decimal */

void dec2dis(int d,int n, char *ter); /* decimal to binary or ternary */

void random_permutation(int p_max, int *store_permutation); /* generate random permutation of numbers from 0 to p_max-1*/

void random_attractors(int n, int *no_of_attractors, int *store_attractors);

void bubbleSort(int *numbers, int array_size); /* sorting */

void select_predictors(int n,int pred_max, int fixed,int *predictors, int *number_predictors);

void select_predictors_from_file(int n,int pred_max, int fixed, int *predictors, int *number_predictors);

int fill_singleton_attractors(int *predictors, int *number_predictors, int *attractors, int no_attractors, int n, int *Truth_Table, int *store_permutation);

int Fill_single_connection(int *Truth_Table, int n, int A, int B, int *predictors, int *number_predictors);

void fill_table(int n, int *Truth_Table);

void Generate_NS(int *predictors, int *number_predictors, int n, int *Truth_Table, int *NS);

int random_connections_between_attractors(int *attractors, int no_attractors, int n, int *store_permutation);

int mark_attractors(int *attractors, int no_attractors, int *store_permutation, int *markers);

//void write_basin_sizes(char* output_directory, int *attractors, int no_attractors, int m,int states, int *markers,int *store_basin, int bn);

void write_basin_sizes( char* folder, int *attractors, int no_attractors, int m,int states, int *markers,int *store_basin, int bn);

int fill_cyclic_attractors(int *predictors, int *number_predictors, int *attractors, int no_attractors, int n, int *Truth_Table, int *store_permutation);

int Check_for_Cycles(int n, int *NS, int *attractors, int no_attractors);

int Check_for_Cycles_less_restricted(int n, int *NS, int *attractors, int no_attractors);

int Check_for_Max_Level(int n, int *NS, int *store_basin);

int Check_for_Min_Level(int n, int *NS);


/* function to generate random numbers */
void   RandomInitialise(int,int);
double RandomUniform(void);
int    RandomInt(int,int);
double RandomDouble(double,double);

