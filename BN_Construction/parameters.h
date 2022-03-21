/*
Header file containing the inclusion of standard header files and
defining of constant parameters for the algorithm

Program Written by : Ranadip Pal
Last Modified      : 15th July, 2005
*/



#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>    // For time()
#include <sys/types.h>
//#include <cstdlib>  // For srand() and rand()

/* define the constants */

#define Min_Attractors 1
#define Max_Attractors 6

/* if fixed =0 then any number of predictors between 1 and Max_Predictors can be selected
for each gene, if Fixed=1 then number of predictors for each gene is Max_Predictors */
#define Fixed 0

/*maximum number of predictors*/
#define Max_Predictors 3

/* if attractors are not prespecified, the number of different
random attractor sets to be tried*/
#define Attractor_Tries  1000

/* number of different random predictor sets should be tried */
#define Predictor_Tries  2500

/* number of random filling of truth tables if attractors rae feasible to
the predictor set */
#define TruthTable_Tries 1000

/* this is 2 for binary and 3 for ternary, the implementation for
ternary hasn't been completely done yet */
#define Discretization   2

/* maximum length of cycle allowed, for singleton the cycle length is one */
#define Max_Cycle_Length 1

/* minimum and maximum levels of the generated boolean networks */
#define Min_Level 2
#define Max_Level 10


/*the number of genes */
//#define Number_of_Genes 5

//#define ATTRACTORS_INFO "attractor_info.txt"

//#define ATTRACTORS "attractor_sets.txt"

#define PREDICTOR_FILE "probable_predictors.txt"



//#define SEL_PREDICTORS select_predictors_from_file
#define SEL_PREDICTORS select_predictors


/* Function1 should be fill_singleton_attractors  if the attractors are
assumed to be only singletons or else it should be set to random_connectivity_of_attractors
where cycles are allowed between attractors */
#define Function1 fill_singleton_attractors
//#define Function1 fill_cyclic_attractors



/* Cycle_Checking should be set to Check_for_Cycles if unwanted attractors are to be avoided
or else if other attractors than the specified ones are allowed, then set Cycle_Checking to
Check_for_Cycles_less_restricted */

//#define Cycle_Checking Check_for_Cycles_less_restricted
#define Cycle_Checking Check_for_Cycles


/* constants for the random function */
#define FALSE 0
#define TRUE 1
