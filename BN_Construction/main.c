#include "parameters.h"
#include "common_usage.h"
#include "common_usage.c"
#include <string.h>
#include <stdio.h>

int main(int argc, char* argv[])
{
  char* att_info = NULL;
  char* attractors_file = NULL;
  char* output_directory = NULL;
	char* num_genes = NULL;

  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i], "--att_info_file"))  att_info = argv[++i];
    else if (!strcmp(argv[i], "--attractors_file"))  attractors_file = argv[++i];
    else if (!strcmp(argv[i], "--output_directory")) output_directory = argv[++i];
		else if (!strcmp(argv[i], "--num_genes")) num_genes = argv[++i];
	}
	int Number_of_Genes=atoi(num_genes);
  char s[100];
  char file_bn[100];
  double ij,kl;

  int no_of_attractor_sets,total_size_att,*store_attractor_set_sizes,*store_attractors;
	int i,j,tem1,tem2,tem3,p_t,t_t,a_t, *predictors, *number_predictors ,*Truth_Table, *store_Truth_Table, no_attractors,n, *NS, no_of_states;
	int *attractors,*store_permutation,counter1,counter2,counter3, next_step_flag;
	int start_point, d_sets,m, *store_basin,*markers;
	time_t t1,t2;
	FILE *fp1,*fp2,*fp3;

     (void) time(&t1);
     srand((long) t1);
     ij = ( 31328 * rand()) / (RAND_MAX +1);
	 kl= ( 30081 * rand()) / (RAND_MAX +1);
	 /* Initialize the random generator by seeds using system random generator rand */
     RandomInitialise(ij,kl);

    tem1 = pow(Discretization, Max_Predictors);
    no_of_states = pow(Discretization, Number_of_Genes);

	fp2=fopen(att_info,"r");
    fscanf(fp2,"%d",&no_of_attractor_sets);

	store_attractor_set_sizes = malloc(sizeof(int) * no_of_attractor_sets);

	for(i=0;i<no_of_attractor_sets;i++)
      fscanf(fp2,"%d",store_attractor_set_sizes+i);

    fscanf(fp2,"%d",&total_size_att);
    fclose(fp2);

    store_attractors = malloc(sizeof(int) * total_size_att);
	fp2=fopen(attractors_file,"r");
    for(i=0;i<total_size_att;i++)
      fscanf(fp2,"%d",store_attractors+i);

     fclose(fp2);

    start_point=0;


for(d_sets=0;d_sets<no_of_attractor_sets;d_sets++)
{
    time_t mytime = time(NULL);
    char * time_str = ctime(&mytime);
    time_str[strlen(time_str)-1] = '\0';
    //printf("Current Time : %s\n", time_str);
   // printf(" The Attractor set number %d\n", (d_sets+1));

    	counter1=0;
	counter2=0;
	counter3=0;



		attractors = malloc(sizeof(int)*(*(store_attractor_set_sizes + d_sets)));


		for(i=start_point;i< (start_point+*(store_attractor_set_sizes + d_sets));i++)
    		*(attractors+i-start_point) = *(store_attractors+i);



		start_point = start_point + *(store_attractor_set_sizes + d_sets);
		no_attractors = *(store_attractor_set_sizes + d_sets);



    for(a_t=0;a_t<Attractor_Tries;a_t++)
	{
        //printf("%d\t",a_t);
	   next_step_flag=1;
	   store_permutation = malloc(sizeof (int) * no_attractors);


	   if(Max_Cycle_Length >1)
           next_step_flag = random_connections_between_attractors(attractors,no_attractors,Number_of_Genes,store_permutation);



	  if(next_step_flag==1)
	  {

		 for(p_t=0;p_t<Predictor_Tries; p_t++)
		 {

			predictors = malloc (sizeof (int) * Number_of_Genes * Max_Predictors);

			number_predictors = malloc (sizeof(int) * Number_of_Genes);

			/* predictors are selected, 3rd input Fixed=0 signifies any number of predictors between 1 and
			pred_max can be selected for each gene */
			SEL_PREDICTORS(Number_of_Genes, Max_Predictors, Fixed, predictors, number_predictors);



			Truth_Table = malloc(sizeof (int) * Number_of_Genes * tem1);

			for(i=0;i<tem1;i++)
			for(j=0;j<Number_of_Genes;j++)
			*(Truth_Table + Number_of_Genes*i + j) = -10;


			 if( Function1(predictors, number_predictors, attractors, no_attractors, Number_of_Genes, Truth_Table,store_permutation)==1)
			 {

			   store_Truth_Table = malloc(sizeof (int) * Number_of_Genes * tem1);
			   for(i=0;i<tem1;i++)
			   for(j=0;j<Number_of_Genes;j++)
			   *(store_Truth_Table + Number_of_Genes*i + j) = *(Truth_Table + Number_of_Genes*i + j);

				for(t_t=0;t_t<TruthTable_Tries ; t_t ++)
				{

				   counter1++;

					for(i=0;i<tem1;i++)
					 for(j=0;j<Number_of_Genes;j++)
					*(Truth_Table + Number_of_Genes*i + j) = *(store_Truth_Table + Number_of_Genes*i + j);

				   fill_table( Number_of_Genes, Truth_Table);

				   NS = (int *) malloc (sizeof (int) * no_of_states);
				   Generate_NS(predictors, number_predictors, Number_of_Genes, Truth_Table, NS);


				   if( Cycle_Checking(Number_of_Genes, NS, attractors, no_attractors) != 1)
				   {
					  counter2++;
						store_basin = calloc(sizeof(int),no_of_states);

						if( (Check_for_Max_Level(Number_of_Genes, NS, store_basin)==1) && (Check_for_Min_Level(Number_of_Genes, NS)==1) )
						 {
							counter3++;
                            //printf("The BN number %d\n", counter3);

							markers =  calloc(sizeof(int),no_of_states);
							m= mark_attractors(attractors, no_attractors, store_permutation, markers);

						//write_basin_sizes(output_directory, attractors, no_attractors,m, no_of_states, markers,store_basin, d_sets);
						write_basin_sizes(output_directory,  attractors, no_attractors,m,no_of_states, markers,store_basin, d_sets);
			 //  printf("dfsdgadgadf\n");

						free(markers);

						/*fp1 = fopen("bns.txt","a");*/


						char file_out[100];
				                sprintf(file_out, output_directory);
						//printf(output_directory);
						//strcat(time_out, "/time.txt");
						sprintf(file_bn, "/bns_%d.txt",d_sets);
						//printf(file_bn);
						sprintf(s, strcat(file_out,file_bn));
						//sprintf(s,"bns_%d.txt",d_sets);

   						//printf("xdsdgadfgsdfhdsfdfsdgadgadf\n");

	                            fp1=fopen(s,"a");
									 for(i=0;i<tem1;i++)
									 {
											 for(j=0;j<Number_of_Genes;j++)
											 {
												 fprintf(fp1,"%d\t",*(Truth_Table + Number_of_Genes*i + j));
											 }
										fprintf(fp1,"\n");
									 }

									tem2=0;
									tem3=-10;
										for(i=0;i<Number_of_Genes;i++)
										{

											for(j=0;j<*(number_predictors+i)-tem2;j++)
												fprintf(fp1,"%d\t",(*(predictors+j+tem2))+1);


											for(j= *(number_predictors+i)-tem2; j < Max_Predictors; j++)
												 fprintf(fp1,"%d\t",tem3);

											fprintf(fp1,"\n");
											tem2= *(number_predictors+i);

										}




									 fprintf(fp1,"\n");



							  fclose(fp1);
						//free kkk;
                             if(counter3>100)
							 {
                               a_t=Attractor_Tries+1;
						       p_t=Predictor_Tries+1;
							   t_t=TruthTable_Tries+1;
							 }

						} // check for level if ends

						free(store_basin);
				   }// check for cycles if ends

				   free(NS);
				}// truth table tries ends


			free(store_Truth_Table);

		 }// attractor feasibility if ends

		 free(Truth_Table);
		 free(number_predictors );
		 free(predictors);
	}// predictor tries ends

    }

	free(store_permutation);

}


free(attractors);
//printf(time_out);
//fp1 = fopen(time_out,"a");
//fp1 = fopen("time.txt","a");
//fprintf(fp1,"%d\t%d\t%d\t%d\t%d\t%d\n",Number_of_Genes, Max_Predictors,counter1,counter2,counter3,d_sets);
//fclose(fp1);
}

//(void) time(&t2);
//fp1 = fopen(time_out,"a");
//fprintf(fp1,"time taken = %lf\n\n\n",difftime(t2,t1));


//fclose(fp1);

    return 0;

}
