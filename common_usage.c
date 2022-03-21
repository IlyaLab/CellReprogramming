/*

Declaration of functions required by the algorithm

Program Written by : Ranadip Pal
Last Modified      : 15th July, 2005
*/


#include "common_usage.h"

/* convert binary or ternary states (x) to decimal value, n is the number of bits */
int dis2dec(char *x, int n){
	int i, decvalue;
	decvalue = 0;

    for (i=0; i<n; i++)
		decvalue = Discretization * decvalue + x[i] - 48;

	return decvalue;
}


/* convert decimal integer (d)  to binary or ternary (ter), n is the number of bits */
void dec2dis(int d, int n, char *ter)
{
	int i;
	for (i=0; i < n; i++){
		ter[n-i-1] = d % Discretization + '0';
		d = d / Discretization;
	}
}


/* generate random permutation of numbers from 0 to p_max-1 and store in 'store_permutation' */
void random_permutation(int p_max, int *store_permutation)
{
	int i,j, k;


	for(i=0;i<p_max;i++)
		*(store_permutation +i)=i;

	for(i=0;i< p_max;i++)
	{

      k = RandomInt(0,(p_max-i)-1);

	  j= *(store_permutation + k);
      *(store_permutation + k) = *(store_permutation + p_max-i-1);
      *(store_permutation + p_max-i-1) = j;
	}

}


void random_attractors(int n, int *no_of_attractors, int *store_attractors)
{

	int i,j,k,*temp1,no_of_states;
	no_of_states=pow(Discretization,n);

	temp1=calloc(sizeof(int),no_of_states);
	if( (*no_of_attractors) == 0)
         *no_of_attractors = RandomInt(Min_Attractors,Max_Attractors);

    for(i=0;i<*no_of_attractors; i++)
	{
		k = RandomInt(0,no_of_states-1);
		while(*(temp1+k)!=0)
           k = RandomInt(0,no_of_states-1);
		*(temp1+k)=1;
	    *(store_attractors + i) = k;
	}
	free(temp1);
}




/* sorting an array ( numbers) of size 'array_size' using bubble sort */
void bubbleSort(int *numbers, int array_size)
{
	int finished = 0, i, j, temp1;

	for (i = (array_size - 1); i >= 0; i--)
	{
		if(finished)
			break;
		finished = 1;
		for (j = 1; j <= i; j++)
		{
			if (*(numbers + j-1) > *(numbers + j))
			{
				finished = 0;
				temp1 = *(numbers + j-1);
				*(numbers + j-1) = *(numbers + j);
				*(numbers + j) = temp1;
			}
		}
	}
}



/* Select Predictors for each gene, if the integer 'fixed' is '1' then all the genes has pred_max number of
predictors and if 'fixed' is '0' then the number of predictors for each gene varies between one and pred_max,
n is the number of genes, predictors stores the predictors as an array and number of predictors stores the
number of predictors for each gene
n= Number of Genes,
pred_max = Maximum number of predictors for each gene
predictors = memory where predictors are stored as a single array and hence
number_predictors gives the number of predictors for each gene such that the
predictors array can be demarcated. for example if predictors= 2 3 1 0 2 and
number_predictors = 2 3, then the predictors for first gene is 2 and 3 and predictors
for 2nd gene is 1,0 and 2.
*/

void select_predictors(int n, int pred_max, int fixed, int *predictors, int *number_predictors)
{
	int i,j,m,a[n],k;

	if(fixed == 1)
	{
		m=0;
		for(i=0;i<n;i++)
		{
				*(number_predictors+i) = m + pred_max;

					/* initialization of an array to keep track of predictors already selected for a gene*/
				   for(j=0;j<n;j++)
					   a[j]=0;

				   for(j=m;j< *(number_predictors+i);j++)
				   {

						k = RandomInt(0,n-1);

						/* while loop to avoid repeated predictors for a gene*/
						while(a[k]!=0)
							k= RandomInt(0,n-1);

						*(predictors + j)= k;
						a[k]=1;
					}
				  m= *(number_predictors+i);
		  }
	}
	else
	{
		  m=0;
          for(i=0;i<n;i++)
		  {
				 /* variable number of predictors selected for each gene */
			     *(number_predictors+i) = m + RandomInt(1,pred_max );


					for(j=0;j<n;j++)
					   a[j]=0;

					for(j=m;j< *(number_predictors+i);j++)
					{
                          k= RandomInt(0,n-1);
						  /* while loop to avoid repeated predictors for a gene*/
						  while(a[k]!=0)
							 k= RandomInt(0,n-1);
						   *(predictors + j)= k;

						   a[k]=1;
					}

					m = *(number_predictors+i);
		  }
	}

}





void select_predictors_from_file(int n,int pred_max, int fixed, int *predictors, int *number_predictors)
{
	int i,j,m,a[n],k, *temp_number_predictors, *temp_predictors,sum_p, tmp1,tmp2,tmp_max;
	FILE *fp;

	fp=fopen(PREDICTOR_FILE,"r");
	temp_number_predictors = malloc(sizeof(int) * n);

	for(i=0;i<n;i++)
		fscanf(fp,"%d",temp_number_predictors+i);


	fscanf(fp,"%d",&sum_p);
    temp_predictors = malloc(sizeof(int) * sum_p);

	for(i=0;i<sum_p;i++)
		fscanf(fp,"%d",temp_predictors+i);

    fclose(fp);

	if(fixed == 1)
	{
		m=0;
		tmp2=0;
		for(i=0;i<n;i++)
		{
				*(number_predictors+i) = m + pred_max;

				tmp1= *(temp_number_predictors+i);
					/* initialization of an array to keep track of predictors already selected for a gene*/
				   for(j=0;j<tmp1;j++)
					   a[j]=0;

				   for(j = m;j< *(number_predictors+i);j++)
				   {

						k = RandomInt(0,tmp1-1);

						/* while loop to avoid repeated predictors for a gene*/
						while(a[k]!=0)
							k= RandomInt(0,tmp1-1);

						*(predictors + j)= *(temp_predictors+tmp2+k);
						a[k]=1;
					}
				  m= *(number_predictors+i);
				  tmp2=tmp2+ tmp1;
		  }
	}
	else
	{
		  m=0;
		  tmp2=0;
          for(i=0;i<n;i++)
		  {

			tmp1= *(temp_number_predictors+i);
			   /* if number of probable predictors are less than pred_max, then the number of predictors
			   selected is between 1 and number of probable predictors */
			     if(tmp1 >= pred_max )
					 tmp_max = pred_max;
				 else
					 tmp_max= tmp1;

				 /* variable number of predictors selected for each gene */
			     *(number_predictors+i) = m + RandomInt(1,tmp_max);


					for(j=0;j<tmp1;j++)
					   a[j]=0;

					for(j=m;j< *(number_predictors+i);j++)
					{
                          k= RandomInt(0,tmp1-1);

						  /* while loop to avoid repeated predictors for a gene*/
						  while(a[k]!=0)
							 k= RandomInt(0,tmp1-1);

						   *(predictors + j)= *(temp_predictors + tmp2 + k);
						   a[k]=1;
					}

					m = *(number_predictors+i);
					tmp2=tmp2+tmp1;
		  }
	}

}



/* this function returns '-1' if the singleton attractors are not compatible with the predictors and returns '1' if compatible.
In case of compatibility, the entries of the 'Truth_Table' corresponding to the attractors are filled.

Predictors and number_predictors are same as in last function,
attractors is an array containing the attractors and no_attractors gives the number of attractors.
n= Number of Genes, Truth_Table is a memory containing the Truth_Table.

store_predictors is redundant here, it's kept only to have the same inputs as the fill_cyclic_attractors function

*/

int fill_singleton_attractors(int *predictors, int *number_predictors, int *attractors, int no_attractors, int n, int *Truth_Table, int *store_permutation)
{
  char *discrete1, *temp1;
  int i,j,k,tem1,tem2,tem3;
  for(i=0;i<no_attractors;i++)
		*(store_permutation +i)=i;
  for(i=0;i<no_attractors;i++)
  {	      //printf("%d\t",*(attractors+i));
		  discrete1 = (char *) malloc (sizeof (char) * n);
		  /* convert the decimal value of the attractor to binary or ternary value and store in discrete1*/
		  dec2dis(*(attractors+i),n,discrete1);
					  tem1=0;
					  for(j=0;j<n;j++)
					  {
						   tem3=*(number_predictors+j) - tem1;
						   temp1 = (char *) malloc(sizeof (char) * tem3);

								 for(k = tem1; k < *(number_predictors+j);k++)
								 {
                                      tem2=k-tem1;
                                 //printf("%d\t", *(predictors + k));
                                // printf("%d\t",*(discrete1));
									  *(temp1 + tem2)= *(discrete1 + *(predictors + k));
								  }
                             //  printf("%s\n","---");
							tem2=dis2dec(temp1,tem3);
							free(temp1);
								  if(*(Truth_Table + tem2*n + j) == -10 || *(Truth_Table + tem2*n + j)== *(discrete1 + j) -48)
								  {

									 *(Truth_Table + tem2*n + j) = *(discrete1 + j) -48 ;

								  }
								  else
									  return -1;
						tem1 = *(number_predictors+j);
					  }
             //  printf("%s\n","attractor düzeyi");

		  free(discrete1);
	  }
  return 1;
}



/* try to check whether a connection from A to B is possible with the current predictors and if possible
fill the 'Truth_Table' corresponding to that transition and return 1, else return -1 if that transition is
not feasible with the current predictors */

int Fill_single_connection(int *Truth_Table, int n, int A, int B ,int *predictors, int *number_predictors)
{

  char *discrete1,*discrete2, *temp1;
  int j,k,tem1,tem2,tem3;


		  discrete1 = (char *) malloc (sizeof (char) * n);
          discrete2 = (char *) malloc (sizeof (char) * n);
		  /* convert the decimal value of the states A and B to binary or ternary value and store in discrete1 and discrete2*/
		  dec2dis(A,n,discrete1);
          dec2dis(B,n,discrete2);

					  tem1=0;
					  for(j=0;j<n;j++)
					  {
						   tem3=*(number_predictors+j) - tem1;
						   temp1 = (char *) malloc(sizeof (char) * tem3);

								 for(k = tem1; k < *(number_predictors+j);k++)
								 {
									  tem2=k-tem1;
									  *(temp1 + tem2)= *(discrete1 + *(predictors + k));
								  }


							tem2=dis2dec(temp1,tem3);
							free(temp1);

								  if(*(Truth_Table + tem2*n + j) == -10 || *(Truth_Table + tem2*n + j)== *(discrete2 + j) - 48)
								  {

									 *(Truth_Table + tem2*n + j) = *(discrete2 + j) - 48 ;
								  }
								  else
								  {
									  free(discrete1);
		                              free(discrete2);
									  return -1;
								  }

						tem1 = *(number_predictors+j);

					  }

		  free(discrete1);
		  free(discrete2);


  return 1;
}







/* It randomly fills the remaining entries of the Truth_Table
n= number of genes */
void fill_table(int n, int *Truth_Table)
{
   int maxlength,i,j;


    maxlength = pow(Discretization, Max_Predictors);

	for(i=0;i<maxlength;i++)
	for(j=0;j<n;j++)
	{
		if(*(Truth_Table + n*i + j) == -10)
           *(Truth_Table + n*i + j)= RandomInt(0,Discretization-1);


	}
}


/* generate the transition Array and store in NS
n= number of genes */
void Generate_NS(int *predictors, int *number_predictors, int n, int *Truth_Table, int *NS)
{

	int i,j,k,no_of_states, tem1,tem2,tem3;
    char *discrete1, *temp1,*temp2;

	no_of_states = pow(Discretization, n);

			for(i=0;i<no_of_states;i++)
			{
			  discrete1 = (char *) malloc (sizeof (char) * n);
			  dec2dis(i,n,discrete1);

				temp2 = (char *) malloc(sizeof (char) * n);
				tem1=0;
						for(j=0;j<n;j++)
						{
								tem3=*(number_predictors+j) - tem1;
								temp1 = (char *) malloc(sizeof (char) * tem3);
										for(k = tem1; k < *(number_predictors+j);k++)
										{
										tem2=k-tem1;
										*(temp1 + tem2)= *(discrete1 + *(predictors + k));
										}
								tem2=dis2dec(temp1,tem3);
								free(temp1);
								*(temp2+j)= *(Truth_Table + tem2*n + j)+48 ;
								tem1 = *(number_predictors+j);
						}
			 *(NS + i) = dis2dec(temp2,n);
            //printf("%d\t",*(NS+i));
			 free(temp2);
			 free(discrete1);

			}


}


/* generate a random connection between the attractors. The random connection is generated by linking the attractors
to a random permutation of the attractors, the check for the maximum cycle length of the attractors is  done in this function */
int random_connections_between_attractors(int *attractors, int no_attractors, int n, int *store_permutation)
{
   int i,j,k,*temp1;



   random_permutation(no_attractors, store_permutation);




   for(i=0;i<no_attractors;i++)
   {
       temp1= calloc(sizeof(int),no_attractors);
       j=i;
	   k=0;
	   while(*(temp1+j)==0)
	   {
              *(temp1+j)=1;
			  j=*(store_permutation+j);
              k=k+1;
	   }

	   free(temp1);



	   if(k > Max_Cycle_Length)
		  return -1;


   }

	    return 1;
}



/* Mark the attractors such that calculating basin sizes for the attractors become easy */
int mark_attractors(int *attractors, int no_attractors, int *store_permutation, int *markers)
{
   int i,j,*temp1,m;

   m=0;

   for(i=0;i<no_attractors;i++)
   {
       temp1= calloc(sizeof(int),no_attractors);
       j=i;
	   m=m+1;
	   while(*(temp1+j)==0)
	   {
              *(temp1+j)=1;
			  *(markers+*(attractors+j))=m;
			  j=*(store_permutation+j);
	   }

	   free(temp1);

   }
return m;

}


void write_basin_sizes(char* folder, int *attractors, int no_attractors, int m, int states, int *markers, int *store_basin, int bn)
{
    int i,j,k, *temp1,*temp2,*temp3;
      FILE *fp;
      char s[100];
      char file_out[100];
      sprintf(file_out, folder);

      float t1;
      char file_basinsize[100];


	sprintf(file_basinsize, "/basinsize_%d.txt",bn);
	sprintf(s, strcat(file_out, file_basinsize));
	fp=fopen(s,"a");

    temp1=calloc(sizeof(int),states);
    temp2= calloc(sizeof(int),m);
    temp3= calloc(sizeof(int),m);

	for(i=0;i<states;i++)
	{
		if(*(markers+i)>0)
			*(temp3+ *(markers+i) -1) = *(temp3+ *(markers+i) -1)+1;

		*(temp2+ *(markers + *(store_basin+i))-1) = *(temp2+ *(markers + *(store_basin+i))-1) + 1;
	}

for(i=0;i<no_attractors;i++)
	fprintf(fp,"%d\t",*(attractors+i));
fprintf(fp,"\n");

for(i=0;i<no_attractors;i++)
{
	t1= (float) *(temp2 + *(markers+*(attractors+i))-1)/(*(temp3+*(markers+*(attractors+i))-1));
	fprintf(fp,"%f\t",t1);
}
fprintf(fp,"\n\n");

fclose(fp);
free(temp1);
free(temp2);
free(temp3);

}



/* The connections between the attractors checked for feasiblity with the current predictor set and
if feasible fill the truth table corresponding to those random connections and return 1 else return -1.
 */
int fill_cyclic_attractors(int *predictors, int *number_predictors, int *attractors, int no_attractors, int n, int *Truth_Table, int *store_permutation)
{
   int i;

       for(i=0;i< no_attractors; i++)
	   {
		  if( Fill_single_connection(Truth_Table, n, *(attractors + i), *(attractors + *(store_permutation+i)), predictors, number_predictors )== -1 )
			  return -1;
	   }

	    return 1;
}


/* check whether cycles of size more than Max_Cycle_Length is present or not
if cycle is present return 1 else return 0
if the network contains attractors other than the specified attractors, then it is rejected
*/
int Check_for_Cycles(int n, int *NS, int *attractors, int no_attractors)
{
	int i,j,k,m, no_of_states, *temp1, *temp2, *temp3;
	no_of_states= pow(Discretization,n);

	temp3= calloc(sizeof (int), no_of_states );

    temp2= calloc(sizeof (int) ,no_of_states);

    for(j=0;j<no_attractors;j++)
					*(temp3 + *(attractors+j))=1;


	if(Max_Cycle_Length >1)
	{

				for(i=0;i<no_of_states;i++)
				{
		   			   if(*(temp2+i)==0)
						{
							   k=0;
							   temp1= malloc(sizeof (int) * no_of_states);

		   					   for(j=0;j<no_of_states;j++)
							   *(temp1+j)=0;

							   m=i;
							   k++;
							   *(temp1+m)= k;


								   while( (*(temp1 + *(NS+m))==0) & (*(temp2 + *(NS+m))==0) )
								   {

									   k++;
									   m = *(NS+m);
									   *(temp1+m)= k;
								   }



							   if (*(temp2 + *(NS+m))==0)
							   {
									   if (  ( (  *(temp1+m) -  *(temp1+ *(NS+m)) ) > (Max_Cycle_Length-1) ) | (*(temp3 + *(NS+m))!=1))
									   {

											   free(temp1);
											   free(temp2);
											   free(temp3);

										   return 1;
									   }
							   }

							   for(j=0;j<no_of_states;j++)
							   {
								if(*(temp1+j)!=0)
									*(temp2+j)=1;
							   }

       					  free(temp1);
						}
				}

	}
	else  // seperate code for singleton attractors only
	 {

					for(i=0; i< no_of_states; i++)
					  {
							   if(*(temp2+i)==0)
								{
									k=0;
									temp1= malloc(sizeof (int) * no_of_states);

		   								for(j=0;j<no_of_states;j++)
										*(temp1+j)=0;


									   m=i;
									   k++;
									   *(temp1+m)= k;

										   while((*(temp1 + *(NS+m))==0) & (*(temp2 + *(NS+m))==0))
										   {

											   k++;
											   m = *(NS+m);
											   *(temp1+m)= k;

										   }


											   if( (*(temp2 + *(NS+m))==0 ) & (*(temp1+m)- *(temp1+ *(NS+m)) >= Max_Cycle_Length ))
											   {   //free(temp1);
											       //free(temp2);
												   //free(temp3);
												   return 1;
											   }
											   if( (*(temp1+m)- *(temp1+ *(NS+m)) == 0) & (*(temp3 + *(NS+m))!=1))
											   {   //free(temp1);
												   //free(temp2);
											       //free(temp3);
												   return 1;
											   }

											   for(j=0;j<no_of_states;j++)
											   {
													if(*(temp1+j)!=0)
														*(temp2+j)=1;
											   }

								   free(temp1);

								}

					  }

		}

    free(temp3);
    free(temp2);
	return 0;
}







/* check whether cycles of size more than Max_Cycle_Length is present or not
if cycle is present return 1 else return 0
if the network contains attractors other than the specified attractors, then also it can be accepted
*/
int Check_for_Cycles_less_restricted(int n, int *NS, int *attractors, int no_attractors)
{
	int i,j,k,m, no_of_states, *temp1, *temp2;

	no_of_states= pow(Discretization,n);
    temp2= calloc(sizeof (int) , no_of_states);

	if(Max_Cycle_Length >1)
	{
				for(i=0;i<no_of_states;i++)
				{
		   			   if(*(temp2+i)==0)
						{
						   k=0;
						   temp1= calloc(sizeof (int) , no_of_states);

						   m=i;
						   k++;
						   *(temp1+m)= k;

						   while( (*(temp1 + *(NS+m))==0) & (*(temp2 + *(NS+m))==0) )
						   {
							   k++;
							   m = *(NS+m);
							   *(temp1+m)= k;
						   }

						   if(  (*(temp2 + *(NS+m))==0) & ((*(temp1+m)- *(temp1+ *(NS+m))) > (Max_Cycle_Length-1)))
						   {
									   free(temp1);
									   free(temp2);
									   return 1;
						   }

						   for(j=0;j<no_of_states;j++)
						   {
							if(*(temp1+j)!=0)
								*(temp2+j)=1;
						   }
       					   free(temp1);
						}

				}
	}
	else  // seperate code for singleton attractors only
	 {


					for(i=0; i< no_of_states; i++)
					  {
							   if(*(temp2+i)==0)
								{
									k=0;
									temp1= malloc(sizeof (int) * no_of_states);

		   							   for(j=0;j<no_of_states;j++)
									   *(temp1+j)=0;

									   m=i;
									   k++;
									   *(temp1+m)= k;

										   while((*(temp1 + *(NS+m))==0) & (*(temp2 + *(NS+m))==0))
										   {

											   k++;
											   m = *(NS+m);
											   *(temp1+m)= k;

										   }


											   if( (*(temp2 + *(NS+m))==0 ) & (*(temp1+m)- *(temp1+ *(NS+m)) >= Max_Cycle_Length ))
											   {   free(temp1);
											       free(temp2);
												   return 1;
											   }


											   for(j=0;j<no_of_states;j++)
											   {
													if(*(temp1+j)!=0)
														*(temp2+j)=1;
											   }

								   free(temp1);

								}

					  }
		}


	free(temp2);
	return 0;
}







/* returns -1 if maximum level of the boolean network with transition array NS is more than
system contstant Max_Level or else return 1 */

int Check_for_Max_Level(int n, int *NS, int *store_basin)
{
	int i,j, no_of_states, *temp1, *temp2, *temp3;
	no_of_states= pow(Discretization,n);


							   temp1= calloc(sizeof (int) , no_of_states);
                               temp2= calloc(sizeof (int) ,no_of_states);
                               temp3= calloc(sizeof (int) ,no_of_states);

            for(i=0;i<no_of_states;i++)
				{
						   *(temp1+i)= *(NS+i) ;
						   *(temp2+i)= *(NS+i) ;
				}


			for(j=0;j<Max_Level-1;j++)
				{
					for(i=0;i<no_of_states;i++)
					{

							   *(temp2+i)= *(NS+ *(temp1+i)) ;
					}
					for(i=0;i<no_of_states;i++)
					{

							   *(temp1+i)= *(temp2+i) ;
					}

				}

			for(i=0;i<no_of_states;i++)
				{

						   *(temp2+i)= *(NS+ *(temp1+i)) ;
						   *(store_basin +i) = *(temp2+i);
						   *(temp3 + *(temp2+i))=1;
				}

            for(i=0;i<no_of_states;i++)
				{


						  if( *(temp3 + *(temp1+i)) !=1)
						  {
							  free(temp1);
							  free(temp2);
							  free(temp3);
							  return -1;
						  }
				}





    free(temp1);
    free(temp2);
	free(temp3);
	return 1;
}



/* returns 1 if minimum level of the boolean network with transition array NS is more than
system contstant Min_Level or else return -1 */

int Check_for_Min_Level(int n, int *NS)
{
	int i,j, no_of_states, *temp1, *temp2, *temp3;
	no_of_states= pow(Discretization,n);


							   temp1= calloc(sizeof (int) , no_of_states);
                               temp2= calloc(sizeof (int) ,no_of_states);
                               temp3= calloc(sizeof (int) ,no_of_states);

            for(i=0;i<no_of_states;i++)
				{
						   *(temp1+i)= *(NS+i) ;
						   *(temp2+i)= *(NS+i) ;
				}


			for(j=0;j<Min_Level-1;j++)
				{
					for(i=0;i<no_of_states;i++)
					{

							   *(temp2+i)= *(NS+ *(temp1+i)) ;
					}
					for(i=0;i<no_of_states;i++)
					{

							   *(temp1+i)= *(temp2+i) ;
					}

				}

			for(i=0;i<no_of_states;i++)
				{

						   *(temp2+i)= *(NS+ *(temp1+i)) ;
						   *(temp3 + *(temp2+i))=1;
				}

            for(i=0;i<no_of_states;i++)
				{


						  if( *(temp3 + *(temp1+i)) == 0 )
						  {
							   free(temp1);
							   free(temp2);
							   free(temp3);
							   return 1;
						  }
				}





    free(temp1);
    free(temp2);
	free(temp3);
	return -1;
}






/* functions to generate random numbers */

/* Globals */
double u[97],c,cd,cm;
int i97,j97;
int test = FALSE;

/*
   This Random Number Generator is based on the algorithm in a FORTRAN
   version published by George Marsaglia and Arif Zaman, Florida State
   University;

   This is the initialization routine for the random number generator.
   NOTE: The seed variables can have values between:    0 <= IJ <= 31328
                                                        0 <= KL <= 30081
   The random number sequences created by these two seeds are of sufficient
   length to complete an entire calculation with. For example, if sveral
   different groups are working on different parts of the same calculation,
   each group could be assigned its own IJ seed. This would leave each group
   with 30000 choices for the second seed. That is to say, this random
   number generator can create 900 million different subsequences -- with
   each subsequence having a length of approximately 10^30.
*/

void RandomInitialise(int ij,int kl)
{
   double s,t;
   int ii,i,j,k,l,jj,m;

   /*
      Handle the seed range errors
         First random number seed must be between 0 and 31328
         Second seed must have a value between 0 and 30081
   */
   if (ij < 0 || ij > 31328 || kl < 0 || kl > 30081) {
		ij = 1802;
		kl = 9373;
   }

   i = (ij / 177) % 177 + 2;
   j = (ij % 177)       + 2;
   k = (kl / 169) % 178 + 1;
   l = (kl % 169);

   for (ii=0; ii<97; ii++) {
      s = 0.0;
      t = 0.5;
      for (jj=0; jj<24; jj++) {
         m = (((i * j) % 179) * k) % 179;
         i = j;
         j = k;
         k = m;
         l = (53 * l + 1) % 169;
         if (((l * m % 64)) >= 32)
            s += t;
         t *= 0.5;
      }
      u[ii] = s;
   }

   c    = 362436.0 / 16777216.0;
   cd   = 7654321.0 / 16777216.0;
   cm   = 16777213.0 / 16777216.0;
   i97  = 97;
   j97  = 33;
   test = TRUE;
}


/*
   This is the random number generator proposed by George Marsaglia in
   Florida State University Report: FSU-SCRI-87-50
*/

double RandomUniform(void)
{
   double uni;

   /* Make sure the initialisation routine has been called */
   if (!test)
   	RandomInitialise(1802,9373);

   uni = u[i97-1] - u[j97-1];
   if (uni <= 0.0)
      uni++;
   u[i97-1] = uni;
   i97--;
   if (i97 == 0)
      i97 = 97;
   j97--;
   if (j97 == 0)
      j97 = 97;
   c -= cd;
   if (c < 0.0)
      c += cm;
   uni -= c;
   if (uni < 0.0)
      uni++;

   return(uni);
}


/*
   Return random integer within a range, lower -> upper INCLUSIVE
*/
int RandomInt(int lower,int upper)
{
   return((int)(RandomUniform() * (upper - lower + 1)) + lower);
}

/*
   Return random float within a range, lower -> upper
*/
double RandomDouble(double lower,double upper)
{
   return((upper - lower) * RandomUniform() + lower);
}
