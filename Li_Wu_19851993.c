/******************************************************************************/
// Crann: Detecting Adaptive evolution in Protein Coding DNA sequences.
//
// Copyright 2000 2001 2002 2003 Chris Creevey   
//
//
//	This file is part of Crann
//
//  Crann is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA 
//
//  To Reference this program please use: 
//  For the Software:
//		Creevey, C.J. and McInerney, J.O. (2003) CRANN: detecting adaptive evolution in protein-coding DNA sequences. Bioinformatics 19(13): 1726.
//  For the Algorithm:
//		Creevey, C.J. and McInerney, J.O. (2002) An algorithm for detecting directional and non-directional positive selection, neutrality and negative selection in protein coding DNA sequences. Gene 300(1-2):43-51
//
//
/*****************************************************************************/

#include "externals.h"
#include "Li_Wu_19851993.h"
#include "evolve.h"

/* This program is called as the first part of the Li-Wu algorithim, with a pointer to the two sequences to be checked */
/* and pointers to the values of L0, L1, and L2. Also passed is the value of code which specifies what genetic         */
/* code the sequence is specified in. The program will go through both sequences and count the numbber of 0, 2 & 4     */
/* fold degenerate sites in either sequence. It will then calculate the average of each type between the two sequences */
/* and assign those values L0, L1 & L2 accordingly.																	   */


void count_degen(int *seq1,int *seq2,float *L,int code,int startw,int endw)
	{
	int i=startw, j=0, seq1f[3] = {0,0,0}, seq2f[3] = {0,0,0};
	
	while(i <= endw && seq1[i] != 193 && seq2[i] != 193)
		{
		if(seq1[i] != 64 && seq2[i] != 64 && genetic_codes[code][seq1[i]] != 0 && genetic_codes[code][seq2[i]] != 0 && !deletion[i])  /* If that codon hasn't been marked for complete deletion in any analysis */
			{
			for(j=0; j<3; j++)
				{
				switch(degenerate_sites[code][seq1[i]][j])
					{
					case 0:
						seq1f[0]++;
						break;
					case 2:
						seq1f[1]++;
						break;
					case 4:
						seq1f[2]++;
						break;
					default:
						break;
					}
				switch(degenerate_sites[code][seq2[i]][j])
					{
					case 0:
						seq2f[0]++;
						break;
					case 2:
						seq2f[1]++;
						break;
					case 4:
						seq2f[2]++;
						break;
					default:
						break;
					} 
				}
			}
		i++;
		}
	for(i=0; i<3; i++)
		{
		L[i] = (seq1f[i]+seq2f[i])/2;
		}
	}	
	
	
/* This function is passed with variables P[i], Q[i], L[i], (where i = 0,1,2). It reads through the sequences comparing */
/* them for nucliotide differences. If there is a difference it calls the function transition? to determine whether the */
/* difference is transitional ar transversional, based on that it counts the number of transitions and transversions    */
/* in both sequences that occur at 0,2 and 4 fold degenerate sites (s[0],s[1],s[3] representing 0,2 & 4 fold degenerate */
/* site differences which are tranisitions and v[0], v[1], v[2] representing 0, 2 & 4 fold degenerate site differences  */
/* which are transversions). It then calcultes the values of P[0],P[1] & P[2], which are equal to s[i]/l[i] (i= 0,1,2,) */
/* and represent the proportion of transitional differences of both sequences out of the total number of 0, 2, or 4 fold*/
/* degenerate sites. It also calculates Q[0], Q[1] & Q[2], which is the same as P[i] except for transversions.          */  


void trans_diff(int *seq1,int *seq2, float *P, float *Q, float L[3], int startw, int endw)
	{
	int i=startw, j=0; 
	float s[3] = {0,0,0}, v[3] = {0,0,0};
	
	while(i <= endw && seq1[i] != 193 && seq2[i] != 193)
		{
		if(seq1[i] != 64 && seq2[i] != 64 && genetic_codes[code][seq1[i]] != 0 && genetic_codes[code][seq2[i]] != 0 && !deletion[i]) /* pairwise deletion of gaps and stop codons and If the codon has been marked for complete deletion */
			{	
			for(j=0; j<3; j++)
				{
				switch(is_transitional(codons[seq1[i]][j], codons[seq2[i]][j]))
					{
					case 1:
						
						if(degenerate_sites[code][seq1[i]][j] <= degenerate_sites[code][seq2[i]][j])
							{	
							s[(degenerate_sites[code][seq1[i]][j])/2]++;
							}
						else{ s[(degenerate_sites[code][seq2[i]][j])/2]++; }
						
						break;
					case 0:	
					
						if(degenerate_sites[code][seq1[i]][j] <= degenerate_sites[code][seq2[i]][j])
							{
							v[(degenerate_sites[code][seq1[i]][j])/2]++;
							}
						else{ v[(degenerate_sites[code][seq2[i]][j])/2]++; }
						
						break;
					default:
						break;
					}
				}
			}
		i++;	
		}
	for(i=0; i<3; i++)
		{
		if(L[i] != 0)/**************CC************/
			{
			P[i] = s[i]/L[i];
			Q[i] = v[i]/L[i];
			}
		}
	}
	
	
	
	
/* this function when given two nucliotides will determine whether or not the difference between the two is a transition */
/* or a transversion                                                                                                     */

int is_transitional(char base1, char base2)
	{
	int transition = 0;
	
	switch(toupper(base1))
		{
		case 'U':
			if(base2 == 'C') transition = 1;
			else if(base2 == 'U') transition = 2;
			break;
		case 'C':
			if(base2 == 'U') transition = 1;
			else if(base2 == 'C') transition = 2;
			break;
		case 'A':
			if(base2 == 'G') transition = 1;
			else if(base2 == 'A') transition = 2;
			break;
		case 'G':
			if(base2 == 'A') transition = 1;
			else if(base2 == 'G') transition = 2;
			break;
		default:
			transition = 2;
			break;
		}
	return(transition);
	}
	


/* This function uses Kimura's two-parameter method to estimate the numbers of transitional (A[0], A[1], A[2]) and  */
/* transversional (B[0], B[1], B[2]) substitutions per ith type of site, where i is 0 2 and 4 fold degenerate sites */

void  kimura_method(float *A, float *B, float *a, float *b, float *c, float *varA, float *varB, float *K, float P[3], float Q[3], float L[3] )
	{
	int i=0;
	
	for(i=0; i<3; i++)
		{
		a[i] = 1/(1- (2 * P[i]) - Q[i]);
		b[i] = 1/(1- (2 * Q[i]));
		c[i] = (a[i] - b[i])/2;
		/* next section: mean number of transitional and transversional substitutions per type of site */
		
		A[i] = ((log(a[i]))/2) - ((log(b[i]))/4);
		B[i] = (log(b[i]))/2;
		/* next section: variance */
		
		varA[i] = (((pow(a[i], 2) * P[i]) + (pow(c[i], 2) * Q[i])) - (pow(((a[i] * P[i])+(c[i] * Q[i])), 2))/L[i]);
		varB[i] = (((pow(b[i], 2)) * (Q[i] * (1 - Q[i])))/L[i]);
		
		/*next section: total number of substitutions per ith type of site */
		
		K[i] = A[i]+B[i];
		}
	}
	
	
/* This function is the method outlined in Li & Wu 1985 for estimating the number of substitutions per synonymous site (Ks) */
/* and the number of substitutions per non synonymous site (Ka).															*/									

void li_wu_85(float L[3], float A[3], float B[3], float *K, float *Ks, float *Ka)
	{
	
	/* next section: mean */
	
	*Ks = ((L[1] * A[1])+(L[2] * K[2]))/((L[1]/3)+L[2]);
/*	if(!(*Ks > 0) || !(*Ks < 5)) *Ks = -1; 
*/	
	*Ka = ((L[1] * B[1])+(L[0] * K[0]))/(((2 * L[1])/3)+L[0]);
/*	if(!(*Ka > 0) || !(*Ka < 5)) *Ka = -1;
*/	
	/* next section: variance (not programmed in) */ /*
	if(!(*varKs > 0) || !(*varKs < 5)) *varKs = -1;
	if(!(*varKa > 0) || !(*varKa < 5)) *varKa = -1;
*/	
	}
	
	
/* This function is the method outlined in Li 1993 which revised the formula of the 1985 paper, estimating the number of */
/* substitutions per synonymous site (Ks) and the number of substitutions per non synonymous site (Ka).					 */

void li_wu_93(float L[3], float A[3], float B[3], float varA[3], float varB[3],float a[3], float b[3], float Q[3], float P[3], float c[3], float *Ks, float *Ka,float *varKs, float *varKa)
	{
	
	/* next section: mean */
	
	*Ks = (((L[1] * A[1])+(L[2] * A[2]))/(L[1]+L[2])) + B[2];
/*	if(!(*Ks >= 0) || !(*Ks < 5)) *Ks = -1; 
*/	
	*Ka = A[0] + (((L[0] * B[0])+(L[1] * B[1]))/(L[0]+L[1]));
/*	if(!(*Ka >= 0) || !(*Ka < 5)) *Ka = -1;
*/	/* next section: variance */
	
	*varKs = ((((pow(L[1],2) * varA[1])+(pow(L[2],2) * varA[2]))/pow((L[1]+L[2]),2)) + varB[2]) - ((b[2] * Q[2]) * (((2 * (a[2] * P[2]))-(c[2] * (1-Q[2])))))/(L[1]+L[2]);
/*	if(!(*varKs >= 0) || !(*varKs < 5)) *varKs = -1;
*/	*varKa = (varA[0] + (((pow(L[0],2) * B[0])+(pow(L[1],2) * varB[1]))/(pow((L[0]+L[1]),2))) - ((b[0] * Q[0])*(((2 * (a[0] * P[0]))-(c[0]*(1-Q[0])))/(L[0]+L[1]))));
/*	if(!(*varKa >= 0) || !(*varKa < 5)) *varKa = -1;
*/	
	}
	 
	 



/* this is the main function for the Li-Wu method, which defines the main variables in use, and calls the other functions */

void Li_Wu(void)
	{
	float L[3] = {0,0,0}, P[3] = {0,0,0}, Q[3] = {0,0,0}, A[3] = {0,0,0}, B[3] = {0,0,0}, a[3] = {0,0,0}, b[3] = {0,0,0},
	c[3] = {0,0,0}, varA[3] = {0,0,0}, varB[3] = {0,0,0}, K[3] = {0,0,0}, Ks = 0, Ka = 0, varKs = 0, varKa = 0; 

	int i = 0, j = 0;
	char choice = '\0';
	
	struct sequence *seq1 = start;
	struct sequence *seq2 = start->next;
	struct synon *new = '\0';
	
	clear_results();	
	if(gen_opt[3] == 1){ startw = 0; endw = (start->length/3 -1); }	
	while(seq1->next)
		{		
		if(seq1->tag)
			{			
	/* change1 		
			new = (struct synon *)malloc(sizeof(li_wu_result));  */
			new = malloc(sizeof(li_wu_result));
				if(!new)
					{
					printf("\n\t Out of memory\n");
					clean_exit();
					}
			
			new->seq_num = seq1->seq_num;
			
			new->Ks = malloc((num_of_seqs - seq1->seq_num)*(sizeof(int)));
			if(!new->Ks)
				{
				printf("\n\t Out of memory\n");
				clean_exit();
				}

			new->varKs = malloc((num_of_seqs - seq1->seq_num)*(sizeof(int)));
			if(!new->varKs)
				{
				printf("\n\t Out of memory\n");
				clean_exit();
				}

			new->Ka = malloc((num_of_seqs - seq1->seq_num)*(sizeof(int)));
			if(!new->Ka)
				{
				printf("\n\t Out of memory\n");
				clean_exit();
				}

			new->varKa = malloc((num_of_seqs - seq1->seq_num)*(sizeof(int)));
			if(!new->varKa)
				{
				printf("\n\t Out of memory\n");
				clean_exit();
				}

			/* Assign 0 values to all the arrays in the results list */
			
			for(i=0; i<(num_of_seqs - seq1->seq_num); i++)
				{
				new->Ks[i] = 0;
				new->varKs[i] = 0;
				new->Ka[i] = 0;
				new->varKa[i] = 0;
				}
		
			/* now assign the pointers to the previous and next sequence */
			if (li_wu_end == '\0' || li_wu_start == '\0') 	/* If this is a new list */
				{
				li_wu_end = new;
				li_wu_start = new;
				}		
			else (li_wu_end)->next = new;
			new->next = '\0';
			new->previous = li_wu_end;
			li_wu_end = new;			
			
			seq2 = seq1->next;
			while(seq2 != '\0')
				{
				if(seq2->tag)
					{
					
					for(i=0; i<3; i++)
						{
						L[i] = P[i] = Q[i] = A[i] = B[i] = a[i] = b[i] = c[i] = varA[i] = varB[i] = 0;
						K[i] = 0;
						}
						Ks = 0; Ka = 0; varKs = 0; varKa = 0; 
					
				
					count_degen(seq1->bases, seq2->bases, L, code, startw, endw);
				
					trans_diff(seq1->bases, seq2->bases, P, Q, L, startw, endw);
					
					kimura_method(A, B, a, b, c, varA, varB, K, P, Q, L);
					
					switch(gen_opt[4])  /* gen_opt 4 specifies this choice **/
						{
						case 1:
							li_wu_85(L, A, B, K, &Ks, &Ka);
							break;
						default:
							li_wu_93(L, A, B, varA, varB, a, b, c, Q, P, &Ks, &Ka, &varKs, &varKa);
							break;
						}
					new->Ks[seq2->seq_num - (seq1->seq_num + 1)] = Ks;
					new->varKs[seq2->seq_num - (seq1->seq_num + 1)] = varKs;						
					new->Ka[seq2->seq_num - (seq1->seq_num + 1)] = Ka;
					new->varKa[seq2->seq_num - (seq1->seq_num + 1)] = varKa;	
 					
					}
				seq2 = seq2->next;
				}
			}
		seq1 = seq1->next;
		}
	}


	
/* This function is passed two pointers to sequences which will be analysed with the Li93 method, Returned to the calling
	function is an array holding the results for those sequences results[0] = Ka, results[1] = Ks, results[2] = Varka, results[3] = VarKs */  
void li_wu_complete(int *seq1, int *seq2, float *results)
	{
	
	float L[3] = {0,0,0}, P[3] = {0,0,0}, Q[3] = {0,0,0}, A[3] = {0,0,0}, B[3] = {0,0,0}, a[3] = {0,0,0}, b[3] = {0,0,0},
	c[3] = {0,0,0}, varA[3] = {0,0,0}, varB[3] = {0,0,0}, K[3] = {0,0,0}, Ks = 0, Ka = 0, varKs = 0, varKa = 0; 
	
	
	count_degen(seq1, seq2, L, code, startw, endw);
				
	trans_diff(seq1, seq2, P, Q, L, startw, endw);
					
	kimura_method(A, B, a, b, c, varA, varB, K, P, Q, L);
					
	li_wu_93(L, A, B, varA, varB, a, b, c, Q, P, &Ks, &Ka, &varKs, &varKa);

	results[0] = Ka;
	results[1] = Ks;
	results[2] = varKa;
	results[3] = varKs;
	
	}









void which_sequences(void)
	{
	int choice1 = 0;
	int choice2 = '\0';
	struct sequence *position = start;
	
	tag_all();

		printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
		printf("\n\n\tWould you like to compare all sequences in memory,(1) (default)");
		printf("\n\t Or would you like to select those sequences not to be included (2) ?");
		choice1 = getint("\n\n\tPlease select 1 or 2", 1,2,1);
		switch(choice1)
			{
			case 1:
				break;
			default:
				printf("\nsequence numbers, and names as follows");
				do {
					position = start;
					do{
						if(position->tag) printf("\nNo: %d Name: %s ", ((position->seq_num)+1), (position->name));
						position = position->next;
						}while(position != '\0');		
					choice2 = getint("\n\nPlease select the number of a sequence to be omitted from the analysis, and enter 0 when finished\n", 0, num_of_seqs, 0);
					if(choice2 == 0) choice2 = '\0';
					else
						{
						omit_sequence((choice2 - 1));
						untagged++;
						}
					}while(choice2 != '\0');
				break;

			}

	}
	

void omit_sequence(int choice)
	{
	struct sequence *position = start;
	
	printf("\n\nUntagging the selected sequence no: %d", choice);
	while(position)
		{
		if(position->seq_num == choice)
			{
			printf("\nsequence number: %d", (position->seq_num)+1);
			printf("\nUntagging the selected sequence");
			position->tag = FALSE;
			break;
			}
		position = position->next;
		}
	}


void tag_all(void)
	{
	struct sequence *position = start;
	
	while(position)
		{
		position->tag = TRUE;
		position = position->next;
		}
	untagged = 0;
	}
	
	
	
	
void window_size(int *startw, int *endw)
	{

	
	printf("\n\n\tPlease select the desired start postition for analysis of the sequences (in codons)");
	*startw = getint("\n\t(Press return for the begining of the sequences)", 0, ((start->length)/3) -1, 0);
	printf("\n choice = %d", *startw);
	printf("\n\n\tPlease select the desired end postition for analysis of the sequences (in codons)");
	*endw = getint("\n\t(Press return for the end of the sequences)", 0, ((start->length)/3) -1, ((start->length)/3 - 1));
	printf("\n choice = %d", *endw);
	}

void shiftnwin_size(int *window, int *shift)
	{
	
	printf("\n\n\tPlease select the window size for analysis (in codons)");
	*window = getint("\n\t( or Press return for 10% of the total length)", 1, ((start->length)/3) -1, (((start->length)/3) -1)/10);
	printf("\n choice = %d", *window);
	printf("\n\n\tPlease select the shift size (in codons)");
	*shift = getint("\n\t( or Press return for 50% of window size )", 1, *window, *window/2);
	printf("\n choice = %d", *shift);
	}


void display_Li_Wu(void)
	{
	int i = 0;
	struct synon *position = li_wu_start;
	
	printf("\nDs values\n");
	while(position)
		{
		printf("Sequence num: %d ", position->seq_num);
		for(i=0; i<((num_of_seqs - position->seq_num) - untagged); i++)
			{
			if(fmod(i, 30) == 0) printf("\n");
			printf("%f ", position->Ks[i]);
			}
		position = position->next;
		printf("\n\n");
		}
	position = li_wu_start;
	printf("\nDn values\n");
	while(position)
		{
		printf("Sequence num: %d ", position->seq_num);
		for(i=0; i<((num_of_seqs - position->seq_num) - untagged); i++)
			{
			if(fmod(i, 30) == 0) printf("\n");
			printf("%f ", position->Ka[i]);
			}
		position = position->next;
		printf("\n\n");
		}
		
	position = li_wu_start;
	printf("\nDn\\Ds values\n");
	while(position)
		{
		printf("Sequence num: %d ", position->seq_num);
		for(i=0; i<((num_of_seqs - position->seq_num) - untagged); i++)
			{
			if(fmod(i, 30) == 0) printf("\n");
			printf("%f ", (position->Ka[i])/(position->Ks[i]));
			}
		position = position->next;
		printf("\n\n");
		}
	}










/* this is the main function for the Li-Wu method, which defines the main variables in use, and calls the other functions */

void Li_Wu_movwin(void)
	{
	float L[3] = {0,0,0}, P[3] = {0,0,0}, Q[3] = {0,0,0}, A[3] = {0,0,0}, B[3] = {0,0,0}, a[3] = {0,0,0}, b[3] = {0,0,0},
	c[3] = {0,0,0}, varA[3] = {0,0,0}, varB[3] = {0,0,0}, K[3] = {0,0,0}, Ks = 0, Ka = 0, varKs = 0, varKa = 0; 



	int i = 0, j = 0, startw = 0, endw = 0, moves = 0, window = 0, shift = 0;
	int choice = 0;
	
	struct sequence *seq1 = start;
	struct sequence *seq2 = start->next;
	struct synon *new = '\0';
	

	tag_all();
	clear_results();	
	which_sequences();
	shiftnwin_size(&window, &shift);

	

		printf("\n\tWhich method would you like to use? ");
		printf("\n\t1 = 1985 method ");
		printf("\n\t2 = 1993 method ");
		choice = getint("\n\n\tPlease choose either 1 or 2", 1,2,2);
						
		switch(choice)
			{
			case 1:
				break;
			default:
				break;

			}

			
	
	while(seq1->next)
		{
		seq2 = seq1->next;		
		if(seq1->tag)
			{		
			while(seq2)
				{
				if(seq2->tag)
					{
					/* change3
					new = (struct synon *)malloc(sizeof(li_wu_result));  */
					new = malloc(sizeof(li_wu_result));
					if(!new)
						{
						printf("\n\tOut of memory\n");
						clean_exit();
						}
					
					new->seq_num = seq2->seq_num;
					new->Ks = malloc(((start->length/3 -1)/shift)*sizeof(int));
					if(!new->Ks)
						{
						printf("\n\t Out of memory\n");
						clean_exit();
						}

					new->varKs = malloc(((start->length/3 -1)/shift)*sizeof(int));
					if(!new->varKs)
						{		
						printf("\n\t Out of memory\n");
						clean_exit();
						}

					new->Ka = malloc(((start->length/3 -1)/shift)*sizeof(int));
					if(!new->Ka)
						{
						printf("\n\t Out of memory\n");
						clean_exit();
						}

					new->varKa = malloc(((start->length/3 -1)/shift)*sizeof(int));
					if(!new->varKa)
						{
						printf("\n\t Out of memory\n");
						clean_exit();
						}

/*
					new->Ks = calloc(((start->length/3 -1)), sizeof(int));
					new->Ka = calloc(((start->length/3 -1)/shift) + 2, sizeof(int));
					new->varKs = calloc(((start->length/3 -1)/shift) + 2, sizeof(int));
					new->varKa = calloc(((start->length/3 -1)/shift) + 2, sizeof(int));
					
*/					
					/* assign pointers to the prevoius and next sequence */
					
					if(li_wu_start == '\0' || li_wu_end == '\0')   /* if this is a new list */
						{
						li_wu_end = new;
						li_wu_start = new;
						}
					else li_wu_end->next = new;
					
					new->next = '\0';
					new->previous = li_wu_end;
					li_wu_end = new;
					
					
					startw = 0;
					endw = startw + window;
					moves = 0;
					
					while( (moves*shift) + window <= (start->length/3) -1 )
						{
						for(i=0; i<3; i++)
							{
						    L[i] = 0; P[i] = 0; Q[i] = 0; A[i] = 0; B[i] = 0; a[i] = 0; b[i] = 0;
							c[i] = 0; varA[i] = 0; varB[i] = 0; K[i] = 0;
							}
						 Ks = 0; Ka = 0; varKs = 0; varKa = 0;

						startw = moves * shift;
						endw = (startw + window) -1;
						
						count_degen(seq1->bases, seq2->bases, L, code, startw, endw);
					
						trans_diff(seq1->bases, seq2->bases, P, Q, L, startw, endw);
						
						kimura_method(A, B, a, b, c, varA, varB, K, P, Q, L);
						
						switch(toupper(choice))
							{
							case 'A':
								li_wu_85(L, A, B, K, &Ks, &Ka);
								break;
							default:
								li_wu_93(L, A, B, varA, varB, a, b, c, Q, P, &Ks, &Ka, &varKs, &varKa);
								break;
							}
							
						
						new->Ks[moves] = Ks;
						new->Ka[moves] = Ka;
						new->varKs[moves] = varKs;
						new->varKa[moves] = varKa;
						
						moves++;
						}
					new->Ks[moves] = -2;
					new->Ka[moves] = -2;
					new->varKs[moves] = -2;
					new->varKa[moves] = -2;	
					}
				seq2 = seq2->next;
				}
			}
		seq1 = seq1->next;
		}
	display_movwin(window, shift);
	}
	
			
		


void display_movwin(int window, int shift)
	{
	int i = 0, start = 0;
	float total = 0;
	struct synon *position = li_wu_start;
	
	fprintf(outfile, "\n Accumulative Ds values\n");
	start = 0;
	while(li_wu_start->Ks[i] != -2)
		{
		position = li_wu_start;
		total = 0;
		while(position != '\0')
			{
		/*	fprintf(outfile, "%f \t ", position->Ks[i]);   */
			total += position->Ks[i];
			position = position->next;
			}
		i++;
		fprintf(outfile, "%d-%d\t%f\n",start,(start+window), total);
		start += shift;
		}
	i = 0;
	position = li_wu_start;
	fprintf(outfile, "\nAccumulative Dn values\n");
	start = 0;
	while(li_wu_start->Ks[i] != -2)
		{
		position = li_wu_start;
		total = 0;
		while(position != '\0')
			{
		/*	fprintf(outfile, "%f \t ", position->Ka[i]); */
			total += position->Ka[i];
			position = position->next;
			}
		i++;
		fprintf(outfile, "%d-%d\t%f\n",start,(start+window), total);
		start += shift;

		}
	i = 0;
	position = li_wu_start;
/*	fprintf(outfile, "\nKs/Ka values\n");
	while(li_wu_start->Ks[i] != -2)
		{
		position = li_wu_start;
		total = 0;
		while(position != '\0')
			{
		/*	if(position->Ka[i] != 0) fprintf(outfile, "%f \t ", position->Ks[i]/position->Ka[i]); */
/*			if(position->Ka[i] != 0) total += position->Ks[i]/position->Ka[i];
			position = position->next;
			}
		i++;
		fprintf(outfile, "%f\n", total);
		}
	i = 0;
	position = li_wu_start;
	fprintf(outfile, "\nKa/Ks values\n");
	while(li_wu_start->Ks[i] != -2)
		{
		position = li_wu_start;
		total = 0;
		while(position != '\0')
			{
		/*	if(position->Ks[i] != 0) fprintf(outfile, "%f \t ", position->Ka[i]/position->Ks[i]); */
	/*		if(position->Ks[i] != 0)  total += position->Ka[i]/position->Ks[i];
		/*	else fprintf(outfile, "0.000000 \t");  */
	/*		position = position->next;
			}
		i++;
		fprintf(outfile, "%f\n", total);
		}
*/

	i = 0;
	position = li_wu_start;
	fprintf(outfile, "\nAccumulative varDs values\n");
	start = 0;
	while(li_wu_start->Ks[i] != -2)
		{
		position = li_wu_start;
		total = 0;
		while(position != '\0')
			{
		/*	fprintf(outfile, "%f \t ", position->varKs[i]); */
			total += position->varKs[i];
			position = position->next;
			}
		i++;
		fprintf(outfile, "%d-%d\t%f\n",start,(start+window), total);
		start += shift;
		}
	i = 0;
	position = li_wu_start;
	fprintf(outfile, "\nAccumulative varDn values\n");
	start = 0;
	while(li_wu_start->Ks[i] != -2)
		{
		position = li_wu_start;
		total = 0;
		while(position != '\0')
			{
		/*	fprintf(outfile, "%f \t ", position->varKa[i]); */
			total += position->varKa[i];
			position = position->next;
			}
		i++;
		fprintf(outfile, "%d-%d\t%f\n",start,(start+window), total);
		start += shift;
		}
		if(fflush(outfile) == 0) printf("\n\n\tTable results written to file\n");
	}

	



						
					
			



	
