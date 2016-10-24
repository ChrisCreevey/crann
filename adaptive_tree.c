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
#include "adaptive_tree.h"
#include "evolve.h"
#include "linked tree.h"

/* this is the program to run the algorithm devised by Kreitman and McDonald using trees */


void adaptive_tree(int **ratio)
	{
	int i = 0, j=0, count = 0, mutation = FALSE, x = 0, y = 0, number = 0, gaps = 0, last = FALSE, done = FALSE;
	struct sequence *seq1 = start, *position = NULL;
	float ***subst_matrix;   /* holds the likelihood of nucliotide substitutions at each fold site */
	char nuc = '\0';

	subst_matrix = malloc(3*sizeof(float**));
	if(subst_matrix == NULL)
		{
		printf("Out of memory\n");
		clean_exit();
		}
		for(y=0; y<3; y++)
			{		
			subst_matrix[y] = malloc(5*sizeof(float*));
			if(subst_matrix[y] == NULL)
				{
				printf("Out of memory\n");
				clean_exit();
				}
			for(x=0; x<5; x++)
				{
				subst_matrix[y][x] = malloc(5*sizeof(float));
				if(subst_matrix[y][x] == NULL)
					{
					printf("Out of memory\n");
					clean_exit();
					}
				}
			}
	for(i=0; i<3; i++)
		{
		for(x=0; x<5; x++)
			for(y=0; y<5; y++)
				subst_matrix[i][x][y] = 0;
		}
	substitution_matrix(subst_matrix);  /* calculate the substitution matrix for the given data */
	
/*	printf("after making tree\n");
		number = 0;
	if(tree_top != NULL)
		{
		check_tree(tree_top, &number);  *//* make sure the tree has been created correctly */	
/*		}
	else printf("tree not defined\n");
	printf("press return to continue - adaptive tree");
	getchar();
*/

	define_outgroup(); /* define the outgroup and give direction to the tree */
/*	printf("after define outgroup");
		number = 0;
	if(tree_top != NULL)
		{
		check_tree(tree_top, &number);  *//* make sure the tree has been created correctly */	
/*		}
	else printf("tree not defined\n");
	printf("press return to continue - adaptive tree");
	getchar();
*/

	
	done = prune_tree(tree_top);
	while(done == TRUE)  /* This checks to make sure that there are no obsolete nodes on the tree */
		done = prune_tree(tree_top);

/*	printf("Finished pruning tree (in Adaptive_tree) \n");
	number = 0;
	if(tree_top != NULL)
		{
		check_tree(tree_top, &number);  *//* make sure the tree has been created correctly */	
/*		}
	else printf("tree not defined\n");
	printf("press return to continue - adaptive tree");
	getchar();
*/

	i = 0;
	for(i=0; i<(start->length)/3; i++)
		{
		for(j=0; j<3; j++)
	 		{
			
			ancestral_codon(i, subst_matrix);	
	/*		ancestral_nuc(i, j, subst_matrix); *//* calculate the ancestral nucliotide for this position on each node */									
			
										
			}
		for(j=0; j<3; j++)
			{
			nuc = 's'; mutation = FALSE;
			is_mutation(tree_top, &mutation, &nuc, i, j); /* we need to know if a mutation has occured in any sequence at this point */
			if(mutation)
				{
				calculate_gaprun(tree_top, i);   /* claculate gapruns at every position on the tree for i */
				count = 0;
				is_fixed(tree_top, i, j, ratio, &count, last);  /* if a mutation has occured then calculate is_fixed otherwise don't */	
				last = TRUE;
				}
			else
				{ 
				last = FALSE;
				if(toupper(nuc) != 'X')/* this is in case every sequence in the alignment has a gap here, while unlikely, without this it would think that a gaprun  */
					increment_gaprun(tree_top);  /* increment every gaprun variable since there is a gap here at every sequence */
				else  
					reset_gaprun(tree_top);  /* see explaination with function */
				}
			}
		}

	if(seq1->bases[i] == 193) close_ances(tree_top, i);
	
	
	

/*	number = 0;
	if(tree_top != NULL)
		{
		check_tree(tree_top, &number);  *//* make sure the tree has been created correctly */
/*		printf("press return to continue - adaptive tree2    ---- after counting\n");
		getchar();	
		}
	else printf("tree not defined\n");
*/	
	
	if(ances_file != NULL)
		{
		fclose(ances_file);
		}
	if((ances_file = fopen("ancestor.out", "w")) == NULL)		/* check to see if the file is there */
		{
		printf("\n\n\tCannot open the  file, ancestor.out\n\t ancestor file not written\n");
		}
	else
		{
		count = 0;
		output_ancestors(tree_top, &count);
		fflush(ances_file);
		fclose(ances_file);
		ances_file = NULL;
		printf("\n\n\n\n\n\n\n\nAncestors and current taxa written to file ancestor.out \n");
		printf("Phylogenetic tree with significances written to result-tree.ph\n");
		}
	
	/* free up allocated memory */
	if(subst_matrix != NULL)
		{
		for(i=0; i<3; i++)
			{
			for(j=0; j<5; j++)
				free(subst_matrix[i][j]);
			free(subst_matrix[i]);
			}
		free(subst_matrix);
		}
	
	}	
	
/* for each i we check if there is a gap at any position, and increment the gaprun variable for that taxa or ancestor
	if there is no gap at a particular position gaprun is assigned to 0.   Doing this at this stage means we don't risk counting 
	the same position on the tree more than once, which may happen when the tree is traversed from within the program */
void calculate_gaprun(struct node *position, int i)
	{
	if(position->node1 != NULL) calculate_gaprun(position->node1, i);
	if(position->node2 != NULL) calculate_gaprun(position->node2, i);
	
	if(position->seq_num1 != NULL)
		{	
		if((position->seq_num1)->bases[i] == 64)
			(position->seq_num1)->gaprun++;
		else
			(position->seq_num1)->gaprun = 0;
		}
	if(position->seq_num2 != NULL)
		{	
		if((position->seq_num2)->bases[i] == 64)
			(position->seq_num2)->gaprun++;
		else
			(position->seq_num2)->gaprun = 0;
		}
	if(position->ances_seq[i] == 64)
		position->gaprun++;
	else
		position->gaprun = 0;

	}
	
	
	




/* This function is called when we come across a section of the sequences where there is no change. It resets all the gaprun variables to 0
	so the program doesn't mistakenly think that the the last difference is part of this new one. */
void reset_gaprun(struct node *position)
	{
	if(position->node1 != NULL) reset_gaprun(position->node1);
	if(position->node2 != NULL) reset_gaprun(position->node2);
	
	if(position->seq_num1 != NULL)
		(position->seq_num1)->gaprun = 0;
		
	if(position->seq_num2 != NULL)
		(position->seq_num2)->gaprun = 0;
		
	position->gaprun = 0;
	
	}	


/* in the unlikely event that there is a gap at this position in every sequence, this increments gaprun at every position forthe tree */
void increment_gaprun(struct node *position)
	{
	if(position->node1 != NULL) reset_gaprun(position->node1);
	if(position->node2 != NULL) reset_gaprun(position->node2);
	
	if(position->seq_num1 != NULL)
		(position->seq_num1)->gaprun++;
		
	if(position->seq_num2 != NULL)
		(position->seq_num2)->gaprun++;
		
	position->gaprun++;
	
	}

	
	
/* Is gap travels down through the tree checking if any of the mutations are a gap, if it is it reports this 
	and depending on whether the previous mutation was next to it AND contained gap this mutation is not counted in the 
	McD and K algorithm, since a row of gaps on a tree would generally be from a single insertion or deletion,
	and not individual events as the program would otherwise assume */
	
int is_gap(int i, int j, struct node *position)
	{
	int ans = FALSE;
	
	if(position->seq_num1 != NULL)
		if((position->seq_num1)->bases[i] == 64)
			ans = TRUE;
	if(position->seq_num2 != NULL)
		if((position->seq_num2)->bases[i] == 64)
			ans = TRUE;
	
	
	if(position->node1 != NULL && ans == FALSE) ans = is_gap(i,j, position->node1);
	if(position->node2 != NULL && ans == FALSE) ans = is_gap(i,j, position->node2);


	return(ans);
	}



/* This function assigns the value 193 to the end of the ancestral sequence to signify that that is the end */
void close_ances(struct node *position, int i)
	{
	
	if(position->node1 != NULL) close_ances(position->node1, i);
	if(position->node2 != NULL) close_ances(position->node2, i);
	
	position->ances_seq[i] = 193;
	
	}


/* this function outputs the ancestral sequences to a file along with the original sequences - all in fasta format */
void output_ancestors(struct node *position, int *count)
	{
	int i = 0, j = 0;
	
	
	if(position->node1 != NULL && position->seq_num1 == NULL) output_ancestors(position->node1, count);
	if(position->node2 != NULL && position->seq_num2 == NULL) output_ancestors(position->node2, count);
	
	
		
		fprintf( ances_file, ">Node%d Ancestor sequence\n", *count); 

		for(i=0; i< (start->length)/3; i++)
			{
			for(j=0; j<3; j++)
				{
				fprintf(ances_file, "%c", codons[position->ances_seq[i]][j]);
				}
			if(fmod(i+1, 20) == 0 ) fprintf(ances_file, "\n");
			}
		fprintf(ances_file, "\n");
		

		if(position->seq_num1 != NULL)
			{
			i = 0;
			fprintf( ances_file, "> %s\n", (position->seq_num1)->name );
			for(i=0; i< (start->length)/3; i++)
				{
				for(j=0; j<3; j++)
					{
					fprintf(ances_file, "%c", codons[(position->seq_num1)->bases[i]][j]);
					}
				if(fmod(i+1, 20) == 0) fprintf(ances_file, "\n");
				}
			fprintf(ances_file, "\n");
			}

		if(position->seq_num2 != NULL)
			{
			i = 0;
			fprintf( ances_file, "> %s\n", (position->seq_num2)->name );
			for(i=0; i< (start->length)/3; i++)
				{
				for(j=0; j<3; j++)
					{
					fprintf(ances_file, "%c", codons[(position->seq_num2)->bases[i]][j]);
					}
				if(fmod(i+1, 20) == 0) fprintf(ances_file, "\n");
				}
			fprintf(ances_file, "\n");
			}	

	*count = *count +1;
	}
		
		
		
		
/* This function travels down through the tree to see if there are any mutations at this position (if there is is_fixed
	will be called to see if for each node the mutation is fixed or not ) */
	
void is_mutation(struct node *position,int *mutation, char *nuc, int i, int j)
	{
	

	if(position->seq_num1 != NULL)
		{
		if(*nuc == '\0') {*nuc = codons[(position->seq_num1)->bases[i]][j];
			}
		else
			{
			if(*nuc != codons[(position->seq_num1)->bases[i]][j])
				{
				*mutation = TRUE;
				}
			}
		}
	if(position->seq_num2 != NULL)
		{
		if(*nuc == '\0') { *nuc = codons[(position->seq_num2)->bases[i]][j];
			}
		else
			{
			if(*nuc != codons[(position->seq_num2)->bases[i]][j])
				{
				*mutation = TRUE;
				}
			}
		}
	

	if(*mutation == FALSE)
		if(position->node1 != NULL)
			is_mutation(position->node1, mutation, nuc, i, j);
	
	if(*mutation == FALSE)
		if(position->node2 != NULL)
			is_mutation(position->node2, mutation, nuc, i, j);
			
	} 



/* this function is called with the following variables passed to it, i = codon number we are at, j the position of the codon we are at */
/* tree_description a pointer to the tree array created earlier, and diff, which is the number of the nucliotide (0 = U, 1 = C, 2 = A, 3 = G)*/
/*which occurs as the most common difference in that position. It uses the protocol defined by Mc Donald and Kreitman to determine whether a mutation is fixed or */
/* polymorphic, given the phylogenetic tree array. It increments the array ratio which is passed to it.                 		  */

void is_fixed(struct node *position, int i, int j, int **ratio, int *count, int last)
	{
	char nuc = 'a';
	int fixed = TRUE, synon = TRUE, value = 0, polymorphic = FALSE;
	int fixednew = TRUE, repl = 0, silent = 0;	
	
	if(position->node1 != NULL) is_fixed(position->node1, i, j, ratio, count, last);


	if(position->node2 != NULL) is_fixed(position->node2, i, j, ratio, count, last);
	
	if(position->prev != NULL)
		{
		
		
		if(genetic_codes[code][position->ances_seq[i]] != 0)	 /* It won't count it if the ancestral codon is a stop */
			{
			snapshot(position, ratio, &nuc, count, i, j);
			}	
		
				
		*count = *count + 1;
		}
	}			

/* This function implements an idea that the measure of whether a population is under selection or not
	is a function of the complete timeline from the snapshot in time represented by the internal branch, to 
	the moment the sequences in the clade were sampled. The question you have to ask is, whether an change occured
	int the enviornment at this point in time from which point on the sequences in the clade have been trying to optimise
	their survival rate in the enviornment?  Adaption doesn't occur at one timepoint, but it is a progressive thing that
	occurs after a change in the enviornment.
*/

void snapshot(struct node *position, int **ratio, char *nuc, int *count, int i, int j)
	{
	int fixed = TRUE, silent = FALSE, value = 0;

	/* The ratio is divided as follows: there is a ratio for each node in the tree. that ratio is */
	/* 0: repl/fix	1: repl/poly	2: synon/fix	3: synon/poly								*/
	/* for each character the appropriate ratio is incremented according to which node it is at */
	
	/* The variable Graphs is used to hold the results of the histograms at each node */

	/* Check for fixations within the clade, based on a change in the ancestors */
	if(codons[(position->prev)->ances_seq[i]][j] != codons[position->ances_seq[i]][j])
		{
		check_fixed_in(position, &fixed, &codons[position->ances_seq[i]][j], i, j);
			
		value = (position->prev)->ances_seq[i] - transform_base(codons[(position->prev)->ances_seq[i]][j], j) + transform_base(codons[position->ances_seq[i]][j], j);
		if(value > 64) value = 64; /* if there is a non standard charater tha whole codon is treated as a gap, and given the value of 64 to represent that */			
		if(genetic_codes[code][value] != genetic_codes[code][(position->prev)->ances_seq[i]])			
			silent = FALSE;
		else
			silent = TRUE;
			
				
		if(position->gaprun < 2 && ((position->prev)->gaprun < 2))  /* if this ancestor of the previous has not got a run of gaps here */
			{

			
			if(fixed == TRUE)
				{
			
				if(silent == TRUE)
					{
					ratio[*count][2] = ratio[*count][2] + 1;
					graphs[*count][2][i]++;
			/*		position->graph[2][i] = position->graph[2][i] + 1;
			*/		}
				else 
					{
					ratio[*count][0] = ratio[*count][0] + 1;
					graphs[*count][0][i]++;
			/*		position->graph[0][i] = position->graph[0][i] + 1;
			*/		}
				}
			else
				{
				if(silent == TRUE)
					{
					ratio[*count][3] = ratio[*count][3] + 1;
					graphs[*count][3][i]++;
			/*		position->graph[3][i] = position->graph[3][i] + 1;
			*/		}
				else
					{
					ratio[*count][1] = ratio[*count][1] + 1;
					graphs[*count][1][i]++;
			/*		position->graph[1][i] = position->graph[1][i] + 1;
			*/		}
				}
			}
		
		}

	/* count polymorphisms */		
	if(position->seq_num1 != NULL)
		{
		
		if(codons[position->ances_seq[i]][j] != codons[(position->seq_num1)->bases[i]][j])
			{
						
			if((position->seq_num1)->gaprun < 2 && position->gaprun < 2) /* if the ancestor or the taxa is not in a run of gaps */
				{
				value = position->ances_seq[i] - transform_base(codons[position->ances_seq[i]][j], j) + transform_base(codons[(position->seq_num1)->bases[i]][j], j);
				if(value > 64) value = 64;  /* if there is a non standard charater tha whole codon is treated as a gap, and given the value of 64 to rpreent that */
				if(genetic_codes[code][value] != genetic_codes[code][position->ances_seq[i]])
					{
					ratio[*count][1] = ratio[*count][1] + 1;
					graphs[*count][1][i]++;
				/*	position->graph[1][i] = position->graph[1][i] + 1;
				*/	}
				else
					{
					ratio[*count][3] = ratio[*count][3] + 1;
					graphs[*count][3][i]++;
				/*	position->graph[3][i] = position->graph[3][i] + 1;
				*/	}
				}
			}
		}

	if(position->seq_num2 != NULL)
		{
		
		if(codons[position->ances_seq[i]][j] != codons[(position->seq_num2)->bases[i]][j])
			{
			
				
			if((position->seq_num2)->gaprun < 2 && position->gaprun < 2)  /* if the ancestor or the taxa is not part of a run of gaps */
				{
				value = position->ances_seq[i] - transform_base(codons[position->ances_seq[i]][j], j) + transform_base(codons[(position->seq_num2)->bases[i]][j], j);
				if(value > 64) value = 64;
				if(genetic_codes[code][value] != genetic_codes[code][position->ances_seq[i]])
					{
					ratio[*count][1] = ratio[*count][1] + 1;
					graphs[*count][1][i]++;
				/*	position->graph[1][i] = position->graph[1][i] + 1;
				*/	}
				else
					{
					ratio[*count][3] = ratio[*count][3] + 1;
					graphs[*count][3][i]++;
	/*				position->graph[3][i] = position->graph[3][i] + 1;
	*/				}
				}
			}
		}
	
		 

	if(position->node1 != NULL) snapshot(position->node1, ratio, nuc, count, i, j);
	if(position->node2 != NULL) snapshot(position->node2, ratio, nuc, count, i, j);
	
	}








/* This function will descend through the tree counting every change that occurs within the clade defined by the 
	internal node we first call the funtion with. during counting it also checks to see if the change wass synonymous
	or replacement, This function is not called if there is a fixed mutation within the clade as this would mean that 
	there are no other changes within the tree.  (20/9/00) */
	


void count_polymorphisms(struct node *position, int *repl, int *silent, int i, int j)
	{
	int value = 0;

	if(position->node1 != NULL) count_polymorphisms(position->node1, repl, silent, i, j);
	if(position->node2 != NULL) count_polymorphisms(position->node2, repl, silent, i ,j);
	
	if(position->node1 != NULL)
		{
		
		if(codons[position->ances_seq[i]][j] != codons[(position->node1)->ances_seq[i]][j])
			{
			value = position->ances_seq[i] - transform_base(codons[position->ances_seq[i]][j], j) + transform_base(codons[(position->node1)->ances_seq[i]][j], j);

			if(genetic_codes[code][value] != genetic_codes[code][position->ances_seq[i]])
				*repl = *repl +1;
			else
				*silent = *silent +1;
			}
		}
		
	if(position->node2 != NULL)
		{
		
		if(codons[position->ances_seq[i]][j] != codons[(position->node2)->ances_seq[i]][j])
			{
			value = position->ances_seq[i] - transform_base(codons[position->ances_seq[i]][j], j) + transform_base(codons[(position->node2)->ances_seq[i]][j], j);

			if(genetic_codes[code][value] != genetic_codes[code][position->ances_seq[i]])
				*repl = *repl +1;
			else
				*silent = *silent +1;
			}
		}
		
	if(position->seq_num1 != NULL)
		{
		
		if(codons[position->ances_seq[i]][j] != codons[(position->seq_num1)->bases[i]][j])
			{
			value = position->ances_seq[i] - transform_base(codons[position->ances_seq[i]][j], j) + transform_base(codons[(position->seq_num1)->bases[i]][j], j);

			if(genetic_codes[code][value] != genetic_codes[code][position->ances_seq[i]])
				*repl = *repl +1;
			else
				*silent = *silent +1;
			}
		}

	if(position->seq_num2 != NULL)
		{
		
		if(codons[position->ances_seq[i]][j] != codons[(position->seq_num2)->bases[i]][j])
			{
			value = position->ances_seq[i] - transform_base(codons[position->ances_seq[i]][j], j) + transform_base(codons[(position->seq_num2)->bases[i]][j], j);

			if(genetic_codes[code][value] != genetic_codes[code][position->ances_seq[i]])
				*repl = *repl +1;
			else
				*silent = *silent +1;
			}
		}
		
	}			


/* This function implements the new rules for determining whether a polymorphism is replacement or not,
	This travels down the tree from the internal branch specified, firstly looking for differences in any ancestral sequences
	it passes from the ancestor preceeding it. The algorithm increments the ratios (depending on synon or not) and then travels
	down the clade, checking any other ancestors it comes across
	and it also checks any taxa it comes across against the preceeding ancestor for differences and increments the ratios until
	the end of the tree is found.
	
	This method has an advantage in thata it will only count a mutation when it occurs in the tree, and not count for every occurance of the
	difference that decendants may contain. It also ensures that every difference is counted. 
	
*/

void find_polymorphisms(struct node *position, int **ratio, int **ratio1, int *count, int i, int j)
	{
	int value = 0;
	
	
	if(codons[position->ances_seq[i]][j] != codons[(position->prev)->ances_seq[i]][j])
		{
		value = (position->prev)->ances_seq[i] - transform_base(codons[(position->prev)->ances_seq[i]][j], j) + transform_base(codons[position->ances_seq[i]][j], j);
			
		if(genetic_codes[code][value] != genetic_codes[code][(position->prev)->ances_seq[i]])
			{
			ratio[*count][1] = ratio1[*count][1] + 1;
			}
		else
			{
			ratio[*count][3] = ratio1[*count][3] + 1;
			}
	
		}	
	
	if(position->seq_num1 != NULL)
		if(codons[(position->seq_num1)->bases[i]][j] != codons[position->ances_seq[i]][j])
		{
		value = position->ances_seq[i] - transform_base(codons[position->ances_seq[i]][j], j) + transform_base(codons[(position->seq_num1)->bases[i]][j], j);
			
		if(genetic_codes[code][value] != genetic_codes[code][position->ances_seq[i]])
			{
			ratio1[*count][1] = ratio[*count][1] + 1;
			ratio[*count][1] = ratio1[*count][1] + 1;
			}
		else
			{
			ratio1[*count][3] = ratio[*count][3] + 1;
			ratio[*count][3] = ratio1[*count][3] + 1;
			}
	
		}
	
	if(position->seq_num2 != NULL)
		if(codons[(position->seq_num2)->bases[i]][j] != codons[position->ances_seq[i]][j])
		{
		value = position->ances_seq[i] - transform_base(codons[position->ances_seq[i]][j], j) + transform_base(codons[(position->seq_num2)->bases[i]][j], j);
			
		if(genetic_codes[code][value] != genetic_codes[code][position->ances_seq[i]])
			{
			ratio1[*count][1] = ratio[*count][1] + 1;
			ratio[*count][1] = ratio1[*count][1] + 1;
			}
		else
			{
			ratio1[*count][3] = ratio[*count][3] + 1;
			ratio[*count][3] = ratio1[*count][3] + 1;
			}
	
		}
	
	if(position->node1 != NULL) find_polymorphisms(position->node1, ratio, ratio1, count, i, j);
	if(position->node2 != NULL) find_polymorphisms(position->node2, ratio, ratio1, count, i, j);
	

	}
	
	
	
	
	

/* This funcction implements a new rule for whether a mutation is fixed or not. The rule is as follows:
	If the path from the internal branch we are at to the root (not counting the root as we can't resolve
	an ambiguity at the root) doesn't contain an ancestor or *nucliotide* the same as the ancestor at the 
	internal branch we are at, then the mutation is fixed at that branch. - The mutation must be fixed in
	that clade also- */
	
void check_fixed_ances_out(struct node *position, int *fixed, char *nuc, int i, int j)
	{
	
	
	if(codons[position->ances_seq[i]][j] == *nuc) *fixed = FALSE;
	if(position->seq_num1 != NULL)
		if(codons[(position->seq_num1)->bases[i]][j] == *nuc) *fixed = FALSE;
	if(position->seq_num2 != NULL)
		if(codons[(position->seq_num2)->bases[i]][j] == *nuc) *fixed = FALSE;
		
	if(position->prev != NULL)
		{	
		if(*fixed) check_fixed_ances_out(position->prev, fixed, nuc, i, j);
		}
	}	
		
		




/* Called by is_fixed to travel down the tree from the node we are at to check if the nucliotide was fixed in that branch */

void check_fixed_in(struct node *position, int *fixed, char *nuc, int i, int j)
	{
	
	if(position->seq_num1 != NULL)
		if(*nuc != codons[(position->seq_num1)->bases[i]][j]) *fixed = FALSE;
	
	if(position->seq_num2 != NULL)
		if(*nuc != codons[(position->seq_num2)->bases[i]][j]) *fixed = FALSE;
			
		

	if(position->node1 != NULL && fixed) check_fixed_in(position->node1, fixed, nuc, i, j);
	if(position->node2 != NULL && fixed) check_fixed_in(position->node2, fixed, nuc, i, j);
	
	}
	
/* This checks to make sure that the mutation doesn't occur in any other position on the tree othe than that specified
	by place (and children) */
		
void check_fixed_out(struct node *position, struct node *place, int *fixed, char *nuc, int i, int j	)
	{
	
	if(position != place)
		{
		
		if(position->seq_num1 != NULL)
			if(*nuc == codons[(position->seq_num1)->bases[i]][j]) *fixed = FALSE;
			
		if(position->seq_num2 != NULL)
			if(*nuc == codons[(position->seq_num2)->bases[i]][j]) *fixed = FALSE;
			
		if(position->node1 != NULL && fixed) check_fixed_out(position->node1, place, fixed, nuc, i, j);
		if(position->node2 != NULL && fixed) check_fixed_out(position->node2, place, fixed, nuc, i, j);
		
		}
	}

/* This function travels through the tree, printing out the names of the sequences to check that the tree is being built properly */
void check_tree(struct node *position, int *count)
	{
	
	printf("(");
	if(position->node1 != NULL)
		{
		 check_tree(position->node1, count);
		 printf(",");
		}
	if(position->seq_num1 != NULL)
		{
		printf("%s,", position->seq_num1->name);
		}
	if(position->node2 != NULL) check_tree(position->node2, count);
	if(position->seq_num2 != NULL)
		{
		printf("%s", position->seq_num2->name);
		}
	printf(")");

	
	*count = *count+1;
	
/*	if(position->node1 != NULL) printf("node1 is assigned! %d\n", *count);
	if(position->seq_num1 != NULL)
		{
		printf("seqnum1\n%s %d\n", (position->seq_num1)->name, *count);
		printf("%s %d\n", (position->seq_num1)->nickname, *count);
		printf("sequence no:%d %d\n", (position->seq_num1)->seq_num, *count);
		printf("TAG:%d %d\n", (position->seq_num1)->tag, *count);
		printf("Outgroup: %d %d\n", (position->seq_num1)->outgroup, *count);
		printf("Length: %d %d\n", (position->seq_num1)->length, *count);
		printf("Numofstopcodons%d %d\n", (position->seq_num1)->numofstpcodons, *count);
		}
	if(position->node2 != NULL) printf("node2 is assigned! %d\n", *count);
	if(position->seq_num2 != NULL) 
		{
		printf("seqnum2\n%s %d\n", (position->seq_num2)->name, *count);
		printf("%s %d\n", (position->seq_num2)->nickname, *count);
		printf("sequence no:%d %d\n", (position->seq_num2)->seq_num, *count);
		printf("TAG:%d %d\n", (position->seq_num2)->tag, *count);
		printf("Outgroup: %d %d\n", (position->seq_num2)->outgroup, *count);
		printf("Length: %d %d\n", (position->seq_num2)->length, *count);
		printf("Numofstopcodons%d %d\n", (position->seq_num2)->numofstpcodons, *count);
		}
	if(position->seq_num1 == NULL && position->seq_num2 == NULL) printf("-\n");
*/
	
	}
	
	 			

/* function starts the checking process by calling assign_ances_up, and assign_ances_down which are both recursive and calls assign_codon_num, to calculate the ancestral sequence */
void ancestral_nuc(int i, int j, float *** subst_matrix)
	{
	struct node *position = NULL;
	
	position = tree_top;
	
	assign_ances_up(i, j, position);

	assign_root_ances(i,j);  /* this tries to solve ambiguities at  the root */

	assign_ances_down(position, subst_matrix);

	assign_codon_num(i, j, position);
	}
	
/* This function assigns the codon number to the ancestral sequence data in each node, depending on the position of the current nucliotide */
void assign_codon_num(int i, int j, struct node *position)
	{
	
	if(position->node1 != NULL) assign_codon_num(i, j, position->node1);
	if(position->node2 != NULL) assign_codon_num(i, j, position->node2);

	if( j == 0) position->ances_seq[i] = 0;
	position->ances_seq[i] += transform_base(position->ancestor[0], j);
	
/*	if(position == tree_top && i == 0) printf("ances tree_top j = %d - %d >> ances = %c\n", j, position->ances_seq[i], position->ancestor[0]);
	if(position == tree_top && i == 0) printf("ances next j = %d  - %d >> ances = %c\n", j, (position->node1)->ances_seq[i], (position->node1)->ancestor[0]);
*/	if(position->ances_seq[i] > 64) position->ances_seq[i] = 64;
	
	}



/* This function travels up through the tree assigning the ancestral nucliotides */	
void assign_ances_up(int i, int j, struct node *position)
	{
	char tmp[5] = {'\0','\0','\0','\0','\0'},tmp1[5] = {'\0','\0','\0','\0','\0'}, *pointer = NULL;
	int l = 0, k = 0;

	if(position->node1 != NULL) assign_ances_up(i, j, position->node1);
	if(position->node2 != NULL) assign_ances_up(i, j, position->node2);
	
	for(l=0; l<5; l++) position->ancestor[l] = '\0';
	l = 0;
	/*First assign those nucliotides which are at this node to the ancestors */
	if(position->seq_num1 != NULL)
		{
		position->ancestor[0] = codons[(position->seq_num1)->bases[i]][j];
		position->ancestor[1] = '\0';
		}
	if(position->seq_num2 != NULL && position->seq_num1 == NULL)
		{
		position->ancestor[0] = codons[(position->seq_num2)->bases[i]][j];
		position->ancestor[1] = '\0';
		}
	if(position->seq_num2 != NULL && position->seq_num1 != NULL)
		{
		if(codons[(position->seq_num1)->bases[i]][j] == codons[(position->seq_num2)->bases[i]][j])
			{	}
		else
			{
			position->ancestor[1] = codons[(position->seq_num2)->bases[i]][j];
			position->ancestor[2] = '\0';
			}
		}
		
	/* Second, check any children to see what the ancestor was assigned for those */	
		
	if(position->node1 != NULL) 
		{
		strcpy(tmp, (position->node1)->ancestor);
		while((pointer = strpbrk(tmp, position->ancestor)) != NULL)
			{
			pointer[0] = 'q';
			}
		if((pointer = strchr(tmp, 'q')) == NULL)
			{
			strcat(position->ancestor, tmp);
			}
		else
			{
			l = 0;
			for(k=0; k<5; k++)
				if(tmp[k] == 'q')
					{
					position->ancestor[l] = (position->node1)->ancestor[k];
					l++;
					}
			position->ancestor[l] = '\0';
			}
		}	

	if(position->node2 != NULL) 
		{
		strcpy(tmp, (position->node2)->ancestor);
		while((pointer = strpbrk(tmp, position->ancestor)) != NULL)
			{
			pointer[0] = 'q';
			}
		if((pointer = strchr(tmp, 'q')) == NULL)
			{
			strcat(position->ancestor, tmp);
			}
		else
			{
			l = 0;
			for(k=0; k<5; k++)
				if(tmp[k] == 'q')
					{
					position->ancestor[l] = (position->node2)->ancestor[k];
					l++;
					}
			position->ancestor[l] = '\0';
			}
		}	

	}

/* This function simply solves abiguities at the root, before assign_ances_down is called, so that there can never be an ambiguity in the tree. */
/* It solves ambiguities by assigning the ancesral root to that of the outgroup */
/* this however still leaves the chance for ambiguous sites, with an arbitrary descision being made if there is ambiguity a the top of the */
/* clade that defines the outgroup: For these purposes, you can only garantee unambiguous descisions if you only define ONE out group */
/* However this is not a great way of solving this, the other way would be to define the outgroup in the middle of the outgroup, so that the */
/* tree would still be rooted correctly, and the separation between the outgroup and the rest would be taken away from the root, where the ambiguities */
/* may lie. */
 

void assign_root_ances(int i, int j)
	{
	int x = 0, found = FALSE;
	char c = '\0';

	while(tree_top->ancestor[x] != '\0') x++;
		if(x > 1)
			{
			/* if the tree is rooted about a single sequence */
			if(tree_top->seq_num1 != NULL)  
				{
				if((tree_top->seq_num1)->outgroup == TRUE)
					{
					tree_top->ancestor[0] = codons[(tree_top->seq_num1)->bases[i]][j];
					found = TRUE;
					}
				}
			if(tree_top->seq_num2 != NULL)
				{
				if((tree_top->seq_num2)->outgroup == TRUE)
					{
					tree_top->ancestor[0] = codons[(tree_top->seq_num2)->bases[i]][j];
					found = TRUE;
					}
				}
			/* if the tree is rooted about an internal node */	
			if(!found)
				{
				if(tree_top->node1 != NULL)
					{
					if((look_for_ances(tree_top->node1)) != FALSE)
						{
						found = TRUE;
						tree_top->ancestor[0] = (tree_top->node1)->ancestor[0];
						} 
					}
				if(!found)
					if(tree_top->node2 != NULL)
						{
						if((look_for_ances(tree_top->node2)) != FALSE)
							{
							found = TRUE;
							tree_top->ancestor[0] = (tree_top->node2)->ancestor[0];
							}
						}
				}
			tree_top->ancestor[1] = '\0';
			}
	}
	
											
/* this is a recursive function which checks if a particular clade from the root is the outgroup*/
/* it returns TRUE if the clade is in the outgroup */
int look_for_ances(struct node *position)
	{
	char out = FALSE;
	
	if(position->seq_num1 != NULL)
			{
			if((position->seq_num1)->outgroup == TRUE)
				{
				out = TRUE;
				}
			}
	else
		{
		if(position->seq_num2 != NULL)
			{
			if((position->seq_num2)->outgroup == TRUE)
				{
				out = TRUE;
				}
			}
		else
			{
			if(position->node1 != NULL)
				out = look_for_ances(position->node1);
					/* if neither of the sequences are defined in this node, then both of nodes must be */
					/* we only need ths function to check one, since there has to be a sequence somewhere down that clade */
				  /* we don't need to check any more possible situations, since we only need to see 1 sequence */
			}
					
			
		}
	return(out);
	}




/* this function travels down the tree solving any ambiguities for ancestors, by checking previous branches
	or if necessary using a substitution matrix to solve it...... it doesn't solve for an ambiguity at the 
	root though      */
void assign_ances_down(struct node *position, float *** subst_matrix)
	{
	int x = 0, i = 0, j = 0;
	float total = 0;
	char tmp[5] = {'\0','\0','\0','\0','\0'}, *pointer = NULL, nuc = '\0';

	while(position->ancestor[x] != '\0') x++;
	
	if(x == 1)
		{
		/* First check the left part of the tree */
		if(position->node1 != NULL)
			{
			x = 0;
			while((position->node1)->ancestor[x] != '\0') x++;

			if(x > 1)
				{

				strcpy(tmp, (position->node1)->ancestor);
				while((pointer = strpbrk(tmp, position->ancestor)) != NULL)
					{
					pointer[0] = 'q';
					}
				if((pointer = strchr(tmp, 'q')) != NULL)
					{
					strcpy((position->node1)->ancestor, position->ancestor);
					}
				else   /* If there is an ambiguity which needs the substitution matrix */
					{

					x= 0;
					while(what[x] != position->ancestor[0]) x++;
					i = 0; nuc = '\0'; total = 0;
					while((position->node1)->ancestor[i] != '\0')
						{
						for(j=0; j<5; j++)
							{
							if((position->node1)->ancestor[i] == what[j])
								{
								if((subst_matrix[0][x][j] + subst_matrix[1][x][j] + subst_matrix[2][x][j]) > total)
									{
									total = (subst_matrix[0][x][j] + subst_matrix[1][x][j] + subst_matrix[2][x][j]);  /* this calculation doesn't need the  values split up into their ?-fold sites */
									nuc = what[j];
									}
								}
							}
						i++;
						}
					(position->node1)->ancestor[0] = nuc;
					(position->node1)->ancestor[1] = '\0';
					}
				}
			}
		/* Now check the right part of the tree */		
		if(position->node2 != NULL)
			{
			x = 0;
			while((position->node2)->ancestor[x] != '\0') x++;
			
			if(x > 1)
				{
				strcpy(tmp, (position->node2)->ancestor);
				while((pointer = strpbrk(tmp, position->ancestor)) != NULL)
					{
					pointer[0] = 'q';
					}
				if((pointer = strchr(tmp, 'q')) != NULL)
					{
					strcpy((position->node2)->ancestor, position->ancestor);
					}
				else   /* If there is an ambiguity which needs the substitution matrix */
					{
					x= 0;
					while(what[x] != position->ancestor[0]) x++;
					i = 0; nuc = '\0'; total = 0;
					while((position->node2)->ancestor[i] != '\0')
						{
						for(j=0; j<5; j++)
							{
							if((position->node2)->ancestor[i] == what[j])
								{
								if((subst_matrix[0][x][j] + subst_matrix[1][x][j] + subst_matrix[2][x][j]) > total)
									{
									total = (subst_matrix[0][x][j] + subst_matrix[1][x][j] + subst_matrix[2][x][j]);  /* this calculation doesn't need the  values split up into their ?-fold sites */
									nuc = what[j];
									}
								}
							}
						i++;
						}
					(position->node2)->ancestor[0] = nuc;
					(position->node2)->ancestor[1] = '\0';
					}
				}
			}

		}
	if(position->node1 != NULL) assign_ances_down(position->node1, subst_matrix);
	if(position->node2 != NULL) assign_ances_down(position->node2, subst_matrix);
	
	}
					
/* This function checks for nan values in the distances */					
int check_distances(void)
	{
	struct synon *position = li_wu_start;
	int i = 0, Kavalue = FALSE, Ksvalue = FALSE;
	float Kamax = 0, Ksmax = 0;
	
	
	/* travel through the results linked list and calculate all Ka and Ks */
	/* If any value is incalculable, then we set it to twice the maximum value anywhere else in the distance matrix */
	
	
	while(position)  /* Find all occurances of Nans in the distances, and mark them with a -1*/
		{
		for(i=0; i<((num_of_seqs - (position->seq_num + 1)) - untagged); i++)
			{
			
				if(position->Ka[i] < 100 && position->Ka[i] > -100 )
					{
					if(position->Ka[i] > Kamax) Kamax = position->Ka[i];
					}
				else
					{
					Kavalue = TRUE;
					position->Ka[i] = -1;
					}
					
			
			
				if(position->Ks[i] < 100 && position->Ks[i] > -100 )
					{
					if(position->Ks[i] > Ksmax) Ksmax = position->Ks[i];
					}
				else 
					{
					Ksvalue = TRUE;
					position->Ks[i] = -1;
					}		
				
			}
			
		position = position->next;
		}
	
	position = li_wu_start;
	if(Kavalue == TRUE || Ksvalue == TRUE)
		{
		while(position)  /* replace any -1 in the distance matrix with twice the max value from that matrix */
			{
			for(i=0; i<((num_of_seqs - (position->seq_num + 1)) - untagged); i++)
				{
				if(position->Ka[i] == -1) position->Ka[i] = 2*Kamax;
				if(position->Ks[i] == -1) position->Ks[i] = 2*Ksmax;
				}
			position = position->next;
			}
		}
				
	
	/* print results to the main output file */
	if(Kavalue == TRUE || Ksvalue == TRUE)
		{
		fprintf(outfile, "\n\nNaN values were detected in the");
		if(Kavalue == TRUE) fprintf(outfile, " Dn");
		if(Kavalue == TRUE && Ksvalue == TRUE) fprintf(outfile, " &");
		if(Ksvalue == TRUE) fprintf(outfile, " Ds");
		fprintf(outfile, " distance matrix, These values have been set to twice the largest distance in the distance matrix\n");
		}
	if(Kavalue == TRUE || Ksvalue == TRUE) Kavalue = TRUE;
	return(Kavalue);
	}
					
					



/* This function calculates looks at all the Ka and Ks values calculated and calculates all the possible Ka/Ks. 
   It then tallies the number of results that are greater, less than and equal to 1. These tallies are then printed
   to the main output file as part of the results.
*/


void tally_distances(void)
	{
	struct synon *position = li_wu_start;
	int greater = 0, less = 0, equal = 0, i = 0;
	
	/* travel through the results linked list and calculate all Ka/Ks */
	
	while(position)
		{
		for(i=0; i<((num_of_seqs - (position->seq_num + 1)) - untagged); i++)
			{
			if(position->Ks[i] == 0)
				{ 
				if(position->Ka[i] == 0) less++;
					/* 0/0 is assumed to be 0 since there has been no change */
				else
					{
					if(position->Ka[i] - position->varKa[i] > 0) greater++; /* If the variance doesn't span 0, so Ka is definitly greater than 0 */
					else less++;
					}
				}
			else
				{
				if(position->Ka[i]/position->Ks[i] > 1) greater++;
				if(position->Ka[i]/position->Ks[i] < 1) less++;
				if(position->Ka[i]/position->Ks[i] == 1) equal++;	
				}
			}
		position = position->next;
		}
	
	/* print results to the main output file */

	
	fprintf(outfile, "\n\n Totals for Dn/Ds:\n\tGreater than 1:\t%d\n\tLess than 1:\t%d\n\tEqual to 1: \t%d\n\n", greater, less, equal);
	
	}
					
				
					

/* this function takes the results produced by the Li Wu method and puts them into an array which can be read by the neighbour 
   joining tree function. This array is an num_of_seqsxnum_of_seqs array which is only filled in in the bottom half.
   The function first prints the three tables to the distance.out file and then it assigns the values to the distance array
*/



void allocate_distances(int k)
	{
	struct synon *position = li_wu_start;
	struct sequence *place = start;
	int i = 0, j = 0, last = FALSE, l = 0;
	float tmp = 0, largest = 2;

	/* assign distances array */
	
	
	distances = malloc((num_of_seqs-untagged)*sizeof(float*));
	if(!distances)
		{
		printf("\nError: Out of memory\n");
		clean_exit();
		}
	for(i=0; i<num_of_seqs-untagged; i++)
		{
		distances[i] = malloc((num_of_seqs-untagged)*sizeof(float));
		if(!distances[i])
			{
			printf("\nError: Out of memory\n");
			clean_exit();
			}
		}

	for(i=0; i<num_of_seqs-untagged; i++)
		for(j=0; j<num_of_seqs-untagged; j++)
			distances[i][j] = 0;

	

	
	for(l=0; l<4; l++)
		{
	
		 /* initialise the array tmat */
		for(i=0; i<num_of_seqs-untagged; i++)
			for(j=0; j<num_of_seqs-untagged; j++)
				distances[i][j] = 0;
		if(distance_written == TRUE) l = 4;
		if(l == 3)
			{
			last = TRUE;
			distance_written = TRUE;
			l = k;
			}
		

		/* assign the results of distances from li_wu to the array tmat */		
		j = 0;
		position = li_wu_start;
		while(position)
			{
			for(i=0; i<((num_of_seqs - (position->seq_num + 1)) - untagged); i++)
				{
				switch(l)
					{
					case 0:
						distances[j+i+1][j] = position->Ka[i];
				 		break;
				 	case 1:
				 		distances[j+i+1][j] = position->Ks[i];
				 		break;
				 	default:
				 		if(position->Ks[i] != 0) distances[j+i+1][j] = (position->Ka[i]/position->Ks[i]);
				 		else 
				 			{
				 			if(position->Ka[i] == 0) distances[j+i+1][j] = 0; /* if Ks/Ka = 0/0 then equal to 0 */
				 			else
				 				{	
				 				if(position->Ka[i] - position->varKa[i] > 0) distances[j+i+1][j] = -2; /* tag those results that are significantly larger than 0 */
								else distances[j+i+1][j] = 0;
								}
				 			}
				 		break;
				 	}
				 	
				}
			j++;
			position = position->next;
			}	
	

		largest = 2;
		for(i=0; i<num_of_seqs-untagged; i++)
			for(j=i+1; j<num_of_seqs-untagged; j++)
				{
				distances[i][j] = distances[j][i];
				if(distances[i][j] > largest) largest = distances[i][j];
				}
		
		/* the previous sections when they encountered a situation where Ka was some +ve value and Ks was 0, they 
		   checked the variance of the result and minused it from the value. If the result was greater than 0 then
		   the variance didn't include 0 and Ka was definitly not equal to 0. So this meant that Ka/Ks >> 1 and so 
		   the result was tagged with a -2. This next section then searches for all occurances of -2, and assign 
		   them to twice the value of the largest calculated before this. This gives a good representation of the 
		   huge Ka/Ks value that the result would represent.
		*/ 
		for(i=0; i<num_of_seqs-untagged; i++)
			for(j=0; j<num_of_seqs-untagged; j++)
				if(distances[i][j] == -2) distances[i][j] = largest*2;	
	

		if(!distance_written)
			{

			switch(l)
				{
				case 0:
					/* Ks values */
					if(dist == NULL)
						if((dist = fopen("Dn.dis", "w")) == NULL)		/* check to see if the file can be opened/created */
							printf("\n\n\tCannot open the output file, named Dn.dis,\n Distance matrix not written\n\n");				

					fprintf(dist,"%d\n", num_of_seqs - untagged);		
					break;
				case 1:
					/* Ka values */
					if(dist == NULL)
						if((dist = fopen("Ds.dis", "w")) == NULL)		/* check to see if the file can be opened/created */
							printf("\n\n\tCannot open the output file, named Ds.dis,\n Distance matrix not written\n\n");				

					fprintf(dist,"%d\n", num_of_seqs - untagged);		
					break;
				default:
					/* Ka/ks values */
					if(dist == NULL)
						if((dist = fopen("DnDs.dis", "w")) == NULL)		/* check to see if the file can be opened/created */
							printf("\n\n\tCannot open the output file, named DnDs.dis,\n Distance matrix not written\n\n");				

					fprintf(dist,"%d\n", num_of_seqs - untagged);		
					break;	
				}	

			place = start;
			for(i=0; i<num_of_seqs-untagged; i++)
				{
				fprintf(dist, "%-10.10s", place->name);
				for(j=0; j<num_of_seqs-untagged; j++)
					{
					if(distances[i][j] >= 0 && distances[i][j] < 100)
						fprintf(dist, " %f", distances[i][j]);
					else 
						{
						if(distances[i][j] < 0) 
							{
							distances[i][j] = 0;
							fprintf(dist, " %f", distances[i][j]);
							}
						else
							fprintf(dist, " NaN  ");
						}
					}
 				fprintf(dist, "\n");
				place = place->next;
				}
			fclose(dist);
			dist = NULL;
			}
		if(last) l = 4;
		}

	

	/* free up the memory used by li_wu_results */

		while(li_wu_start)
			{
			position = li_wu_start->next;
			free(li_wu_start->Ks);
			free(li_wu_start->varKs);
			free(li_wu_start->Ka);
			free(li_wu_start->varKa);
			free(li_wu_start);
			li_wu_start = position;
			}

	printf("\n\nDistances written to Dn.dis, Ds.dis & DnDs.dis\n");
	}





void McDonald_Kreitman(void)
/* 
   Calculate a tree using the distances in the num_of_seqs*num_of_seqs array d_mat.
*/
{	int i=0, j=0, k = 0, method = 0, exit = FALSE, num = 0, error = FALSE, number = 0, count =0, title = TRUE, ans = FALSE, bracket = FALSE;
	char choice = '\0', overflow = '\0', *pointer = NULL, last = '\0', c = '\0';
	float **pvalue = NULL, **pvalue2 = NULL, *gChi = NULL, *gChi2 = NULL, **subs = NULL, **ratio2 = NULL, **tmp_ratio = NULL;
	float ** standard_tree = NULL;
	struct sequence *position = start;
	struct node *previous = NULL;
	int  **ratio = NULL; /* used to record the ratio calculated by the McDonald-Kreitman method. */
	
	FILE *nested = NULL, *subsfile = NULL; /* used for reading in a tree in nest parentheses format if desired */



	distance_written = FALSE;
	if(num_of_seqs<4) 
		{
		printf("\nAlignment has only %d sequences\n",num_of_seqs);
		return;
		}
			
	/* assign standard_tree */
	standard_tree   = malloc((num_of_seqs-untagged)*sizeof(float*));
	if(!standard_tree)
		{
		printf("ERROR: Out of memory\n");
		clean_exit();
		}

	for(i=0; i<num_of_seqs-untagged; i++)
		{	
		standard_tree[i]  = malloc((num_of_seqs-untagged)*sizeof(float));
		if(!standard_tree[i])
			{
			printf("ERROR: Out of memory\n");
			clean_exit();
			}
		}
	
	for(i=0; i<num_of_seqs-untagged; i++)
		for(j=0; j<num_of_seqs - untagged; j++)
			standard_tree[i][j] = 0;		

	/* assign graphs */
	graphs = malloc(((num_of_seqs - untagged)*sizeof(int **)));
	if(!graphs)
		{
		printf("ERROR: Out of memory\n");
		clean_exit();
		}
	for(i=0; i<num_of_seqs - untagged; i++)
		{
		graphs[i] = malloc(4*sizeof(int *));
		if(!graphs[i])
			{
			printf("ERROR: Out of memory\n");
			clean_exit();
			}
		for(j=0; j<4; j++)
			{
			graphs[i][j] = malloc((start->length/3)*sizeof(int));
			if(!graphs[i][j])
				{
				printf("ERROR: Out of memory\n");
				clean_exit();
				}
			}
		}
	
	for(i=0; i<num_of_seqs - untagged; i++)
		for(j=0; j<4; j++)
			for(k=0; k<start->length/3; k++)
				graphs[i][j][k] = 0;
	
		
	/* assign ratio;  the array is divided up as follows: */
	ratio = malloc(((num_of_seqs - untagged )*sizeof(int *)));
	if(!ratio)
		{
		printf("ERROR: Out of memory\n");
		clean_exit();
		}
	ratio2 = malloc(((num_of_seqs - untagged )*sizeof(float *)));
	if(!ratio2)
		{
		printf("ERROR: Out of memory\n");
		clean_exit();
		}	
	tmp_ratio = malloc(((num_of_seqs - untagged )*sizeof(float *))); /* the input for the g-test needed to be changed to a float, so instead of changing every mention of ratio
																		to refer to a float, this array will be assigned to the values that are in ratio, and will be passed in its
																		stead. This leaves the original as it is, and meets the requirements for the function G_test */ 
	if(!tmp_ratio)
		{
		printf("ERROR: Out of memory\n");
		clean_exit();
		}
	for(i=0; i<(num_of_seqs - untagged ); i++)
		{
		ratio[i] = malloc(4*sizeof(int));
		if(!ratio[i])
			{
			printf("ERROR: Out of memory\n");
			clean_exit();
			}
		ratio2[i] = malloc(4*sizeof(float));
		if(!ratio2[i])
			{
			printf("ERROR: Out of memory\n");
			clean_exit();
			}
		tmp_ratio[i] = malloc(4*sizeof(float));
		if(!tmp_ratio[i])
			{
			printf("ERROR: Out of memory\n");
			clean_exit();
			}
		}
	for(i=0; i<(num_of_seqs - untagged ); i++)
		for(j=0; j<4; j++)
			{
			ratio[i][j] = 0;	
			ratio2[i][j] = 0;
			tmp_ratio[i][j] = 0;
			}

	/* Assign gChi: this is used to stroe the chi vaules that are used to calculate the pvalues at one degree of freedom, it is being 
		remembered so that if can be printed out for simulation testing (we can then look over 1000 runs of neutral data, the distribution of chi)
		*/
	gChi = malloc(((num_of_seqs - untagged )*sizeof(float)));
	if(!gChi)
		{
		printf("ERROR: Out of memory\n");
		clean_exit();
		}		
	for(j=0; j<num_of_seqs - untagged; j++)
		{
		gChi[j] = 0;
		}
		
	gChi2 = malloc(((num_of_seqs - untagged )*sizeof(float)));
	if(!gChi2)
		{
		printf("ERROR: Out of memory\n");
		clean_exit();
		}		
	for(j=0; j<num_of_seqs - untagged; j++)
		{
		gChi2[j] = 0;
		}	
	
			
	/* assign pvalue */
	
	pvalue = malloc(((num_of_seqs - untagged )*sizeof(float *)));
	if(!pvalue)
		{
		printf("ERROR: Out of memory\n");
		clean_exit();
		}
	pvalue2 = malloc(((num_of_seqs - untagged )*sizeof(float *)));
	if(!pvalue2)
		{
		printf("ERROR: Out of memory\n");
		clean_exit();
		}

	for(i=0; i<(num_of_seqs - untagged ); i++)
		{
		pvalue[i] = malloc((3*sizeof(float)));
		if(!pvalue[i])
			{
			printf("ERROR: Out of memory\n");
			clean_exit();
			}
		pvalue2[i] = malloc((3*sizeof(float)));
		if(!pvalue2[i])
			{
			printf("ERROR: Out of memory\n");
			clean_exit();
			}
		}

	for(i=0; i<(num_of_seqs - untagged ); i++)
		for(j=0; j<3; j++)
			{
			pvalue[i][j] = 0;     /* pvalue[x][0] will be maxP, pvalue[x][1] will be minP */
			pvalue2[i][j] = 0;    /* this is to be used to check if the number of silent and replacement substitutions deviates from what is expected by the neutral rate */
			}
	
	/* Assign subs array (to be used in the function subs_inclade) */
	subs = malloc(((num_of_seqs - untagged)*sizeof(float *)));
	if(!subs)
		{
		printf("ERROR: out of memory\n");
		clean_exit();
		}
	
	for(i=0; i<(num_of_seqs - untagged); i++)
		{
		subs[i] = malloc((2*sizeof(float)));
		if(!subs[i])
			{
			printf("ERROR: out of memory\n");
			clean_exit();
			}
		}
	for(i=0; i<(num_of_seqs - untagged); i++)
		{
		subs[i][0] = 0;  
		subs[i][1] = 0;
		}

		
	/* create tree array */
	fprintf(outfile,"--------------------------------------------------------\n");
	fprintf(outfile, "\nResults from input file: %s\n", filename);
	if(gen_opt[5] == 0)
		{
		if(!distance_written) tally_distances();
		
		allocate_distances(gen_opt[6]); 					/* allocate tmat */
				
		
		n_joining_tree(standard_tree);				/* calculate the phylogenetic tree */
		
		linked_tree(standard_tree); /* This creates the binary linked tree that represents the phylogentic tree */
		ans = val_tree(tree_top);
		if(ans == TRUE) printf("HAD TO FIX TREE, PLEASE CHECK THAT THE TREE CONTAINS ALL TAXA\n");
		
/*		printf("After Linked_tree\n");
		number = 0;
		if(tree_top != NULL)
			{
			check_tree(tree_top, &number);  *//* make sure the tree has been created correctly */	
/*			}
		else printf("tree not defined\n");
		printf("press return to continue - adaptive tree");
		getchar();

*/	

		}





	/********print results to output file*******/	
	fprintf(outfile, "\nReults of ");
	if(gen_opt[5] == 0) fprintf(outfile, "neighbour joining and ");
	fprintf(outfile, "Creevey, McInerney method for the file named %s\n", filename);

	if(gen_opt[5] == 0)
		{
		fprintf(outfile,"Distances used:  ");
		switch(gen_opt[6])
			{
			case 1:
				fprintf(outfile, "Ds values");		
				break;
			case 0:
				fprintf(outfile, "Dn values");
				break;
			default:
				fprintf(outfile, "Rate of evolution (Dn/Ds)");
				break;
			}

		
		
		/*
		Print out the distance values that built the tree in a tree array
		position = start;
		while(position && !position->tag) position = position->next;
	
		fprintf(outfile, "\t\t                                       ");
		for(i=0; i<num_of_seqs-untagged; i++) fprintf(outfile, "\t%d ", (i+1));
		fprintf(outfile, "\n" );
		for(i=0; i<num_of_seqs-untagged; i++)
			{
			fprintf(outfile, "%d %-41.41s", i+1, position->name);
			for(j=0; j<num_of_seqs-untagged; j++)
				fprintf(outfile, " \t%f", standard_tree[i][j]);
			fprintf(outfile, "\n");
			position = position->next;
			while(position && !position->tag) position = position->next;

			}
		*/

		/* print out the tree using + and - in a tree array */

		position = start;
		while(position && !position->tag) position = position->next;

		fprintf(outfile, "\n\n\n\n\t\t                                       ");
		for(i=0; i<num_of_seqs-untagged; i++) fprintf(outfile, "\t%d ", i+1);
		fprintf(outfile, "\n" );
		for(i=0; i<num_of_seqs-untagged; i++)
			{
			fprintf(outfile, "%d %-41.41s", i+1, position->name);
			for(j=0; j<num_of_seqs-untagged; j++)
				{
				if(standard_tree[i][j] == 0) fprintf(outfile, " \t-");
				else fprintf(outfile, " \t+");
				}
			fprintf(outfile, "\n");
			position = position->next;
			while(position && !position->tag) position = position->next;	

			}
		fprintf(outfile, "\n\n\n\n");
		}

	if(!error)
		{
/*		number = 0;
		if(tree_top != NULL)
			{
			printf("\n\n\n");
			check_tree(tree_top, &number);  *//* make sure the tree has been created correctly */
/*			printf("press return to continue - after input_treefile\n");
			getchar();	
			}
		else printf("tree not defined\n");
*/	
		adaptive_tree(ratio);/*perform McDonald and Kreitmans method */
		
		/*  calculate the pvalue of the ratio using the G test & Fishers Exact test */
		/* first assign everything in the array tmp_ratio to whats in ratio. Tmp_ratio is in float while ratio is in int, and so we can pass the right type the g_test() without changing the array ratio */
		
		for(i=0; i<num_of_seqs - untagged - 2; i++)
			for(j=0; j<4; j++)
				{
				tmp_ratio[i][j] = ratio[i][j];
				}
		g_test(tmp_ratio, pvalue, gChi);



/*********************** New SECTION */

		/* next count the number of replacements or silent sites in each clade defined gy each internal branch */
		count = 0;
		subs_inclade(tree_top, &count, subs); 
		
		
		/* Subs[0] contains the number of replacement sites subs[1] contains the number of silent sites */
		/* I'm going to use the 2d array called ratio2 to  hold the expected and observed to be tested using the gtest */
		/* The observed are the totals of replacemetns and silents as calculated by the Md&K method and will be stored in ratio2[i][0] and ratio2[i][1] respectively */
		/* The expected replacements and silents will be calculated using the numbers stored in the subs array and stored in ratio2[i][2] and ratio[i][3] respectively*/
		/* The pvalue for this result will be stored in the array pvaule2 */
		
		for(i=0; i<(num_of_seqs - untagged); i++)
			{
			 /* Fill in the observed values */
			 ratio2[i][0] = ratio[i][0] + ratio[i][1];   /* total up all the replacements in this calde */
			 ratio2[i][1] = ratio[i][2] + ratio[i][3];   /* total up all the silents in this clade */
			 
			 /* fill in the expected values */
			if(subs[i][0] + subs[i][1] != 0)  /* so we don't divide by 0 by accident */
				 ratio2[i][2] = (subs[i][0]/(subs[i][0] + subs[i][1])) * (ratio2[i][0] + ratio2[i][1]);    /* calculate the expected replacements ---   */
				/* this equation is calculaing this by basing it on the ratio of repl/silent as from the expected but reducing it to the number of substitutions in the observed */
			else
				ratio2[i][2] = 0;
				
			 ratio2[i][3] = (ratio2[i][0] + ratio2[i][1]) - ratio2[i][2];    /* calculate the expected silents, by minusing the number of replacements caluculated above from the total number of replacemetns in the observed */
			}
			
		/* we can now put these values into the g-test to test for independance */
		g_test(ratio2, pvalue2, gChi2);
		
			
			
		/* Print out the results ......... This is to a new file temporarily called substitutions.out  */
		bracket = FALSE;
		
		if((subsfile = fopen("substitutions.out", "w")) == NULL)		/* check to see if the file can be opened/created */
			{
			printf("\n\n\tCannot open the  file, substitutions.out\n");
			}
		else
			{
			fprintf(subsfile,"The calculated values of replacement and silent substitutions at each internal branch and associated pvalues \n\n");
			fprintf(subsfile,"Node\tObserved\t\tExpected\t\tTotal sites\t\tChanges per site\t\tpvalue (Observed vs expected)\n");
			fprintf(subsfile,"\tReplacement\tSilent\tReplacement\tSilent\tReplacement\tSilent\tReplacement\tSilent\n");
			for(i=0; i<(num_of_seqs - untagged) - 2; i++)
				{
				fprintf(subsfile,"%d\t%-10.0f\t%-10.0f\t%-10.3f\t%-10.3f\t%-10.1f\t%-10.1f\t%-10f\t%-10f", i, ratio2[i][0], ratio2[i][1], ratio2[i][2], ratio2[i][3], subs[i][0], subs[i][1], (ratio2[i][0]/subs[i][0]), (ratio2[i][1]/subs[i][1]));

				fprintf(subsfile, "\tG = %f ", gChi2[i]);
				
				if(pvalue2[i][1] < 0.05 && pvalue2[i][2] > 0.05 && pvalue2[i][2] < 1)
					{
					fprintf(subsfile,"\t(Gtest:%.3f > pvalue > %.3f)", pvalue2[i][0], pvalue2[i][1]);
					bracket = TRUE;
					}
				else
					fprintf(subsfile,"\t Gtest:%.3f > pvalue > %.3f ", pvalue2[i][0], pvalue2[i][1]);

				if(pvalue2[i][2] > 0 && pvalue2[i][2] < 1)
					fprintf(subsfile,"\tfishers: p = %f", pvalue2[i][2]);
				else
					fprintf(subsfile,"\tfishers: incalculable");
					
				if(pvalue2[i][1] < 0.05 || (pvalue2[i][2] <= 0.05 && pvalue2[i][2] > 0)) fprintf(subsfile," *\n");
				else fprintf(subsfile,"\n");

				}
			
			if(bracket) fprintf(subsfile, "\n\n\n Legend:\n\tThose significant Gtest results that are in brackets are statistically unreliable due to the size of the results,\n\tIn these cases the result calculated from the Fishers exact should be taken as the 'true' p value"); 
			}
		fflush(subsfile);
		fclose(subsfile);

/********************* end of NEW SECTION */


		/* Print out the graph file */
		
		if((graphfile = fopen("seqgraph.out", "w")) == NULL)		/* check to see if the file can be opened/created */
			{
			printf("\n\n\tCannot open the  file, seqgraph.out\n");
			}
		else
			{
			for(i=0; i<num_of_seqs - untagged - 2; i++)
				{
					fprintf(graphfile, "Node %d\t", i);
					fprintf(graphfile, "repl/invar\t");
					fprintf(graphfile, "repl/var\t");
					fprintf(graphfile, "synon/invar\t");
					fprintf(graphfile, "synon/var\t");
				}
			fprintf(graphfile, "\n");
			
				
				for(k=0; k<start->length/3; k++)
					{
					for(i=0; i<num_of_seqs - untagged - 2; i++)
						{
						for(j=0; j<4; j++)
							{
						 	fprintf(graphfile, "\t%d", graphs[i][j][k]);
							}
						fprintf(graphfile, "\t");
						}
					fprintf(graphfile, "\n");
					}
					
		
			fclose(graphfile);
			}
		 tree_pairwise_distances(tree_top); /* calculate the pairwise distances for the tree */

		/* go through the sequences and get rid of any parentheses in the names, replace with "-"  */
		
		position = start;
		while(position)
			{
			while((pointer = strchr(position->nickname, '(')) != NULL) pointer[0] = '_';
			while((pointer = strchr(position->nickname, ')')) != NULL) pointer[0] = '_';
			while((pointer = strchr(position->name, '(')) != NULL) pointer[0] = '_';
			while((pointer = strchr(position->name, ')')) != NULL) pointer[0] = '_';
			while((pointer = strchr(position->name, '-')) != NULL) pointer[0] = '_';
			while((pointer = strchr(position->name, '-')) != NULL) pointer[0] = '_';
			while((pointer = strchr(position->nickname, '-')) != NULL) pointer[0] = '_';
			while((pointer = strchr(position->nickname, '-')) != NULL) pointer[0] = '_';
			while((pointer = strchr(position->nickname, ',')) != NULL) pointer[0] = '_';
			while((pointer = strchr(position->name, ',')) != NULL) pointer[0] = '_';
			while((pointer = strchr(position->nickname, ':')) != NULL) pointer[0] = '_';
			while((pointer = strchr(position->name, ':')) != NULL) pointer[0] = '_';

			position = position->next;
			}		

		fflush(outfile);
		/* Write tree file to output */
		
	
		if(outtree == NULL)
			if((outtree = fopen("result-tree.ph", "w")) == NULL)		/* check to see if the file can be opened/created */
				printf("\n\n\tCannot open the output file, named result-tree.ph,\n Tree file not written\n\n");	
		fprintf(outtree,"[Results from input file: %s\nCreevey, McInerney method: internal nodes labeled indicate where the ratios differed significantly]\n", filename);
		write_tree1(tree_top, &last, &num, pvalue);
		fprintf(outtree, ";\n");
		last = '\0'; num = 0;
		
		yadf = fopen("yadf.out", "w");
		assign_node_nums(tree_top, 0);
		fprintf(yadf, "Which branch of the tree?\tDn\tDs\tDn/Ds\tvar Dn\tvar Ds\n");
		print_tree_pair_dist(tree_top);
		fclose(yadf);

	/*  This has been commented out since the file yadf.out describes much better what the second tree used to....

		fprintf(outtree,"\n[Dn/Ds ratios between neighbouring nodes, labels are where Ka/Ks > 1. for either the above branch (1) or the lower branch (2)]\n");
		write_tree2(tree_top, &last, &num);
		fprintf(outtree, ";\n");
	*/	

		last = '\0'; num = 0;
		fprintf(outtree,"\n[This tree has every node labeled just for reference, there are no results in the tree]\n");
		write_tree3(tree_top, &last, &num, pvalue);
		fprintf(outtree, ";\n\n");
		  
		
		bracket = FALSE;	
		
		
					
		/* Print the results of the Mc Donald & Kreitman test */
		fprintf(outfile, "\n\n\nResults of Creevey, McInerney Test for each internal branch\n\n");
		fprintf(outfile, "no\trepl:invar  \trepl:var \tsynon:invar \tsynon:var\n");
		for(i=0; i<(num_of_seqs - untagged - 2); i++)
			{
			fprintf(outfile, "branch %d", i);
			for(j=0; j<4; j++) 
				{
				fprintf(outfile, "\t%-10d ", ratio[i][j]);
				}
				
			fprintf(outfile, "\tG = %f ", gChi[i]); /* printf out the value of chi for this branch */
			
			if(pvalue[i][2] > 0.05 && pvalue[i][2] < 100 && pvalue[i][1] < 0.05) fprintf(outfile, "\t(");
			else fprintf(outfile, "\t ");
			
			fprintf(outfile, "Gtest:%f > pvalue > %f", pvalue[i][0], pvalue[i][1]);
			
			if(pvalue[i][2] > 0.05 && pvalue[i][2] < 100 && pvalue[i][1] < 0.05)
				{
				fprintf(outfile, ")  ");
				bracket = TRUE;
				}
			else fprintf(outfile, "   ");
			
			if(pvalue[i][2] > 0 && pvalue[i][2] < 100) fprintf(outfile, "\tFishers: p = %f ", pvalue[i][2]);
			else fprintf(outfile, "\tFishers: incalculable ");
									 
			
			if(pvalue[i][1] < 0.050 ||( pvalue[i][2] <= 0.05 && pvalue[i][2] != 0))
				{
				if(ratio[i][2] != 0 && ratio[i][3] != 0)
					{ 
					if((float)ratio[i][0]/(float)ratio[i][2] > (float)ratio[i][1]/(float)ratio[i][3]) fprintf(outfile, "*\n");
					else fprintf(outfile, "%%\n");
					}
				else
					{
					if(ratio[i][2] == 0 && ratio[i][3] != 0)
						{
						if((float)ratio[i][0] > (float)ratio[i][1]/(float)ratio[i][3]) fprintf(outfile, "*\n");
						else fprintf(outfile, "%%\n");
						}
					if(ratio[i][2] != 0 && ratio[i][3] == 0)
						{
						if((float)ratio[i][0]/(float)ratio[i][2] > (float)ratio[i][1]) fprintf(outfile, "*\n");
						else fprintf(outfile, "%%\n");
						}
					if(ratio[i][2] == 0 && ratio[i][3] == 0)
						{
						if((float)ratio[i][0] > (float)ratio[i][1]) fprintf(outfile, "*\n");
						else fprintf(outfile, "%%\n");
						}
					}
				}
			else fprintf(outfile, "\n");
			
	
			}
		if(bracket) fprintf(outfile, "\n\n\n Legend:\n\tThose significant Gtest results that are in brackets are statistically unreliable due to the size of the results,\n\tIn these cases the result calculated from the Fishers exact should be taken as the 'true' p value"); 

		fflush(outfile);

		}
	else printf("Error in treefile, the tree is not properly formatted\n");

	/* free up allocated memory */
	count = 0;
	dismantle(tree_top, &count);
	tree_top = NULL;
	
	if(pvalue != NULL) 
		{
		for(i=0; i<num_of_seqs-untagged; i++) free(pvalue[i]);
		free(pvalue);
		}
	
	if(standard_tree != NULL) 
		{
		for(i=0; i<num_of_seqs-untagged; i++) free(standard_tree[i]);
		free(standard_tree);
		}
	
		if(ratio != NULL) 
		{
		for(i=0; i<num_of_seqs-untagged; i++) free(ratio[i]);
		free(ratio);
		}
	
	if(pvalue2 != NULL) 
		{
		for(i=0; i<num_of_seqs-untagged; i++) free(pvalue2[i]);
		free(pvalue2);
		}
	
	if(ratio2 != NULL) 
		{
		for(i=0; i<num_of_seqs-untagged; i++) free(ratio2[i]);
		free(ratio2);
		}
	
	if(subs != NULL)
		{
		for(i=0; i<num_of_seqs-untagged; i++) free(subs[i]);
		free(subs);
		}
	
	if(tmp_ratio != NULL)
		{
		for(i=0; i<num_of_seqs - untagged; i++) free(tmp_ratio[i]);
		free(tmp_ratio);
		}

}


int tree_choice(void)
	{
	int error = FALSE, exit = FALSE, method = 1, number = 0;
	char c = '\0', overflow = '\0', ans = '\0';
	struct node *previous = NULL;
	
	
	do{
		do{
			error = FALSE;
			printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
			printf("\n\n\tPhylogenetic tree:\n\t\tWould you like to:\n\t\t1 = read in a tree file in nested parentheses format?\n\t\t2 = create a tree using the neighbour joining algorithm?\n\n");
			printf("\tPlease choose 1 or 2");
			ans = getletter("\tPlease choose 1 or 2");
			
			switch(ans)
				{
				case '1':
					method = 1;
					break;
				case '2':
					method = 2;
					break;
				default:
					ans = '\0';
					break;
				}
			}while(ans == '\0');
	

		if(method == 1)
			{
			if(parenthesis == NULL)
				{
				do{
					printf("\n\n\tName of the tree file in nested parentheses format: ");
					nestname[0] = '\0'; 
					scanf("%s%c", string1, string2);
					strcpy(nestname, string1);
					if((parenthesis = fopen(nestname, "r")) == NULL)		/* check to see if the file is there */
						{
						printf("\n\n\tCannot open the tree file, named %s\n\n", filename);
						method = getint("\n\n\t\nPress 1 to try again\nOr press 2 to build a neighbour joining tree", 1, 2, 1);
						if(method == 2) exit = TRUE;
						}

					}while(parenthesis == NULL && !exit);
				}
			if(method == 1)	
				{
				error = check_treefile();
				c = getc(parenthesis);
				if(!error) error = input_tree(&c, previous, 0);
				if(error) printf("\nError in the tree file, the tree is not formatted properly\n"); 
				fclose(parenthesis);
				parenthesis = NULL;	
				
			/*	number = 0;
				if(tree_top != NULL)
					{
					check_tree(tree_top, &number);  *//* make sure the tree has been created correctly */
			/*		printf("press return to continue (Just inputted the tree from the file)");
					getchar();	
					}
				else printf("tree not defined\n");
			*/	
				prune_tree(tree_top);
			/*	number = 0;
				if(tree_top != NULL)
					{
					check_tree(tree_top, &number);  *//* make sure the tree has been created correctly */
			/*		printf("press return to continue (Just inputted the tree from the file, after pruning)");
					getchar();	
					}
				else printf("tree not defined\n");

			*/				
				}
			}
		}while(error == TRUE && method == 1);


	return(method);
	
	}



 /*	This function applies the neighbour joining tree method described by Saitou and Nei in 1987, and described in Molecular 
	Systematics 2nd Ed. by David Hillis.
 	This function needs passed to it a double pointer <passed by reference>  This will contain the results when finished
 	 It needs access to a global array distances calculated earlier <an Num_of_seqsXNum_of_seqs array 
	symetrical about the diagonal. The diagonal contains 0's.> The double pointer will have the tree file written into 
	it at the end of this function <this array is allocated memory in this function>.     
 */

/*new*/
void n_joining_tree(float **tree)
	{
	int i = 0, j = 0, smli = 0, smlj = 0, x = 0, y = 0, k =0, nodes = 0, **seqs = NULL, **tmp2 = NULL;
	float tot = 0, sml_dist = 0, assigned = 0, **d_mat = NULL, **tmp = NULL, second = 0, **net_div = NULL, total_div = 0, dr = 0, d2r = 0, to = 0;
	
	
	nodes = num_of_seqs-untagged; /* nodes is the number of terminal nodes */
	
	net_div = malloc((nodes)*(sizeof(float*)));
	if(net_div == NULL) 
		{
		printf("Out of memory\n");
		clean_exit();
		}
	for(i=0; i<nodes; i++) 
		{
		net_div[i] = malloc((2)*(sizeof(float)));
		if(net_div[i] == NULL) 
			{
			printf("Out of memory\n");
			clean_exit();
			}
		}
	
	/* d_mat is the distance matrix passed to the function */ 
	d_mat = malloc(nodes*(sizeof(float *)));
	if(d_mat == NULL) 
		{
		printf("Out of memory\n");
		clean_exit();
		}
	for(i=0; i<nodes; i++)	
		{
		d_mat[i] = malloc((nodes)*(sizeof(float)));
		if(d_mat[i] == NULL) 
			{
			printf("Out of memory\n");
			clean_exit();
			}
		}
	
	for(i=0; i<nodes; i++)
		for(j=0; j<nodes; j++)
			d_mat[i][j] = distances[i][j];    /* assign the array equal to the distance array paeed to the fuction */
	
	seqs = malloc(nodes*(sizeof(int*)));
	if(seqs == NULL) 
		{
		printf("Out of memory\n");
		clean_exit();
		}	
	for(i=0; i<nodes; i++) 
		{
		seqs[i] = malloc(2*(sizeof(int)));   /* Use the first row of the array to signify which sequence is in each column of the matrices */ 	
		if(seqs[i] == NULL) 
			{
			printf("Out of memory\n");
			clean_exit();
			}
		}

	for(i=0; i<nodes; i++)
		{
		seqs[i][0] = i;
		seqs[i][1] = FALSE;   /* stores whether or not the column is a node */
		}

	/* tmp is a matrix used during calculation */
	tmp = malloc((nodes-1)*(sizeof(float *)));

	if(tmp == NULL) 
		{
		printf("Out of memory\n");
		clean_exit();
		}

	for(i=0; i<(nodes-1); i++)
		{
		tmp[i] = malloc((nodes-1)*(sizeof(float)));
		if(tmp[i] == NULL) 
			{
			printf("Out of memory\n");
			clean_exit();
			}
		}	
/* tmp2 is used to transfer information between iterations */
		tmp2 = malloc((nodes-1)*sizeof(int*));
		if(tmp2 == NULL) 
			{
			printf("Out of memory\n");
			clean_exit();
			}
		for(i=0; i<nodes-1; i++) 
			{
			tmp2[i] = malloc(2*(sizeof(int)));
			if(tmp2[i] == NULL) 
				{
				printf("Out of memory\n");
				clean_exit();
				}
			}


	
		/* 1: Given a matrix of pairwise distances (t_mat), for each terminal node i calculate 
			  its net divergence (net_div[i][x]) from all other taxa, by summing them */
		  





	while(nodes > 2)  
		{
		

		for(i=0; i<(num_of_seqs-untagged -1); i++)
			for(j=0; j<(num_of_seqs-untagged -1); j++)
				tmp[i][j] = 0;  /* initialise arrays */
		

		for(i=0; i<nodes; i++)
			{
			net_div[i][0] = 0; /* this will hold the net divergence */
			net_div[i][1] = 0; /* this will hold the net divergence divided by n-2 for later calculations */
			for(j=0; j<nodes; j++)
				{
				net_div[i][0] += d_mat[i][j];
				net_div[i][1] = net_div[i][0]/(nodes - 2);
				}
			}
		
		
		/* 2: calculate a rate-corrected distance in which the elements are defined by:
				 d_matij - (net_divi - net_divj)/(nodes-2)
			  for all i and j when j>i. Then remember the i and j value of the minimum 
			  value for later.*/
		

		sml_dist = 1000;
		total_div = 0;
		tot = 0;
		for(i=0; i<nodes; i++)
			for(j=0; j<i; j++)
				total_div += d_mat[i][j];
		
		for(i=0; i<nodes; i++)
			{
			for(j=0; j<i; j++)
				{
				
				d2r = (net_div[i][0] + net_div[j][0] -(2*d_mat[i][j]));
				dr = total_div - d_mat[i][j] - d2r;
				to = d2r +((nodes -2) * d_mat[i][j]) + (dr * 2);
				to = to/(2*(nodes - 2));
				
				if(to < sml_dist)
					{
					sml_dist = to;
					smli = i;
					smlj = j;
					}	
				}
			}
				











/*		sml_dist = 1000;
		
		for(i=0; i<nodes; i++)
			{
			for(j=(i+1); j<nodes; j++)
				{
				total_div = d_mat[i][j] - (net_div[i][0] + net_div[j][0])/(nodes - 2);
				if(total_div < sml_dist)
					{
					sml_dist = total_div;
					smli = i;
					smlj = j;
					}
					
				}
			}
		
*/		
		/* 3: Assign the line to the tree file */
		
		/* We ned to take four cases into account:
			A) smli & smlj both are NOT nodes
			B) smli is a node and smlj is NOT
			C) smlj is a node and smli is NOT
			D) BOTH smli & smlj are nodes.
		   For these reasons the following section is divided into 4 parts using if-else statements
		   each section dealing with the individual case.
		 */
		
		if(nodes > 3) /* There are num_of_seqs -3 internal branches in a tree */
			{
			
	/*		tree[num_of_seqs-untagged - nodes][seqs[smli][0]] = (d_mat[smli][smlj])/2 + (net_div[smli][1] - net_div[smlj][1])/2;
			if(tree[num_of_seqs-untagged - nodes][seqs[smli][0]] == 0) tree[num_of_seqs-untagged - nodes][seqs[smli][0]] = 0.000005;  *//* set a minimum distance to 0.000005 */
			
	/* A */
			
			if(!seqs[smlj][1] && !seqs[smli][1])  /* if smlj and smli are not nodes */
				{
				tree[num_of_seqs-untagged - nodes][seqs[smli][0]] = (d_mat[smli][smlj])/2 + (net_div[smli][1] - net_div[smlj][1])/2;
				if(tree[num_of_seqs-untagged - nodes][seqs[smli][0]] == 0) tree[num_of_seqs-untagged - nodes][seqs[smli][0]] = 0.000005;  /* set a minimum distance to 0.000005 */

				tree[num_of_seqs-untagged - nodes][seqs[smlj][0]] = (d_mat[smli][smlj] - tree[num_of_seqs-untagged - nodes][seqs[smli][0]]);
				if(tree[num_of_seqs-untagged - nodes][seqs[smlj][0]] == 0) tree[num_of_seqs-untagged - nodes][seqs[smlj][0]] = 0.000005; /* set a minimum distance to 0.000005 */
				}
			
			else
				{
	/* B */			
				if(!seqs[smlj][1] && seqs[smli][1]) 
					{
					tree[num_of_seqs-untagged - nodes][seqs[smlj][0]] = (d_mat[smli][smlj])/2 + (net_div[smli][1] - net_div[smlj][1])/2;
					if(tree[num_of_seqs-untagged - nodes][seqs[smlj][0]] == 0) tree[num_of_seqs-untagged - nodes][seqs[smlj][0]] = 0.000005;  /* set a minimum distance to 0.000005 */					

					/* go the the line on the tree file (specified by seqs) where the node was last used and assign the distances to all seqs on that line before */				
					for(j=0; j<num_of_seqs-untagged; j++)
						{
						if(tree[seqs[smli][0]][j] != 0)	/* if we have defined this node within a species previously */
							{
							tree[num_of_seqs-untagged - nodes][j] = (d_mat[smli][smlj] - tree[num_of_seqs-untagged - nodes][seqs[smlj][0]]); /* assign the distance of the previous node to this node to the sequence on this line */
							if(tree[num_of_seqs-untagged - nodes][j] == 0) tree[num_of_seqs-untagged - nodes][j] = 0.000005;  /* set a minimum distance to 0.000005 */

							}
						}
					}
				else
	/* C */			{
					if(!seqs[smli][1] && seqs[smlj][1]) 
						{
						tree[num_of_seqs-untagged - nodes][seqs[smli][0]] = (d_mat[smli][smlj])/2 + (net_div[smli][1] - net_div[smlj][1])/2;
						if(tree[num_of_seqs-untagged - nodes][seqs[smli][0]] == 0) tree[num_of_seqs-untagged - nodes][seqs[smli][0]] = 0.000005;  /* set a minimum distance to 0.000005 */					

						/* go the the line on the tree file (specified by seqs) where the node was last used and assign the distances to all seqs on that line before */				
						for(j=0; j<num_of_seqs-untagged; j++)
							{
							if(tree[seqs[smlj][0]][j] != 0)	/* if we have defined this node within a species previously */
								{
								tree[num_of_seqs-untagged - nodes][j] = (d_mat[smli][smlj] - tree[num_of_seqs-untagged - nodes][seqs[smli][0]]); /* assign the distance of the previous node to this node to the sequence on this line */
								if(tree[num_of_seqs-untagged - nodes][j] == 0) tree[num_of_seqs-untagged - nodes][j] = 0.000005;  /* set a minimum distance to 0.000005 */
		
								}
							}
						}
					else
	/* D */				{
						
						/* go the the line on the tree file (specified by seqs) where the smlj was last used and assign the distances to all seqs on that line before */				
						for(j=0; j<num_of_seqs-untagged; j++)
							{
							if(tree[seqs[smlj][0]][j] != 0)	/* if we have defined this node within a species previously */
								{
								tree[num_of_seqs-untagged - nodes][j] = (d_mat[smli][smlj] - tree[num_of_seqs-untagged - nodes][seqs[smli][0]]); /* assign the distance of the previous node to this node to the sequence on this line */
								if(tree[num_of_seqs-untagged - nodes][j] == 0) tree[num_of_seqs-untagged - nodes][j] = 0.000005;  /* set a minimum distance to 0.000005 */

								}
							}
						/* go the the line on the tree file (specified by seqs) where the smli was last used and assign the distances to all seqs on that line before */				
						for(j=0; j<num_of_seqs-untagged; j++)
							{
							if(tree[seqs[smli][0]][j] != 0)	/* if we have defined this node within a species previously */
								{
								tree[num_of_seqs-untagged - nodes][j] = (d_mat[smli][smlj] - tree[num_of_seqs-untagged - nodes][seqs[smlj][0]]); /* assign the distance of the previous node to this node to the sequence on this line */
								if(tree[num_of_seqs-untagged - nodes][j] == 0) tree[num_of_seqs-untagged - nodes][j] = 0.000005;  /* set a minimum distance to 0.000005 */

								}
							}
					
						}
					}
				}
			}							
				
		else  /* we don't want it to write a line for the last two joinings, they are paired together since there can be only n-3 internal branches */
			{
			second = d_mat[smli][smlj]/2 + (net_div[smli][1] - net_div[smlj][1])/2;
			}
		
		
		/* 4: Defines the distance from the current node to the remaining terminal nodes */
		x = y = 0;
		for(i=0; i<nodes; i++)
			{
			if(i != smli && i != smlj) 
				{
				y = 0;
				for(j=0; j<nodes; j++)
					{
					if(i != smli && i != smlj && j != smli && j !=  smlj)
						{
						tmp[x][y] = d_mat[i][j];
						y++;
						}
					}
				x++;	
				}
			}
		
		y = 0;
		for(i=0; i<nodes; i++)
			{
			if(i != smli && i != smlj)
				{
				tmp[nodes-2][y] = tmp[y][nodes-2] = ((d_mat[smli][i] + d_mat[smlj][i]) - d_mat[smli][smlj])/2;
				y++;
				}
			}
		tmp[nodes-2][nodes-2] = 0;
		
		
		/* 5: get rid of the old matrices and assign  new ones equal in size to the new nodes sides */
		
		
		/* initialise tmp2*/
		
		for(i=0; i<num_of_seqs - untagged - 1; i++)
			{
			tmp2[i][0] = i;
			tmp2[i][1] = FALSE;
			}
	
		/* transfer information for the seqs array */
		
		j = 0;
		for(i=0; i<nodes; i++)
			if(i != smli && i != smlj)
				{
				tmp2[j][0] = seqs[i][0];
				tmp2[j][1] = seqs[i][1];
				j++;
				}
		tmp2[nodes-2][0] = num_of_seqs-untagged - nodes;  /* tells the tree line number that the node was last used */
		tmp2[nodes-2][1] = TRUE;  /* this is now a node */
		
		/* initialise seqs*/
		
		for(i=0; i<num_of_seqs - untagged - 1; i++)
			{
			seqs[i][0] = i;
			seqs[i][1] = FALSE;
			}
				

		for(i=0; i<nodes-1; i++)
			{
			seqs[i][0] = tmp2[i][0];
			seqs[i][1] = tmp2[i][1];
			}
			
		
		/* initialise dmat array */
		
		for(i=0; i<num_of_seqs - untagged - 1; i++)
			for(j=0; j<num_of_seqs - untagged - 1; j++)
				d_mat[i][j] = 0;



		/* transfer the information from the d_mat array */
			
		for(i=0; i<nodes-1; i++)
			for(j=0; j<nodes-1; j++)
				d_mat[i][j] = tmp[i][j];
				
	
		nodes--;
		
		}

	/* if an outgroup exsists, then write an extra line on the tree to check the outgroup as a clade */
			
	for(i=0; i<num_of_seqs-untagged; i++) 
		{
		assigned = FALSE;
		for(j=0; j<num_of_seqs-3-untagged; j++)
			if(tree[j][i] != 0) assigned = TRUE;
		if(!assigned)
			{
			if(second == 0) tree[num_of_seqs-3-untagged][i] = 0.000005;
			else tree[num_of_seqs-3-untagged][i] = second;
			}
		}
					
	
	/* free up used memory before exiting */
		
	for(i=0; i<num_of_seqs - untagged; i++) free(net_div[i]);
	free(net_div);
	net_div = NULL;

	for(i=0; i<num_of_seqs - untagged; i++) free(d_mat[i]);
	free(d_mat);
	d_mat = NULL;
	
	for(i=0; i<(num_of_seqs - untagged); i++) free(distances[i]);
	free(distances);
	distances = NULL;
	
	for(i=0; i<num_of_seqs - untagged; i++) free(seqs[i]);
	free(seqs);
	seqs = NULL;

	for(i=0; i<num_of_seqs - untagged -1; i++) free(tmp[i]);
	free(tmp);
	tmp = NULL;

	for(i=0; i<num_of_seqs - untagged -1; i++) free(tmp2[i]);
	free(tmp2);
	tmp2 = NULL;
	
	
	}
/* end new */	


/* This function performs the g test for independance on a 2 by 2 array which contain 2 ratios */

void g_test(float **ratio, float **pvalue, float *gChi )
	{
	float totf = 0, totp = 0, totr = 0, tots = 0, fulltot = 0, a = 0, b = 0, c = 0, d = 0, expected[4], g = 0, value1 = 0, value2 = 0, value3 = 0, value4 = 0, williams = 0;
	float a2 = 0, b2 = 0, c2 = 0, d2 = 0, first = -2, antilog = 0;
	float chi[2][10] ={0.99,   0.95,   0.90,  0.50, 0.20, 0.10,  0.05,  0.025,  0.01,  0.005,
					   0.0002, 0.0039, 0.016,  0.45, 1.64, 2.706, 3.841, 5.024,  6.635, 7.87 };

	int i, j;
	
	for(j=0; j<num_of_seqs - untagged - 2; j++)  /* This needs to be completed for every branch of the tree */
		{

		/* initialise variables */		
		totf = 0; totp = 0; totr = 0; tots = 0; fulltot = 0; g = 0;
		a = 0; b = 0; c = 0; d = 0;
		first = -2;
		for(i=0; i<4; i++)
			expected[i] = 0;

		a = ratio[j][0];
		b = ratio[j][1];
		c = ratio[j][2];
		d = ratio[j][3];
		
		a2 = ratio[j][0];
		b2 = ratio[j][1];
		c2 = ratio[j][2];
		d2 = ratio[j][3];
			
		totf = ratio[j][0] + ratio[j][2];
		totp = ratio[j][1] + ratio[j][3];
		totr = ratio[j][0] + ratio[j][1];
		tots = ratio[j][2] + ratio[j][3];
		fulltot = totf + totp;




		/* Since Fishers exact test cannot handle big numbers there is no point using all that computation time to try to calculate it if the numbers are big...... and I really takes a long time if the numberrs are big */
		/* So what we'll do is not calculate the fishers value if the total number of all four values is greater than 200..... Fishers wouldn't be able to calculate more than 80 anyway */	
		
		if(a+b+c+d <= 200)
			{	

			/* fishers exact test */
	/*		printf("totf %f, totp %f, totr %f, tots %f\n", totf, totp, totr, tots);
			printf("a %f, b %f, c %f, d %f\n", a, b, c, d);
			printf("totf! %f, totp! %f, totr! %f, tots! %f, fulltot! %f\n", fctrl(totf), fctrl(totp), fctrl(totr), fctrl(tots), fctrl(fulltot));
	*/		
		
			value1 = log10(fctrl(totf)) + log10(fctrl(totp)) + log10(fctrl(totr)) + log10(fctrl(tots)) - log10(fctrl(fulltot));
	/*		printf("value1 = %f\t", value1);
	*/		value3 = 0;
		
			do
				{
				value2 = log10(fctrl(a)) + log10(fctrl(b)) + log10(fctrl(c)) + log10(fctrl(d));
	/*			printf("value2 %f\t", value2);
	*/				
				antilog = pow(10, (value1 - value2));
	/*			printf("antilog = %f\t", antilog); 
	*/			value3 = value3 + antilog ;
	/*			printf("value3 = %f\n", value3);
	*/			if(first == -2) first = value3;
		
				if((a*d - b*c) < 0)
					{
					a = a -1;
					d = d -1;
					b = b +1;
					c = c +1;
					}
				else
					{
					b = b -1;
					c = c -1;
					a = a +1;
					d = d +1;
					}
	/*			printf("value3 = %f\n", value3);		
	*/		
		
				}while(a >= 0 && b >= 0 && c >= 0 && d >= 0);
			
			if((a2*d2 - b2*c2) > 0)	
				{
				if(a2 > d2)
					{
					d2 = 0;
					a2 = a2 - d2;
					b2 = b2 + d2;
					c2 = c2 + d2;
					}
				else
					{
					a2 = 0;
					d2 = d2 - a2;
					b2 = b2 + a2;
					c2 = c2 + a2;
					}
				}
			else
				{
				if(b2 > c2)
					{
					c2 = 0;
					b2 = b2 - c2;
					a2 = a2 + c2;
					d2 = d2 + c2;
					}
				else
					{
					b2 = 0;
					c2 = c2 - b2;
					a2 = a2 + b2;
					d2 = d2 + b2;
					}
				}
			
			do
				{
				value2 = log10(fctrl(a2)) + log10(fctrl(b2)) + log10(fctrl(c2)) + log10(fctrl(d2));
	/*			printf("2value2 %f\n", value2);
	*/				
				antilog = pow(10, (value1 - value2));
	/*			printf("antilog %f\n", antilog);	
	*/			if(antilog < first) value3 = value3 + antilog;
			
		
				if((a*d - b*c) < 0)
					{
					a2 = a2 -1;
					d2 = d2 -1;
					b2 = b2 +1;
					c2 = c2 +1;
					}
				else
					{
					b2 = b2 -1;
					c2 = c2 -1;
					a2 = a2 +1;
					d2 = d2 +1;
					}
		
		
				}while(antilog < first && a2 > -1 && b2 > -1 && c2 > -1 && d2 > -1);
			
			pvalue[j][2] = value3;  /* position 2 will hold the fishers exact result */
			
			}
		
			/*  G test */
		
		a = ratio[j][0];
		b = ratio[j][1];
		c = ratio[j][2];
		d = ratio[j][3];
		
	
		value1 = (a*logE(a))+(b*logE(b))+(c*logE(c))+(d*logE(d));
		value2 = (totf*logE(totf))+(totp*logE(totp))+(totr*logE(totr))+(tots*logE(tots));
		value3 = (fulltot*logE(fulltot));
		value4 = 2*(value1 - value2 + value3);
		
	
	
		williams = 1 + ((((fulltot/totf)+(fulltot/totp)- 1)*((fulltot/totr)+(fulltot/tots) - 1))/(6 * fulltot));
		
		g = value4/williams;
		
		
		gChi[j] = g;  /* recording the value of chi */
		i = 0;
		while(g > chi[1][i])
			{
			if(i == 9) break;
			else i++;
			}
		
		switch(i)
			{
			case 0:

				pvalue[j][0] = 1.0;  /* max p-value */
				pvalue[j][1] = chi[0][i]; /*min p-value */
				break;
			case 9:
				pvalue[j][0] = chi[0][i];
				pvalue[j][1] = 0.0;
				break;
			default:
				pvalue[j][0] = chi[0][i-1];
				pvalue[j][1] = chi[0][i];
				break;
		
		
			}
		}
		
		
	}

/* calculates the natural log of value passed, and if the value is 0, returns 0 */
	
float logE(float value)
	{
	if(value > 0)
		value = log(value);
	else value = 0;	
	return(value);
	}


/* calculates factorial of given number */
float fctrl(float value)
	{
	float tmp , i;
	
	tmp = value;
	if(value > 1)
		for(i=value -1; i>=1; i--) tmp = tmp *  i;
	else
		tmp = 1;
	return(tmp);
	
	}

	
	
	
/*	Function: output_tree	Created 11/4/00    by Chris Creevey
	This function outputs the nj tree in nested parenthses format to the output file njtree.out.
	
	Procedure:
		We have a 2 dimensional dynamically allocated array (called build) to hold the different parts of the tree as we are building it.
		The array has to be (num_of_seqs - untagged - 2) X (length) of type char. (ie they are strings).
		As we travel down through the tree array (built using nj method) for every pairing we give them a seperate line on the 
		array, unless it pairs with something previously paired. Ifit is pairing with something previously paired, then 
		whichever is nearest the top of the array, is used as the base on which everything is added. When all the internal branches
		are traversed, if we still have several lines, they are paired together in parentheses on one line. This 'string' is then
		written to the file "njtree.out".
		
	note:
		When joining in parentheses there are three different cases:
		1)	First joining of two taxa.
		2)	A single taxa joining an internal branch.
		3)	Two internal branches joining.
		
		To keep track of these situations, an array is defined (called joined) which will represent the sequences on the tree array. 
		If the taxa have been part of some branch before (TRUE or FALSE), this will record it, along with the line on the 'build' array that it was 
		last included in a join.
		So on each line of the tree array check the status of whether the taxa had been included before or not. (ie TRUE or FALSE).
		
		So for case:
		1)	taxa will read FALSE, since they have never been joined before, so join them on a new line.
		2)	All of the taxa except for 1 will be TRUE so look for where you have both TRUE and FALSE on one line, and add the
			new taxa to the exsisting internal branch.
		3)	All taxa will read TRUE for this case, so join both on which ever is on the line closest to the top of the array.

*/



void output_tree(float **tree)
	{
	char **build = NULL, *temp = NULL, *pointer = NULL;
	int i=0, j=0, k=0, joins = 0, paired = FALSE, notpaired = FALSE, first = TRUE, second = TRUE, third = TRUE, place = 0, **joined = NULL, which = 0, reference = 0, count = 0;
	struct sequence *position = NULL;
	FILE *test = NULL;
	


	/***declaration and initialisation **/	
	build = malloc((num_of_seqs - untagged - 2)*(sizeof(char*)));
	if(!build)
			{
			printf("\n\t Out of memory\n");
			clean_exit();
			}
	for(i=0; i<(num_of_seqs - untagged - 2); i++)
		{
		 build[i] = malloc((100 * num_of_seqs)*(sizeof(char)));
		if(!build[i])
			{
			printf("\n\t Out of memory\n");
			clean_exit();
			}
		}
	for(i=0; i<num_of_seqs - untagged - 2; i++)
		for(j=0; j<(10 * num_of_seqs); j++)
			build[i][j] = '\0';
	
	
	joined = malloc((num_of_seqs - untagged)*(sizeof(int*)));
	if(!joined)
			{
			printf("\n\t Out of memory\n");
			clean_exit();
			}
	for(i=0; i<num_of_seqs - untagged; i++)
		{
		joined[i] = malloc((2)*(sizeof(int)));
		if(!joined[i])
			{
			printf("\n\t Out of memory\n");
			clean_exit();
			}
		}
			
	for(i=0; i<num_of_seqs - untagged; i++)
		{
		joined[i][0] = FALSE; /* this records if a taxa has been joined before */
		joined[i][1] = -1; /* this records where the last joining took place */
		}
	
	temp = malloc((100 * num_of_seqs)*(sizeof(char)));  /* used to manipualte strings */
	if(!temp)
			{
			printf("\n\t Out of memory\n");
			clean_exit();
			}
	for(j=0; j<(100 * num_of_seqs); j++)
		temp[j] = '\0';
		
	/* go through the sequences and get rid of any parentheses in the names, replace with "-"  */
	
	position = start;
	while(position)
		{
		while((pointer = strchr(position->nickname, '(')) != NULL) pointer[0] = '_';
		while((pointer = strchr(position->nickname, ')')) != NULL) pointer[0] = '_';
		while((pointer = strchr(position->name, '(')) != NULL) pointer[0] = '_';
		while((pointer = strchr(position->name, ')')) != NULL) pointer[0] = '_';
		while((pointer = strchr(position->name, '-')) != NULL) pointer[0] = '_';
		while((pointer = strchr(position->name, '-')) != NULL) pointer[0] = '_';
		while((pointer = strchr(position->nickname, '-')) != NULL) pointer[0] = '_';
		while((pointer = strchr(position->nickname, '-')) != NULL) pointer[0] = '_';
		position = position->next;
		}		


	/** start of the algorithm **/
	
	for(i=0; i<num_of_seqs - untagged - 2; i++)
		{
		first = TRUE;
		second = TRUE;
		third = TRUE;
		paired = notpaired = FALSE;
		which = 0;
		for(j=0; j<num_of_seqs - untagged; j++)
			{
			
			if(tree[i][j] != 0)
				{
				if(joined[j][0] == FALSE) notpaired = TRUE;
				else paired = TRUE;
				}
			}
		
		position = start;
		while(position && !position->tag) position = position->next;

		for(j=0; j<(num_of_seqs - untagged); j++)
			{
	
			if(tree[i][j] != 0)
				{
				if(notpaired && !paired)  /* if this is a first joining of two taxa */
					{
					if(first)
						{
						place = joins;
						strcpy(build[place], "(");
						strcat(build[place], position->nickname);
						strcat(build[place], ",");
						joined[j][0] = TRUE;
						joined[j][1] = place;
						first = FALSE;
						joins++;
						}
					else
						{
						strcat(build[place], position->nickname);
						strcat(build[place], ")");
						joined[j][0] = TRUE;
						joined[j][1] = place;
						}
					}
				
				if(notpaired && paired)  /* if one taxa is joining an internal branch */
					{
					if(joined[j][0] && which == 0)   /* if we meet the internal branch first */
						{
						strcpy(temp, "(");
						place = joined[j][1];
						strcat(temp, build[place]);
						strcpy(build[place], temp);
						which = 1;
						}
					else
						{	
						if(!joined[j][0] && which == 0)  /* if we meet the taxa first */
							{
							strcpy(temp, ",");
							strcat(temp, position->nickname);
							strcat(temp, ")");
							place = j;
							joined[j][0] = TRUE;
							which = 2;
							}
						else	
							{
							if(joined[j][0] && which == 2)  /* we have already met the taxa and this is the first time we meet the branch */
								{
								strcat(build[joined[j][1]], temp);
								strcpy(temp, build[joined[j][1]]);
								strcpy(build[joined[j][1]], "(");
								strcat(build[joined[j][1]], temp);
								joined[place][1] = joined[j][1];
								which = 3;
								}
							else
								{
								if(!joined[j][0] && which == 1)   /* we have already met the branch and now we meet the taxa */
									{
									strcat(build[place], ",");
									strcat(build[place], position->nickname);
									strcat(build[place], ")" );
									joined[j][0] = TRUE;
									joined[j][1] = place;
									which = 3;
									}
								}
							}
						}
					}
				if(!notpaired && paired)  /* if we are joining two internal branches */
					{
					if(first)
						{
						strcpy(temp, "(");
						place = joined[j][1];
						strcat(temp, build[place]);
						strcpy(build[place], temp);
						strcat(build[place], ",");
						first = FALSE;
						}
					else
						if(joined[j][1] != place && second)
							{
							strcat(build[place], build[joined[j][1]]);
							strcat(build[place], ")");
							second = FALSE;
							reference = joined[j][1];
							if(place < reference)
								{
								strcpy(build[reference], "\0" );
								for(k=0; k<num_of_seqs - untagged; k++) 
									if(joined[k][1] == reference && joined[k][0]) joined[k][1] = place;  /*reference the taxon to the correct line in the build array */
								}
							else
								{
								strcpy(build[joined[j][1]], build[place]);
								strcpy(build[place], "\0");
								for(k=0; k<num_of_seqs - untagged; k++) 
									if(joined[k][1] == place && joined[k][0]) joined[k][1] = joined[j][1];  /*reference the taxon to the correct line in the build array */
								}
							}
					}
				}
				
			position = position->next;
			while(position && !position->tag) position = position->next;
			}

		}
	
	
	
	/* now join those seperate internal branches together */	
	first = second = TRUE;
	
	/* check to make sure that there are more than 1 internal branches left */	
	for(i=0; i<num_of_seqs - untagged -2; i++)
		{
		if(strcmp(build[i], "\0") != 0)
			{
			if(!first)
				second = FALSE;
			else first = FALSE;
			}
		}				
	if(!second) /* If there is more than one internal branch left to join */
		{
		first = TRUE;
		for(i=0; i<num_of_seqs - untagged -2; i++)
			{
			if(strcmp(build[i], "\0") != 0)
				{
				if(first) 
					{
					strcpy(temp, "(" );
					j = 0;
					while(build[i][j] != '\0')   /* check to make sure that the last grouping wasn't a single taxa */
						{
						if(build[i][j] == '(') count++;
						if(build[i][j] == ')') count--;
						j++;
						}
					if(count != 0)    /* if it was a single taxa, then modify the line to delete the first ( and , to make the tree valid */
						{
						pointer = strchr(build[i], '('); pointer[0] = ' ';
						pointer = strchr(build[i], ','); pointer[0] = ' ';
						}
					strcat(temp, build[i]);
					strcpy(build[i], temp);
					place = i;
					first = FALSE;					
					}
				else
					{
					j = 0;
					while(build[i][j] != '\0')   /* check to make sure that the last grouping wasn't a single taxa */
						{
						if(build[i][j] == '(') count++;
						if(build[i][j] == ')') count--;
						j++;
						}
					if(count != 0)     /* if it was a single taxa, then modify the line to delete the first ( and , to make the tree valid */
						{
						pointer = strchr(build[i], '('); pointer[0] = ' ';
						pointer = strchr(build[i], ','); pointer[0] = ' ';
						}

					strcat(build[place], ",");
					strcat(build[place], build[i]);
					strcpy(build[i], "\0");
					}
				}
			}
		strcat(build[place], ");" );
		}

	/* Write to output file */

	if(parenthesis == NULL)
			if((parenthesis = fopen("njtree.out", "w")) == NULL)		/* check to see if the file can be opened/created */
			printf("\n\n\tCannot open the output file, named njtree.out,\n Tree file not written\n\n");	

	fprintf(parenthesis, "%s", build[place]);
		
	fflush(parenthesis);  
/*	fclose(parenthesis);
	parenthesis = NULL;
*/
		
	/* free up memory used  */
	
	for(i=0; i<(num_of_seqs - untagged - 2); i++) free(build[i]);
	free(build);
		
	for(i=0; i<num_of_seqs - untagged; i++) free(joined[i]);
	free(joined);
	
	free(temp);
		
	
	
		
	}			



						


/* This function compares the sequences pairwise and counts the number of occurances of 'mutations' between the sequences
	In reality it counts how many times a certain nucliotide is replaced by another, or not as the case may be. The results
	of this count are then divided by the total number of pairwise comparisons made, to give a percentile result. This is 
	outputted on the array 'matrix' that is passed by reference to the function, this also includes gaps in its calculation */ 

void substitution_matrix(float ***ratio)
	{
	struct sequence *position1 = start, *position2 = start->next;
	int i = 0, j = 0, x = 0, y = 0, z = 0, total = 0, nuc_count[3][5] = {0,0,0,0,0 , 0,0,0,0,0 , 0,0,0,0,0};
	float tmp = 0;
	for(i=0; i<3; i++)
		for(j=0; j<5; j++)
			for(x=0; x<5; x++)
					ratio[i][j][x] = 0;
	
	
	while(position1 != NULL)
		{
		position2 = position1->next;
		while(!position1->tag && position1 != NULL) 
			{
			position1 = position1->next;
			position2 = position1->next;
			}
		while(position2 != NULL)
			{
			while(!position2->tag && position2 != NULL)
				{
				position2 = position2->next;
				}
			i = 0, j = 0;
			while(position1->bases[i] != 193)
				{
				for(j = 0; j<3; j++)
					{
					total++;
					for(x=0; x<5; x++)
						{
						if(what[x] == codons[position1->bases[i]][j]) break;
						}
					for(y=0; y<5; y++)
						{
						if(what[y] == codons[position2->bases[i]][j]) break;
						}
					nuc_count[(degenerate_sites[code][position1->bases[i]][j])/2][x]++;  /* count how many times we check a change from this codon at either 0 2 or 4 fold degen sites */
					nuc_count[(degenerate_sites[code][position2->bases[i]][j])/2][y]++;	
					ratio[(degenerate_sites[code][position1->bases[i]][j])/2][x][y]++;
					ratio[(degenerate_sites[code][position2->bases[i]][j])/2][y][x]++;
					
						
					}
				i++;
				}
			position2 = position2->next;
			}
		position1 = position1->next;
		}
	fprintf(outfile, "Nucleotide substitution matrices\n");	
	for(z=0; z<3; z++)
		{
		fprintf(outfile, "\n%d-fold degenerate sites\n", z*2);
		for(x=0; x<5; x++)
			{
			for(y=0; y<5; y++)
				{
			
				if(nuc_count[z][x] != 0)
					{
					ratio[z][x][y] = ratio[z][x][y]/(float)nuc_count[z][x];  /* Total * 2 because we count the change to and from each of the ?-fold degenerate sites, so we have twice as many comparisons */
					}
				else ratio[z][y][x] = ratio[z][x][y] = 0;
		
				if(y==0)
					{
					
					fprintf(outfile, "From: %c ", what[x]);
					}
				fprintf(outfile, "%f ", ratio[z][x][y]);
				
				}
			fprintf(outfile, "\n");
			}
		}
	}			





/* This is the function through which the user can define the outgroup on the tree. 
	It then calls check_outgroup to validate the choice	, and if it is valid, then it calls
	tree_directionality, which does the work */				
					
void define_outgroup(void)
	{

	int choice2 = '\0', valid = FALSE, quit = FALSE, count = 0;
	struct sequence *position = NULL;


	/* Initialise all sequences to not being part of the outgroup */
	position = start;
	while(position != NULL)
		{
		position->outgroup = FALSE;
		position = position->next;
		}

	/* SELECT THOSE SEQUENCES WHICH ARE TO BE PART OF THE OUTGROUP */
	do{
		count = 0;
		position = start;
		printf("\n\nIn order to give the tree directionality, please define those sequences which form the outgroup");
		printf("\n\nSequence numbers, and names as follows:");
		do {
			position = start;
			do{
				if(position->tag && !position->outgroup) printf("\nNo: %d Name: %s ", (position->seq_num)+1, (position->name));
				position = position->next;
				}while(position != NULL);		
			choice2 = getint("\n\nPlease select the number of a sequence belonging to the outgroup, and press return\nEnter 0 when finished\n", 0, num_of_seqs, 0);
			if(choice2 == 0) choice2 = '\0';
			else
				{
				position = start;
				while(position != NULL && (position->seq_num != (choice2 - 1))) position = position->next;
				if(position != NULL)
					{
					position->outgroup = TRUE;
					count++;
					}
				}
			}while(choice2 != '\0');
		
		if(check_outgroup(count) == TRUE)
			{
			 valid = TRUE;
			
			}
		else
			{
			position = start;
			while(position)
				{
				position->outgroup = FALSE;
				position = position->next;
				}
			printf("\n\nThe chosen sequences are not valid outgroups");
			choice2 = getint("\n\tPress 1 to try again\n\tPress 2 to choose the first sequence in the file\n\tOr press 3 to choose the last sequence in the file\n", 1, 3, 1);
			if(choice2 == 2)
				{
				quit = TRUE;
				start->outgroup = TRUE;
				check_outgroup(1);
				}
			if(choice2 == 3)
				{
				quit = TRUE;
				last->outgroup = TRUE;
				check_outgroup(1);
				}
			}
		}while(!valid && !quit); 
	
	}
	
	


void assign_codon_up(int codon, struct node *position)
	{
	int tmp[66], other[66];  /* There are 64 codons plus a gap codon, so to make them into strings we add a'\0' at the end so the arrays need to be 66 in length */
	int l = 0, k=0, x=0, found = FALSE;
	
	for(l=0; l<66; l++)	  /* initialise the arrays */
		{
		tmp[l] = FALSE;
		other[l] = FALSE;
		} 
	/* Travel down through the tree to the bottom */	
	if(position->node1 != NULL) assign_codon_up(codon, position->node1);
	if(position->node2 != NULL) assign_codon_up(codon, position->node2);
	
	/* initialise the present node before assigning the ancestors */
	for(l=0; l<66; l++) position->codon_ances[l] = FALSE;
	l=0;
	
	/*First assign those nucliotides which are at this node to the ancestors */
	
	if(position->seq_num1 != NULL)
		{
		position->codon_ances[(position->seq_num1)->bases[codon]] = TRUE;
		}
		
	if(position->seq_num2 != NULL)
		{
		position->codon_ances[(position->seq_num2)->bases[codon]] = TRUE;
		}
		
		
		
		/* Second, check any children to see what the ancestor was assigned for those */	

	if(position->node1 != NULL) 
		{
		
		/* copy the child into the tmp */
		for(l=0; l<66; l++)
			{
			tmp[l] = (position->node1)->codon_ances[l];
			}
	
		found = FALSE;
		/* if any of the states in the children are the same as in the present node then those are the only ones from that child that are kept at this position */
		for(l=0; l<66; l++)
			{
			if(position->codon_ances[l] == TRUE && tmp[l] == TRUE)
				{
				other[l] = TRUE; /* other will hold all the codons that are the same */
				found = TRUE;
				
				}
			}
		
		if(found == TRUE) /* if there are the same codons in the child as in the present ancestor */
			{
			for(l=0; l<66; l++)
				{
				position->codon_ances[l] = other[l];  /* assign all other to the present position */
				other[l] = FALSE;  /* initialise other as we go */
				}
			}
		else   /* there are no common codons, so we add all codons in this child to the ancestor here */
			{
			for(l=0; l<66; l++)
				{
				if(tmp[l] == TRUE)
					position->codon_ances[l] = tmp[l];
				}
			}	
			
		}




	found = FALSE;
	for(l=0; l<66; l++)
		{
		other[l] = '\0';
		tmp[l] = '\0';
		}

	if(position->node2 != NULL) 
		{
		
		/* copy the child into the tmp */
		for(l=0; l<66; l++)
			{
			tmp[l] = (position->node2)->codon_ances[l];
			}


		/* if any of the states in the children are the same as in the present node then those are the only ones from that child that are kept at this position */
		
		
		for(l=0; l<66; l++)
			{
			if(position->codon_ances[l]== TRUE && tmp[l] == TRUE)
				{
				other[l] = tmp[l]; /* other will hold all the codons that are the same */
				found = TRUE;
				}
			}
		
		if(found == TRUE) /* if there are the same codons in the child as in the present ancestor */
			{
			for(l=0; l<66; l++)
				{
				
				position->codon_ances[l] = other[l];  /* assign all other to the present position */
				other[l] = FALSE;  /* initialise other as we go */
				}
			}
		else   /* there are no common codons, so we add all child codons calculated to the ancestor here */
			{
			
			for(l=0; l<66; l++)
				{
				if(tmp[l] == TRUE)
					position->codon_ances[l] = tmp[l];
				}
			}	
			
		}

	
	}
	
	
	
	




/* This function simply solves abiguities at the root, before assign_ances_down is called, so that there can never be an ambiguity in the tree. */
/* It solves ambiguities by assigning the ancesral root to that of the outgroup */
/* this however still leaves the chance for ambiguous sites, with an arbitrary descision being made if there is ambiguity a the top of the */
/* clade that defines the outgroup: For these purposes, you can only garantee unambiguous descisions if you only define ONE out group */
/* However this is not a great way of solving this, the other way would be to define the outgroup in the middle of the outgroup, so that the */
/* tree would still be rooted correctly, and the separation between the outgroup and the rest would be taken away from the root, where the ambiguities */
/* may lie. */
 

void assign_root_codon(int codon)
	{
	int x = 0, y = 0, found = FALSE;
	char c = '\0';

	for(y=0; y<66; y++)
		{
		if(tree_top->codon_ances[y] == TRUE) x++;
		}
		if(x > 1)
			{
			/* if the tree is rooted about a single sequence */
			if(tree_top->seq_num1 != NULL)  
				{
				if((tree_top->seq_num1)->outgroup == TRUE)
					{
					for(y=0; y<66; y++) tree_top->codon_ances[y] = FALSE;
					tree_top->codon_ances[(tree_top->seq_num1)->bases[codon]] = TRUE;
					found = TRUE;
					}
				}
			if(tree_top->seq_num2 != NULL)
				{
				if((tree_top->seq_num2)->outgroup == TRUE)
					{
					for(y=0; y<66; y++) tree_top->codon_ances[y] = FALSE;
					tree_top->codon_ances[(tree_top->seq_num2)->bases[codon]] = TRUE;
					found = TRUE;
					}
				}


			/* if the tree is rooted about an internal node */	
			if(!found)
				{
				if(tree_top->node1 != NULL)
					{
					if((look_for_ances(tree_top->node1)) != FALSE)
						{
						found = TRUE;
						for(y=0; y<66; y++) tree_top->codon_ances[y] = FALSE;
						for(y=0; y<66; y++)
							{
							if((tree_top->node1)->codon_ances[y] == TRUE)
								{
								tree_top->codon_ances[y] = TRUE;
								break;
								}
							}
							
						} 
					}
				if(!found)
					if(tree_top->node2 != NULL)
						{
						if((look_for_ances(tree_top->node2)) != FALSE)
							{
							found = TRUE;
							for(y=0; y<66; y++) tree_top->codon_ances[y] = FALSE;
							for(y=0; y<66; y++)
								{
								if((tree_top->node2)->codon_ances[y] == TRUE)
									{
									tree_top->codon_ances[y] = TRUE;
									break;
									}
								}
								
							}
						}
				}
			}
	}










/* this function travels down the tree solving any ambiguities for ancestors, by checking previous branches
	or if necessary using a substitution matrix to solve it...... it doesn't solve for an ambiguity at the 
	root though      */
void assign_codons_down(int codon, struct node *position, float ***subst_matrix)
	{
	int x = 0, y = 0, i = 0, j = 0, num = 0, l = 0, k = 0, highest_codon = 0, found = FALSE;
	float total = 0, highest_total = 0;
	char tmp[66], other[66];
	
	for(i=0; i<66; i++)
		{
		tmp[i] = FALSE;
		other[i] = FALSE;
		}

	for(i=0; i<66; i++)
		{
		if(position->codon_ances[i] == TRUE) x++;
		}
	



	if(x == 1)
		{
		/* First check the left part of the tree */
		if(position->node1 != NULL)
			{
			x = 0;
			for(i=0; i<66; i++)
				{
				if((position->node1)->codon_ances[i] == TRUE) x++;
				}

			if(x > 1)
				{
				/* copy ancestor into tmp */
				for(i=0; i<66; i++)
					{
					tmp[i] = (position->node1)->codon_ances[i];
					}
					
				/* Identify any codons which are the same in the present position and the child */
				for(i=0; i<66; i++)
					{
					
					if(position->codon_ances[i] == TRUE && tmp[i] == TRUE)
						{
						other[i] = TRUE; /* other will hold all the codons that are the same */
						found = TRUE;
						
						}
					}
		
				if(found == TRUE) /* if there are the same codons in the child as in the present ancestor */
					{
					for(i=0; i<66; i++)
						{
						(position->node1)->codon_ances[i] = other[i];  /* assign all in other to the childs ancestor since there can be only one ancestor, then there can be only one TRUE in other at this time */
						other[i] = FALSE;  /* initialise other as we go */
						}
					}
				else   /* if there is an ambiguity, the next level is to look at the amino acids coded*/
					{
					
					
					/* find the codon number of the ancestor of the present position for comparison in the next seection */
					for(i=0; i<66; i++) if(position->codon_ances[i] == TRUE) break; /* i is now the codon number of the present ancestor */
					
					for(j=0; j<66; j++)
						{
						if(tmp[j] == TRUE) /* if that codon is one of the childs ancestors */
							if(genetic_codes[code][position->codon_ances[i]] == genetic_codes[code][tmp[j]]) /* if the two codons code for the same amino acid */ 
								{
								other[j] = TRUE;  /* other will hold all of those codons that are ancestors and code for the same amino acid */						
								}
						}
					
					
					/* count how many codons are in other */
					l=0;
					for(j=0; j<66; j++) if(other[j] == TRUE) l++;


					if(l != 1)  /* if its not equal to 1 then either we have multiples that are different or multiples that are the same */
						{
						if(l == 0)   /* if none code for the same amino acid */ 
							{
							for(j=0; j<66; j++) /* this is to make sure that regardless of whether we found loads of the same Amino acid or none, all candidates are in the array other */
								{
								other[j] = tmp[j];
								}
							}
					
						total = y = k = 0;
						highest_total = 0;
						for(y=0; y<66; y++)
							if(other[y] == TRUE)
								{
								for(j=0; j<3; j++)
									{
									x= 0;
									while(what[x] != codons[position->codon_ances[i]][j]) x++;  /* i points at the codon in the present ancestor, and j points to the position in the codon */
									/* now x is equal the nucleotide at position j at the present position on the tree */
									k = 0;
									while(what[k] != codons[y][i]) k++;   /* now k points to the nucleotide in position i on the codon in other[l] */
								
									total += subst_matrix[degenerate_sites[code][position->codon_ances[i]][j]][x][k];  /* keep a running total of the likelihood of the changes for this codon */
									}
								if(total > highest_total)
									{
									highest_total = total;
									highest_codon = y;
									}
								total = 0;
								}
								
						}
					else  /* If there is only one answer after the amino acid check, then we make that the answer */
						{
						for(j=0; j<66; j++)
							{
							if(other[j] == TRUE)
								highest_codon = j;
							}
						}
						
					for(i=0; i<66; i++)
						(position->node1)->codon_ances[i] = FALSE;
						
					(position->node1)->codon_ances[highest_codon] = TRUE;  /* Since the substitution matrix was based on the alignment, and represents the chances of any nucleotide
																			changing to another, (or not changing) based on this dataset, when we add up the likelihoods,the codon 
																			with the highest likelihood is chosen as the descendant of this node. This means that there is a pressure to stay the same 
																			If the dataset shows that the nucleotide is likely to change, then it will force the final decision to represent this */
															
					}
				}
			}



		/* initialise before the second part of the algorithm */
		for(i=0; i<66; i++)
			{
			tmp[i] = '\0';
			other[i] = '\0';
			}

		found = FALSE;






		/* next check the right part of the tree */
		if(position->node2 != NULL)
			{
			x = 0;
			for(i=0; i<66; i++)
				{
				if((position->node2)->codon_ances[i] == TRUE) x++;
				}

			if(x > 1)
				{
				/* copy ancestor into tmp */
				for(i=0; i<66; i++)
					{
					tmp[i] = (position->node2)->codon_ances[i];
					}
					
				/* Identify any codons which are the same in the present position and the child */
				for(i=0; i<66; i++)
					{
					
					if(position->codon_ances[i] == TRUE && tmp[i] == TRUE)
						{
						other[i] = TRUE; /* other will hold all the codons that are the same */
						found = TRUE;
						
						}
					}
		
				if(found == TRUE) /* if there are the same codons in the child as in the present ancestor */
					{
					for(i=0; i<66; i++)
						{
						(position->node2)->codon_ances[i] = other[i];  /* assign all in other to the childs ancestor since there can be only one ancestor, then there can be only one TRUE in other at this time */
						other[i] = FALSE;  /* initialise other as we go */
						}
					}
				else   /* if there is an ambiguity, the next level is to look at the amino acids coded*/
					{
					
					
					/* find the codon number of the ancestor of the present position for comparison in the next seection */
					for(i=0; i<66; i++) if(position->codon_ances[i] == TRUE) break; /* i is now the codon number of the present ancestor */
					
					for(j=0; j<66; j++)
						{
						if(tmp[j] == TRUE) /* if that codon is one of the childs ancestors */
							if(genetic_codes[code][position->codon_ances[i]] == genetic_codes[code][tmp[j]]) /* if the two codons code for the same amino acid */ 
								{
								other[j] = TRUE;  /* other will hold all of those codons that are ancestors and code for the same amino acid */						
								}
						}
					
					
					/* count how many codons are in other */
					l=0;
					for(j=0; j<66; j++) if(other[j] == TRUE) l++;


					if(l != 1)  /* if its not equal to 1 then either we have multiples that are different or multiples that are the same */
						{
						if(l == 0)   /* if none code for the same amino acid */ 
							{
							for(j=0; j<66; j++) /* this is to make sure that regardless of whether we found loads of the same Amino acid or none, all candidates are in the array other */
								{
								other[j] = tmp[j];
								}
							}
					
						total = y = k = 0;
						highest_total = 0;
						for(y=0; y<66; y++)
							if(other[y] == TRUE)
								{
								for(j=0; j<3; j++)
									{
									x= 0;
									while(what[x] != codons[position->codon_ances[i]][j]) x++;  /* i points at the codon in the present ancestor, and j points to the position in the codon */
									/* now x is equal the nucleotide at position j at the present position on the tree */
									k = 0;
									while(what[k] != codons[y][i]) k++;   /* now k points to the nucleotide in position i on the codon in other[l] */
								
									total += subst_matrix[degenerate_sites[code][position->codon_ances[i]][j]][x][k];  /* keep a running total of the likelihood of the changes for this codon */
									}
								if(total > highest_total)
									{
									highest_total = total;
									highest_codon = y;
									}
								total = 0;
								}
								
						}
					else  /* If there is only one answer after the amino acid check, then we make that the answer */
						{
						for(j=0; j<66; j++)
							{
							if(other[j] == TRUE)
								highest_codon = j;
							}
						}
						
					for(i=0; i<66; i++)
						(position->node2)->codon_ances[i] = FALSE;
						
					(position->node2)->codon_ances[highest_codon] = TRUE;  /* Since the substitution matrix was based on the alignment, and represents the chances of any nucleotide
																			changing to another, (or not changing) based on this dataset, when we add up the likelihoods,the codon 
																			with the highest likelihood is chosen as the descendant of this node. This means that there is a pressure to stay the same 
																			If the dataset shows that the nucleotide is likely to change, then it will force the final decision to represent this */
															
					}
				}
			}

		}
	if(position->node1 != NULL) assign_codons_down(codon, position->node1, subst_matrix);
	if(position->node2 != NULL) assign_codons_down(codon, position->node2, subst_matrix);
	
	}
	
	
	
	
	
	
	/* This function assigns the codon number to the ancestral sequence data in each node, depending on the position of the current nucliotide */
void assign_codon(int codon, struct node *position)
	{
	int i = 0;
	
	if(position->node1 != NULL) assign_codon(codon, position->node1);
	if(position->node2 != NULL) assign_codon(codon, position->node2);

	for(i=0; i<66; i++)
		if(position->codon_ances[i] == TRUE) position->ances_seq[codon] = i;  /* This assigns the chosen codon to the ancestral sequence overall */
	
	}

/* function starts the checking process by calling assign_ances_up, and assign_ances_down which are both recursive and calls assign_codon_num, to calculate the ancestral sequence */
void ancestral_codon(int codon, float ***subst_matrix)
	{
	struct node *position = NULL;
	
	position = tree_top;
	
	assign_codon_up(codon, position);

	assign_root_codon(codon);  /* this tries to solve ambiguities at  the root */

	assign_codons_down(codon, position, subst_matrix);
	
	assign_codon(codon, position);
	
	}

			
						
	
	
	
	
	
						
	
