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
#include "linked tree.h"
#include "evolve.h"
#include "adaptive_tree.h"
#include "Li_Wu_19851993.h"


void linked_tree(float **tree_array)
	{
	int *done_col, row = 0, col = 0, taxa = 0, one = 0, i = 0, seq = 0, j =0, joining = FALSE;
	
	struct node **build = NULL, *position = NULL;
	
	struct sequence *list_pos = start;
	
	/* start at the top of the tree_array and join together the nodes */ 
	
	/* This keeps track of those taxa we have already joined */
	done_col = malloc((num_of_seqs - untagged)*sizeof(int));
	if(!done_col)
		{
		printf("\n Out of memory \n");
		clean_exit();
		}
	for(i=0; i<num_of_seqs - untagged; i++)
		done_col[i] = -1;
		
	
	/* this keeps the diferent parts of the tree during the building process */
	build = malloc((num_of_seqs - untagged - 2) * sizeof(struct node *));
	if(!build)
		{
		printf("\n Out of memory \n");
		clean_exit();
		} 
	for(i=0; i<num_of_seqs - untagged - 2; i++)
		build[i] = NULL;
	
	
	for(row = 0; row< num_of_seqs - untagged; row++) /* While we are not at the end of the array */
		{
		joining = FALSE;
		for(col = 0; col < num_of_seqs - untagged; col++)
			{
			if(tree_array[row][col] != 0) 
				joining = TRUE;						/* The last two lines of the array may comtain nothing, so we don't want a node created in this case, so this checks if there is a joining at this row if there is not, then that row is skipped */
			}
		if(joining)
			{
		
			position = malloc(sizeof(node_type));    /* Each line of the tree arrray represents a node, (the joining of two neighbours), so for every line we create a new node */
			if(!position)
				{
				printf("\n Out of memory \n");
				clean_exit();
				}
			/* initialise all parts ofthe node */		
			position->seq_num1 = NULL;
			position->seq_num2 = NULL;
			position->node1 = NULL;
			position->node2 = NULL;
			position->prev = NULL;
			position->ancestor[0] = '\0';
			position->ances_seq = malloc(((start->length)/3 +1 )*sizeof(int));
			if(!position->ances_seq)
				{
				printf("\n Out of memory \n");
				clean_exit();
				}
			for(j=0; j<(start->length)/3 +1; j++) position->ances_seq[j] = 0;
			
			position->gaprun = 0;
		


			taxa = 1;    /* keeps track of how many taxa we have joined (never more than 2) */
			one = -1;
			for(col = 0; col < num_of_seqs - untagged; col ++)
				{
				if(done_col[col] == -1 && tree_array[row][col] != 0)  /* if we haven't joined this taxa yet and itis identified for joining at this node */
					{
					if(taxa == 1)  /* signifies we have joined one taxa at this node (max of two) */
						{
						list_pos = start;
						seq = 0;
						while(!list_pos->tag && list_pos != NULL) list_pos = list_pos->next;   /* as long as the sequence is tagged and we are not at the end of the list, go to the next sequence in the list */
						while(list_pos != NULL && seq != col)         /* While we haven't found the sequence identified by the position on the tree array */
							{
						 	list_pos = list_pos->next;
						 	seq++;
							while(!list_pos->tag && list_pos != NULL) list_pos = list_pos->next;
				 			}
						if(list_pos != NULL)
							{
							position->seq_num1 = list_pos;  /* each node on the tree will point to the actual sequence it refences in the linked list */
							taxa =  2;
							done_col[col] = row;   /* assign the number of the row to that position in the done_col array, so that we keep track of the last time we joined this taxa */
							}                      /* This also means we know what clade this taxa belongs to, so we know to join those clades together */
						}
					else
						{	
						if(taxa == 2)              /* if we are joining the second taxa to this node */
							{
							list_pos = start;
							seq = 0;	
							while(!list_pos->tag && list_pos != NULL) list_pos = list_pos->next;
							while(list_pos != NULL && seq != col)
								{
							 	list_pos = list_pos->next;
							 	seq++;
								while(!list_pos->tag && list_pos != NULL) list_pos = list_pos->next;
					 			}
							if(list_pos != NULL)
								{
								position->seq_num2 = list_pos;
								taxa = 3;
								done_col[col] = row;
								}
							}
						}
					}
				else     /*  if we identify this sequence as being part of this node (within the clade ) BUT previously joined */
					{
					if(done_col[col] != -1 && tree_array[row][col] != 0)
						{
						if(taxa == 1)
							{
							position->node1 = build[done_col[col]]; /* node1 of this node is assigned to the joining last made with this clade, we know that this joining includes the clade with this taxa in it, and something else */
							build[done_col[col]]->prev = position;
							build[done_col[col]] = NULL;
							taxa = 2;
							one = done_col[col];
							done_col[col] = row;
							}
						else
							{	
							if(taxa == 2)
								{
								if(done_col[col] == one) done_col[col] = row;
								else
									{
									position->node2 = build[done_col[col]];
									build[done_col[col]]->prev = position;
									build[done_col[col]] = NULL;
									taxa = 3;
									done_col[col] = row;
									}
								}
							else
								{
								if(taxa == 3)
									{
									done_col[col] = row;
									}
								}
							}
						}
					}
				}
			build[row] = position;
			position = NULL;
			}
		}



	/* We are finished putting the clades together as described by the tree array, now we need to join those clades that are left seperate by the array */

	taxa = 1;
	for(i=0; i<num_of_seqs - untagged - 2; i++)
		{
		if(build[i] != NULL)
			{
			if(taxa == 1)
				{
				tree_top = build[i];
				build[i] = NULL;
				taxa = 2;
				}
			else
				{	
				if(taxa != 1)    /* if there is a single sequence at the end to be joined */
					{
					if(((build[i]->seq_num1 != NULL &&  build[i]->seq_num2 == NULL)  || (build[i]->seq_num1 == NULL && build[i]->seq_num2 != NULL)) && (build[i]->node1 == NULL && build[i]->node2 == NULL))
						{
						if(build[i]->seq_num1 != NULL)
							{
							build[i]->node2 = tree_top;   /* With only one taxa we can make this new node the tree top with the free noe pointing to the rest of the tree. */
							build[i]->prev = NULL;           
							tree_top->prev = build[i];
							tree_top = build[i];
							build[i] = NULL;
							}
						else
							{
							build[i]->node1 = tree_top;
							build[i]->prev = NULL;
							tree_top->prev = build[i];
							tree_top = build[i];
							build[i] = NULL;
							}
						}
					else
						{
						position = malloc(sizeof(node_type));
						if(!position)
							{
							printf("Out of memory\n");
							clean_exit();
							}
					
						position->seq_num1 = NULL;
						position->seq_num2 = NULL;
						position->node1 = build[i];
						position->node2 = tree_top;
						position->prev = NULL;
						position->ancestor[0] = '\0';
						position->ances_seq = malloc((((start->length)/3) +1) * sizeof(int));
						if(!position->ances_seq)
							{
							printf("\n Out of memory \n");
							clean_exit();
							}
						for(j=0; j<(((start->length)/3) +1); j++) position->ances_seq[j] = 0;
	
						position->gaprun = 0;
	

						build[i]->prev = position;
						tree_top->prev = position;
						
					
						build[i] = NULL;
						
						tree_top = position;
					
						position = NULL;
						}											
					}
				}
			}
		}
	

	if(done_col != NULL)
		{
		free(done_col);	
		done_col = NULL;
		}
	
	if(build != NULL)
		{
		for(i=0; i<(num_of_seqs - untagged - 2); i++) build[i] = NULL;
		free(build);
		build = NULL;
		}
	}	
	
	

/* I apoligise for this function in advance.... There is a problem in some platforms with the tree, in that for some reason the seq_num1 value
	is not assigned to NULL when the the two nodes are assigned, The problem I cannot find, which makes me think its an internal problem
	with the code being implemented on the particular environment. Since I cannot find anything to fix, I decided to come with an ad-hoc solution.
	This function will travel through the tree checking to see if there is a case of the two nodes being assigned, AND one of the seqNUM being assigned.
	it will then assign the pointer to NULL. By doing this I am assuming ALOT, namely that this is the only problem that can occur, since it could actually assign a 
	valid seqnum pointer to null, if for some reason one of the node pointers got assigned wrongly.
*/
int val_tree(struct node *position)
	{
	int result = FALSE, returned = FALSE;
	
	if(position->node1 != NULL) result = val_tree(position->node1);
	if(result == TRUE) returned = TRUE;
	if(position->node2 != NULL) result = val_tree(position->node2);
	if(result == TRUE) returned = TRUE;


	if(position->node1 != NULL && position->node2 != NULL)
		{
		if(position->seq_num1 != NULL)
			{
			position->seq_num1 = NULL;
			returned = TRUE;
			}
		if(position->seq_num2 != NULL)
			{
			position->seq_num2 = NULL;
			returned = TRUE;
			}
		}
	
	return(returned);
	
	}
			
				



/* This recursive function travels through the tree depth first to search for the first occurance of an outgroup
	which is present on a leaf */

struct node * find_outgroup (struct node *position)
	{
	struct node *now = NULL;


	if((position->seq_num1) != NULL)
		if((position->seq_num1)->outgroup == TRUE) now = position;
	if((position->seq_num2) != NULL)
		if((position->seq_num2)->outgroup == TRUE) now = position;
	
	/* we want to find the first part of the outgroup that we come across */	
	if(position->node1 != NULL && now == NULL) now = find_outgroup(position->node1);
	if(position->node2 != NULL && now == NULL) now = find_outgroup(position->node2);
	

	return(now);
	}
	

/* this goes down the tree from the position found and counts all the  outgroups found
	along with the total number of taxa found.  */ 
void travel_down(struct node *position, int *found, int *total)
	{
	
	if(position->seq_num1 != NULL)
		{
		if(outgroup[(position->seq_num1)->seq_num]) found++;
		total++;
		}
	if(position->seq_num2 != NULL)
		{
		if(outgroup[(position->seq_num2)->seq_num]) found++;
		total++;
		}
	
	if(position->node1 != NULL)
		travel_down(position->node1, found, total);
		
	if(position->node2 != NULL)
		travel_down(position->node2, found, total);
		
	}
		

	
/* this will check if the outgroup is divided on either side of the tree_top which makes it very difficult to
	check if the outgroup is valid. if this is divided the function will reroot the tree so that the
	outgroup is all on one side of the tree */	
int check_divided(struct node *position)
	{
	struct node *left = NULL, *right = NULL;
	int answer = FALSE;
	
	if(position->node1 != NULL)
		left = find_outgroup(position->node1);
		
	if(position->node2 != NULL)
		right = find_outgroup(position->node2);
		
	
	if(left != NULL && right != NULL)
		{
		if(left->node1 != NULL || left->node2 != NULL)
			{
			prepare_tree(left, 1);
			answer = check_divided(tree_top);
			}
		if(answer == FALSE)
			{
			if(right->node1 != NULL || right->node2 != NULL)
				{
				prepare_tree(right, 1);
				answer = check_divided(tree_top);
				}
			}
		}
	else answer = TRUE;
	
	
	return(answer);
	}

/* this function is called in the times the tree needs to be rerooted about a position */
/* It does nothing but blindly introduces a new node before the position specified, and
   then reroots the tree about that new node */
   
void do_rerooting(struct node *position)
	{		
	struct node *newnode = NULL;
	int j = 0;
	
	
	newnode = malloc(sizeof(node_type));
	if(!newnode)
		{
		printf("ERROR\nOut of memory\n");
		clean_exit();
		}

	newnode->seq_num1 = NULL;
	newnode->seq_num2 = NULL;
	newnode->node1 = NULL;
	newnode->node2 = NULL;
	newnode->prev = NULL;
	newnode->ancestor[0] = '\0';
	newnode->ances_seq = malloc(((start->length)/3) +1 * sizeof(int));
	if(!newnode->ances_seq)
		{
		printf("ERROR\nOut of memory\n");
		clean_exit();
		}
	else
		for(j=0; j<(((start->length)/3) +1); j++) position->ances_seq[j] = 0;

	/* put in the newnode before the present position */
			
	newnode->node2 = position;
	newnode->prev = position->prev;
	position->prev = newnode;
	if(newnode->prev != NULL)
		{
		if((newnode->prev)->node1 == position)
			(newnode->prev)->node1 = newnode;
		else (newnode->prev)->node2 = newnode;	
		}
	else tree_top = newnode;
	position = newnode;		
	/* reroot the tree about the present position */
	if(newnode != tree_top)	
		{		
		reroot_tree(position);
				
										
		if(position->node1 == NULL)
			position->node1 = position->prev;
	   	else
  		 	position->node2 = position->prev;
	 	position->prev = NULL;
		}
	}



/* Because we have to introduce a new node if we try to reroot the tree on a node with two daughter branches
	this function does this if necessary, and then uses this new node to reroot the tree using the function
	reroot tree. It then cleans up the process by assigning the pointers for the top of the tree */

int prepare_tree(struct node *position, int count)
	{
	struct node *newnode = NULL, *tmp = NULL;
	int j = 0, valid = TRUE;
	
	
	/***********  allocate new node *****************/
	newnode = malloc(sizeof(node_type));

	newnode->seq_num1 = NULL;
	newnode->seq_num2 = NULL;
	newnode->node1 = NULL;
	newnode->node2 = NULL;
	newnode->prev = NULL;
	newnode->ancestor[0] = '\0';
	newnode->ances_seq = malloc((((start->length)/3) +1) * sizeof(int));
	if(!newnode->ances_seq)
		{
		printf("ERROR\nOut of memory\n");
		clean_exit();
		}
	else
		for(j=0; j<(((start->length)/3 )+1); j++) position->ances_seq[j] = 0;
	
	
	/********** start preparing ****************/
	newnode->gaprun = 0;
	if(count == 1)  /* if there is only 1 taxa defined as the out group */
		{
		if(position->seq_num1 != NULL && position->seq_num2 != NULL)  /* if the outgroup has been paired with a taxa */
			{
			if((position->seq_num1)->outgroup) /* if seq_num1 points to the outgroup */
				{
				newnode->seq_num1 = position->seq_num1;
				position->seq_num1 = NULL;
				position->node1 = newnode;
				}
			else
				{
				newnode->seq_num2 = position->seq_num2;
				position->seq_num2 = NULL;
				position->node2 = newnode;
				}
			newnode->prev = position;
			position = newnode;
			reroot_tree(position);
			if(position->seq_num1 == NULL)
				position->node1 = position->prev;
   			else
   				position->node2 = position->prev;
  			position->prev = NULL;
			}
			
		else  /* if the outgroup is on its own */
			{
			newnode->node2 = position;
			newnode->prev = position->prev;
			position->prev = newnode;
			
			if((newnode->prev)->node1 == position)
				{
				(newnode->prev)->node1 = newnode;
				}
			else
				{
				(newnode->prev)->node2 = newnode;	
				}
				reroot_tree(newnode);
			if(newnode->node1 == NULL)
				newnode->node1 = newnode->prev;
   			else
   				newnode->node2 = newnode->prev;
  			newnode->prev = NULL;

			
			if(position->node1 != NULL)
				{
				tree_top->node2 = position->node1;
				(position->node1)->prev = tree_top;
				position->node1 = tree_top;
				tree_top->prev = position;
				tree_top = position;
				tree_top->prev = NULL;
				}
			else
				{
				tree_top->node2 = position->node2;
				(position->node2)->prev = tree_top;
				position->node2 = tree_top;
				tree_top->prev = position;
				tree_top = position;
				tree_top->prev = NULL;
				}
			}	
		}	
	else /* if there are more than one taxa defined as the out groups */
		{
		/* first check to see if the taxa are valid. */
		tmp = find_non_outgroup(position);
		if(tmp == NULL)
			{
			if(count_taxa(position) == count)
				{
				/* put in the newnode before the present position */
			
				newnode->node2 = position;
				newnode->prev = position->prev;
				position->prev = newnode;
				if((newnode->prev)->node1 == position)
					(newnode->prev)->node1 = newnode;
				else (newnode->prev)->node2 = newnode;	
			
				/* reroot the tree about the present position */
				reroot_tree(newnode);
										
				if(newnode->node1 == NULL)
					newnode->node1 = newnode->prev;
   				else
   					newnode->node2 = newnode->prev;
  				newnode->prev = NULL;

			
				}
			else
			 	valid = FALSE;
			}
		else
			valid = FALSE;
			
		}
	return(valid);
	}



 /* This simply counts how many taxa there are from the position on the tree that it is called from to the bottom */

 int count_taxa(struct node *position)
 	{
 	int i = 0;
 	
 	if(position->node1 != NULL) i += count_taxa(position->node1);
 	if(position->node2 != NULL) i += count_taxa(position->node2);
 	
 	if(position->seq_num1 != NULL) i++;
 	if(position->seq_num2 != NULL) i++;

	return(i);
	}
	



/* This function will reeroot the tree to the node specified by the input
	To complete the rerooting you need to add the following lines after whenver you call this function */
/* if(position->node1 == '\0')
	position->node1 = position->prev;
   else
   	position->node2 = position->prev;
   position->prev = '\0';
   
  without these lines the tree will not be fully rerooted */
void reroot_tree(struct node *position)
	{
	
	if((position->prev)->node1 == position)
		(position->prev)->node1 = (position->prev)->prev;
	if((position->prev)->node2 == position)
		(position->prev)->node2 = (position->prev)->prev;

	
	if((position->prev)->prev != NULL)
		{
		reroot_tree(position->prev);
		}

	(position->prev)->prev = position;
	tree_top = position;
	
	}
	
/* This function looks for occurances where we have rerooted a tree and in doing so we had to add another
	node. Since this node is then obsolete if we reroot the tree again, this function deletes that node, 
	which would contain no information for the rerooted tree */
int prune_tree(struct node *position)
	{
	struct node *place = NULL;
	int done = FALSE;
	
	if(position->seq_num1 == NULL && position->seq_num2 == NULL)
		{
		if(position->node1 == NULL || position->node2 == NULL)
			{
			if(position->node1 == NULL && position->node2 == NULL)   /* If the node is at the end of branch (leaf) */
				{
				if((position->prev)->node1 == position)
					(position->prev)->node1 = NULL;
				else (position->prev)->node2 = NULL;
				
				place = position;
				position = place->prev;
				if(place->ances_seq != NULL) free(place->ances_seq);
				place->ances_seq = NULL;
				if(place != NULL) free(place);
				place = NULL;
				done = TRUE;
				}
			else  /* if the node appears in the middle of a branch */
				{
				if(position->node1 != NULL)
					{
					place = position->node1;
					(position->node1)->prev = position->prev;
					}
				else
					{
					place = position->node2;
					(position->node2)->prev = position->prev;
					}
				
				if((position->prev)->node1 == position)
					(position->prev)->node1 = place;
				else (position->prev)->node2 = place;
				
				place = position;
				position = place->prev;
				if(place->ances_seq != NULL) free(place->ances_seq);
				place->ances_seq = NULL;
				if(place != NULL) free(place);
				place = NULL;
				done = TRUE;
				}
			}
		}
	/* Ths next section checks to see if there is a single taxa on its own at the bottom of a branch, 
		This breaks the rules about branchs, so the taxa is included in the node previous to it */
	if(position->node1 == NULL && position->node2 == NULL)
		if(position->seq_num1 == NULL || position->seq_num2 == NULL)
			{
			if((position->prev)->node1 == position)
				{
				if(position->seq_num1 != NULL)
					(position->prev)->seq_num1 = position->seq_num1;
				else (position->prev)->seq_num1 = position->seq_num2;
				(position->prev)->node1 = NULL;
				}
			else
				{
				if(position->seq_num1 != NULL)
					(position->prev)->seq_num2 = position->seq_num1;
				else (position->prev)->seq_num2 = position->seq_num2;
				(position->prev)->node2 = NULL;
				}
			
			place = position;
			position = place->prev;
			if(place->ances_seq != NULL) free(place->ances_seq);
			place->ances_seq = NULL;
			if(place != NULL) free(place);
			place = NULL;
			done = TRUE;
			}
			
	
	
	
	if(position->node1 != NULL)
		if(prune_tree(position->node1) == TRUE) done = TRUE;
		
	if(position->node2 != NULL)
		if(prune_tree(position->node2) == TRUE) done = TRUE;
		
	return(done);	
	}
	
/* This checks the validity of the selected outgroup */	
	
int check_outgroup(int count)
	{
	struct node *position = NULL;
	int answer = FALSE, proceed = TRUE;
		
	/* because of the nature of the storage of the tree, if the tree is originally rooted in the middle
	   of the outgroups, it will not be possible to verify the selection of the outgroups. The best way to
	   stop this from happening is to firstly make sure that the tree is rooted somewhere away from the out groups */
	
	
	if((position = find_non_outgroup(tree_top)) != NULL)
		{
		if(position != tree_top)
			prepare_tree(position, 1);
		
/*		number = 0;
		if(tree_top != '\0')
			{
			check_tree(tree_top, &number);  *//* make sure the tree has been created correctly */
/*			printf("\npress return to continue (rerooted about an ingroup)");
			getchar();	
			}
		else printf("tree not defined\n");
*/
		while(prune_tree(tree_top)){}
/*		
		number = 0;
		if(tree_top != '\0')
			{
			check_tree(tree_top, &number);  *//* make sure the tree has been created correctly */
/*			printf("\npress return to continue (rerooted about an ingroup)");
			getchar();	
			}
		else printf("tree not defined\n");
*/
		if(count == 1) /* If there is only 1 outgroup, then we need to make sure that the tree is not already rooted about it before proceeding */
			{
			if(tree_top->seq_num1 != NULL && tree_top->seq_num2 == NULL)
				{
				if(tree_top->seq_num1->outgroup == TRUE)
					proceed = FALSE;
					answer = TRUE;
				}
			if(tree_top->seq_num2 != NULL && tree_top->seq_num1 == NULL)
				{
				if(tree_top->seq_num2->outgroup == TRUE)
					proceed = FALSE;
					answer = TRUE;
				}
			}
				
		if(proceed)
			{
			if((position = find_outgroup(tree_top)) != NULL) /* Find a part of the out group */
				{
				if(prepare_tree(position, count) == TRUE)
					 answer = TRUE;		
				else answer = FALSE;
				}
				
			else
				{
				answer = FALSE;	
				printf("Cannot find an outgroup\n");
				}
			}
		}
	else
		{
		answer = FALSE;
		printf("Cannot find an ingroup\n");
		}
	return(answer);
	
	}
	
/* this function tries to find a non-outgroup taxa below where it is called from */
		
struct node * find_non_outgroup(struct node *position)
	{
	struct node *found = NULL;
	

	if(position->seq_num1 != NULL)
		{
		if(!(position->seq_num1)->outgroup)
			found = position;	
		}
	if(found == NULL && position->seq_num2 != NULL)
		{
		if(!(position->seq_num2)->outgroup)
			found = position;
		}
		
	if(found == NULL && position->node1 != NULL)	found = find_non_outgroup(position->node1);
	if(found == NULL && position->node2 != NULL)    found = find_non_outgroup(position->node2);
	
	return(found);
	
	}
	
	
/* This function dismantles the binary tree */

void dismantle(struct node *position, int *count)
	{
	
	struct node *place = NULL;
	
	if(position->node1 != NULL) dismantle(position->node1, count);
	if(position->node2 != NULL) dismantle(position->node2, count);
	
	if(position->ances_seq != NULL) free(position->ances_seq); 
	position->ances_seq = NULL;
	 
		
	place = position;
	position = place->prev;
	if(place != NULL) free(place);
	place = NULL;
	
	*count = *count + 1;
	}
	
/*
void write_graph(struct node *position, int *count, int i, int title)
	{
	
	if(position->node1 != '\0') write_graph(position->node1, count, i, title);
	if(position->node2 != '\0') write_graph(position->node2, count, i, title);

	if(position != tree_top)
		{

			if(title ) {
				
				
				fprintf(graphfile, "Node %d\t", *count);
				fprintf(graphfile, "repl/fix\t");
				fprintf(graphfile, "repl/poly\t");
				fprintf(graphfile, "synon/fix\t");
				fprintf(graphfile, "synon/poly\t");
				 }
				
		else  fprintf(graphfile, "\t%d", position->graph[0][i] );
		*count = *count + 1;
		}
	}
	

*/
	
	
/* this writes a tree file from the tree in memory, This outputs the tree with the results from the Maynooth methods 1 & 2 */

void write_tree1(struct node *position, char *last, int *count, float **pvalue)
	{
	int i = 0;
	
	if(*last == ')') fprintf(outtree, ",");
	fprintf(outtree, "(");
	if(position->seq_num1 != NULL)
		{
		fprintf(outtree, "%s,", (position->seq_num1)->name);
		i++;
		}
	if(position->seq_num2 != NULL)
		{
		fprintf(outtree, "%s", (position->seq_num2)->name);
		if(i == 0) fprintf(outtree, ",");
		}

	*last = ',';

	if(position->node1 != NULL) write_tree1(position->node1, last, count, pvalue);
	if(position->node2 != NULL) write_tree1(position->node2, last, count, pvalue);
	fprintf(outtree, ")");
	if(position->prev != NULL)
		{
		if(pvalue[*count][1] < 0.050 || (pvalue[*count][2] <= 0.05 && pvalue[*count][2] != 0)) fprintf(outtree, "%d ",*count);
		}
	*count = *count + 1;
	*last = ')';
	
	}

/* This writes a tree to the tree file, including the results from the pairwise analysis of the tree (incuding ancestors) */
void write_tree2(struct node *position, char *last, int *count)
	{
	int i = 0;
	float li = 0;
	
	if(*last == ')') fprintf(outtree, ",");
	fprintf(outtree, "(");
	if(position->seq_num1 != NULL)
		{
		fprintf(outtree, "%s,", (position->seq_num1)->name);
		i++;
		}
	if(position->seq_num2 != NULL)
		{
		fprintf(outtree, "%s", (position->seq_num2)->name);
		if(i == 0) fprintf(outtree, ",");
		}

	*last = ',';

	if(position->node1 != NULL) write_tree2(position->node1, last, count);
	if(position->node2 != NULL) write_tree2(position->node2, last, count);
	fprintf(outtree, ")");
	if(position->prev != NULL)
		{
		li = (position->li93_1[0]/position->li93_1[1]);
		if(position->li93_1[1] != 0 && li > 1 ) fprintf(outtree, "1_%f", li);
		li = (position->li93_2[0]/position->li93_2[1]);
		if(position->li93_2[1] != 0 && li > 1 ) fprintf(outtree, "_2_%f", li);
		}
	*count = *count + 1;
	*last = ')';
	
	}

/* this writes a tree file from the tree in memory, Every internal branch is labeled so we know the number of all the branches */

void write_tree3(struct node *position, char *last, int *count, float **pvalue)
	{
	int i = 0;
	
	if(*last == ')') fprintf(outtree, ",");
	fprintf(outtree, "(");
	if(position->seq_num1 != NULL)
		{
		fprintf(outtree, "%s,", (position->seq_num1)->name);
		i++;
		}
	if(position->seq_num2 != NULL)
		{
		fprintf(outtree, "%s", (position->seq_num2)->name);
		if(i == 0) fprintf(outtree, ",");
		}

	*last = ',';

	if(position->node1 != NULL) write_tree3(position->node1, last, count, pvalue);
	if(position->node2 != NULL) write_tree3(position->node2, last, count, pvalue);
	fprintf(outtree, ")");
	if(position->prev != NULL)
		{
		fprintf(outtree, "%d ",*count);
		}
	*count = *count + 1;
	*last = ')';
	
	}

	
/* This checks to see if the treefile is properly formatted */

int check_treefile(void)
	{
	char c = '\0';
	int left = 0, right = 0, error = FALSE;

	c = getc(parenthesis);
	while(c == ' ' || c == '\r' || c == '\n' || c =='\t') c = getc(parenthesis);
	if(c != '(') error = TRUE;
	else
		{
		left++;
		while((c = getc(parenthesis)) != ';' && !feof(parenthesis) && left >= right)
			{
			switch(c)
				{
				case '(':
					left++;
					break;
				case ')':
					right++;
					break;
				default:
					break;
				}
			}
		if(c != ';' || left != right) error = TRUE;
		
		}
	rewind(parenthesis);
	return(error);
	}


/* This function implements the heirarcial structure of a tree which is given in nested parenthesis format in 
	a given file */
	
int input_tree(char *c, struct node *previous, int num)
	{
	char temp[13] = {'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0'};
	int i = 0, end = FALSE, error = FALSE, count = 0, score = 0, highest = 0;
	struct sequence *seq = start, *same = NULL;
	struct node *extra = NULL, *position  = NULL;

	position = malloc(sizeof(node_type));
	if(previous == NULL) tree_top = position;
				
	position->seq_num1 = NULL;
	position->seq_num2 = NULL;
	position->node1 = NULL;
	position->node2 = NULL;
	position->prev = previous;
	position->ancestor[0] = '\0';
	position->ances_seq = malloc((((start->length)/3 )+1) * sizeof(int));
	if(!position->ances_seq)
		{
		printf("ERROR: Out of memory\n");
		clean_exit();
		}
	for(i=0; i<((start->length)/3) +1; i++) position->ances_seq[i] = 0;

	position->gaprun = 0;
	

	
	if(previous != NULL)
		{
		if(num == 0)
			(position->prev)->node1 = position;
		else
                        (position->prev)->node2 = position;
		}
	previous = position;	
	*c = getc(parenthesis);
	
	while(!error && *c != ')' && *c != ';')
		{
	
		switch(*c) 
			{
			case '(':
                                if(count > 1)
                                    {
                                    extra = malloc(sizeof(node_type));
				
                                    extra->seq_num1 = NULL;
                                    extra->seq_num2 = NULL;
                                    extra->node1 = NULL;
                                    extra->node2 = NULL;
                                    extra->prev = NULL;
                                    extra->ancestor[0] = '\0';
                                    extra->ances_seq = malloc((((start->length)/3) +1 )* sizeof(int));
                                    if(!extra->ances_seq)
                                            {
                                            printf("ERROR: Out of memory\n");
                                            clean_exit();
                                            }
                                    for(i=0; i<((start->length)/3) +1; i++) extra->ances_seq[i] = 0;
						
                                    extra->gaprun = 0;
						

                                    extra->node1 = position;
                                    extra->prev = position->prev;
                                    position->prev = extra;
                                    if(extra->prev == NULL) tree_top = extra;
                                    position = extra;
                                    previous = position;
                                    }

				error = input_tree( c, previous, count);
				count++;
				*c = getc(parenthesis);
				break;
		
			case ',':
				*c = getc(parenthesis);
				break;
			
			case ')':
				break;
			
			case ';':
				break;
		
			default:
				for(i=0; i<13; i++) temp[i] = '\0';
				i = 0;
				end = FALSE;
				do{
					if(*c != ' ' && !end && i < 12)
						{
						temp[i] = *c;
						i++;
						}
					else end = TRUE;
				
					}while((*c = getc(parenthesis)) != ')' && *c != ',' && *c != '(');
				
				/* This next section tries to identify which of the sequences is being refered to by the name on the tree at this position.
				   This is done by giving a similarity score to every sequence name and the largest score is assigned to the tree.
				   This works by adding 1 to the total for every character in the right position in the sequence name */ 		
				seq = start;
				same = NULL;
				highest = 0;
				while(seq != NULL)
					{
					score = 0;
					for(i=0; i<13; i++)
						{
						if(temp[i] == seq->nickname[i]) score++;
						if(temp[i] == '\0') i = 13;
						}
					if(score > highest)
						{
						same = seq;
						highest = score;
						}
					seq = seq->next;
					}
				
				seq = same;
				if(seq == NULL) error = TRUE;
				else
					{			
					if(count == 0)	position->seq_num1 = seq;
					if(count == 1) position->seq_num2 = seq;
					if(count > 1)
						{
						extra = malloc(sizeof(node_type));
				
						extra->seq_num1 = NULL;
						extra->seq_num2 = NULL;
						extra->node1 = NULL;
						extra->node2 = NULL;
						extra->prev = NULL;
						extra->ancestor[0] = '\0';
						extra->ances_seq = malloc((((start->length)/3) +1 )* sizeof(int));
						if(!extra->ances_seq)
							{
							printf("ERROR: Out of memory\n");
							clean_exit();
							}
						for(i=0; i<((start->length)/3) +1; i++) extra->ances_seq[i] = 0;
						
						extra->gaprun = 0;
						

						extra->seq_num1 = seq;
						extra->node2 = position;
                                                extra->prev = position->prev;
						position->prev = extra;
                                                
						if(extra->prev == NULL) tree_top = extra;
                                                position = extra;
                                                previous = position;
						}

					count++;
					}
				
				break;				
		
			}
		
		}
		
	return(error);
	
	}
	
	
void tree_pairwise_distances(struct node *position)
	{
	
	if(position->node1 != NULL) tree_pairwise_distances(position->node1);
	if(position->node2 != NULL) tree_pairwise_distances(position->node2);
	
	
	if(position->prev != NULL)
		{
		if(position->node1 != NULL)
			{
			li_wu_complete(position->ances_seq, (position->node1)->ances_seq, position->li93_1);
			
			}
		if(position->node2 != NULL)
			{
			li_wu_complete(position->ances_seq, (position->node2)->ances_seq, position->li93_2);
			
			}
		if(position->seq_num1 != NULL)
			{
			li_wu_complete((position->seq_num1)->bases, position->ances_seq, position->li93_1);
			
			}
		if(position->seq_num2 != NULL)	
			{
			li_wu_complete((position->seq_num2)->bases, position->ances_seq, position->li93_2);
			
			}
		}
		
	}



/* this goes through the tree, assigning the node numbers */	
int assign_node_nums(struct node * position, int num)
	{
	
	if(position->node1 != NULL)
		num = assign_node_nums(position->node1, num);
	if(position->node2 != NULL)
		num = assign_node_nums(position->node2, num);
		
	position->nodenum = num;
	num++;
	
	return(num);
	}


/* prints the pairwise comparisons of on the tree to the yadf file */
void print_tree_pair_dist(struct node * position)
	{
	
	if(position->node1 != NULL)
		print_tree_pair_dist(position->node1);
	if(position->node2 != NULL)
		print_tree_pair_dist(position->node2);
		
		
	if(position->seq_num1 !=NULL)
		fprintf(yadf, "node%d to %s\t",  position->nodenum, (position->seq_num1)->name);
	else
		fprintf(yadf, "node%d to node%d\t", position->nodenum, (position->node1)->nodenum ); 
	
	
	fprintf(yadf, "%f\t%f\t%f\t%f\t%f\n", position->li93_1[0], position->li93_1[1], (position->li93_1[0]/position->li93_1[1]), position->li93_1[2], position->li93_1[3]);
	
	if(position->seq_num2 !=NULL)
		fprintf(yadf, "node%d to %s\t",  position->nodenum, (position->seq_num2)->name);
	else
		fprintf(yadf, "node%d to node%d\t", position->nodenum, (position->node2)->nodenum ); 
	
	
	fprintf(yadf, "%f\t%f\t%f\t%f\t%f\n", position->li93_2[0], position->li93_2[1], (position->li93_2[0]/position->li93_2[1]), position->li93_2[2], position->li93_2[3]);


	}





/* This function travels down into the tree and at each internal node calcaulates the number
	of Replacement and silent sites in the clade defined by that internal branch.
	Each of the sequences, whether they are reconstructed ancestral sequences or real taxa are evaluated.
	At each nucleotide position, we determine whether it is a 0, 2 or 4-fold degenerate site (the matrices for this
	have already been defined in definitions.h. A 0-fold site means a substitution here is alway a replplacement and
	 so one of these sites gets 1 added to the replacement count.
	A 2-fold site means that half the time its a replacement, and the other half its silent, so one of these means that the 
	silent count gets 0.5 added and the replacement count gets 0.5 added.
	A 4-fold site means that every substitution is silent and so one of these sites gets 1 added to the silent count.
	This will give us a count of the number of opportunities for change per replacement or silent site within the clade defined
	by the branch we are examining.
	This function is the first part which travels down the tree totalling for each branch the nuber of replacements and silents
	but to count how many there are in the clade we call count_subs which travels down through the clade counting the sites. 
	<<< Think of it like a nested for loop, subs_include is the outer loop which calls count_subs - the inner loop-.
	These functions are recursive and optimise the data structures implemented in the software >>>
	*/
void subs_inclade(struct node * position, int *count, float **subs)
	{
	float replacements = 0, silents = 0;
	
	if(position->node1 != NULL && position->seq_num1 == NULL) 
		subs_inclade(position->node1, count, subs);
	
	if(position->node2 != NULL && position->seq_num2 == NULL)
		subs_inclade(position->node2, count, subs);
		
	
	count_subs(position, &replacements, &silents);
	
	subs[*count][0] += replacements;
	
	subs[*count][1] += silents;
	
	*count = *count + 1;
	
	}
	





/* This is called by subs_include. For a description see the comments preceding that function */
void count_subs(struct node * position, float *replacements, float *silents)
	{
	int i=0, j =0;
	
	
	
	if(position->node1 != NULL && position->seq_num1 == NULL) /* If node_1 points to a daughter node */
		count_subs(position->node1, replacements, silents);
	else   /* count the sites in the first sequence */
		{
		i=0;
		for(i=0; i<(start->length)/3; i++)
			{
			for(j=0; j<3; j++)
				{
				switch(degenerate_sites[code][(position->seq_num1)->bases[i]][j])
					{
					case 0:
						*replacements = *replacements +1;
						break;
					case 2:
						*replacements = *replacements + (0.5);
						*silents = *silents + (0.5);
						break;
					default:
						*silents = *silents + 1;
						break;
					}
				}
			}
		}
			
			
			
			
			
	
	if(position->node2 != NULL) /* if node_2 points to a daughter node */
		count_subs(position->node2, replacements, silents);
	else /* count the sites in the second sequence */
		{
		i=0;
		for(i=0; i<(start->length)/3; i++)
			{
			for(j=0; j<3; j++)
				{
				switch(degenerate_sites[code][(position->seq_num2)->bases[i]][j])
					{
					case 0:
						*replacements = *replacements +1;
						break;
					case 2:
						*replacements = *replacements + (1/2);
						*silents = *silents + (1/2);
						break;
					default:
						*silents = *silents + 1;
						break;
					}
				}
			}
		}
	
	
	
	
	/* Lastly count the number of sites in the ancestral sequence */
	i=0;
	for(i=0; i<(start->length)/3; i++)
		{
		for(j=0; j<3; j++)
			{
			switch(degenerate_sites[code][position->ances_seq[i]][j])
				{
				case 0:
					*replacements = *replacements +1;
					break;
				case 2:
					*replacements = *replacements + (1/2);
					*silents = *silents + (1/2);
					break;
				default:
					*silents = *silents + 1;
					break;
				}
			}
		}
		
		
		
	}
		
		
