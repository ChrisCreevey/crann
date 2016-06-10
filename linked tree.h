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

/* found in linked tree.c  */
void linked_tree(float **tree_array);
int val_tree(struct node *position);
struct node * find_outgroup (struct node *position);
void travel_down(struct node *position, int *found, int *total);
int check_divided(struct node *position);
void do_rerooting(struct node *position);
int prepare_tree(struct node *position, int count);
int count_taxa(struct node *position);
void reroot_tree(struct node *position);
int prune_tree(struct node *position);
int check_outgroup(int count);
struct node * find_non_outgroup(struct node *position);
void dismantle(struct node *position, int *count);
void write_tree1(struct node *position, char *last, int *count, float **pvalue);
void write_tree2(struct node *position, char *last, int *count);
void write_tree3(struct node *position, char *last, int *count, float **pvalue);
int check_treefile(void);
int input_tree(char *c, struct node *previous, int num);
void tree_pairwise_distances(struct node *position);
int assign_node_nums(struct node * position, int num);
void print_tree_pair_dist(struct node * position);
void subs_inclade(struct node * position, int *count, float **subs);
void count_subs(struct node * position, float *replacements, float *silents);
