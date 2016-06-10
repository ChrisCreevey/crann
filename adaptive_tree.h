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

/* found in adaptive_tree.c*/
void adaptive_tree(int **ratio);
void calculate_gaprun(struct node *position, int i);
void reset_gaprun(struct node *position);
void increment_gaprun(struct node *position);
int is_gap(int i, int j, struct node *position);
void close_ances(struct node *position, int i);
void output_ancestors(struct node *position, int *count);
void is_mutation(struct node *position,int *mutation, char *nuc, int i, int j);
void is_fixed(struct node *position, int i, int j, int **ratio, int *count, int last);
void snapshot(struct node *position, int **ratio, char *nuc, int *count, int i, int j);
void count_polymorphisms(struct node *position, int *repl, int *silent, int i, int j);
void find_polymorphisms(struct node *position, int **ratio, int **ratio1, int *count, int i, int j);
void check_fixed_ances_out(struct node *position, int *fixed, char *nuc, int i, int j);
void check_fixed_in(struct node *position, int *fixed, char *nuc, int i, int j);
void check_fixed_out(struct node *position, struct node *place, int *fixed, char *nuc, int i, int j	);
void check_tree(struct node *position, int *count);
void ancestral_nuc(int i, int j, float *** subst_matrix);
void assign_codon_num(int i, int j, struct node *position);
void assign_ances_up(int i, int j, struct node *position);
void assign_root_ances(int i, int j);
int look_for_ances(struct node *position);
void assign_ances_down(struct node *position, float *** subst_matrix);
int check_distances(void);
void tally_distances(void);
void allocate_distances(int k);
void McDonald_Kreitman(void);
int tree_choice(void);
void n_joining_tree(float **tree);
void g_test(float **ratio, float **pvalue, float * gChi );
float logE(float value);
float fctrl(float value);
void output_tree(float **tree);
void substitution_matrix(float ***ratio);
void define_outgroup(void);
void assign_codon_up(int codon, struct node *position);
void assign_root_codon(int codon);
void assign_codons_down(int codon, struct node *position, float ***subst_matrix);
void assign_codon(int codon, struct node *position);
void ancestral_codon(int codon, float ***subst_matrix);
