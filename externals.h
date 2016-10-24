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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
/*#include <console.h>*/

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif


extern const int maxnamlen;

extern struct sequence{

	int *bases;					/* The sequence which is stored in dynamically allocted memory  */
	char name[61];	/* the name of the sequence */
	char nickname[20];			/* the shortened name of the sequence */
	int seq_num;				/* the number of the sequence in the file */
	int tag;					/* A tag for the seguence */
	int outgroup;				/* This tells whether the sequence is part of the outgroup or not */
	int length;					/* the length of the sequence */
	int numofstpcodons;			/* A number representing how many stop codons were in the sequence */
	int *stopcodons;			/* A dynamically allocated array which will store the posotions of stop codons in the sequence */ 
	int gaprun;					/* This keeps track of how many gaps there are in a row when we are going through the sequence in the McD&K alorithm */
	struct sequence *next;		/* pointer to the next sequence */
	struct sequence *previous;	/* pointer to the previous sequence */

	} list_entry;

extern struct sequence *start, *last;

extern struct synon{
	
	int seq_num;  			/* The number of the sequence that all the rest of the sequences were compared to */
	float *Ks;				/* The result of Ks (synonymous) against the other sequences */
	float *varKs;				/* The result of varKs (variance of Ks ) against the other sequences */
	float *Ka;				/* The result of Ka (non-synonymous) against the other sequences */
	float *varKa;				/* The result of varKa (variance of Ka ) against the other sequences */
	struct synon *next;		/* pointer to the results against the next sequence */
	struct synon *previous;	/* pointer to the results against the previous sequence */	
	
	} li_wu_result; 
	
extern struct synon *li_wu_start; /* Pointer to the first set of results */
extern struct synon *li_wu_end;  /* pointer to the last set of results */

extern struct node{

	struct sequence *seq_num1;			/* holds the number  of the first sequence in included in this node (if any) */
	struct sequence *seq_num2;			/* holds the number  of the second sequence in included in this node (if any) */
	int *ances_seq;			/* This holds the actual ancestral sequence as its being calculated */
	char ancestor[5];		/* holds the ancestor of the node */
	int codon_ances[66];    /* Used to calculate teh ancestral codon state in an attempt to make the reconstruction more reliable */
	struct node *node1;		/* points to the first node joined here (if any) */
	struct node *node2;		/* points to the second node joined here (if any) */
	struct node *prev;		/* points to the previous node always except if its the top of the tree */
	float li93_1[4];          /* holds results of the li_wu analysis of the left branch based on the tree - 0 = Ka, 1 = Ks, 2 = varKa, 3 = varKs */
	float li93_2[4];          /* holds results of the li_wu analysis of the right branch based on the tree - 0 = Ka, 1 = Ks, 2 = varKa, 3 = varKs */
	int gaprun;
	int nodenum;
	
	} node_type;
	

extern struct node *tree_top;     /* pointer to the top of the tree */


extern FILE *file, *outfile, *dist, *parenthesis, *ances_file, *outtree, *graphfile, *yadf;
extern char filename[36], outfilename[36], nestname[36], string1[1000], string2[1000];
extern int code, num_of_seqs, untagged, distance_written, startw, endw, gen_opt[7], ***graphs; /* distance_written is used to tell if the distance.out file was already written */
extern float **distances; /* used to store distances when computing a neighbour joining tree */

extern int transform_values [5][3];
extern char amino_acids[22][3];
extern char codons[65][3];
extern char what[5];
extern int genetic_codes[13][65];
extern int degenerate_sites[13][65][3];
extern char nonstandchars[100]; 
extern int *deletion;
extern int *outgroup;

