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

#define STD_CODON_NUM 1000			/* standard length of the gene in codons */
#define maxnamlen 1000			/* maximum length of the gene name */
#define progname "Evolve 0.6 by Chris Creevey"

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
/*#include <console.h> */ /* Uncomment for mac */

struct sequence{

	int *bases;					/* The sequence which is stored in dynamically allocted memory  */
	char name[maxnamlen + 1 ];	/* the name of the sequence */
	char nickname[20];			/* the shortened name of the sequence */
	int seq_num;				/* the number of the sequence in the file */
	int tag;					/* A tag for the sequence */
	int outgroup;				/* This tells whether the sequence is part of the outgroup or not */
	int length;					/* the length of the sequence */
	int numofstpcodons;			/* A number representing how many stop codons were in the sequence */
	int *stopcodons;			/* A dynamically allocated array which will store the posotions of stop codons in the sequence */ 
	int gaprun;					/* This keeps track of how many gaps there are in a row when we are going through the sequence in the McD&K alorithm */
	struct sequence *next;		/* pointer to the next sequence */
	struct sequence *previous;	/* pointer to the previous sequence */

	} list_entry;


struct sequence *start = NULL;  /*pointer to first sequence in the list */
struct sequence *last = NULL; /*point to last sequece in the list */


struct synon{
	
	int seq_num;  			/* The number of the sequence that all the rest of the sequences were compared to */
	float *Ks;				/* The result of Ks (synonymous) against the other sequences */
	float *varKs;				/* The result of varKs (variance of Ks ) against the other sequences */
	float *Ka;				/* The result of Ka (non-synonymous) against the other sequences */
	float *varKa;				/* The result of varKa (variance of Ka ) against the other sequences */
	struct synon *next;		/* pointer to the results against the next sequence */
	struct synon *previous;	/* pointer to the results against the previous sequence */	
	
	} li_wu_result; 
	
struct synon *li_wu_start = NULL; /* Pointer to the first set of results */
struct synon *li_wu_end = NULL;  /* pointer to the last set of results */


struct node{

	struct sequence *seq_num1;			/* Points to the sequence in the linked list in included in this node (if any) */
	struct sequence *seq_num2;			/* Points to the the sequence in the linked list in included in this node (if any) */
	int *ances_seq;			/* Holds the full ancestral sequence as it is calculated for each node */
	char ancestor[5];		/* holds the ancestor of the node */
	int codon_ances[66];    /* Used to calculate teh ancestral codon state in an attempt to make the reconstruction more reliable */
	struct node *node1;		/* points to the first node joined here (if any) */
	struct node *node2;		/* points to the second node joined here (if any) */
	struct node *prev;		/* points to the previous node always except if its the top of the tree */
	float li93_1[4];          /* holds results of the li_wu analysis of the left branch based on the tree - 0 = Ka, 1 = Ks, 2 = varKa, 3 = varKs */
	float li93_2[4];          /* holds results of the li_wu analysis of the right branch based on the tree - 0 = Ka, 1 = Ks, 2 = varKa, 3 = varKs */
	int gaprun;				/* Keeps track of runs of gaps so that they are not counted as individual polymorphisms. */
	int nodenum;
	} node_type;
	
struct node *tree_top = NULL;     /* pointer to the top of the tree */


FILE *file = NULL, *outfile = NULL, *dist = NULL, *parenthesis = NULL, *ances_file = NULL, *outtree = NULL, *usagefile = NULL, *graphfile = NULL, *yadf = NULL;
char filename[1000], outfilename[1000], nestname[1000], string1[1000], string2[1000];
int code = 0, num_of_seqs = 0, untagged = 0, distance_written = FALSE, startw = 0, endw = 0, gen_opt[7] = {0, 0, 1, 1, 0, 0, 0}, ***graphs = NULL; /* distance_written is used to tell if the distance.out file was already written */
float **distances = NULL; /* used to store distances when computing a neighbour joining tree */


/* This array represents the values we give the nucliotides based on their position in the codon */
/* The code will add up the three values for a given codon. Going downwards, the array           */
/* represents 't/u', 'c', 'a', 'g', other (gap values). The values are ordered so that the codons*/
/* differing only in third position are numerically close to each other.						 */ 
int transform_values [5][3]  = {0,0,0,   
								16,4,1,
								32,8,2,
								48,12,3,
								64,64,64};
								


/* This array represents the amino acids as they are defined in the genetic_codes array. The Amino Acids are numbered from 0 to 21 in the following order; */
/* Stop = X, Phe = F, Trp = W, Tyr = Y, His = H, Met = M, Leu = L, Ile = I, Val = V, Pro = P, Cys = C, Ala = A, Gly = G, Thr = T, Ser = S, Gln = Q, Asn = N, Lys = K, Arg = R, Glu = E, Asp = D, Gap = -.										   */								
char amino_acids[22] = {
  'X', 'F', 'W', 'Y', 'H', 'M', 'L', 'I', 'V', 'P', 'C', 'A', 'G', 'T', 'S', 'Q', 'N', 'K', 'R', 'E', 'D', '-'};

/* this array contains the amino acid categories as defined in "Prediction of protein secondary structure and active sites using the alignment of homologous sequences",
    M.J> Zvelebil, G.J. Barton, W.R. Taylor, and Sternberg, (J. Mol. Biol., 195, p957 - 1987)
    Wherever a particular amino acid belongs to a category, there is a 1, otherwise a 0 (or TRUE and FALSE), The order of the amino acids is the same as the array
    amino_acids and so stop and gaps are included even though they don't belong to any category, this is to ensure compatibility
    There are nine categories listed and are in the following order:
    hydrophobic, positive, negative, polar, charged, small, tiny, aromatic, aliphatic */
int AA_categories[22][9] = {
	0,0,0,0,0,0,0,0,0,  /*X */
	1,0,0,0,0,0,0,1,0,	/*F */
	1,0,0,1,0,0,0,1,0,	/*W */
	1,0,0,1,0,0,0,1,0,	/*Y */
	1,1,0,1,1,0,0,1,0,	/*H */
	1,0,0,0,0,0,0,0,0,	/*M */
	1,0,0,0,0,0,0,0,1,	/*L */
	1,0,0,0,0,0,0,0,1,	/*I */
	1,0,0,0,0,1,0,0,1,	/*V */
	0,0,0,0,0,1,0,0,0,	/*P */
	1,0,0,0,0,1,0,0,0,	/*C */
	1,0,0,0,0,1,1,0,0,	/*A */
	1,0,0,0,0,1,1,0,0,	/*G */
	1,0,0,1,0,1,0,0,0,	/*T */
	0,0,0,1,0,1,1,0,0,	/*S */
	0,0,0,1,0,0,0,0,0,	/*Q */
	0,0,0,1,0,1,0,0,0,	/*N */
	1,1,0,1,1,0,0,0,0,	/*K */
	0,1,0,1,1,0,0,0,0,	/*R */ 
	0,0,1,1,1,0,0,0,0,	/*E */
	0,0,1,1,1,1,0,0,0,	/*D */
	0,0,0,0,0,0,0,0,0  }; /*- */
    


/* This array gives the actual make up of the codons in nucliotides, compared to their codon number, as they are stored in memory */
char codons[65][3] = {
'U','U','U', 'U','U','C', 'U','U','A', 'U','U','G', 'U','C','U', 'U','C','C', 'U','C','A', 'U','C','G', 'U','A','U', 'U','A','C', 'U','A','A', 'U','A','G', 'U','G','U', 'U','G','C', 'U','G','A', 'U','G','G',
'C','U','U', 'C','U','C', 'C','U','A', 'C','U','G', 'C','C','U', 'C','C','C', 'C','C','A', 'C','C','G', 'C','A','U', 'C','A','C', 'C','A','A', 'C','A','G', 'C','G','U', 'C','G','C', 'C','G','A', 'C','G','G',
'A','U','U', 'A','U','C', 'A','U','A', 'A','U','G', 'A','C','U', 'A','C','C', 'A','C','A', 'A','C','G', 'A','A','U', 'A','A','C', 'A','A','A', 'A','A','G', 'A','G','U', 'A','G','C', 'A','G','A', 'A','G','G',
'G','U','U', 'G','U','C', 'G','U','A', 'G','U','G', 'G','C','U', 'G','C','C', 'G','C','A', 'G','C','G', 'G','A','U', 'G','A','C', 'G','A','A', 'G','A','G', 'G','G','U', 'G','G','C', 'G','G','A', 'G','G','G',
'X','X','X' };

 char what[5] = {'U', 'C', 'A', 'G', 'X'}; 		     


/* this array represents the various standard genetic codes. The numbers refer to the array amino_acids, which contains the three letter abreviations of all the amino acids */
int genetic_codes[13][65] =
  {1,1,6,6,14,14,14,14,3,3,0,0,10,10,0,2,6,6,6,6,9,9,9,9,4,4,15,15,18,18,18,18,7,7,7,5,13,13,13,13,16,16,17,17,14,14,18,18,8,8,8,8,11,11,11,11,20,20,19,19,12,12,12,12, 21,		/* Universal genetic code */	
  1,1,6,6,14,14,14,14,3,3,0,0,10,10,2,2,6,6,6,6,9,9,9,9,4,4,15,15,18,18,18,18,7,7,5,5,13,13,13,13,16,16,17,17,14,14,0,0,8,8,8,8,11,11,11,11,20,20,19,19,12,12,12,12, 21,		/* Vertebrate standard mitochondrial genetic code */
  1,1,6,6,14,14,14,14,3,3,0,0,10,10,2,2,13,13,13,13,9,9,9,9,4,4,15,15,18,18,18,18,7,7,5,5,13,13,13,13,16,16,17,17,14,14,18,18,8,8,8,8,11,11,11,11,20,20,19,19,12,12,12,12, 21,	/* Yeast mitochondrial genetic code */
  1,1,6,6,14,14,14,14,3,3,0,0,10,10,2,2,6,6,6,6,9,9,9,9,4,4,15,15,18,18,18,18,7,7,7,5,13,13,13,13,16,16,17,17,14,14,18,18,8,8,8,8,11,11,11,11,20,20,19,19,12,12,12,12, 21,		/* Mycoplasma/Spiroplasma/Mold/Protozoan/Coelenterate genetic code */
  1,1,6,6,14,14,14,14,3,3,0,0,10,10,2,2,6,6,6,6,9,9,9,9,4,4,15,15,18,18,18,18,7,7,5,5,13,13,13,13,16,16,17,17,14,14,14,14,8,8,8,8,11,11,11,11,20,20,19,19,12,12,12,12, 21,		/* Invertebrate mitochondrial genetic code */
  1,1,6,6,14,14,14,14,3,3,15,15,10,10,0,2,6,6,6,6,9,9,9,9,4,4,15,15,18,18,18,18,7,7,7,5,13,13,13,13,16,16,17,17,14,14,18,18,8,8,8,8,11,11,11,11,20,20,19,19,12,12,12,12, 21,	/* Ciliate genetic code */
  1,1,6,6,14,14,14,14,3,3,0,0,10,10,2,2,6,6,6,6,9,9,9,9,4,4,15,15,18,18,18,18,7,7,7,5,13,13,13,13,16,16,16,17,14,14,14,14,8,8,8,8,11,11,11,11,20,20,19,19,12,12,12,12, 21,		/* Echinoderm mitochondrial genetic code */
  1,1,6,6,14,14,14,14,3,3,0,0,10,10,10,2,6,6,6,6,9,9,9,9,4,4,15,15,18,18,18,18,7,7,7,5,13,13,13,13,16,16,17,17,14,14,18,18,8,8,8,8,11,11,11,11,20,20,19,19,12,12,12,12, 21,		/* Euplotid genetic code */
  1,1,6,6,14,14,14,14,3,3,0,0,10,10,0,2,6,6,6,6,9,9,9,9,4,4,15,15,18,18,18,18,7,7,7,5,13,13,13,13,16,16,17,17,14,14,18,18,8,8,8,8,11,11,11,11,20,20,19,19,12,12,12,12, 21,		/* Bacterial genetic code (same as universal) */
  1,1,6,6,14,14,14,14,3,3,0,0,10,10,0,2,6,6,6,14,9,9,9,9,4,4,15,15,18,18,18,18,7,7,7,5,13,13,13,13,16,16,17,17,14,14,18,18,8,8,8,8,11,11,11,11,20,20,19,19,12,12,12,12, 21,		/* Alternative Yeast Nuclear genetic code */
  1,1,6,6,14,14,14,14,3,3,0,0,10,10,2,2,6,6,6,6,9,9,9,9,4,4,15,15,18,18,18,18,7,7,5,5,13,13,13,13,16,16,17,17,14,14,12,12,8,8,8,8,11,11,11,11,20,20,19,19,12,12,12,12, 21,		/* Ascidian mitochondrial genetic code */
  1,1,6,6,14,14,14,14,3,3,3,0,10,10,2,2,6,6,6,6,9,9,9,9,4,4,15,15,18,18,18,18,7,7,7,5,13,13,13,13,16,16,16,17,14,14,14,14,8,8,8,8,11,11,11,11,20,20,19,19,12,12,12,12, 21,		/* Flatworm mitochondrial genetic code */
  1,1,6,6,14,14,14,14,3,3,0,15,10,10,0,2,6,6,6,6,9,9,9,9,4,4,15,15,18,18,18,18,7,7,7,5,13,13,13,13,16,16,17,17,14,14,18,18,8,8,8,8,11,11,11,11,20,20,19,19,12,12,12,12, 21};	/* Blepharisma Nuclear genetic code */




/* These are the calculated degenerate sites for each of the genetic codes. Every position of a stop codon was considered a 0-fold site, as was every position of a gap */
int degenerate_sites[13][65][3] = 
{0,0,2, 0,0,2, 2,0,2, 2,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,0, 0,0,0, 0,0,2, 0,0,2, 0,0,0, 0,0,0, 0,0,4, 0,0,4, 2,0,4, 2,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 2,0,4, 2,0,4, 
 0,0,2, 0,0,2, 0,0,2, 0,0,0, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 2,0,4, 2,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,0,

 0,0,2, 0,0,2, 2,0,2, 2,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,0, 0,0,0, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 2,0,4, 2,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 
0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,0, 0,0,0, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,0,

0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,0, 0,0,0, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 2,0,4, 2,0,4, 
0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 2,0,2, 2,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,0,

0,0,2, 0,0,2, 2,0,2, 2,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,0, 0,0,0, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 2,0,4, 2,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 2,0,4, 2,0,4, 
0,0,2, 0,0,2, 0,0,2, 0,0,0, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 2,0,2, 2,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,0,

0,0,2, 0,0,2, 2,0,2, 2,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,0, 0,0,0, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 2,0,4, 2,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 
0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,0,

0,0,2, 0,0,2, 2,0,2, 2,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 2,0,2, 2,0,2, 0,0,2, 0,0,2, 0,0,0, 0,0,0, 0,0,4, 0,0,4, 2,0,4, 2,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 2,0,2, 2,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4,
0,0,2, 0,0,2, 0,0,2, 0,0,0, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,0,

0,0,2, 0,0,2, 2,0,2, 2,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,0, 0,0,0, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 2,0,4, 2,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 
0,0,2, 0,0,2, 0,0,2, 0,0,0, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,0, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,0,

0,0,2, 0,0,2, 2,0,2, 2,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,0, 0,0,0, 0,0,2, 0,0,2, 0,0,2, 0,0,0, 0,0,4, 0,0,4, 2,0,4, 2,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 2,0,4, 2,0,4, 
0,0,2, 0,0,2, 0,0,2, 0,0,0, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 2,0,2, 2,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,0,

0,0,2, 0,0,2, 2,0,2, 2,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,0, 0,0,0, 0,0,2, 0,0,2, 0,0,0, 0,0,0, 0,0,4, 0,0,4, 2,0,4, 2,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 2,0,4, 2,0,4, 
0,0,2, 0,0,2, 0,0,2, 0,0,0, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 2,0,2, 2,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,0,

0,0,2, 0,0,2, 2,0,2, 2,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,0, 0,0,0, 0,0,2, 0,0,2, 0,0,0, 0,0,0, 0,0,4, 0,0,4, 2,0,4, 2,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 
0,0,2, 0,0,2, 0,0,2, 0,0,0, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,0,

0,0,2, 0,0,2, 2,0,2, 2,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,0, 0,0,0, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 2,0,4, 2,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 
0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 2,0,2, 2,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 2,0,4, 2,0,4, 0,0,0,

0,0,2, 0,0,2, 2,0,2, 2,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,0, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 2,0,4, 2,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 
0,0,2, 0,0,2, 0,0,2, 0,0,0, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,0, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,0,

0,0,2, 0,0,2, 2,0,2, 2,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,0, 2,0,0, 0,0,2, 0,0,2, 0,0,0, 0,0,0, 0,0,4, 0,0,4, 2,0,4, 2,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 2,0,4, 2,0,4, 
0,0,2, 0,0,2, 0,0,2, 0,0,0, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 2,0,2, 2,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,2, 0,0,2, 0,0,2, 0,0,2, 0,0,4, 0,0,4, 0,0,4, 0,0,4, 0,0,0};

 
 
 
 
 

/*the array nonstandchars stores instances of non standard characters (ie not 'a' 't' 'c' 'g' or 'u') for display after a file has been read in */
char nonstandchars[100]; 

/* the dynamically allocated array deletion specifies whether or not a codon at a position is to be excluded from all the analyses */
int *deletion = NULL;

/* this dynamically allocated array records whether each of the sequences is part of the outgroup or not */
int *outgroup = NULL;








