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

/*The following functions are found in Li_Wu_19851993.c     */
void count_degen(int *seq1,int *seq2,float *L,int code,int startw,int endw);
void trans_diff(int *seq1,int *seq2, float *P, float *Q, float L[3], int startw, int endw);
int is_transitional(char base1, char base2);
void kimura_method(float *A, float *B, float *a, float *b, float *c, float *varA, float *varB, float *K, float P[3], float Q[3], float L[3] );
void li_wu_85(float L[3], float A[3], float B[3], float *K, float *Ks, float *Ka);
void li_wu_93(float L[3], float A[3], float B[3], float varA[3], float varB[3],float a[3], float b[3], float Q[3], float P[3], float c[3], float *Ks, float *Ka,float *varKs, float *varKa);
void Li_Wu(void);
void li_wu_complete(int *seq1, int *seq2, float *results);
void which_sequences(void);
void omit_sequence(int choice);
void tag_all(void);
void window_size(int *startw, int *endw);
void shiftnwin_size(int *window, int *shift);
void display_Li_Wu(void);
void Li_Wu_movwin(void);
void display_movwin(int window, int shift);

