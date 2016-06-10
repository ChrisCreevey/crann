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


/* The following functions are found in evolve.c */
int main(void);
void main_menu(void);
void general_options(void);
void delete_nonstd(int num);
void undelete_nonstd(int num);
void delete_stpcodons(void);
void undelete_stpcodons(void);
void Splash(void);
int check_files(void);
void open_input_file(void);
void open_output_file(void);
void show(void);
void show_tab_dna(void);
int choose_code(void);
int read_file (int fastaformat);
char read_sequence(int seq_num, struct sequence *new);
char getletter(char *instr);
void getstr(char *instr, char *outstr);
void clear_memory(void);
void clear_results(void);
void codon_usage(void);
void show_codons(void);
void show_Amino_Acids(void);
int transform_base(char c, int place);
void nonstandard(char c);
void clean_exit(void);
int getint(char *instr,int minx,int maxx, int def);
char *xgets(char *s);
void checkdata(void);
int found(void);
void count_bases(void);
void summary(void);
char * name_code(int i);
void AAmakeup(void);
