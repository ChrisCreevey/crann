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



/***************************************************************************************************************/																												
/* This program reads in a fasta formatted file which contains multiple sequences along with names etc....     */
/* The sequences are read into a linked list of structures which record the name, length, number and tag of    */
/* the sequence. Each structure conatins a dynamically allocated array, which stores the data representing the */
/* sequence itself. It runs through the input file reading in nucliotide triplets and allocates them a number  */
/* which represents the codon (see the Global array transform_values).										   */
/* The program then asks whether the user specified output file is to be in codon number format, or in amino   */ 
/* acid format. If amino acid is selected, then the program asks in which genetic code the input file is       */
/* written, there are a choice of 13 (see the global array genetic_codes). The amino acids are then outputted  */
/* to the file along with the name, sequence length etc.. for each sequence in the input file.				   */
/*																											   */
/* 					 Programmed by Chris Creevey 13-1-2000		Bioinformatics Lab, NUI Maynooth Ireland.	   */
/***************************************************************************************************************/




#include "definitions.h"
#include "evolve.h"
#include "adaptive_tree.h"
#include "linked tree.h"
#include "Li_Wu_19851993.h"


int main(void)
	{
		
	start = NULL;
	last = NULL;  /* initialise top and bottom pointers of the sequence in memory */
	li_wu_start = NULL;
	li_wu_end = NULL;  /* initialise top and bottom pointers of the li_wu results */

	printf("\n\n\n\n\n\n\n\n");
	main_menu();
	if(file != NULL) fclose(file);
	if(outfile != NULL) fclose(outfile);
	if(dist != NULL) fclose(dist);
	if(outtree != NULL) fclose(outtree);
	outtree = NULL;


	printf("\n\nFinished!!!!\n\n");	
	clean_exit();
	return(0);
	}











void main_menu(void)
	{
	int exit = FALSE, done = FALSE;
	int choice =0;
	do
		{
		choice = 0;
		printf("\n\t\t***************************************************");
		printf("\n\t\t* Main menu                                       *");
		printf("\n\t\t*                 Crann  1.2                      *");
        printf("\n\t\t*                                                 *");
		printf("\n\t\t*                 Copyright (c) 2016 Chris Creevey*");
		printf("\n\t\t***************************************************\n\n");
		
		printf("\t1 = Read new input file            < Current input file =  ");
		if(file == NULL) printf("None>\n");
		else printf("%s>\n", filename);
		printf("\t2 = Open new output file           < Current output file =  ");
		if(outfile == NULL) printf("None>\n");
		else printf("%s>\n", outfilename);	
		printf("\t3 = Output sequences from memory   < Sequences in memory = %d>\n", num_of_seqs);
		printf("\t4 = Perform Creevey, McInerney method\n");
		printf("\t5 = Calculate pairwise distances over the full length of the sequence\n");
		printf("\t6 = Calculate pairwise distances as a moving window analysis\n");
		printf("\t7 = Options\n");
		printf("\t8 = About Crann\n");
		printf("\t9 = Quit program\n");
		printf("\n\n\n");
		if(start == NULL && done == FALSE)
			choice = getint("\n\tPlease select an action (default = 1) ", 1, 9, 1);
		if(start != NULL && outfile == NULL && done == FALSE)
			choice = getint("\n\tPlease select an action (default = 2) ", 1, 9, 2);
		if(start != NULL && outfile != NULL && done != TRUE)
			choice = getint("\n\tPlease select an action (default = 4) ", 1, 9, 4);
		if(done == TRUE)
			choice = getint("\n\tPlease select an action (default = 9) ", 1, 9, 9);
		
		
		switch(choice)
			{			
			case 1:
				open_input_file();
				
				printf("\n\n\n");
				break;				
			case 2:
				open_output_file();
				printf("\n\n\n\n\n\n\n\n\n");
				break;			
			case 3:
				show();
				printf("\n\n\n");
				break;				
			case 4:
				if(check_files())
					{
					Li_Wu(); 
					if(check_distances() == TRUE)
						{
						printf("\nSome sequences were too distantly related from each other to accuratley\ncalculate the Li distances.");
						printf(" In these cases the Li value was set to\ntwice the largest calculated value in the matrix\n");
						}
					McDonald_Kreitman(); 
						
					done = TRUE;
					}
				break;
			case 5:
				if(check_files())
					{
					Li_Wu();
					allocate_distances(2);
					done = TRUE;
					}				
				break;
			case 6:
				if(check_files()) Li_Wu_movwin();
				done = TRUE;
				break;
			case 7:
				if(check_files()) general_options();
				break;
			case 8:
				Splash();			
				break;
			case 9:
				exit = TRUE;
				break;
			default:
				choice = 0;
			}
				
		}while(!exit|| choice == 0);
	}

				

void general_options(void)
	{
	int  choice = 0;
	int count = 0, ans = 0;
	do{
	
		printf("\n\n\n\n\n\n\n\n\t\t***************************************************");
		printf("\n\t\t* General Options Menu                            *");
		printf("\n\t\t*                 Crann  1.1                      *");
		printf("\n\t\t*                                 by Chris Creevey*");
		printf("\n\t\t***************************************************\n\n");
		printf("\t1 = Genetic Code is %s\n", name_code(code) );
	 
		printf("\t2 = Deletion of Non-Standard Characters = ");
		if(gen_opt[0] == 0) printf("Pairwise\n"); else printf("Complete\n");
		printf("\t3 = Deletion of Stop Codons  = ");
		if(gen_opt[1] == 0) printf("Pairwise\n"); else printf("Complete\n");
		printf("\t4 = Analyse all sequences in memory? = ");
		if(gen_opt[2] == 0) printf("No, # omitted = %d\n", untagged); else printf("Yes\n");
		printf("\t5 = Analyse whole sequence length? = ");
		if(gen_opt[3] == 0) printf("No, start = %d, end = %d\n", startw, endw); else printf("Yes\n");

		printf("\t6 = Which Li method to use? = ");
		if(gen_opt[4] == 0) printf("1993\n"); else printf("1985\n");

		printf("\t7 = Input Phylogentic tree? = ");
		if(gen_opt[5] == 0) printf("No\n"); else printf("Yes, File = %s\n", nestname);
		if(gen_opt[5] == 0)
			{
			printf("\t8 = Build a neighbour joining tree with which distances? = ");
				switch(gen_opt[6])
					{
					case 0:
						printf("Dn\n");
						break;
					case 1:
						printf("Ds\n");
						break;
					default:
						printf("Dn/Ds\n");
						break;	
					}
			
			}
		printf("\n\t0 = Return to main menu (default)\n");
	/* Use for Mac version */
/*			printf("\n\nEnter a number to change that option, or press return to continue\n\n\n\n\n\n");
		ans = getch(); */
	/* Use for PC version */
	ans = getint("\n\nEnter a number to change that option, or press return to continue\n", 0, 9, 0 );		
	
		switch(toupper(ans))
			{
			case 1 :
				choose_code();
				if(gen_opt[1] == 1) delete_stpcodons();
				checkdata();
				break;
			case 2:
				if(found() == TRUE)
					{
					gen_opt[0] = 1;
					delete_nonstd(64);
					}
				else 
					{
					gen_opt[0] = 0;
					undelete_nonstd(64);
					}
				break;
			case 3:
				if(found() == TRUE)
					{
					gen_opt[1] = 1;
					delete_stpcodons();
					}
				else
					{
					gen_opt[1] = 0;
					undelete_stpcodons();
					}
				break;
			case 4:
				which_sequences();
				if(untagged == 0) gen_opt[2] = 1; else gen_opt[2] = 0;
				break;
			case 5:
				printf("\n\n\n\n\n\n\n\n\n\n\n\n");
				window_size(&startw, &endw);
				if(startw == 0 && endw == (start->length/3 -1)) gen_opt[3] = 1; else gen_opt[3] = 0;
				break;
			case 6:

					printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\tWhich method would you like to use? ");
					printf("\n\t1 = 1985 method ");
					printf("\n\t2 = 1993 method ");
					choice = getint("\n\n\tPlease choose either 1 or 2", 1,2,2);
						
					switch(choice)
						{
						case 1:
							gen_opt[4] = 1;
							break;
						default:
							gen_opt[4] = 0;
							break;
			
						}

				break;
			case 7:
				if(tree_top != NULL) dismantle(tree_top, &count);
				if(tree_choice() == 1) gen_opt[5] = 1; else gen_opt[5] = 0;
				break;
			case 8:
					/* choose whether to use Ks or Ka values */
				choice = 0;

					printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\tWhich distances would you like to use in the algorithm\n");	
					printf("\t1 = Dn values\n");
					printf("\t2 = Ds values\n");	
					printf("\t3 = Rate of evolution (Dn/Ds)\n");
					choice = getint("\n\tPlease choose either 1, 2 or 3\n", 1,3,1);
				
					switch(toupper(choice))
						{
						case 1:
							gen_opt[6] = 0;
							break;
						case 2:
							gen_opt[6] = 1;
							break;
						default:
							gen_opt[6] = 2;
							break;

						}

				break;
					
			default:
				ans = '\0';
				break;
				
			}
		}while(ans != '\0');
	printf("\n\n\n\n\n\n\n\n");
	}

void delete_nonstd(int num)
	{
	struct sequence *position = start;
	int i = 0;
	
	while(position != NULL)
		{
		i = 0;
		while(position->bases[i] != 193)
			{
			if(position->bases[i] == num) deletion[i] = TRUE;
			i++;
			}
		position = position->next;
		}
	}

void undelete_nonstd(int num)
	{
	struct sequence *position = start;
	int i = 0;
	
	while(position != NULL)
		{
		i = 0;
		while(position->bases[i] != 193)
			{
			if(position->bases[i] == num) deletion[i] = FALSE;
			i++;
			}
		position = position->next;
		}
	}

void delete_stpcodons(void)
	{
	struct sequence *position = start;
	int i = 0;
	
	for(i=0; i<(start->length/3) -1; i++) deletion[i] = FALSE; /* reset the array */
	
	if(gen_opt[0] == 1) delete_nonstd(64); /* if we want the non std chars deleted */
	
	while(position != NULL)
		{
		i = 0;
		while(position->bases[i] != 193)
			{
			if(genetic_codes[code][position->bases[i]] == 0)
				{
				deletion[i] = TRUE;
				}
			i++;
			}
		position = position->next;
		}
	}
	
void undelete_stpcodons(void)
	{
	int i =0;
	
	for(i=0; i<(start->length/3) -1; i++) deletion[i] = FALSE; /* reset the array */
	if(gen_opt[0] == 1) delete_nonstd(64); /* if we want the non std chars deleted */
	}	

void Splash(void)
	{
	char ans[80];
	printf("\n\n");
	printf("/**********************************************************************/\n");
	printf("/*                             Crann 1.2                              */\n");
	printf("/*                                                                    */\n");
	printf("/*                    Written by: Chris Creevey                       */\n");
	printf("/*      Initial release April 2003; Current release October 2016      */\n");
	printf("/*         Address:                                                   */\n");
	printf("/*          Ecological and Evolutionary Genomics Laboratory,          */\n");
	printf("/*          IBERS,                                                    */\n");
	printf("/*          Aberystwyth University,                                   */\n");
	printf("/*          Aberystwyth,                                              */\n");
	printf("/*          SY23 2DA,                                                 */\n");
	printf("/*          United Kingdom.                                           */\n");
	printf("/*                                                                    */\n");
	printf("/*         Email: chris.creevey@gmail.com                             */\n");
	printf("/*         Web:   http://www.creeveylab.org                           */\n");
	printf("/*                                                                    */\n");
	printf("/*                                                                    */\n");
	printf("/*       Any and all suggestions and/or comments are welcome.         */\n");
	printf("/*                                                                    */\n");
	printf("/*  Copyright: Chris Creevey 2003-2016                                */\n");
	printf("/**********************************************************************/\n");
	printf(" Press return to continue: ");
	xgets(ans);
	printf("\n\n\n\n\n\n\n\n");

	}






int check_files(void)
	{
	int pass = TRUE;
	
	if(start == NULL)
		{
		printf("\n\nNo sequences have been read into memory, \nPlease select to read in a new input file from the main menu\n");	
		pass = FALSE;
		}
	else
		{
		if(outfile == NULL) 
			{
			printf("\n\nNo output file has been specified, \nPlease select to open a new output file from the main menu\n");
			pass = FALSE;
			}
		}
	return(pass);
	}




void open_input_file(void)

	{
	char choice;
 	int fastaformat, exit = FALSE, i =0; 
	char overflow = 'a';
	struct sequence *position = start;


	if(file != NULL)
		 {
		 fclose(file);      /* If there has been an input file opened already, close it before opening the new one */
		 file = NULL;
		 }
	clear_memory();
	
	for(i=0; i<100; i++) nonstandchars[i] = '\0'; /* initialise nonstandchars array before reading new file */
	
	do{
		do{
			fastaformat = TRUE;
			printf("\n\n\tName of the fasta formatted file ");
			scanf("%s%c", string1, &overflow);
			filename[0] = '\0'; strcpy(filename, string1); 
			if((file = fopen(filename, "r")) == NULL)		/* check to see if the file is there */
				{
				printf("\n\n\tCannot open the sequence file, named %s\n\n", filename);
				choice = getletter("\n\n\tPress 1 to return to the Main Menu\n\tOR press 2 to try again");
				if(toupper(choice) == '1') exit = TRUE;
				}

			}while(file == NULL && !exit);
		if(!exit)
			{
			printf("opened %s.....\n", filename);
			
			fastaformat = read_file(fastaformat);
			if(!fastaformat)
				{
				choice = getletter("\n\n\tPress 1 to return to the Main Menu\n\tOr press 2 to try a different file");
				if(toupper(choice) == '1') exit = TRUE;
				}
			}
  		}while(!fastaformat && !exit );


  	if(!exit) 
  		{
  		checkdata();
/*	 	codon_usage(); 
	 	AAmakeup();
	count_bases();		
*/	
		}
	
		

	}
	
	
	
void open_output_file(void)
	{
  	
  	char choice, overflow = 'a';
  	int exit = FALSE;
  	
  	if(outfile != NULL)
  		{
  		fclose(outfile); /* If there has been an output file opened already, close it before opening the new one */	
		outfile = NULL;
		}
	
  	do{		
		printf("\n\n\tName of the output file: ");
		scanf("%s%c", string1, &overflow);
		outfilename[0] = '\0'; strcpy(outfilename, string1); 
		if((outfile = fopen(outfilename, "w")) == NULL)		/* check to see if the file can be opened/created */
			{
			printf("\n\n\tCannot open the output file, named %s\n\n", outfilename);
			choice = getletter("\n\n\tPress 1 to return to the Main Menu\n\nOR press 2 to try again: ");
			if(toupper(choice == '1')) exit = TRUE;
			}
		}while(outfile == NULL && !exit);
		
	}
	
	
	



void show(void)
	{

	int choice1 = 0; 

	
	if(check_files())
		{	

			choice1 = getint("\n\n\n Please select the output file format:\n\t1 = Codon numbers\n\t2 = Amino acid format\n\t3 = Dna sequences in tab delimited format (in columns)\n\t4 = Return to the main menu\n\n\t please select 1, 2, 3 or 4: ", 1,4,4);
				
			
			switch(choice1)
				{
				case 1:
					show_codons();
					break;
				case 2:
					show_Amino_Acids();
					break;
				case 3:
					show_tab_dna();
					break;	
				default:
					break;
				}

		}	
		
	}



void show_tab_dna(void)
	{
	struct sequence *seq = start;
	int i = 0, j = 0;
	
	while(seq != NULL)
		{
		fprintf(outfile, "%s\t", seq->name);
		seq = seq->next;
		}
	fprintf(outfile, "\n");
	seq = start;
	while(seq->bases[i] != 193)
		{
		
		while(seq != NULL)
			{
			for(j=0; j<3; j++)
				{
				fprintf(outfile, "%c", codons[genetic_codes[code][seq->bases[i]]][j]);
				}
			fprintf(outfile, "\t");
			seq = seq->next;
			}
		fprintf(outfile, "\n");
		i++;
		seq = start;
		}
	}
			
			
			
			

int choose_code(void)
	{
	int i = '\0';
	char choice = '\0';

			printf("\n\n\n\n\n\n\n\n\n\n\n\n");
			printf("\n\nPlease specify which genetic code the input file contains:\n\n");
			printf("\t\t 1 = Universal (default)\n");
			printf("\t\t 2 = Vertebrate standard mitochondrial\n");
			printf("\t\t 3 = Yeast mitochondrial\n");	
			printf("\t\t 4 = Mycoplasma/Spiroplasma/Mold/Protozoan/Coelenterate\n");
			printf("\t\t 5 = Invertebrate mitochondrial\n");
			printf("\t\t 6 = Ciliate\n");
			printf("\t\t 7 = Echinoderm mitochondrial\n");
			printf("\t\t 8 = Euplotid\n");	
			printf("\t\t 9 = Bacterial (same as universal)\n");
			printf("\t\t 10= Alternative Yeast Nuclear\n");
			printf("\t\t 11= Ascidian mitochondrial\n");
			printf("\t\t 12= Flatworm mitochondrial\n");		
			printf("\t\t 13= Blepharisma Nuclear\n");
			
			i = getint("\n\n\t please select one of 1 to 13 from above: ", 1,13,1);
	

		code = i -1;
		return(i);
	}
	
					


int read_file(int fastaformat)
	{
	
	struct sequence *new = NULL;
	int j;
	char c;
	while(!feof(file) && ((c= getc(file)) == ' ' || c == '\t' || c == '\n')); /* skip past invisible characters */

	if(c == '>')
		{
		fastaformat = TRUE;
		j = 0;
		do{
				new = malloc(sizeof(list_entry));
				if(new == NULL)
					{
					printf("\n\t Out of memory\n");
					clean_exit();
					}
			
			c = read_sequence(j, new);
			j++;
			
			}while((c == '>') && (feof(file) != 1));
		}
			
	else
		{
		/* if the first character in the file is not a '>' then  we assume the format is wrong */
 		printf("\n\n\tThis does not seem to be a fasta formatted file\n\n");
		fastaformat = FALSE;
		}
	num_of_seqs = j;	
	return(fastaformat);
	
	}
			
			
			
			
			
			

char read_sequence(int seq_num, struct sequence *new)
	{
	
	char c = '\0';
	int i = 0, j = 0, place = 0, value = 0, memory_allocations = 0;
	


	/* assign the sequence number */
	new->seq_num = seq_num;
	

	/* initialise the tag (equal to TRUE) */
	new->tag = TRUE;
	
	
	/* All sequences are assumed to not be in the outgroup in the begining */
	new->outgroup = FALSE;

	/* read in the name of the sequence */
	j = 0;
	c = getc(file);
	while(c == ' ') c = getc(file);
	do{				
	 
	 	new->name[j] = c;
		if(j < maxnamlen) j++;
		

		}while(feof(file) == 0 && (c = getc(file)) != '\n' && c != '\r');
	new->name[j] = '\0';				/* append a '\0' terminator */

	/* create the nickname */
	for(i=0; i<19; i++)
		{
		if(new->name[i] != '\0' && new->name[i] != ' ')
			{
			new->nickname[i] = new->name[i];
			}
		else
			{
			new->nickname[i] = '\0';
			i=20;
			}
		}
	new->nickname[i] = '\0';
	
	/* read in the sequence */
    i = 0;
    j = 0;
    
    new->bases = malloc((STD_CODON_NUM * (sizeof(int))));   /* Allocate the memory to store the sequence in codons */
	memory_allocations = 1;
	c = getc(file);
	do{	

		place = 0;
		value = 0;
		do{			
			if(c != '\n' && c != '\r' && feof(file) == 0 && c != ' ')  /* if not an end of line or end of file */
				{
				value = value + transform_base(c, place);      /* call transform_base to calculate the codon number */
				place++;
				j++;			
				}
				
		}while(((c= getc(file)) != '>') && (place < 3) && (feof(file) == 0));

		if(value > 64) value = 64;     /* if there is a gap, then  the codon value is assigned to 64 */

		if(i >= ((STD_CODON_NUM * memory_allocations) - 1))  /* if the sequence has exceeded the length of STD_CODON_NUM * the number of memory allocations */
			{
			memory_allocations++;
			new->bases = realloc(new->bases, ((STD_CODON_NUM * (sizeof(int))) * memory_allocations));    /* reallocate the memory space to the size required to fit the STD_CODON_NUM * new memory allocations */
			if(new->bases == NULL)
				{
				printf("\n\t Out of memory\n");
				clean_exit();
				}
			}
		if(c != '>' && place == 3)   /* the check to make sure c is not a space is incase there are any spaces at the end of the file */
			{		
			new->bases[i] = value;
			i++;
			}	
		
	}while((c != '>') && (feof(file) == 0));	
	new->bases[i] = 193;	/* 193 is the terminator value for the sequence */
	

	/* assign the length of the sequence */			
	new->length = j;
	
	/* initialise the numofstpcodons value */
	new->numofstpcodons = 0;
	

	/* initialise the gaprun value (Used in McD&K algorithm) */
	new->gaprun = 0;
	
	/* initialise the stopcodon array */
	
	new->stopcodons = (malloc(1 * (sizeof (int))));
	if(new->stopcodons == NULL)
			{
			printf("\n\t Out of memory\n");
			clean_exit();
			}

	/* now assign the pointers to the previous and next sequence */
	if (last == NULL || start == NULL) 	/* If this is a new list */
		{
		last = new;
		start = new;
		}		
	else (last)->next = new;
	new->next = NULL;
	new->previous = last;
	last = new;

	return(c);


	}


char getletter(char *instr)
	{
	char outchar = '\0', overflow = '\0';
	printf("%s: ", instr);
	
		outchar = getc(stdin);  /* Use for Mac version */
	
/*	scanf("%1c%c", &outchar, &overflow); */ /* Use for PC version */
	return(outchar);
	}


void getstr(char *instr, char *outstr)
	{
	printf("%s: ", instr);
	xgets(outstr);
	}
	

void clear_memory(void)
	{
	int i;
	struct sequence *new = NULL;
	struct synon *place = NULL;

	if(start != NULL)
		{
		while(start != NULL) 
			{
			 new = start->next;
			 free(start->bases);
			 free(start); 
			 start = new;
			}
		start = NULL;
		last = NULL;
		}
	if(li_wu_start != NULL)
		{
		while(li_wu_start != NULL)
			{
			place = li_wu_start->next;
			free(li_wu_start->Ks);
			free(li_wu_start->varKs);
			free(li_wu_start->Ka);
			free(li_wu_start->varKa);
			free(li_wu_start);
			li_wu_start = place;
			}
		li_wu_start = NULL;
		li_wu_end = NULL;
		}
	if(deletion != NULL)
		{
		free(deletion);
		deletion = NULL;
		}
	num_of_seqs = 0;
	if(outgroup != NULL) 
		{
		free(outgroup);
		outgroup = NULL;
		}

	if(distances != NULL)
		{
		for(i=0; i<num_of_seqs-untagged; i++)
			free(distances[i]);
		free(distances);
		}
	distances = NULL;
	
	if(file != NULL)
		{
		fclose(file);
		file = NULL;
		}
	}
	
void clear_results(void)
	{

	struct synon *place = NULL;
	if(li_wu_start != NULL)
		{
		while(li_wu_start != NULL)
			{
			place = li_wu_start->next;
			free(li_wu_start->Ks);
			free(li_wu_start->varKs);
			free(li_wu_start->Ka);
			free(li_wu_start->varKa);
			free(li_wu_start);
			li_wu_start = place;
		}
		li_wu_start = NULL;
		li_wu_end = NULL;
		}
	}

/* This file counts the codon usage data for the dataset entered and writes it to the file codonUsage.out */

void codon_usage(void)
	{
	struct sequence *position = start;
	int i;
	float usage[64], total = 0;
	
	for(i=0; i<64; i++) usage[i] = 0;
	
	while(position != NULL)
		{
		i = 0;
		while(position->bases[i] != 193) 
			{
			if(position->bases[i] < 64)
				{
				usage[position->bases[i]]++;
				total++;
				}
			i++;
			}
		position = position->next;
		}
		
	if(usagefile != NULL)
 		{
  		fclose(usagefile); /* If there has been an output file opened already, close it before opening the new one */	
		usagefile = NULL;
		}
	
	if((usagefile = fopen("codonUsage.out", "w")) == NULL)
		{
		printf("Could not open cononUsage.out\n codon usage data not written\n");
		}
	else
		{
	
		for(i=0; i<64; i++)
			{
			fprintf(usagefile, "%f ", usage[i]/total);
			if(fmod(i+1,4) == 0 && i != 0) fprintf(usagefile, "\n");
			}
		printf("Codon usage data written to codonUsage.out");	
		
		}
	}



void show_codons(void)
	{
	int i;
	struct sequence *info = start;
	
	while(info != NULL) {		/* while we're not at the end of the linked list */
			i = 0;
			fprintf(outfile, "\n\nname:......%s\n", info->name);
			fprintf(outfile, "sequence length: %d\n", info->length);
			fprintf(outfile, "sequence number: %d\n", info->seq_num);
			fprintf(outfile, "codon numbers as follows:\n");
				while(info->bases[i] != 193)  /* while the codon number is not the terminator value */
					{
					fprintf(outfile, "%d ", info->bases[i]);
					i++;
					if(fmod(i, 30) == 0) fprintf(outfile, "\n");
					}
			
			info = info->next; /* get next sequence */
			
		}
	}

void show_Amino_Acids(void)
	{
	int i;
	struct sequence *info = start;
	
	while(info != NULL) {		/* while we're not at the end of the linked list */
			i = 0;
			fprintf(outfile, ">%s\n", info->name);
	/*		fprintf(outfile, "sequence length: %d\n", info->length);
			fprintf(outfile, "sequence number: %d\n", info->seq_num);
			fprintf(outfile, "Amino Acids as follows:\n (note that 'XXX' maens stop, & '---' means space or gap)\n");
	*/			while(info->bases[i] != 193)  /* while the codon number is not the terminator value */
					{					
					fprintf(outfile, "%c", amino_acids[genetic_codes[code][info->bases[i]]]);

					i++;
					if(fmod(i, 60) == 0) fprintf(outfile, "\n");
					}
			fprintf(outfile, "\n");
			info = info->next; /* get next sequence */
			
		}
	}
	





int transform_base(char c, int place)
	{
	int value;

	value = 0;	
	switch(toupper(c))
		{
		case 'T':
		case 'U':
			value = transform_values [0][place]; 
			break;
		case 'C':
			value = transform_values [1][place];
			break;
		case 'A':
			value = transform_values [2][place];
			break;
		case 'G':
			value = transform_values[3][place];
			break;
		default:
			nonstandard(c);
			value = transform_values[4][place];
		}
	return(value);
	
	}
	
			
void nonstandard(char c)
	{
	int i = 0;
	
	while(nonstandchars[i] != '\0' && nonstandchars[i] != c && i < 100) i++;
			
	if(nonstandchars[i] == '\0')
		{
		nonstandchars[i] = c;
		nonstandchars[i+1] = '\0';
		}
	}
	
	
void clean_exit(void)
	{
	clear_memory();
	exit(1);
	}

int getint(char *instr,int minx,int maxx, int def)
{
        int ret, status;
        char line[3];
        while(TRUE) {
                fprintf(stdout,"%s (%d..%d)    [%d]: ",instr,(int)minx,(int)maxx,(int)def);
        	    xgets(line);
                status=sscanf(line,"%d", &ret);
                if(status == EOF || status == 0) return def;
                if(ret>maxx) {
                        fprintf(stdout,"ERROR: Max. value=%d\n\n",(int)maxx);
                        continue;
                }
                if(ret<minx) {
                        fprintf(stdout,"ERROR: Min. value=%d\n\n",(int)minx);
                        continue;
                }
                break;
        }		
        return ret;
}


/* A simple version of the standard gets() library function which is unsafe on some systems......... code taken from C: the complete reference
	by Herbert Schildt
*/
char *xgets(char *s)
	{
	char ch, *p = NULL;
	int t = 0;
	
	p = s;  /* gets returns a pointer to s */
	
	for(t=0; t<80; ++t) {
		ch = getchar();
		
		switch(ch) {
			case '\n':
				s[t] = '\0'; /* terminate the string */
				return(p);
				break;
			case '\r':
				s[t] = '\0'; /* terminate the string */
				return(p);
				break;
			case '\b':
				if(t>0) t--;
				break;
			default:
				s[t] = ch;
				break;
			}
	}
	s[79] = '\0';
	return p;
}


void checkdata(void)
	{
	int i = 0, question1 = FALSE, question2 = FALSE, length = 0, exit = FALSE;
	struct sequence *position = start;
	
	
	
	length = start->length;

	if(deletion != NULL) free(deletion);
	deletion = malloc(length * (sizeof(int))); /* Make the array deletion the size of the sequences */
	if(deletion == NULL) clean_exit();
	for(i =0; i<length; i++) deletion[i] = FALSE;
	
	

		/* check length of sequences */

	while(position != NULL && !exit)
		{
		if(position->length != length)
			{
			printf("\n\n\n\n\n\nERROR: Sequences are not the same length\n");
			printf("\n\tClearing sequences from memory\n");
			clear_memory();
			printf("\n\tPlease only input files with aligned sequences\n");
			file = NULL;
			exit = TRUE;
			}
		position = position->next;
		}
		
		/* end of check length */
		
	position = start;

	while(position != NULL && !exit)
		{
		i = 0;
		while(position->bases[i] != 193 && !exit)
			{
		/* start of section to check for gaps/non-standard characters */
			if(position->bases[i] == 64)
				{
				if(!question1)
					{
					printf("\n\n\n\n\nWarning: Non-standard characters and/or gaps were encountered\n");
					question1 = TRUE;
					}
				}

		/* end of section to check for gaps */
		
		
			i++;
			}
		position = position->next;
		}
		
	/* Start of section to record all stop codons in the sequence (except the last) */
        position = start;
        while(position != NULL)
            {
            if(position->stopcodons != NULL) free(position->stopcodons);
            position->numofstpcodons = 0;
            position->stopcodons = malloc(1*sizeof(int));
            position = position->next;
            }
        
        
        
	 position = start;
	 	while(position != NULL)
			{
			i = 0;
			while(position->bases[i] != 193)
				{
				if(genetic_codes[code][position->bases[i]] == 0 && i != ((start->length)/3) - 1)
					{
					
					position->numofstpcodons++;
					position->stopcodons = realloc(position->stopcodons, ((position->numofstpcodons) * (sizeof(int))));
					position->stopcodons[position->numofstpcodons -1] = i;
					}
				i++;
				}
			position = position->next;
			}
	/* End of stop codon check */
	summary();
	}



/* this function is called by checkdata to ask the user their preferences when dealing with gaps/non-standard characters */	
int found(void)
	{
	int choice = 0;
	int complete = FALSE;
	
		printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
		printf("\n\n\tIn these cases would you prefer to:");
		printf("\n\t1 = Use pairwise deletion");
		printf("\n\t2 = Use complete deletion");
		choice = getint("\n\n\tPlease select from option 1 or 2", 1, 2, 1);
		
		switch(choice)
			{
			case 1:
				complete = FALSE;
				break;
			case 2:
				complete = TRUE;
				break;
			default:
				choice = '\0';
				break;
			}
		
	return(complete);
	}	





void count_bases(void)
	{
	struct sequence *position = NULL;
	int i = 0, j = 0;
	int AA[22];
	FILE *count = NULL;
	
	if((count = fopen("basecount.out", "w")) == NULL)		/* check to see if the file is there */
		{
		printf("\n\n\tCannot open the sequence file, named basecount.out\n\n");
		}
	else
		{
		i =0;
		for(j=0; j<22; j++)
			fprintf(count, "%c\t", amino_acids[j]);

		fprintf(count, "\n");
		position = start;
		while( position->bases[i] != 193)
			{
			
			for(j=0; j<22; j++) AA[j] = 0;
			
			while(position != NULL)
				{
			
				AA[genetic_codes[code][position->bases[i]]]++;
				position = position->next;
				}
			
			for(j=0; j<22; j++)
				{
				fprintf(count, "%f\t", (float)AA[j]/(float)num_of_seqs);
				}
			fprintf(count, "\n");
			position = start;
			i++;
			}
		}	
	fflush(count);
	fclose(count);	

}



/* this function prints a summary of warnings and/or choices made during data input */	
void summary(void)
	{
	int i = 0, done = FALSE, stop = FALSE;
	struct sequence *position = start;
	char ans[80];
	
	printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\t//****************************************************//\n");
	printf("\t                  Input data summary                  \n");
	printf("\t                                                    \n");
	printf("\t Number of sequences in memory = %d\n", num_of_seqs);
	printf("\t genetic code = %s ", name_code(code));

	if(stop)
		{
		printf("\n\n\t data input was halted due to a data format error:  \n");
		printf("\t The sequences were all not the same length         \n");
		printf("\t                                                    \n");
		}
	else
		{
		i = 0;
		if(nonstandchars[0] != '\0')
			{
			printf("\n\n\t The following non-standard characters were found;  \n\t");
			while(nonstandchars[i] != '\0' && i < 100)
				{
				printf("'%c' ", nonstandchars[i]);
				i++;
				}				

			printf("\n");
			printf("\t Non standard charaters will be treated as follows: \n");
			printf("\t Any codons containing nonstandard charaters will be\n");
			if(gen_opt[0] == 1)
				{
				printf("\t excluded through complete deletion of that codon in\n");
				printf("\t all sequences.                                     \n");
				}
			else
				{				
				printf("\t excluded through pairwise deletion of that codon in\n");
				printf("\t comparisons with other sequences                   \n");
				}
			}
		else printf("\n\n\t Number of Non-standard characters = 0\n");
		position = start;
		if(position != NULL) 
			{
			while(position != NULL)
				{
				if(position->numofstpcodons != 0)
					{
					if(!done)
						{
						stop = TRUE;
						printf("\n\n\t Stop codons were found in the following sequences\n");
						done = TRUE;
						}
					printf("\n\t");
					printf("%d in '%s'\n",position->numofstpcodons, position->name);
					printf("\tAt positions:");
					for( i=0; i<position->numofstpcodons; i++)
						{
						printf("'%d' ", position->stopcodons[i]);
						}
		
					}
				position = position->next;
				}
			if(done)
			{
				if(gen_opt[1] == 1) printf("\n\n\t Stop codons will be excluded using complete deletion\n");
				else printf("\n\n\t Stop codons will be excluded using pairwise deletion\n");
			}
			printf("\n");			

			}	
		if(!stop) printf("\t Number of Stop codons = 0\n");
		}
	printf("\n");
	printf("\n");
	printf("\t//****************************************************//\n");
	
	/* Use for PC vesion*/
	printf("\tPress return to continue\n\n\n\n\n\n\n\n");
	xgets(ans);

	

	/* Use for Mac version */
/*	printf("\tPress any key to continue\n");
	ans = getch();
	printf("\n\n\n\n\n\n\n");
*/	}
		


/* This function is called in summary to return the name of the genetic code chosen by the user */		
char * name_code(int i)
	{
	switch(i)
		{
		case 0:		
			return("Universal genetic code\0");
			break;
		case 1:
			return("Vertebrate standard mitochondrial genetic code\0");
			break;
		case 2: 
			return("Yeast mitochondrial genetic code\0");
			break;
		case 3:
			return("Mycoplasma/Spiroplasma/Mold/Protozoan/Coelenterate genetic code\0");
			break;
		case 4:
			return("Invertebrate mitochondrial genetic code\0");
			break;
		case 5:
			return("Ciliate genetic code");
			break;
		case 6:
			return("Echinoderm mitochondrial genetic code\0");
			break;
		case 7:
			return("Euplotid genetic code\0");
			break;
		case 8:
			return("Bacterial genetic code (same as universal)\0");
			break;
		case 9:
			return("Alternative Yeast Nuclear genetic code\0");
			break;
		case 10:
			return("Alternative Yeast Nuclear genetic code\0");
			break;
		case 11:
			return("Flatworm mitochondrial genetic code\0");
			break;
		case 12:
			return("Blepharisma Nuclear genetic code\0");
			break;
		default:
			return("User defined exit\0");
			break;
		}
	}	 
			
/* This function is added to give more details as to the make up of the sequences, Added 18/6/01 it is not impotant to the function of the core program and
will just give an extra file detailing the percentage makeup of the different Amino acids. This can be useful if selection favoured a certain number of a particular 
Amino acid, and not just particular positions on a sequence, commissioned by Mary o connell. */

void AAmakeup(void)
	{
	FILE *AAfile = NULL;
	struct sequence *position = start;
	int AA[22] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, i = 0, j=0, cat = 0;
	

	if((AAfile = fopen("AAmakeup.out", "w")) == NULL)
		{
		printf("unable to open AAmakeup.out\n Amino acid makeup deatils not written..... continuing program.\n");
		}
	else
		{
		for(cat = 0; cat<10; cat++)
			{
			position = start;
			switch (cat)
				{
				case 0:
					fprintf(AAfile, "\n\nhydrophobic Amino acids\n");
					break;
				case 1:
					fprintf(AAfile, "\n\npositive Amino acids\n");
					break;
				case 2:
					fprintf(AAfile, "\n\nnegative Amino acids\n");
					break;
				case 3:
					fprintf(AAfile, "\n\npolar Amino acids\n");
					break;
				case 4:
					fprintf(AAfile, "\n\ncharged Amino acids\n");
					break;
				case 5:
					fprintf(AAfile, "\n\nsmall Amino acids\n");
					break;
				case 6:
					fprintf(AAfile, "\n\ntiny Amino acids\n");
					break;
				case 7:
					fprintf(AAfile, "\n\naromatic Amino acids\n");
					break;
				case 8:
					fprintf(AAfile, "\n\naliphatic Amino acids\n");
					break;
				default:
					fprintf(AAfile, "\n\nAll Amino acids\n");
					break;
				}	

			for(i=0; i<22; i++)
				if(AA_categories[i][cat] == 1 || cat == 9) fprintf(AAfile, "\t%c", amino_acids[i]);
			fprintf(AAfile, "\n");	

			while(position != NULL)
				{
			
				for(i=0; i<22; i++) AA[i] = 0;
				i = 0;
				while(position->bases[i] != 193) /* as long as we havn't reached the end of the sequence */
					{
					AA[genetic_codes[code][position->bases[i]]]++;
					i++;
					}
		
				fprintf(AAfile, "%s", position->name);
				for(j=0; j<22; j++)
					if(AA_categories[j][cat] == 1 || cat == 9) fprintf(AAfile, "\t%f", (float)AA[j]/(float)i);
				fprintf(AAfile, "\n" );
				position = position->next;
				}

			}


		fclose(AAfile);
		}
	}
	
	
			
	
	 				
			
