Crann: Fast heuristic methods of detecting adaptive evolution in protein-coding genes.

CRANN is a software program written in the C programming language that can be used to investigate adaptive evolution in a number of ways. The program requires a set of protein- codingDNAsequences, aligned so that the first residue of the alignment corresponds to the first position of a codon. The program checks to ensure that the alignment has valid codons and gap characters. First of all, the program implements some of the most popular methods of measuring synonymous and non-synonymous distances between a pair of sequences (Li et al., 1985; Li, 1993). 

Phylogenetic hypotheses based upon these distances can also be generated using CRANN, which implements the Neighbor-joining algorithm (Saitou and Nei, 1987). In addition, CRANN can also carry out a moving window analysis using either synonymous (ds) or non-synonymous (dn) dis- tances. This kind of analysis can be useful for locating regions of the protein that appear to be under varying selective pres- sures. CRANN can analyse the entire dataset and compile a cumulative result, or subsets of the dataset can be chosen (perhaps clades of sequences).

The most powerful part of this program is to be found in its ability to detect adaptive evolution along evolutionary lineages. There are two methods implemented in the software—the method described by Messier and Stewart (1997) and the method of Creevey and McInerney (2002). In both cases, hypothetical ancestral sequences are reconstructed at the internal nodes of a given phylogenetic tree.

If you use Crann, cite:
Creevey, C. and J. O. McInerney (2003). CRANN: Detecting adaptive evolution in protein-coding DNA sequences Bioinformatics 19: 1726.
Creevey, C. and J. O. McInerney (2002). An algorithm for detecting directional and non-directional positive selection, neutrality and negative selection in protein coding DNA sequences. Gene 300: 43-51.
