Boxer

A Multiple Sequence Alignment Selection method for Phylogeny Reconstruction

as created for and used in 

Linchangco, G., Foltz, D., Reid, R., Williams, J., Nodzak, C., Kerr, A., Miller, A., Hunter, R., Wilson, N., Nielsen, W., Mah, C., Rouse, G., Wray, G., Janies, D.  The phylogeny of extant starfish (Asteroidea:Echinodermata) including Xyloplax, based on comparative transcriptomics.  Molecular Phylogenetics and Evolution (accepted).

http://www.sciencedirect.com/science/article/pii/S1055790317301653
==========================================================================================================================================

Boxer is a program written in Python that identifies optimal multiple sequence alignments based on alignment occupancy statistics and a user-defined unique taxa threshold.  

It has two inputs:1) a directory with the orignal alignments 2) a directory with alignments that have had spurious sequences removed.

Boxer is a lightweight program that can be run with a single line as shown below:

python /boxer.py -orig original/ -reAl realigned/ --gappyness 95 --taxacutoff 9 -f g95t9 --output g95t9.txt

Where -orig is populated with the path to the directory with the original alignments and -reAl is the path to the directory containing alignments with spurious sequences removed.

The --gappyness flag indicates the percentage of gap characters allowed in given alignment and the --taxacutoff indicates the number of unique taxa that an alignment should meet or exceed. 

The -f flag indicates the path to the output directory in which optimal alignemnts are copied into and the --ouput flag is the name of the output file containing some summary statistics for the boxer run.



