*Note:
I and the original coder of the md-3.x.cc program are no loger affliated with the Baek Group. Please contact them directly for code requests or other inquiries. Thank you!


# Site type discovery
These are the instructions for the use of the site type discovery program as described in Kim et al., <i>Nature Genetics</i> 2016. For details, please refer to the original publication, “General rules for functional microRNA targeting”, doi:10.1038/ng.3694.

# License information
The source code provided is for academic use only. For licensing details regarding the other publically available programs used in the preparation of the input files or elsewhere, please refer to the information provided by each individual program. 

# Pre-requisites 
The site type discovery program was developed in C++, compiled using GCC 4.4.7 or later. The standalone R library is required for conducting the statistical tests. The input files were created using in-house scripts written in Python 3.0 or later. 

# Usage
stdiscovery \<Commands\> \<Fold change list\> \<Site type list\> \<Index start\> \<Index end\> \<Site type pool directory\> \<Output file\> \<Target site cutoff\> \<Bin designation\>

The site type discovery program requires three preprocessed input files and additional input parameters. The descriptions for each of the input files and parameters are detailed below.

\<Commands\>

cst1, default option to conduct the chi-square test and Fisher’s Exact test
  
\<Fold change list\>

The fold change list contains the reference genes sorted according to the fold change values from the microarray data and organized into a predetermined number of bins from highly downregulated to less downregulated. In original study, the reference genes were organized into ten bins. 
-	See sample file, ‘sample_fclist.txt’ for more information.


\<Site type list\>

The site type list is comprised of non-redundant unambiguous site types representing the collection of interactions between each nucleotide from the miRNA and the interacting nucleotide on the mRNA. For an in-depth description regarding the definition of a site type, please refer to the original study. 
In order to compile the site type list, RNAsubopt from the ViennaRNA package (https://www.tbi.univie.ac.at/RNA/index.html) was used to predict the secondary structures between the mature miRNA sequence and the mRNA sequence. The program represents the complex between the two RNA molecules in a dot-bracket notation (example below), which depicts whether each nucleotide is bound (brackets) or not (dots). The site types were determined by converting the dot-bracket notation using the ‘interaction residues’ as described in Kim et al., 2016. 
Example:
-	3’UTR sequence		&nbsp AGAGGUCGGUGCAGACUCGGGGGCCAGUCC
-	miRNA sequence	UAAGGCACGCGGUGAAUGCCA 
-	dot-bracket notation	.(...(((.(((.(.((.((((.....)))).))...)))).)))...)..  
*note due to stochastic sampling, your notation may vary
-	site type		X_X_O_O_B1B1B1O1O_O_W_X_W_O_O_X_X_X_O_X_X_
-	See sample file ‘sample_stlist.txt’ for more information.

RNAsubopt usage example:
RNAsubopt -p 5 \<mRNA sequence\>&\<miRNA sequence\>

The –p option indicates the number of stochastically determined suboptimal structures drawn with probabilities equal to their Boltzmann weights as output. The '&' is used to connect 2 sequences that shall form a complex.

\<Index start\> and \<Index end\>

For multiplexing, the site type list is divided into a user defined number of segments as indicated by the start and end index parameters. 

\<Site type pool directory\>

The site type pool directory contains the compressed binary files that have been designed to optimize the runtime of the site type discovery program. The binary file holds information about the site type in regards to the positions of each miRNA:mRNA interaction (see example below). Additionally, the binary file stores the index of the site type from the site type list described previously. The compression of the files is conducted using the ‘struct’ package in Python.
Example:
-	dot-bracket notation	.(...(((.(((.(.((.((((.....)))).))...)))).)))...)..  
-	site type		X_X_O_O_B1B1B1O1O_O_W_X_W_O_O_X_X_X_O_X_X_
-	site type index		5033336
-	position info		[18, 17, 16, 15, -1, -1, -1, 13, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2]  
*Each position of where there is an miRNA:mRNA interaction is noted while bulges are denoted with -1
-	stored info		(5033336, 18, 17, 16, 15, -1, -1, -1, 13, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2)
-	See sample directory ‘sample_poolfiles’ for more information.

\<Output file\>

Designate the path to the output file. Please see ‘sample_output.txt’ for an example of the output format. 

\<Target site cutoff\> 

Minimum number of target sites required for a site type to be statistically evaluated in the program.

\<Bin designation\>

Decimal value of the binary notation used to determine the P values from which bins would be included in the output file. 
Example (for ten bins):
-	Bin designation value		387
-	Binary notation			110000011
-	The first and the last two bins are reported in the output file, see ‘sample_output.txt’ for more details.
	

#Example Usage
stdiscovery cst1 sample_fc_list.txt sample_st_list.txt 0 20 \<path to site type pool directory\> sample_output.txt
