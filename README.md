

RNAseq analysis pipeline


This snakemake pipeline is designed to take SRA run accession numbers and output both sorted bam alignments and raw read counts for an infinite number of samples. 

To use: 

	<setup environment> - P
		(install snakemake, samtools, HTSeq) 

	Edit the snakemake pipeline file to reflect the locations of your reference file and executables that snakemake will use. 

	Construct a list of run accession IDs from SRA and put them into a file.... 

	Run it! You'll probably want to be on a fairly large machine for this.  

	snakemake -s snakemakefile.py -j 

	the -j flag enables multithreaded operations. 

	You'll want to refer to the snakemake documentation on how to (trivially!) run snakemake efficiently in a cluster environment that requires job submission. 

