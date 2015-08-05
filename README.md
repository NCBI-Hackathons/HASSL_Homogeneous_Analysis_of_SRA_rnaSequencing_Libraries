

RNAseq analysis pipeline


This snakemake pipeline is designed to take SRA run accession numbers and output both sorted bam alignments and raw read counts for an infinite number of samples. 

To use: 

SETUP ENVIRONMENT

We've included a handy setupsnake.py tool to gather references and build the hisat index for you. This will put the GRCh38.p4 reference and annotation files in the lib directory and build the hisat index there as well. 

dependencies: snakemake, hisat-build, gunzip, wget

	snakemake -s setupsnake.py -j 


RUN PIPELINE

Edit the snakemake pipeline file to reflect the locations of your reference files and executables that snakemake will use. The default locationis the lib directory where setupsnake.py puts them. 

Isolate the run_accession IDs from SRA you want to run and put them in a file line by line. The pipeline defaults to look for this input file as'accessions.txt' 

Run it! You'll probably want to be on a fairly large machine for this.  

dependencies: snakemake, hisat (built with SRA support - see note below), samtools, HTSEQ (htseq-count), perl 

	snakemake -s snakemakefile.py -j 

(the -j flag enables multithreaded operations)

You'll want to refer to the snakemake documentation on how to (trivially!) run snakemake efficiently in a cluster environment that requires jobsubmission. 

<note on building HISAT with SRA support> 



