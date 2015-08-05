

# RNAseq analysis pipeline


This snakemake pipeline is designed to take SRA run accession numbers and output both sorted bam alignments and raw read counts for an infinite number of samples. 

To use: 

## SETUP ENVIRONMENT

We've included a handy setupsnake.py tool to gather references and build the
hisat index for you.  This will put the GRCh38.p4 reference and annotation
files in the `lib` directory and build the hisat index there as well.

You will need to install the following before running `setupsnake.py`
* [snakemake](https://bitbucket.org/johanneskoester/snakemake/downloads)
* [hisat](https://github.com/infphilo/hisat/)
* [HTSeq](http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html)
* [Picard](https://broadinstitute.github.io/picard/)
* [samtools](https://github.com/samtools/samtools)
* Other dependencies: gunzip, wget

Then run `setupsnake.py`
	snakemake -s setupsnake.py -j 



## RUN PIPELINE

Edit `cobrasnake.py` file to reflect the locations of your reference files
and executables that snakemake will use.  The default location is the `lib`
directory where `setupsnake.py` puts them. The variables that you need to examine are:
* HISATREF
* GFFFILE
* filename
* HISAT
* PICARD
* HTSEQ
* SAMTOOLS

Isolate the run_accession IDs from SRA you want to run and put them in a
file line by line.  The pipeline defaults to look for this input file
as `~/accessions.txt`

Run it! You'll probably want to be on a fairly large machine for this (16 cpus)

dependencies: snakemake, hisat (built with SRA support - see note below), samtools, HTSEQ (htseq-count), perl 

	snakemake -s snakemakefile.py -j 

(the -j flag enables multithreaded operations)

You'll want to refer to the snakemake documentation on how to (trivially!)
run snakemake efficiently in a cluster environment that requires
jobsubmission.

<note on building HISAT with SRA support> 



