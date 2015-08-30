# HASSL - Homogeneous Analysis of SRA Sequencing Libraries

[![Join the chat at https://gitter.im/DCGenomics/HASSL_Homogeneous_Analysis_of_SRA_rnaSequencing_Libraries](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/DCGenomics/HASSL_Homogeneous_Analysis_of_SRA_rnaSequencing_Libraries?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

HASSL is a snakemake pipeline designed to take RNA-Seq SRA run accession numbers and
output both sorted bam alignments and raw read counts for each run quickly and uniformly. 

## For the lazy (not a hassle if you have all the dependencies installed ;)
```
snakemake -s hassl.py resources -j
echo SRR1200675 > ~/accessions.txt
snakemake -s hassl.py -j --config ACCESSION_FILE='~/accessions.txt'
```

This will produce the following output files:

* SRR1200675.GRCh38.p4.hisat.sorted.bam - HISAT 2-pass mapped BAM file of the reads from SRR1200675
* SRR1200675.hisat.novel.splicesites.txt - HISAT generated file of denovo splicesites
* SRR1200675.hisat.one.log - HISAT log file from the first pass
* SRR1200675.hisat.two.log - HISAT log file form the second pass
* SRR1200675.GRCh38.p4.hisat.crsm - Picard CollectRnaSeqMetrics output file
* SRR1200675.pass - based on the Picard output file, this BAM file passed QC cutoffs
* SRR1200675.GRCh38.p4.featureCounts.counts - featureCounts raw count file


## Setup your environment

HASSL includes a handy resources tool to gather references and build the
hisat index for you.  This will put the GRCh38 reference and annotation
files in the designated reference directory and build the hisat index there as well.

You will need to install the following before running `setupsnake.py`
* [snakemake](https://bitbucket.org/johanneskoester/snakemake/wiki/Documentation#markdown-header-installation)
* [HISAT](https://github.com/infphilo/hisat/) - follow the HISAT directions to compile it with SRA support 
* [featureCounts](http://subread.sourceforge.net/)
* [Picard](https://broadinstitute.github.io/picard/)
* [samtools](https://github.com/samtools/samtools)
* Other dependencies: gunzip, wget

Then run setup:
`snakemake -s hassl.py resources -j `

Be aware of the computational resources required to get these files; unpacking the reference fasta file will take >30G HDD and indexing will need at least >10G RAM. 


## Run the pipeline

Edit `hassl.py` file to reflect the locations of your reference files
and executables that HASSL will use.  The default location at the `/mnt`
directory is where `hassl.py resources` will try put them. Also, you may want to adjust the
threads to be equal to or less than the number of threads on your computer. 
The variables that you need to examine are:

* HISATREF
* GTFFILE
* THREADS
* HISAT
* PICARD
* featureCounts
* SAMTOOLS

Isolate the run_accession IDs from SRA you want to run and put them in a
file line by line.  Do not leave any lines blank.  The pipeline defaults to
look for this input file as `~/accessions.txt`.

Run it! You'll probably want to be on a fairly large machine for this (16 cpus)

dependencies: snakemake, hisat (built with SRA support; see [hisat
manual](https://github.com/infphilo/hisat/blob/master/MANUAL.markdown)),
samtools, featureCounts (from the subread package), perl


You'll want to refer to the snakemake documentation on how to (trivially!)
run snakemake efficiently in a cluster environment that requires
job submission.

