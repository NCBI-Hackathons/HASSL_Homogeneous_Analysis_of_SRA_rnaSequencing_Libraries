# HASSL - Homogeneous Analysis of SRA RNAseq libraries

HASSL is a snakemake pipeline designed to take SRA run accession numbers and
output both sorted bam alignments and raw read counts for an infinite number
of samples.

## For the lazy (not a hassle if you have all the dependencies installed ;)
```
snakemake -s setupsnake.py -j
echo SRR1200675 > ~/accessions.txt
snakemake -s cobrasnake.py -j
```

This will produce the following output files:

* SRR1200675.GRCh38.p4.hisat.sorted.bam - HISAT 2-pass mapped BAM file of the reads from SRR1200675
* SRR1200675.hisat.novel.splicesites.txt - HISAT generated file of denovo splicesites
* SRR1200675.hisat.one.log - HISAT log file from the first pass
* SRR1200675.hisat.two.log - HISAT log file form the second pass
* SRR1200675.GRCh38.p4.hisat.crsm - Picard CollectRnaSeqMetrics output file
* SRR1200675.pass - based on the Picard output file, this BAM file passed QC cutoffs
* SRR1200675.GRCh38.p4.HTSeq.counts - HTSeq raw count file


## Setup your environment

HASSL includes a handy setupsnake.py tool to gather references and build the
hisat index for you.  This will put the GRCh38.p4 reference and annotation
files in the `lib` directory and build the hisat index there as well.

You will need to install the following before running `setupsnake.py`
* [snakemake](https://bitbucket.org/johanneskoester/snakemake/downloads)
* [HISAT](https://github.com/infphilo/hisat/) - follow the HISAT directions to compile it with SRA support 
* [HTSeq](http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html)
* [Picard](https://broadinstitute.github.io/picard/)
* [samtools](https://github.com/samtools/samtools)
* Other dependencies: gunzip, wget

Then run setup:
`snakemake -s setupsnake.py -j `


## Run the pipeline

Edit `cobrasnake.py` file to reflect the locations of your reference files
and executables that HASSL will use.  The default location is the `lib`
directory where `setupsnake.py` puts them.  Also, you may want to adjust the
threads to be equal to or less than the number of threads on your computer. 
The variables that you need to examine are:

* HISATREF
* GFFFILE
* THREADS
* filename
* HISAT
* PICARD
* HTSEQ
* SAMTOOLS

Isolate the run_accession IDs from SRA you want to run and put them in a
file line by line.  Do not leave any lines blank.  The pipeline defaults to
look for this input file as `~/accessions.txt`.

Run it! You'll probably want to be on a fairly large machine for this (16 cpus)

dependencies: snakemake, hisat (built with SRA support; see [hisat
manual](https://github.com/infphilo/hisat/blob/master/MANUAL.markdown)),
samtools, HTSEQ (htseq-count), perl

`snakemake -s snakemakefile.py -j`

(the -j flag enables multithreaded operations)

You'll want to refer to the snakemake documentation on how to (trivially!)
run snakemake efficiently in a cluster environment that requires
jobsubmission.

