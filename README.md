# HASSL - Homogeneous Analysis of SRA Sequencing Libraries

[![Join the chat at https://gitter.im/DCGenomics/HASSL_Homogeneous_Analysis_of_SRA_rnaSequencing_Libraries](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/DCGenomics/HASSL_Homogeneous_Analysis_of_SRA_rnaSequencing_Libraries?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

HASSL is a snakemake pipeline designed to take RNA-Seq SRA run accession numbers and
output both sorted bam alignments and raw read counts for each run quickly and uniformly. 

## For the lazy (not a hassle if you have all the dependencies installed ;)
```
snakemake -s hassl.py resources -j
echo SRR1200675 > accessions.txt
snakemake -s hassl.py -j --config ACCESSION_FILE='accessions.txt'
```

This will produce the following output files (@TODO needs updated for new release)

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
files in a reference directory (default is `/mnt/resources`) and build the hisat index there as well. Make sure you have write access to this directory. If you want to change this directory, change the `REFERENCE_DIR` variable in `hassl.py` to a directory that you have write access to.

In order to setup your environment, you'll need to install two programs and then edit the `HISAT_BUILD` variable in `hassl.py` so it knows where to find it.
* [snakemake](https://bitbucket.org/johanneskoester/snakemake/wiki/Documentation#markdown-header-installation)
* [hisat](https://github.com/infphilo/hisat/) - follow the HISAT directions to compile it with SRA support 
* Other dependencies: gunzip, wget

Then run hassl to setup your environment: `snakemake -s hassl.py resources -j `

Be aware of the computational resources required to get these files; unpacking the reference fasta file will take >30G HDD and indexing will need at least >10G RAM. 


## Run the pipeline

Before running HASSL, you will need to install the following programs:
* [featureCounts](http://subread.sourceforge.net/)
* [Picard](https://broadinstitute.github.io/picard/)
* [samtools rocks] (https://github.com/dnanexus/samtools)
* [samtools](https://github.com/samtools/samtools)

Edit `hassl.py` to reflect the location of the HASSL install directory (`HASSL`), your reference files, and executables that HASSL will use (see the `EXECUTABLE LOCATIONS` section).  The default location at the `/mnt`
directory is where `hassl.py resources` will try put them. Also, you may want to adjust the
threads (`THREADS`) to be equal to or less than the number of threads on your computer.

Isolate the run_accession IDs from SRA you want to run and put them in a
file line by line.  Do not leave any lines blank.  The pipeline defaults to
look for this input file as `accessions.txt`.

Run it! with the command `snakemake -s hassl.py -j`. You'll probably want to be on a fairly large machine for this (16 cpus). Refer to the snakemake documentation on how to (trivially!) run snakemake efficiently in a cluster environment that requires job submission.

HASSL automatically deletes the BAM files to save space. If you want to keep all you bam and bam.bai files, edit `hassl.py` by removing the `temp()` wrapper around the `output` file name in the `index_bam` and `sort_bam` rules.

### Visualizing QC metrics

After you run HASSL on a set of SRA accessions, collate all the Picard and counts output with the following:
`scripts/collate.pl accessions.txt`. This will create two files `qc/collate.qc.tsv` and `counts/collate.counts.tsv`. Then you can create QC plots by running the following `cd qc; Rscript ../scripts/qc_plots.R` (requires ggplot2 and grid R packages). This will create a `qc_histogram.jpg` file with several graphs of important QC metrics.

### Visualizing count histograms of genes
In order to plot the count distributions of high coverage genes, run the following from the `counts` directory: `Rscript ../scripts/densityplots_genes.R` (requires DESeq2 and ggplot2 R packages).

### Dependencies
* [snakemake](https://bitbucket.org/johanneskoester/snakemake/wiki/Documentation#markdown-header-installation)
* [HISAT](https://github.com/infphilo/hisat/) - follow the HISAT directions to compile it with SRA support 
* [featureCounts](http://subread.sourceforge.net/)
* [Picard](https://broadinstitute.github.io/picard/)
* [samtools rocks] (https://github.com/dnanexus/samtools)
* [samtools](https://github.com/samtools/samtools)
* Other dependencies: gunzip, wget, perl




