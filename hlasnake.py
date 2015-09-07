# THIS SCRIPT WAS PRODUCED VIA THE NCBI HACKATHON IN AUGUST 2015 
# WRITTEN BY PAUL CANTALUPO, HIROKO OHMIYA, ALLISSA DILLMAN, AND RUSSELL DURRETT


import os 

#set the number of threads to use for alignment and feature counting 
THREADS=3

# USE ABSOLUTE PATHS!
REFERENCE_DIR="/mnt/resources"
HISAT_REFERENCE_DIR = REFERENCE_DIR + "/hisat_indexes"
HISATREF_BASENAME = "HLA12"   # REPEAT MASKED FASTA
#HISATREF_BASENAME = "Homo_sapiens.GRCh38.dna.toplevel"     # UNMASKED FASTA
FASTA_URL="ftp://ftp.ensembl.org/pub/release-81/fasta/homo_sapiens/dna/" + HISATREF_BASENAME + ".fa.gz"

PICARDFLATFILE_NAME = "GRCh38.77.compatible.ucsc.picard.refflat.txt"
GTFFILE_NAME = "Homo_sapiens.GRCh38.81.gtf"
SPLICEFILE_NAME = "Ensembl.GRCh38.77.splicesites.txt"

HISATREF=HISAT_REFERENCE_DIR + "/" + HISATREF_BASENAME
PICARDFLATFILE=REFERENCE_DIR+ "/" + PICARDFLATFILE_NAME
GTFFILE=REFERENCE_DIR+ "/" + GTFFILE_NAME
SPLICEFILE=REFERENCE_DIR+ "/" + SPLICEFILE_NAME

REFERENCE_BASE_URL="https://s3.amazonaws.com/genomicdata/HASSL"

PICARDFLATFILE_URL=REFERENCE_BASE_URL+ "/" + PICARDFLATFILE_NAME
GTFFILE_URL=REFERENCE_BASE_URL+ "/" + GTFFILE_NAME + ".gz"
SPLICEFILE_URL=REFERENCE_BASE_URL+ "/" + SPLICEFILE_NAME


# EXECUTABLE LOCATIONS (some on path probably)
HASSL=" /home/ubuntu/HASSL"
HISAT=" /home/ubuntu/install/hisat/hisat "
HISAT_BUILD="~/install/hisat/hisat-build"
PICARD=" java -jar /home/ubuntu/install/picard-tools-1.138/picard.jar "
FEATURECOUNTS="/home/ubuntu/install/subread-1.4.6-p4-Linux-x86_64/bin/featureCounts"
SAMTOOLS_ROCKS=" /home/ubuntu/install/samtools_rocksdb/samtools/samtools "
SAMTOOLS=" /home/ubuntu/install/samtools/samtools "
FASTQDUMP=" /home/ubuntu/install/sratoolkit.2.4.2-ubuntu64/bin/fastq-dump "
VDBDUMP=" /home/ubuntu/install/sratoolkit.2.4.2-ubuntu64/bin/vdb-dump -f fastq "

#set the filename of the file with the list of accessions   
try:
  config["ACCESSION_FILE"]
except KeyError: 
  filename = "accessions.txt"
else:
  filename = config["ACCESSION_FILE"]

try: 
  open(filename)
except FileNotFoundError:
  SAMPLES = []
else: 
  SAMPLES = [line.rstrip('\n') for line in open(filename)]


#SAMPLES_FROM_FILE = [line.rstrip('\n') for line in open(filename)]
#SAMPLES = [s for s in SAMPLES_FROM_FILE if s]

try: 
  config["WORKING_DIR"]
except KeyError: 
  WORKING_DIR=os.getcwd()
else: 
  WORKING_DIR=config["WORKING_DIR"]



rule all: 
  input: expand("/home/ubuntu/hla/sams/{sample}.HLA12.hisat.sam", sample=SAMPLES) 




rule clean: 
  shell: "rm -fr log qc bams counts project.featureCounts* "

rule project_counts:
  output: "project.featureCounts"
  input: GTFFILE, expand("bams/{sample}.GRCh38.ens77.hisat.sorted.bam", sample=SAMPLES), GTFFILE
  log: "log/project.featureCounting.log"
  threads: THREADS
  message: "performing overall project feature counting"
  shell: "{FEATURECOUNTS} -T {threads} --primary -F GTF -t exon -g gene_id -a {GTFFILE} -o project.featureCounts bams/*.bam 2>&1 > log/project.featureCounting.log"

rule perform_counting: 
  output: "counts/{sample}.GRCh38.ens77.featureCounts.counts"
  input: "bams/{sample}.GRCh38.ens77.hisat.sorted.bam", GTFFILE
  log: "log/{sample}.featureCounts.log"
  threads: THREADS
  message: "performing featureCounting with {threads} threads on genes in {input}"
  shell: "{FEATURECOUNTS}  -T {threads} --primary -F GTF -t exon -g gene_id -a {GTFFILE} -o counts/{wildcards.sample}.GRCh38.ens77.featureCounts.counts bams/{wildcards.sample}.GRCh38.ens77.hisat.sorted.bam 2> {log}"

rule qc_check: 
  output: touch("log/{sample}.qc_check.done")
  input: "qc/{sample}.GRCh38.ens77.hisat.crsm"
  log: "log/{sample}.perl_qc.log"
  message: "checking quality stats of {input} with perl script"
  shell: " perl {HASSL}/scripts/qc.pl --maplogfile log/{wildcards.sample}.hisat.log --metricsfile qc/{wildcards.sample}.GRCh38.ens77.hisat.crsm --sra {wildcards.sample} 2> {log}"

rule picard_rnaseq_qual: 
  output: "qc/{sample}.GRCh38.ens77.hisat.crsm"
  input: "bams/{sample}.GRCh38.ens77.hisat.sorted.bam", PICARDFLATFILE
  log: "log/{sample}.picard_rnametrics.log"
  message: "running picard rna qc stats on {input}"
  shell: "{PICARD} CollectRnaSeqMetrics REF_FLAT={PICARDFLATFILE} STRAND=NONE INPUT={input[0]} OUTPUT={output} 2> {log}"


rule cp_bam: 
  output: "{OUTDIR}/{sample}.GRCh38.ens77.hisat.sorted.bam"
  input: "bams/{sample}.GRCh38.ens77.hisat.sorted.bam"
  priority: 100 
  log: "log/{sample}.cp_bam.du.log"
  message: "copying bam {input} {output}"
  shell: "cp {input} {output} && du {input} {output} > {log}"


rule index_bam: 
  output: temp("bams/{sample}.GRCh38.ens77.hisat.sorted.bam.bai")
  input: "bams/{sample}.GRCh38.ens77.hisat.sorted.bam"
  message: "indexing bam alignment file {input}"
  shell: " {SAMTOOLS} index {input} {output} "

rule sort_bam:
  output: temp("bams/{sample}.GRCh38.ens77.hisat.sorted.bam")
  input: "bams/{sample}.GRCh38.ens77.hisat.bam"
  threads: THREADS
  message: "sorting {input} to {output}"
  shell: " {SAMTOOLS_ROCKS} sort -@ {threads} {input} bams/{wildcards.sample}.GRCh38.ens77.hisat.sorted "

rule sam_to_bam:
  output: temp("bams/{sample}.GRCh38.ens77.hisat.bam")
  input: "bams/{sample}.GRCh38.ens77.hisat.sam"
  message: "converting sam to bam: {input} to {output}"
  shell: " {SAMTOOLS} view -bS {input} > {output} "


rule hisat_alignment:
  output: temp("bams/{sample}.GRCh38.ens77.hisat.sam")
  input: HISAT_REFERENCE_DIR + "/" + HISATREF_BASENAME + ".rev.2.bt2l" #, "bams/{sample}.GRCh38.ens77.hisat.temp.sam", "splicesites/{sample}.novel.splicesites"
  threads: THREADS
  log: "log/{sample}.hisat.log"
  message: "running hisat alignment on {wildcards.sample} with {threads} threads"
  shell: "  {VDBDUMP}  {wildcards.sample}  |   {HISAT} -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 -x {HISATREF} -p {threads} -U -  -t --known-splicesite-infile {SPLICEFILE} -S bams/{wildcards.sample}.GRCh38.ens77.hisat.sam  2> {log}"   # --novel-splicesite-infile splicesites/{wildcards.sample}.novel.splicesites 



rule hisat_aligned_only_to_sam:
  output: "{TARGET_DIR}/{sample}.HLA12.hisat.sam"
  input: HISAT_REFERENCE_DIR + "/" + HISATREF_BASENAME + ".rev.2.bt2" #, "bams/{sample}.GRCh38.ens77.hisat.temp.sam", "splicesites/{sample}.novel.splicesites"
  threads: THREADS
  log: "log/{sample}.hisat.log"
  message: "running hisat alignment on {wildcards.sample} with {threads} threads"
  shell: "  {VDBDUMP}  {wildcards.sample}  |   {HISAT} --very-sensitive -x {HISATREF} -p {threads} -U -  | samtools view -S - -F 4 > {TARGET_DIR}/{wildcards.sample}.HLA12.hisat.sam  2> {log}"   # --novel-splicesite-infile splicesites/{wildcards.sample}.novel.splicesites 



