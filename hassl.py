# THIS SCRIPT WAS PRODUCED VIA THE NCBI HACKATHON IN AUGUST 2015 
# WRITTEN BY PAUL CANTALUPO, HIROKO OHMIYA, ALLISSA DILLMAN, AND RUSSELL DURRETT


import os 

#set the number of threads to use for alignment and feature counting 
THREADS=12

# USE ABSOLUTE PATHS!
REFERENCE_DIR="/scratch/03505/russd/resources"
HISAT_REFERENCE_DIR = REFERENCE_DIR + "/hisat_indexes"
HISATREF_BASENAME = "Homo_sapiens.GRCh38.dna_rm.toplevel"   # REPEAT MASKED FASTA
#HISATREF_BASENAME = "Homo_sapiens.GRCh38.dna.toplevel"     # UNMASKED FASTA
FASTA_URL="ftp://ftp.ensembl.org/pub/release-81/fasta/homo_sapiens/dna/" + HISATREF_BASENAME + ".fa.gz"

PICARDFLATFILE_NAME = "GRCh38.77.compatible.ucsc.picard.refflat.txt"
GTFFILE_NAME = "Homo_sapiens.GRCh38.81.gtf"
SPLICEFILE_NAME = "Ensembl.GRCh38.77.splicesites.txt"

HISATREF=HISAT_REFERENCE_DIR + "/" + HISATREF_BASENAME
TOPHAT_TRANSCRIPT_INDEX="/scratch/03505/russd/resources/tophat_indexes/Homo_sapiens.GRCh38.dna_rm.toplevel.GTF81.transcriptome.index"

PICARDFLATFILE=REFERENCE_DIR+ "/" + PICARDFLATFILE_NAME
GTFFILE=REFERENCE_DIR+ "/" + GTFFILE_NAME
SPLICEFILE=REFERENCE_DIR+ "/" + SPLICEFILE_NAME

REFERENCE_BASE_URL="https://s3.amazonaws.com/genomicdata/HASSL"

PICARDFLATFILE_URL=REFERENCE_BASE_URL+ "/" + PICARDFLATFILE_NAME
GTFFILE_URL=REFERENCE_BASE_URL+ "/" + GTFFILE_NAME + ".gz"
SPLICEFILE_URL=REFERENCE_BASE_URL+ "/" + SPLICEFILE_NAME


# EXECUTABLE LOCATIONS (some on path probably)
BIN="/home1/03505/russd/bin"
HASSL=" /home1/03505/russd/software/HASSL"
HISAT=" /home1/03505/russd/software//hisat "
HISAT_BUILD="/home1/03505/russd/bin/hisat-build"
PICARD=" java -jar /home/ubuntu/install/picard-tools-1.138/picard.jar "
FEATURECOUNTS="/home/ubuntu/install/subread-1.4.6-p4-Linux-x86_64/bin/featureCounts"
SAMTOOLS_ROCKS=" samtools "
SAMTOOLS=" samtools "
BOWTIE_BUILD= "bowtie2-build "
VDBDUMP=" vdb-dump -f fastq "
TOPHAT=" tophat2 "

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
  input: expand("bams/{sample}.GRCh38.ens77.hisat.sorted.bam", sample=SAMPLES) #, expand("log/{sample}.qc_check.done", sample=SAMPLES)


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

rule tophat_alignment: 
  output: temp("bams/{sample}.GRCh38.ens77.hisat.sam")
  input: "/scratch/03505/russd/{sample}.fastq", HISAT_REFERENCE_DIR + "/" + HISATREF_BASENAME + ".rev.2.bt2" #, "bams/{sample}.GRCh38.ens77.hisat.temp.sam", "splicesites/{sample}.novel.splicesites"
  threads: THREADS
  log: "log/{sample}.tophat2.log"
  message: "running tophat alignment on {wildcards.sample} with {threads} threads"
  shell: " {TOPHAT} -p {threads} -G {GTFFILE} --transcriptome-index {TOPHAT_TRANSCRIPT_INDEX} -z pigz  -o bams/{wildcards.sample}.GRCh38.ens77.hisat.sam  {HISATREF}  {wildcards.sample}.fastq   2> {log}"   

rule vdbdump_reads: 
  output: temp("/scratch/03505/russd/{sample}.fastq")
  message: "downloading reads for {sample}"
  shell: " {VDBDUMP} {wildcards.sample} > /scratch/03505/russd/{wildcards.sample}.fastq "
# rule hisat_alignment:
#   output: temp("bams/{sample}.GRCh38.ens77.hisat.sam")
#   input: HISAT_REFERENCE_DIR + "/" + HISATREF_BASENAME + ".rev.2.bt2l" #, "bams/{sample}.GRCh38.ens77.hisat.temp.sam", "splicesites/{sample}.novel.splicesites"
#   threads: THREADS
#   log: "log/{sample}.hisat.log"
#   message: "running hisat alignment on {wildcards.sample} with {threads} threads"
#   shell: "{HISAT} -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 -x {HISATREF} -p {threads} --sra-acc {wildcards.sample} -t --known-splicesite-infile {SPLICEFILE} -S bams/{wildcards.sample}.GRCh38.ens77.hisat.sam  2> {log}"   # --novel-splicesite-infile splicesites/{wildcards.sample}.novel.splicesites 

# HISAT hangs if you do two pass alignment!???!!! Must Fix 
# rule hisat_alignment_one:
#   output: temp("bams/{sample}.GRCh38.ens77.hisat.temp.sam"), "splicesites/{sample}.novel.splicesites"
#   input: HISAT_REFERENCE_DIR + "/" + HISATREF_BASENAME + ".rev.2.bt2l"
#   threads: THREADS
#   log: "log/{sample}.hisat.log"
#   message: "running primary hisat alignment on {wildcards.sample} with {threads} threads"
#   shell: "{HISAT} -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 -x {HISATREF} -p {threads} --sra-acc {wildcards.sample} -t --known-splicesite-infile {SPLICEFILE} --novel-splicesite-outfile splicesites/{wildcards.sample}.novel.splicesites  -S bams/{wildcards.sample}.GRCh38.ens77.hisat.temp.sam  2> {log}"




#SCRIPT TO GATHER REFERENCES AND ANNOTATION FILES FOR COBRASNAKE RNASEQ PIPELINE

rule resources:
  input:  [PICARDFLATFILE, GTFFILE, SPLICEFILE, HISAT_REFERENCE_DIR+"/"+HISATREF_BASENAME+".rev.2.bt2"]

rule bowtie2_index:
  output: HISAT_REFERENCE_DIR + "/" + HISATREF_BASENAME + ".rev.2.bt2"
  input: REFERENCE_DIR + "/" + HISATREF_BASENAME + ".fa"
  message: "bowtie2 indexing human genome {input}"
  shell: "{BOWTIE_BUILD} {input} {HISAT_REFERENCE_DIR}/{HISATREF_BASENAME}"

rule hisat_index:
  output: HISAT_REFERENCE_DIR + "/" + HISATREF_BASENAME + ".rev.2.bt2"
  input: REFERENCE_DIR + "/" + HISATREF_BASENAME + ".fa"
  message: "hisat-build indexing human genome {input}"
  shell: "{HISAT_BUILD} {input} {HISAT_REFERENCE_DIR}/{HISATREF_BASENAME}"

rule gunzip_reference_fasta:
  output: REFERENCE_DIR + "/" + HISATREF_BASENAME + ".fa"
  input: REFERENCE_DIR + "/" + HISATREF_BASENAME + ".fa.gz"
  message: "extracting human genome fasta {input}"
  shell: "gunzip -c {input} > {output}"
  
rule get_reference_fasta:
  output: temp(REFERENCE_DIR + "/" + HISATREF_BASENAME + ".fa.gz")
  message: "downloading human reference genome from {FASTA_URL}"
  shell: "curl -o {REFERENCE_DIR}/{HISATREF_BASENAME}.fa.gz {FASTA_URL}"


rule get_splicesites:
  output: SPLICEFILE
  message: "downloading splicesites from {SPLICEFILE_URL}"
  shell: "curl -o {REFERENCE_DIR}/{SPLICEFILE_NAME} {SPLICEFILE_URL}"

rule gunzip_gtf: 
  output: GTFFILE
  input: GTFFILE + ".gz"
  message: "gunzipping GTF file {input}"
  shell: "gunzip -c {input} > {output}"

rule get_gtf:
  output: GTFFILE + ".gz"
  message: "downloading GTF from {GTFFILE_URL}"
  shell: "curl -o {REFERENCE_DIR}/{GTFFILE_NAME} {GTFFILE_URL}"
  
  
rule get_refflat:
  output: PICARDFLATFILE
  message: "downloading refflat from {PICARDFLATFILE_URL}"
  shell: "curl -o {REFERENCE_DIR}/{PICARDFLATFILE_NAME} {PICARDFLATFILE_URL}"




