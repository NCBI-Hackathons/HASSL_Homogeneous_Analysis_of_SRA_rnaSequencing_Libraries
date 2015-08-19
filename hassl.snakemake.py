# THIS SCRIPT WAS PRODUCED VIA THE NCBI HACKATHON IN AUGUST 2015 
# WRITTEN BY PAUL CANTALUPO, HIROKO OHMIYA, ALLISSA DILLMAN, AND RUSSELL DURRETT


import os 


# HISATREF="/home/ubuntu/resources/ensembl/hisat_indexes/Homo_sapiens.GRCh38.dna.toplevel"
HISATREF="/mnt/hisat_indexes/Homo_sapiens.GRCh38.dna.toplevel"
PICARDFLATFILE="/home/ubuntu/resources/ensembl/GRCh38.77.compatible.ucsc.picard.refflat.txt"
GTFFILE="/home/ubuntu/resources/ensembl/Ensembl.GRCh38.77.gtf"
SPLICEFILE="/home/ubuntu/resources/ensembl/Ensembl.GRCh38.77.splicesites.txt"
#set the number of threads to use in alignments 
THREADS=5

# EXECUTABLE LOCATIONS (some on path)
HISAT=" /home/ubuntu/install/hisat/hisat "
PICARD=" java -jar /home/ubuntu/install/picard-tools-1.138/picard.jar "
FEATURECOUNTS="/home/ubuntu/install/subread-1.4.6-p4-Linux-x86_64/bin/featureCounts"
SAMTOOLS=" /home/ubuntu/install/samtools_rocksdb/samtools/samtools "


#set the filename of the file with the list of accessions   
try:
  config["ACCESSION_FILE"]
except KeyError: 
  filename = "ACCESSIONS"
else:
  filename = config["ACCESSION_FILE"]
  

SAMPLES_FROM_FILE = [line.rstrip('\n') for line in open(filename)]
SAMPLES = [s for s in SAMPLES_FROM_FILE if s]

if config["WORKING_DIR"]:
  WORKING_DIR=config["WORKING_DIR"]
else: 
  WORKING_DIR=os.getcwd()



rule all: 
  input: expand("{sample}.GRCh38.ens77.featureCounts.counts", sample=SAMPLES)  #, expand("{sample}.qc_check.done", sample=SAMPLES)

# rule all_on_s3: 
#   input: expand("{S3_BUCKET}/{sample}.GRCh38.ens77.featureCounts.counts", sample=SAMPLES, S3_BUCKET=S3_BUCKET)  #, expand("{sample}.qc_check.done", sample=SAMPLES)


rule clean: 
  shell: "rm -fr log qc bams counts  "



rule project_counts:
  output: "project.featureCounts"
  input: GTFFILE, expand("bams/{sample}.GRCh38.ens77.hisat.sorted.bam", sample=SAMPLES)
  log: "log/project.featureCounting.log"
  threads: THREADS 
  message: "performing overall project feature counting"
  shell: "{FEATURECOUNTS} -T {threads} --primary -F GTF -t exon -g gene_id -a {GTFFILE} -o project.featureCounts bams/*.bam 2>&1 > log/project.featureCounting.log"

rule perform_counting: 
  output: "counts/{sample}.GRCh38.ens77.featureCounts.counts"
  input: "bams/{sample}.GRCh38.ens77.hisat.sorted.bam"
  log: "log/{wildcards.sample}.featureCounts.log"
  threads: THREADS
  message: "performing featureCounting with {threads} threads on genes in {input}"
  shell: "{FEATURECOUNTS}  -T {threads} --primary -F GTF -t exon -g gene_id -a {GTFFILE} -o counts/{wildcards.sample}.GRCh38.ens77.featureCounts.counts bams/{wildcards.sample}.GRCh38.ens77.hisat.sorted.bam  2> log/{wildcards.sample}.featureCounts.log"

rule qc_check: 
  output: touch("log/{sample}.qc_check.done")
  input: "qc/{sample}.GRCh38.ens77.hisat.crsm"
  log: "log/{wildcards.sample}.perlqc.log"
  message: "checking quality stats of {input} with perl script"
  shell: " perl qc.pl --maplogfile log/{wildcards.sample}.hisat.log --metricsfile qc/{wildcards.sample}.GRCh38.ens77.hisat.crsm --sra {wildcards.sample} 2> log/{wildcards.sample}.perl_qc.log"

rule picard_rnaseq_qual: 
  output: "qc/{sample}.GRCh38.ens77.hisat.crsm"
  input: "bams/{sample}.GRCh38.ens77.hisat.sorted.bam.bai"
  log: "log/{wildcards.sample}.picard_rnametrics.log"
  message: "running picard rna qc stats on {input}"
  shell: "{PICARD} CollectRnaSeqMetrics REF_FLAT={PICARDFLATFILE} STRAND=NONE INPUT={wildcards.sample}.GRCh38.ens77.hisat.sorted.bam OUTPUT=qc/{output} 2> log/{wildcards.sample}.picard_rnametrics.log"

rule index_bam: 
  output: "bams/{sample}.GRCh38.ens77.hisat.sorted.bam.bai"
  input: "bams/{sample}.GRCh38.ens77.hisat.sorted.bam"
  message: "indexing bam alignment file {input}"
  shell: " {SAMTOOLS} index {input} {output} "

rule sort_bam:
  output: "bams/{sample}.GRCh38.ens77.hisat.sorted.bam"
  input: "bams/{sample}.GRCh38.ens77.hisat.bam"
  threads: THREADS
  message: "sorting {input} to {output}"
  shell: " {SAMTOOLS} sort -@ {threads} {input} bams/{wildcards.sample}.GRCh38.ens77.hisat.sorted "

rule sam_to_bam:
  output: temp("bams/{sample}.GRCh38.ens77.hisat.bam")
  input: "bams/{sample}.GRCh38.ens77.hisat.sam"
  message: "converting sam to bam: {input} to {output}"
  shell: " {SAMTOOLS} view -bS {input} > {output} "

rule hisat_alignment:
  output: temp("bams/{sample}.GRCh38.ens77.hisat.sam")
  input: SPLICEFILE 
  threads: THREADS
  log: "log/{sample}.hisat.alignment.log"
  message: "running second pass hisat alignment on {wildcards.sample} with {threads} threads"
  shell: "{HISAT} -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 -x {HISATREF} -p {threads} --sra-acc {wildcards.sample} -t --known-splicesite-infile {SPLICEFILE} -S bams/{wildcards.sample}.GRCh38.ens77.hisat.sam  2> {log}"



