# THIS SCRIPT WAS PRODUCED VIA THE NCBI HACKATHON IN AUGUST 2015 
# WRITTEN BY PAUL CANTALUPO, HIROKO OHMIYA, ALLISSA DILLMAN, AND RUSSELL DURRETT


import os 

#set the number of threads to use in alignments and sorting
THREADS=3

REFERENCE_BASE_URL="https://s3.amazonaws.com/genomicdata/HASSL"
REFERENCE_DIR="/mnt/hassl/lib"
HISAT_REFERENCE_DIR = REFERENCE_DIR + "/hisat_indexes"
HISATREF_BASENAME = "Homo_sapiens.GRCh38.dna.toplevel"
PICARDFLATFILE_NAME = "GRCh38.77.compatible.ucsc.picard.refflat.txt"
GTFFILE_NAME = "Ensembl.GRCh38.77.gtf"
SPLICEFILE_NAME = "Ensembl.GRCh38.77.splicesites.txt"

#HISATREF=HISAT_REFERENCE_DIR + "/" + HISATREF_BASENAME
HISATREF=REFERENCE_DIR + "/" + HISATREF_BASENAME
PICARDFLATFILE=REFERENCE_DIR+ "/" + PICARDFLATFILE_NAME
GTFFILE=REFERENCE_DIR+ "/" + GTFFILE_NAME
SPLICEFILE=REFERENCE_DIR+ "/" + SPLICEFILE_NAME
PICARDFLATFILE_URL=REFERENCE_BASE_URL+ "/" + PICARDFLATFILE_NAME
GTFFILE_URL=REFERENCE_BASE_URL+ "/" + GTFFILE_NAME
SPLICEFILE_URL=REFERENCE_BASE_URL+ "/" + SPLICEFILE_NAME


# EXECUTABLE LOCATIONS (some on path)
HASSL=" /home/ubuntu/HASSL"
HISAT=" /home/ubuntu/install/hisat/hisat "
HISAT_BUILD="~/install/hisat/hisat-build"
PICARD=" java -jar /home/ubuntu/install/picard-tools-1.138/picard.jar "
FEATURECOUNTS="/home/ubuntu/install/subread-1.4.6-p4-Linux-x86_64/bin/featureCounts"
SAMTOOLS_ROCKS=" /home/ubuntu/install/samtools_rocksdb/samtools/samtools "
SAMTOOLS=" /home/ubuntu/install/samtools/samtools"

#set the filename of the file with the list of accessions   
try:
  config["ACCESSION_FILE"]
except KeyError: 
  filename = "accessions.txt"
else:
  filename = config["ACCESSION_FILE"]


#SAMPLES_FROM_FILE = [line.rstrip('\n') for line in open(filename)]
#SAMPLES = [s for s in SAMPLES_FROM_FILE if s]
SAMPLES = [line.rstrip('\n') for line in open(filename)]

try: 
  config["WORKING_DIR"]
except KeyError: 
  WORKING_DIR=os.getcwd()
else: 
  WORKING_DIR=config["WORKING_DIR"]



rule all: 
  input: expand("counts/{sample}.GRCh38.ens77.featureCounts.counts", sample=SAMPLES) , expand("log/{sample}.qc_check.done", sample=SAMPLES)


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
  shell: "{PICARD} CollectRnaSeqMetrics REF_FLAT={PICARDFLATFILE} STRAND=NONE INPUT={input} OUTPUT={output} 2> {log}"

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
  shell: " {SAMTOOLS_ROCKS} sort -@ {threads} {input} bams/{wildcards.sample}.GRCh38.ens77.hisat.sorted "

rule sam_to_bam:
  output: temp("bams/{sample}.GRCh38.ens77.hisat.bam")
  input: "bams/{sample}.GRCh38.ens77.hisat.sam"
  message: "converting sam to bam: {input} to {output}"
  shell: " {SAMTOOLS} view -bS {input} > {output} "

rule hisat_alignment:
  output: temp("bams/{sample}.GRCh38.ens77.hisat.sam")
  input: HISAT_REFERENCE_DIR + "/" + HISATREF_BASENAME + ".rev.2.bt2"
  threads: THREADS
  log: "log/{sample}.hisat.log"
  message: "running second pass hisat alignment on {wildcards.sample} with {threads} threads"
  shell: "{HISAT} -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 -x {HISATREF} -p {threads} --sra-acc {wildcards.sample} -t -S bams/{wildcards.sample}.GRCh38.ens77.hisat.sam  2> {log}"





#SCRIPT TO GATHER REFERENCES AND ANNOTATION FILES FOR COBRASNAKE RNASEQ PIPELINE

rule resources:
  input:  [PICARDFLATFILE, GTFFILE, SPLICEFILE, "/mnt/hassl/lib/GRCh38.p4.rev.6.bt2"]

    # "{REFERENCE_DIR}/GCF_000001405.30_GRCh38.p4_genomic.gff",
    # "{REFERENCE_DIR}/GCF_000001405.30_GRCh38.p4_genomic.fna",
    # "{REFERENCE_DIR}/GRCh38.p4.1.bt2", "{REFERENCE_DIR}/GRCh38.p4.3.bt2", "{REFERENCE_DIR}/GRCh38.p4.5.bt2", "{REFERENCE_DIR}/GRCh38.p4.rev.1.bt2", "{REFERENCE_DIR}/GRCh38.p4.rev.5.bt2",
    #       "{REFERENCE_DIR}/GRCh38.p4.2.bt2", "{REFERENCE_DIR}/GRCh38.p4.4.bt2", "{REFERENCE_DIR}/GRCh38.p4.6.bt2", "{REFERENCE_DIR}/GRCh38.p4.rev.2.bt2", "{REFERENCE_DIR}/GRCh38.p4.rev.6.bt2"

# Reference genome (h38)
rule index:
  output: "{REFERENCE_DIR}/GRCh38.p4.1.bt2", "{REFERENCE_DIR}/GRCh38.p4.3.bt2", "{REFERENCE_DIR}/GRCh38.p4.5.bt2", "{REFERENCE_DIR}/GRCh38.p4.rev.1.bt2", "{REFERENCE_DIR}/GRCh38.p4.rev.5.bt2",
          "{REFERENCE_DIR}/GRCh38.p4.2.bt2", "{REFERENCE_DIR}/GRCh38.p4.4.bt2", "{REFERENCE_DIR}/GRCh38.p4.6.bt2", "{REFERENCE_DIR}/GRCh38.p4.rev.2.bt2", "{REFERENCE_DIR}/GRCh38.p4.rev.6.bt2"
  input: "{REFERENCE_DIR}/GCF_000001405.30_GRCh38.p4_genomic.fna"
  message: "indexing human genome"
  shell: "{HISAT_BUILD} {input} {REFERENCE_DIR}/GRCh38.p4"

rule gunzipHumRef:
  output: "{REFERENCE_DIR}/GCF_000001405.30_GRCh38.p4_genomic.fna" 
  input: "{REFERENCE_DIR}/GCF_000001405.30_GRCh38.p4_genomic.fna.gz"
  message: "extracting human genome fasta.gz"
  shell: "gunzip -c {input} > {output}"
  
rule get_humanreference:
  output: temp("{REFERENCE_DIR}/GCF_000001405.30_GRCh38.p4_genomic.fna.gz")
  message: "downloading human reference genome"
  shell: "wget -P {REFERENCE_DIR} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000001405.30_GRCh38.p4/GCF_000001405.30_GRCh38.p4_genomic.fna.gz"


rule get_splicesites:
  output: REFERENCE_DIR + "/" + SPLICEFILE_NAME
  message: "downloading splicesites from s3"
  shell: "wget -P {REFERENCE_DIR} {SPLICEFILE_URL}"


rule get_gtf:
  output: REFERENCE_DIR + "/" + GTFFILE_NAME
  message: "downloading GTF from s3"
  shell: "wget -P {REFERENCE_DIR} {REFERENCE_BASE_URL}/{GTFFILE_NAME}"
  
  
rule get_refflat:
  output: REFERENCE_DIR + "/" + PICARDFLATFILE_NAME
  message: "downloading refflat from s3"
  shell: "wget -P {REFERENCE_DIR} {PICARDFLATFILE_URL}"




