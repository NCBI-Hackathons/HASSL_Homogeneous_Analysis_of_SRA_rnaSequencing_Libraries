

REF = 'GRCh38.p4'


FASTAREF='/home/ubuntu/russ/ncbi/GCF_000001405.30_GRCh38.p4_genomic.fna'
HISATREF = "/home/ubuntu/refs/hisat_index/GRCh38.p4"
# STAR50REF = 
# STAR100REF = 


# DATASETS = "SRR1295542".split() 
THREADS = 10 



rule all: 
	input: "SRR959265.GRCh38.p4.hisat.sorted.sra"

sample = 'SRR959265'

# SRA -> PILEUP -> RAW COUNTS OFF NCBI GFF3 


rule kar_sra: 
	output: "{sample}.GRCh38.p4.hisat.sorted.sra"
	input: "{sample}.GRCh38.p4.hisat.sorted.sradir"
	message: "converting sradir {input} to sra archive file {output}"
	shell: "kar -c {output} -d {input}"

rule sort_sra:
	output: "{sample}.GRCh38.p4.hisat.sorted.sradir"
	input: "{sample}.GRCh38.p4.hisat.sradir"
	message: "sorting sradir {input} to {output}" 
	shell: "sra-sort {input} {output}"

rule bam_to_sra:
	output: "{sample}.GRCh38.p4.hisat.sradir"
	input: "{sample}.GRCh38.p4.hisat.bam"
	message: "converting bam to sra: {input} to {output}"
	shell: "bam-load -o {sample}.GRCh38.p4.hisat.sradir -k /home/ubuntu/russ/rnaseq_mapping_hackathon_v002/GCF_000001405.30_GRCh38.p4_genomic.cfg -r {FASTAREF} {input}"

# rule sort_bam:
# 	output: "{sample}.GRCh38.p4.hisat.sorted.bam"
# 	input: "{sample}.GRCh38.p4.hisat.bam"
# 	message: "sorting {intput} to {output}"
# 	shell: "samtools sort {input} {output}.sorted"

rule sam_to_bam:
	output: "{sample}.GRCh38.p4.hisat.bam"
	input: "{sample}.GRCh38.p4.hisat.sam"
	message: "converting sam to bam: {input} to {output}"
	shell: "samtools view -bS {input} > {output}"

rule hisat_alignment:
	output: "{sample}.GRCh38.p4.hisat.sam", "{sample}.hisat.log"
	# input: "{sample}.fastq"
	message: "hisat aligning reads {sample}.fastq to GRCh38.p4 with {THREADS} threads to produce {output}"
	shell: "hisat -x {HISATREF} -p {THREADS} --sra-acc {sample} -S {sample}.GRCh38.p4.hisat.sam 2> {sample}.hisat.log"


# rule star_alignment: 
# 	output: "{sample}.GRCh38.p4.star.sam"
# 	input: "{sample}.fastq"
# 	message: "star aligning reads {sample}.fastq to GRCh38.p4 to produce {output}"
# 	shell: <star code> 






