

REF = 'GRCh38.p4'

HISATREF = "/home/ubuntu/refs/hisat_index/GRCh38.p4"
# STAR50REF = 
# STAR100REF = 


# DATASETS = "SRR1295542".split() 
THREADS = 8 



rule all: 
	input: "SRR1295542.GRCh38.p4.hisat.sorted.bam"



# SRA -> PILEUP -> RAW COUNTS OFF NCBI GFF3 


rule kar_sra: 
	output: "{sample}.GRCh38.p4.{alner}.sorted.sra"
	input: "{sample}.GRCh38.p4.{alner}.sorted.sradir"
	message: "converting sradir {input} to sra archive file {output}"
	shell: "kar -c {output} -d {input}"

rule sort_sra:
	output: "{sample}.GRCh38.p4.{alner}.sorted.sradir"
	input: "{sample}.GRCh38.p4.{alner}.sradir"
	message: "sorting sradir {input} to {output}" 
	shell: "sra-sort {input} {output}"

rule bam_to_sra:
	output: "{sample}.GRCh38.p4.{alner}.sradir"
	input: "{sample}.GRCh38.p4.{alner}.bam"
	message: "converting {input} bam to sra {output}"
	shell: "bam-load -o {sample}.GRCh38.p4.{alner}.sradir -k <config-GI-file> {input}"

rule sort_bam:
	output: "{sample}.GRCh38.p4.{alner}.sorted.bam"
	input: "{sample}.GRCh38.p4.{alner}.bam"
	message: "sorting {intput} to {output}"
	shell: "samtools sort {input} {output}.sorted"

rule sam_to_bam:
	output: "{sample}.GRCh38.p4.{alner}.bam"
	input: "{sample}.GRCh38.p4.{alner}.sam"
	message: "converting sam to bam: {input} to {output}"
	shell: "samtools view -bS {input} > {output}"

rule hisat_alignment:
	output: "{sample}.GRCh38.p4.hisat.sam"
	# input: "{sample}.fastq"
	message: "hisat aligning reads {sample}.fastq to GRCh38.p4 with {THREADS} threads to produce {output}"
	shell: "hisat -x {HISATREF} -p {THREADS} --sra-acc {sample} -S {sample}.GRCh38.p4.hisat.sam"


# rule star_alignment: 
# 	output: "{sample}.GRCh38.p4.star.sam"
# 	input: "{sample}.fastq"
# 	message: "star aligning reads {sample}.fastq to GRCh38.p4 to produce {output}"
# 	shell: <star code> 






