

REF = 'GRCh38.p4'

HISATREF = "/home/ubuntu/russ/ncbi/hisat_index/GRCh38.p4"
# STAR50REF = 
# STAR100REF = 


DATASETS = "SRR1295542".split() 
THREADS = 8 


# rule cufflink_count
# 	output: counts
# 	input: gff3 + sorted.bam 

sample = DATASETS[0] 

rule all: 
	input: "{sample}.{REF}.{alner}.sorted.bam"



# rule index_sam:
# 	output: "{sample}.{REF}.{alner}.sorted.sam.sai"
# 	input: "{sample}.{ref}.{alner}.sorted.sam"
# 	message: "indexing sam {input}"
# 	shell: "{samtools sortsam? }

rule sort_bam:
	output: "{sample}.{REF}.{alner}.sorted.bam"
	input: "{sample}.{REF}.{alner}.bam"
	message: "sorting {intput} to {output}"
	shell: "samtools sort {input} {output}.sorted"

rule sam_to_bam:
	output: "{sample}.{REF}.{alner}.sorted.bam"
	input: "{sample}.{REF}.{alner}.sam"
	message: "converting sam to bam: {input} to {output}"
	shell: "samtools view -bS {input} > {output}"

rule hisat_alignment:
	output: "{sample}.{REF}.hisat.sam"
	input: "{sample}.fastq"
	message: "hisat aligning reads {sample}.fastq to {REF} to produce {output}"
	shell: "hisat -x {HISATREF} -p {THREADS} --sra-acc {sample} -S {sample}.{REF}.hisat.sam"

rule star_alignment: 
	output: "{sample}.{REF}.star.sam"
	input: "{sample}.fastq"
	message: "star aligning reads {sample}.fastq to {REF} to produce {output}"
	shell: <star code> 






