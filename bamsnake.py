

REF = 'GRCh38.p4'


FASTAREF='/home/ubuntu/russ/ncbi/GCF_000001405.30_GRCh38.p4_genomic.fna'
HISATREF = "/home/ubuntu/refs/hisat_index/GRCh38.p4"

# DATASETS = "SRR1295542".split() 
THREADS = 10 

samples = 'SRR959265'.split()

rule all: 
	input: expand("{sample}.transferred", sample=SAMPLES)


	# input: "dSRR959265.GRCh38.p4.hisat.sorted.bam"


# SRA -> PILEUP -> RAW COUNTS OFF NCBI GFF3 
rule transfer_s3: 
	output: touch("{sample}.transferred")
	input: "{sample}.GRCh38.p4.hisat.sorted.bam"
	message: "transferring {input} to S3"
	shell: "s3cmd put {input} s3://ncbi-hackathon-aug/rnamapping/"

rule sort_bam:
	output: "{sample}.GRCh38.p4.hisat.sorted.bam"
	input: "{sample}.GRCh38.p4.hisat.bam"
	message: "sorting {intput} to {output}"
	shell: "samtools sort {input} {output}.sorted"

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






