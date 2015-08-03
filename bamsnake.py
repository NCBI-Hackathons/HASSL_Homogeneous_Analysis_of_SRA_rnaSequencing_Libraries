

REF = 'GRCh38.p4'
FASTAREF='/home/ubuntu/russ/ncbi/GCF_000001405.30_GRCh38.p4_genomic.fna'
HISATREF = "/home/ubuntu/refs/hisat_index/GRCh38.p4"

# DATASETS = "SRR1295542".split() 
THREADS = 10 

SAMPLES = 'SRR959265'.split()

rule all: 
	input: expand("{sample}.hisat1.log", sample=SAMPLES), expand("{sample}.transferred.logs", sample=SAMPLES), expand("{sample}.transferred", sample=SAMPLES), expand("{sample}.transferred.splices", sample=SAMPLES)


	# input: "dSRR959265.GRCh38.p4.hisat.sorted.bam"


# SRA -> PILEUP -> RAW COUNTS OFF NCBI GFF3 

rule transfer_logs_s3: 
	output: touch("{sample}.transferred.logs")
	input:  "{sample}.hisat1.log", "{sample}.hisat2.log"
	message: "transferring {input}'s logs to S3"
	shell: "s3cmd put {sample}.hisat1.log s3://ncbi-hackathon-aug/rnamapping/; s3cmd put {sample}.hisat2.log s3://ncbi-hackathon-aug/rnamapping/"

rule transfer_bam_s3: 
	output: touch("{sample}.transferred")
	input: "{sample}.GRCh38.p4.hisat.sorted.bam"
	message: "transferring {input} to S3"
	shell: "s3cmd put {input} s3://ncbi-hackathon-aug/rnamapping/"

rule transfer_splices_s3:
	output: touch("{sample}.transferred.splices")
	input: "{sample}.hisat.novel.splicesite.txt"
	message: "transferring hisat splice sites {input} to s3"
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


rule hisat_alignment_two:
	output: "{sample}.GRCh38.p4.hisat.sam", "{sample}.hisat2.log"
	input: "{sample}.hisat.novel.splicesite.txt"
	message: "running second pass alignment"
	shell: "hisat -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 -x {HISATREF} -p {THREADS} --sra-acc {sample} --mm -t -S {sample}.GRCh38.p4.hisat.sam --novel-splicesite-infile {sample}.hisat.novel.splicesite.txt 2> {sample}.hisat2.log"


rule hisat_alignment_one:
	output: "{sample}.hisat.novel.splicesite.txt", "{sample}.hisat1.log", "{sample}.GRCh38.p4.hisat_firstpass_tmp.sam"
	message: "hisat aligning reads from {sample} to GRCh38.p4 with {THREADS} threads to produce new splicesites"
	shell: "hisat -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 -x {HISATREF} -p {THREADS} --sra-acc {sample} --mm -t -S {sample}.GRCh38.p4.hisat_firstpass_tmp.sam --novel-splicesite-outfile {sample}.hisat.novel.splicesite.txt 2> {sample}.hisat1.log"


# rule star_alignment: 
# 	output: "{sample}.GRCh38.p4.star.sam"
# 	input: "{sample}.fastq"
# 	message: "star aligning reads {sample}.fastq to GRCh38.p4 to produce {output}"
# 	shell: <star code> 






