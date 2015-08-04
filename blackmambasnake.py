

REF = 'GRCh38.p4'
# FASTAREF='/home/ubuntu/russ/ncbi/GCF_000001405.30_GRCh38.p4_genomic.fna'
HISATREF = "/home/ubuntu/refs/hisat_index/GRCh38.p4"

# DATASETS = "SRR1295542".split() 
# THREADS = 10 

filename = "/home/ubuntu/accessions"
SAMPLES = [line.rstrip('\n') for line in open(filename)]

rule all: 
	input: expand("{sample}.transferred.counts", sample=SAMPLES), expand("{sample}.hisat.novel.splicesites.txt", sample=SAMPLES), expand("{sample}.transferred.log", sample=SAMPLES), expand("{sample}.transferred", sample=SAMPLES), expand("{sample}.transferred.splices", sample=SAMPLES), expand("{sample}.transferred.qual_check", sample=SAMPLES)


# SRA -> PILEUP -> RAW COUNTS OFF NCBI GFF3 


rule transfer_logs_s3: 
	output: touch("{sample}.transferred.log")
	input:  "{sample}.hisat.two.log"
	message: "transferring {input}'s logs to S3"
	shell: "s3cmd put {wildcards.sample}.hisat.one.log s3://ncbi-hackathon-aug/rnamapping/ ; s3cmd put {wildcards.sample}.hisat.two.log s3://ncbi-hackathon-aug/rnamapping/"

rule transfer_bam_s3: 
	output: touch("{sample}.transferred")
	input: "{sample}.GRCh38.p4.hisat.sorted.bam"
	message: "transferring {input} to S3"
	shell: "s3cmd put {input} s3://ncbi-hackathon-aug/rnamapping/"

rule transfer_splices_s3:
	output: touch("{sample}.transferred.splices")
	input: "{sample}.hisat.novel.splicesites.txt"
	message: "transferring hisat splice sites {input} to s3"
	shell: "s3cmd put {input} s3://ncbi-hackathon-aug/rnamapping/"	

rule transfer_qual_check: 
	output: touch("{sample}.transferred.qual_check")
	input: "{sample}.qc_check.done"
	message: "transferring qual_check output to S3"
	shell: "s3cmd put {wildcards.sample}.GRCh38.p4.hisat.crsm s3://ncbi-hackathon-aug/rnamapping/; s3cmd put {wildcards.sample}.pass s3://ncbi-hackathon-aug/rnamapping/; s3cmd put {wildcards.sample}.fail s3://ncbi-hackathon-aug/rnamapping/"

rule transfer_counts: 
	output: touch("{sample}.transferred.counts")
	input: "{sample}.GRCh38.p4.HTSeq.counts"
	message: "transferring {input} counts to s3 "
	shell: "s3cmd put {input} s3://ncbi-hackathon-aug/rnamapping/"

rule perform_counting: 
	output: "{sample}.GRCh38.p4.HTSeq.counts"
	input: "{sample}.GRCh38.p4.hisat.sorted.bam.bai"
	message: "performing counting of reads on genes in {input}"
	shell: "~/HTSeq-0.6.1/build/scripts-2.7/htseq-count -m intersection-nonempty -i gene -s no -f bam {wildcards.sample}.GRCh38.p4.hisat.sorted.bam ~/refs/GCF_000001405.30_GRCh38.p4_genomic.gff"


rule qc_check: 
	output: touch("{sample}.qc_check.done")
	input: "{sample}.GRCh38.p4.hisat.crsm"
	message: "checking quality stats of {input} with perl script"
	shell: "perl qc.pl --maplogfile {wildcards.sample}.hisat.two.log --metricsfile {wildcards.sample}.GRCh38.p4.hisat.crsm --sra {wildcards.sample}"

rule picard_rnaseq_qual: 
	output: "{sample}.GRCh38.p4.hisat.crsm"
	input: "{sample}.GRCh38.p4.hisat.sorted.bam.bai"
	message: "running picard rna qc stats on {input}"
	shell: "java -jar ~/picard-tools-1.137/picard.jar CollectRnaSeqMetrics REF_FLAT=~/refs/refflat/ncbirefflat.txt STRAND=NONE INPUT={wildcards.sample}.GRCh38.p4.hisat.sorted.bam OUTPUT={output}"

rule index_bam: 
	output: "{sample}.GRCh38.p4.hisat.sorted.bam.bai"
	input: "{sample}.GRCh38.p4.hisat.sorted.bam"
	message: "indexing bam alignment file {input}"
	shell: "samtools index {input} {output}"

rule sort_bam:
	output: "{sample}.GRCh38.p4.hisat.sorted.bam"
	input: "{sample}.GRCh38.p4.hisat.bam"
	message: "sorting {input} to {output}"
	shell: "samtools sort {input} {wildcards.sample}.GRCh38.p4.hisat.sorted"

rule sam_to_bam:
	output: "{sample}.GRCh38.p4.hisat.bam"
	input: "{sample}.GRCh38.p4.hisat.sam"
	message: "converting sam to bam: {input} to {output}"
	shell: "samtools view -bS {input} > {output}"


rule hisat_alignment_two:
	output: "{sample}.GRCh38.p4.hisat.sam", "{sample}.hisat.two.log"
	input: "{sample}.hisat.novel.splicesites.txt"
	threads: 10
	message: "running second pass hisat alignment with {threads} threads"
	shell: "hisat -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 -x {HISATREF} -p {threads} --sra-acc {wildcards.sample} -t -S {wildcards.sample}.GRCh38.p4.hisat.sam --novel-splicesite-infile {wildcards.sample}.hisat.novel.splicesite.txt 2> {wildcards.sample}.hisat.two.log"

rule hisat_alignment_one: 
	output: "{sample}.hisat.novel.splicesites.txt", "{sample}.hisat.one.log", "{sample}.GRCh38.p4.hisat.one.sam"
	threads: 10 
	message: "hisat aligning reads from {wildcards.sample} to GRCh38.p4 with {threads} threads to produce splicesites"
	shell: "hisat -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 -x {HISATREF} -p {threads} --sra-acc {wildcards.sample} -t --novel-splicesite-outfile {wildcards.sample}.hisat.novel.splicesites.txt -S {wildcards.sample}.GRCh38.p4.hisat.one.sam  2> {wildcards.sample}.hisat.one.log"




