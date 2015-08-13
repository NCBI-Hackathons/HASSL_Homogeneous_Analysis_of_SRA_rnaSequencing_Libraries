# THIS SCRIPT WAS PRODUCED VIA THE NCBI HACKATHON IN AUGUST 2015 
# WRITTEN (INITIALLY) BY PAUL CANTALUPO, HIROKO OHMIYA, ALLISSA DILLMAN, AND RUSSELL DURRETT


# TO DO - UPDATE TO SUBREAD PACKAGE INSTEAD OF HTSEQ (TOO SLOW)



# FASTAREF='/home/ubuntu/russ/ncbi/GCF_000001405.30_GRCh38.ens77_genomic.fna'
HISATREF = "/resources/GRCh38/indexes/hisat-index/GRCh38.p3.genome"
#HISATREF = "/resources/ensembl/hisat_indexes/Homo_sapiens.GRCh38.dna.toplevel"
#GFFFILE = "/home/ubuntu/refs/GCF_000001405.30_GRCh38.ens77_genomic.gff"

#set the number of threads to use in alignments 
THREADS=12

#set the filename of the file with the list of accessions 
filename = "/home/ubuntu/accessions"


# EXECUTABLE LOCATIONS (some on path)
HISAT=" /home/ubuntu/install/hisat/hisat "
PICARD=" java -jar /home/ubuntu/install/picard-tools-1.138/picard.jar "
FEATURECOUNTS="/home/ubuntu/install/subread-1.4.6-p4-Linux-x86_64/bin/featureCounts"
#HTSEQ=" ~/HTSeq-0.6.1/build/scripts-2.7/htseq-count "

SAMTOOLS=" samtools "

SAMPLES = [line.rstrip('\n') for line in open(filename)]

rule all: 
	input: expand("{sample}.GRCh38.ens77.hisat.sorted.bam.bai", sample=SAMPLES)
#	input: expand("{sample}.GRCh38.ens77.HTSeq.counts", sample=SAMPLES)


rule clean: 
	shell: "rm -fr SRR*"

rule perform_counting: 
	output: "{sample}.GRCh38.ens77.HTSeq.counts"
	input: "{sample}.GRCh38.ens77.hisat.sorted.bam.bai"
	message: "performing counting of reads on genes in {input}"
	shell: "time {HTSEQ} -m intersection-nonempty -i gene -s no -f bam {wildcards.sample}.GRCh38.ens77.hisat.sorted.bam {GFFFILE} > {wildcards.sample}.GRCh38.ens77.HTSeq.counts 2> {wildcards.sample}.HTseq-count.log"

rule qc_check: 
	output: touch("{sample}.qc_check.done")
	input: "{sample}.GRCh38.ens77.hisat.crsm"
	message: "checking quality stats of {input} with perl script"
	shell: "time perl qc.pl --maplogfile {wildcards.sample}.hisat.two.log --metricsfile {wildcards.sample}.GRCh38.ens77.hisat.crsm --sra {wildcards.sample} 2> {wildcards.sample}.perl_qc.log"

rule picard_rnaseq_qual: 
	output: "{sample}.GRCh38.ens77.hisat.crsm"
	input: "{sample}.GRCh38.ens77.hisat.sorted.bam.bai"
	message: "running picard rna qc stats on {input}"
	shell: "time {PICARD} CollectRnaSeqMetrics REF_FLAT=~/refs/refflat/ncbirefflat.txt STRAND=NONE INPUT={wildcards.sample}.GRCh38.ens77.hisat.sorted.bam OUTPUT={output} 2> {wildcards.sample}.picard_qual.log"

rule index_bam: 
	output: "{sample}.GRCh38.ens77.hisat.sorted.bam.bai"
	input: "{sample}.GRCh38.ens77.hisat.sorted.bam"
	message: "indexing bam alignment file {input}"
	shell: "time {SAMTOOLS} index {input} {output} 2> {wildcards.sample}.index_bam.log"

rule sort_bam:
	output: "{sample}.GRCh38.ens77.hisat.sorted.bam"
	input: "{sample}.GRCh38.ens77.hisat.bam"
	message: "sorting {input} to {output}"
	shell: "time {SAMTOOLS} sort {input} {wildcards.sample}.GRCh38.ens77.hisat.sorted 2> {wildcards.sample}.sort_bam.log"

rule sam_to_bam:
	output: temp("{sample}.GRCh38.ens77.hisat.bam")
	input: "{sample}.GRCh38.ens77.hisat.sam"
	message: "converting sam to bam: {input} to {output}"
	shell: "time {SAMTOOLS} view -bS {input} > {output} 2> {wildcards.sample}.sam_to_bam.log"

rule hisat_alignment_two:
	output: temp("{sample}.GRCh38.ens77.hisat.sam"), "{sample}.hisat.two.log"
	input: "{sample}.hisat.novel.splicesites.txt"
	threads: 12
	message: "running second pass hisat alignment with {threads} threads"
	shell: "time {HISAT} -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 -x {HISATREF} -p {threads} --sra-acc {wildcards.sample} -t -S {wildcards.sample}.GRCh38.ens77.hisat.sam --novel-splicesite-infile {wildcards.sample}.hisat.novel.splicesites.txt 2> {wildcards.sample}.hisat.two.log"

rule hisat_alignment_one: 
	output: "{sample}.hisat.novel.splicesites.txt", "{sample}.hisat.one.log", temp("{sample}.GRCh38.ens77.hisat.one.sam")
	threads: 12
	message: "hisat aligning reads from {wildcards.sample} to GRCh38.ens77 with {threads} threads to produce splicesites"
	shell: "time {HISAT} -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 -x {HISATREF} -p {threads} --sra-acc {wildcards.sample} -t --novel-splicesite-outfile {wildcards.sample}.hisat.novel.splicesites.txt -S {wildcards.sample}.GRCh38.ens77.hisat.one.sam  2> {wildcards.sample}.hisat.one.log"




