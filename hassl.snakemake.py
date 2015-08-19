# THIS SCRIPT WAS PRODUCED VIA THE NCBI HACKATHON IN AUGUST 2015 
# WRITTEN (INITIALLY) BY PAUL CANTALUPO, HIROKO OHMIYA, ALLISSA DILLMAN, AND RUSSELL DURRETT


# TO DO - UPDATE TO SUBREAD PACKAGE (for counting) INSTEAD OF HTSEQ (TOO SLOW, SUBREAD SUPER FAST)



#OUTPUT LOCATIONS
#working directory for now 
S3_BUCKET='s3://genomicdata/hassl/'


# FASTAREF='/resources/ensembl/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa'
HISATREF="/home/ubuntu/resources/ensembl/hisat_indexes/Homo_sapiens.GRCh38.dna.toplevel"
# HISATREF="/resources/ensembl/hisat_indexes/Homo_sapiens.GRCh38.dna.toplevel"
#GFFFILE = "/home/ubuntu/refs/GCF_000001405.30_GRCh38.ens77_genomic.gff"
PICARDFLATFILE="/home/ubuntu/resources/ensembl/GRCh38.77.compatible.ucsc.picard.refflat.txt"
GTFFILE="/home/ubuntu/resources/ensembl/Ensembl.GRCh38.77.gtf"
SPLICEFILE="/home/ubuntu/resources/ensembl/Ensembl.GRCh38.77.splicesites.txt"
#set the number of threads to use in alignments 
THREADS=5

#set the filename of the file with the list of accessions   

filename = config["ACCESSION_FILE"]



# EXECUTABLE LOCATIONS (some on path)
HISAT=" /home/ubuntu/install/hisat/hisat "
PICARD=" java -jar /home/ubuntu/install/picard-tools-1.138/picard.jar "
FEATURECOUNTS="/home/ubuntu/install/subread-1.4.6-p4-Linux-x86_64/bin/featureCounts"
#HTSEQ=" ~/HTSeq-0.6.1/build/scripts-2.7/htseq-count "

SAMTOOLS=" /home/ubuntu/install/samtools_rocksdb/samtools/samtools "

SAMPLES_FROM_FILE = [line.rstrip('\n') for line in open(filename)]
SAMPLES = [s for s in SAMPLES_FROM_FILE if s]

rule all: 
  input: expand("{sample}.GRCh38.ens77.featureCounts.counts", sample=SAMPLES)  #, expand("{sample}.qc_check.done", sample=SAMPLES)

rule all_on_s3: 
  input: expand("{S3_BUCKET}/{sample}.GRCh38.ens77.featureCounts.counts", sample=SAMPLES, S3_BUCKET=S3_BUCKET)  #, expand("{sample}.qc_check.done", sample=SAMPLES)

# input: expand("{sample}.GRCh38.ens77.HTSeq.counts", sample=SAMPLES)


rule clean: 
  shell: "rm -fr SRR*; rm log/SRR* "



rule transfer_to_s3: 
  output: "{S3_BUCKET}/{file}"
  input: "{file}"
  log: "log/{file}.xfer.to.s3.log"
  message: "transferring {file} to S3 Bucket {S3_BUCKET}  ->  {S3_BUCKET}/{file}"
  shell: "s3cmd put {wildcards.file} {S3_BUCKET}/{file}"

# rule perform_counting: 
#   output: "{sample}.GRCh38.ens77.HTSeq.counts"
#   input: "{sample}.GRCh38.ens77.hisat.sorted.bam.bai"
#   log: "log/{wildcards.sample}.counting.log"
#   message: "performing counting of reads on genes in {input}"
#   shell: "time {HTSEQ} -m intersection-nonempty -i gene -s no -f bam {wildcards.sample}.GRCh38.ens77.hisat.sorted.bam {GFFFILE} > {wildcards.sample}.GRCh38.ens77.HTSeq.counts 2> {wildcards.sample}.HTseq-count.log"

rule project_counting:
  output: "project.featureCounts"
  input: GTFFILE 
  log: "log/project.featureCounting"
  threads: THREADS 
  message: "performing overall project feature counting"
  shell: "{FEATURECOUNTS} -T {threads} --primary -F GTF -t exon -g gene_id -a {GTFFILE} -o project.featureCounts *.bam 2>&1 > log/project.featureCounting"

rule perform_counting: 
  output: "{sample}.GRCh38.ens77.featureCounts.counts"
  input: "{sample}.GRCh38.ens77.hisat.sorted.bam"
  log: "log/{wildcards.sample}.featureCounts.log"
  threads: THREADS
  message: "performing featureCounting with {threads} threads on genes in {input}"
  shell: " {FEATURECOUNTS}  -T {threads} --primary -F GTF -t exon -g gene_id -a {GTFFILE} -o {wildcards.sample}.GRCh38.ens77.featureCounts.counts {wildcards.sample}.GRCh38.ens77.hisat.sorted.bam  2> {wildcards.sample}.featureCounts.log"

#IF COUNTING THEN JUST REPORT ONE MAX HIT PER READ ?  --primary fixes that? 
# PAIRED END?...  -p  and  -P  

rule qc_check: 
  output: touch("{sample}.qc_check.done")
  input: "{sample}.GRCh38.ens77.hisat.crsm"
  log: "log/{wildcards.sample}.perlqc.log"
  message: "checking quality stats of {input} with perl script"
  shell: " perl qc.pl --maplogfile {wildcards.sample}.hisat.two.log --metricsfile {wildcards.sample}.GRCh38.ens77.hisat.crsm --sra {wildcards.sample} 2> log/{wildcards.sample}.perl_qc.log"

rule picard_rnaseq_qual: 
  output: "{sample}.GRCh38.ens77.hisat.crsm"
  input: "{sample}.GRCh38.ens77.hisat.sorted.bam.bai"
  log: "log/{wildcards.sample}.picard_rnametrics.log"
  message: "running picard rna qc stats on {input}"
  shell: " {PICARD} CollectRnaSeqMetrics REF_FLAT={PICARDFLATFILE} STRAND=NONE INPUT={wildcards.sample}.GRCh38.ens77.hisat.sorted.bam OUTPUT={output} 2> log/{wildcards.sample}.picard_rnametrics.log"

rule index_bam: 
  output: "{sample}.GRCh38.ens77.hisat.sorted.bam.bai"
  input: "{sample}.GRCh38.ens77.hisat.sorted.bam"
  message: "indexing bam alignment file {input}"
  shell: " {SAMTOOLS} index {input} {output} "

rule sort_bam:
  output: "{sample}.GRCh38.ens77.hisat.sorted.bam"
  input: "{sample}.GRCh38.ens77.hisat.bam"
  threads: THREADS
  message: "sorting {input} to {output}"
  shell: " {SAMTOOLS} sort -@ {threads} {input} {wildcards.sample}.GRCh38.ens77.hisat.sorted "

rule sam_to_bam:
  output: temp("{sample}.GRCh38.ens77.hisat.bam")
  input: "{sample}.GRCh38.ens77.hisat.sam"
  message: "converting sam to bam: {input} to {output}"
  shell: " {SAMTOOLS} view -bS {input} > {output} "

rule hisat_alignment:
  output: temp("{sample}.GRCh38.ens77.hisat.sam")
  input: SPLICEFILE 
  threads: THREADS
  log: "log/{sample}.hisat.alignment.log"
  message: "running second pass hisat alignment on {wildcards.sample} with {threads} threads"
  shell: " {HISAT} -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 -x {HISATREF} -p {threads} --sra-acc {wildcards.sample} -t --known-splicesite-infile {SPLICEFILE} -S {wildcards.sample}.GRCh38.ens77.hisat.sam  2> {log}"



