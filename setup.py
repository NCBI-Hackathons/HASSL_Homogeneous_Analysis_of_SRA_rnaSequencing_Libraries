
rule all:
#	input: "lib/ncbirefflat.txt", "lib/GCF_000001405.30_GRCh38.p4_genomic.fna" 
	input: "lib/ncbirefflat.txt",
		"lib/GCF_000001405.30_GRCh38.p4_genomic.gff",
		"lib/GRCh38.p4.1.bt2", "lib/GRCh38.p4.3.bt2", "lib/GRCh38.p4.5.bt2", "lib/GRCh38.p4.rev.1.bt2", "lib/GRCh38.p4.rev.5.bt2",
	        "lib/GRCh38.p4.2.bt2", "lib/GRCh38.p4.4.bt2", "lib/GRCh38.p4.6.bt2", "lib/GRCh38.p4.rev.2.bt2", "lib/GRCh38.p4.rev.6.bt2"


# Reference genome (h38)
rule index:
	output: "lib/GRCh38.p4.1.bt2", "lib/GRCh38.p4.3.bt2", "lib/GRCh38.p4.5.bt2", "lib/GRCh38.p4.rev.1.bt2", "lib/GRCh38.p4.rev.5.bt2",
	        "lib/GRCh38.p4.2.bt2", "lib/GRCh38.p4.4.bt2", "lib/GRCh38.p4.6.bt2", "lib/GRCh38.p4.rev.2.bt2", "lib/GRCh38.p4.rev.6.bt2"
	input: "GCF_000001405.30_GRCh38.p4_genomic.fna"
	message: "indexing human genome"
	shell: "cd lib && hisat-build {input} GRCh38.p4 && cd .."

rule gunzipHumRef:
	output: "lib/GCF_000001405.30_GRCh38.p4_genomic.fna" 
	input: "lib/GCF_000001405.30_GRCh38.p4_genomic.fna.gz"
	message: "extracting human genome fasta.gz"
	shell: "gunzip -c {input} > {output}"
	
rule get_humanreference:
	output: temp("lib/GCF_000001405.30_GRCh38.p4_genomic.fna.gz")
	message: "downloading human reference genome"
	shell: "wget -P lib ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000001405.30_GRCh38.p4/GCF_000001405.30_GRCh38.p4_genomic.fna.gz"



# GFF
rule gunzipgff:
	output: "lib/GCF_000001405.30_GRCh38.p4_genomic.gff"
	input: "lib/GCF_000001405.30_GRCh38.p4_genomic.gff.gz"
	message: "extracting gff.gz"
	shell: "gunzip -c {input} > {output}"

rule get_gff:
	output: temp("lib/GCF_000001405.30_GRCh38.p4_genomic.gff.gz")
	message: "downloading gff from ncbi"
	shell: "wget -P lib ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000001405.30_GRCh38.p4/GCF_000001405.30_GRCh38.p4_genomic.gff.gz"


# REF FLAT
rule convert:
	output: "lib/ncbirefflat.txt"
	input: "lib/refFlat.txt"
	message: "running refflat converter"
	shell: "refflat2ncbi_refflat.pl {input} > {output}"

rule gunzip:
	output: temp("lib/refFlat.txt")
	input: "lib/refFlat.txt.gz"
	message: "extracting refFlat.txt.gz"
	shell: "gunzip -c {input} > {output}"
	
rule get_refflat:
	output: temp("lib/refFlat.txt.gz")
	message: "downloading refflat from ucsc"
	shell: "wget -P lib http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz"

