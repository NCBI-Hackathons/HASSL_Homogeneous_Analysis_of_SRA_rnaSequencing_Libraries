
rule all:
	input: "ncbirefflat.txt", "GCF_000001405.30_GRCh38.p4_genomic.fna"


rule gunzipHumRef:
	output: "GCF_000001405.30_GRCh38.p4_genomic.fna" 
	input: "GCF_000001405.30_GRCh38.p4_genomic.fna.gz"
	message: "extracting human genome fasta.gz"
	shell: "gunzip {input}"
	
rule get_humanreference:
	output: "GCF_000001405.30_GRCh38.p4_genomic.fna.gz"
	message: "downloading human reference genome"
	shell: "wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000001405.30_GRCh38.p4/GCF_000001405.30_GRCh38.p4_genomic.fna.gz"

rule convert:
	output: "ncbirefflat.txt"
	input: "refFlat.txt"
	message: "running refflat converter"
	shell: "refflat2ncbi_refflat.pl {input} > {output} && rm {input}"

rule gunzip:
	output: "refFlat.txt" 
	input: "refFlat.txt.gz"
	message: "extracting refFlat.txt.gz"
	shell: "gunzip {input}"
	
rule get_refflat:
	output: "refFlat.txt.gz"
	message: "downloading refflat from ucsc"
	shell: "wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/{output}"

