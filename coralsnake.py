
snakemake -s coralsnake.py 


REF = 'GRCh38'
alner = 'hisat'

DATASETS = "test test2".split() 

rule index_sam:
	output: "{sample}.{REF}.{alner}.sorted.sam.sai"
	input: "{sample}.{ref}.{alner}.sorted.sam"
	message: "indexing sam {input}"
	shell: {samtools sortsam? }

rule sort_bam:
	output: "{sample}.{REF}.{alner}.sorted.bam"
	input: "{sample}.{REF}.{alner}.bam"
	message: "sorting {intput} to {output}"
	shell: "samtools sort {input} {output}.sorted"

rule sort_sam:
	output: "{sample}.{REF}.{alner}.sorted.sam"
	input: "{sample}.{REF}.{alner}.sam"
	message: "sorting {intput} to {output}"
	shell: "samtools sort {input} {output}"

rule sam_to_bam:
	output: "{sample}.{REF}.{alner}.sorted.bam"
	input: "{sample}.{REF}.{alner}.sam"
	message: "converting sam to bam: {input} to {output}"
	shell: "samtools view -bS {input} > {output}"

rule hisat_alignment:
	output: "{sample}.{REF}.hisat.sam"
	input: "{sample}.fastq"
	message: "hisat aligning reads {sample}.fastq to {REF} to produce {output}"
	shell: <hisat code>

rule fastq_stats: 
	output: dir? 
	input: "{sample}.fastq"
	message: "generating fastq stats on {input}"
	shell: "fastq-stats.py"
	OR
	run: 
		<high perf stats code .python> and potentially QC exclude samples





#### POTENTIAL 


rule annotate:
    output:
        'annotations/{ds}.txt'
    input:
        track=HGTRACK, bam='maphg/{ds}.bam'
    shell:
        # in version 2.16.1 of bedtools this does not result in the correct output (maybe because -wo is ignored)
        'bedtools intersect -wo -bed -abam {input.bam} -b {input.track} | cut -f1-6,8 > {output}'

rule hgtrack:
    output:
        'tracks/Homo_sapiens.GRCh37.67.gtf'
    shell:
        'wget -O {output}.gz ftp://ftp.ensembl.org/pub/release-67/gtf/homo_sapiens/Homo_sapiens.GRCh37.67.gtf.gz; gunzip {output}.gz'


if len(DATASETS) == 0:
	print("No appropriate files were found in the reads/ and readspe/ directories.")
	print("Please copy or symlink .fastq.gz files into the reads/ directory or")
	print("pairs of .1.fastq.gz/.2.fastq.gz files into the readspe/ directory in")
	print("order for this Snakefile to do something.")
	sys.exit(1)



def all_samples_with_data(session=SESSION):
	'''
	Return a set of all samples from the database for which at least
	one FASTQ or bam file exists
	'''
	# TODO take family of patient into account, separate into proper samples after calling.
	units = session.query(Unit).all()
	result = set()
	for unit in get_available_units(units):
		result.add(unit.library.sample.accession)
	result -= IGNORED_SAMPLES
	return result



def 



