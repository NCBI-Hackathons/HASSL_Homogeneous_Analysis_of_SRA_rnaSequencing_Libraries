


STAR --chimericAlignment 




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

