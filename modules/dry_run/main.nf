params.outdir = 'results'

process TEST_GENOME {
	tag "Testing the pipeline on a subset of the real data"
	publishDir params.outdir, mode:'copy'

	input:
	path genome 
	path annotations 
	val testchr

	output:
	path "test/$testchr*"

	shell:
	'''
	zcat < !{annotations} | awk '$1 ~ /!{testchr}/' > !{testchr}.gtf
	zcat < !{genome} > genome.fa
	cdbfasta genome.fa 
	cdbyank -a !{testchr} genome.fa.cidx > !{testchr}.fa
	gzip !{testchr}.gtf
	gzip !{testchr}.fa
	mkdir test
	mv !{testchr}.fa.gz !{testchr}.gtf.gz test
	'''
}

process TEST_FASTQS {
	tag "Subsampling one pair of R1 and R2 fastqs for testing the pipeline"
	publishDir params.outdir, mode:'copy'

	input:
	tuple val(sample_id), path(reads)
	val subsamplefastq

	output:
	path "test/test.R{1,2}.fq.gz"

	shell:
	'''
	R1T=!{reads[0]}
	R2T=!{reads[1]}
	seqtk sample -s100 $R1T !{subsamplefastq} > test.R1.fq 
	seqtk sample -s100 $R2T !{subsamplefastq} > test.R2.fq
	gzip test.R1.fq
	gzip test.R2.fq
	mkdir test
	mv test.R1.fq.gz test.R2.fq.gz test
	'''
}
