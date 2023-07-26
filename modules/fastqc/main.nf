params.outdir = 'results'

process FASTQC {
	tag "FASTQC on $sample_id"
	publishDir params.outdir, mode:'copy'

	input:
	tuple val(sample_id), path(reads)

	output:
	path "fastqc_${sample_id}_logs"

	script:
	"""
	mkdir fastqc_${sample_id}_logs
	fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads} -t ${task.cpus} 
	"""
}

process FASTQC_TEST {
	tag "FASTQC test run"
	label 'secondary'
	publishDir params.outdir, mode:'copy'

	input:
	path(reads)

	output:
	path "fastqc_test_logs"

	script:
	"""
	mkdir fastqc_test_logs
	fastqc -o fastqc_test_logs -f fastq -q ${reads} -t ${task.cpus} 
	"""
}