params.outdir = 'results'

process MULTIQC {
	label 'secondary'
	publishDir params.outdir, mode:'copy'

	input:
	path('*')

	output:
	path "multiqc_report.html"

	script:
	"""
	multiqc .
	"""
}