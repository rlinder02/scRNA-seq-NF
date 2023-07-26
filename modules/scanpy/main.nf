params.outdir = 'results'

process QC_CELLS_TEST {
    echo = true
	tag "QC count matrix"
	label 'big_mem'
    label 'tertiary' 
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	path starsolo

	output:
	path starsolo

	script:

	"""
    qc_count_matrix.py "${starsolo}/Gene/filtered"
    cd output/
    mkdir analyses
    cd ..
    mv data* output/analyses/
    mv figures output/analyses/
    """
	
}