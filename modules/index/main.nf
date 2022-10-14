params.outdir = 'results'

process INDEX {
	tag "Indexing genome"
	label 'big_mem'
	publishDir params.outdir, mode:'copy'

	input:
	path genome
	path annotations
	

	output:
	path "index"

	shell:
	'''
	mkdir index
	zcat < !{genome} > genome.fa
	zcat < !{annotations} > annotations.gtf
	mv genome.fa annotations.gtf index
	echo "Done moving!"
	cd index
	STAR \\
		--runThreadN !{task.cpus} \\
		--runMode genomeGenerate \\
		--genomeDir ./ \\
		--genomeFastaFiles genome.fa \\
		--sjdbGTFfile annotations.gtf \\
		--genomeSAindexNbases 14
	'''
}

process INDEX_TEST {
	tag "Indexing test sample"
	label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	// can't cd in shell script when using AWS; must stay in current directory or else cannot find file
	path genomeanns

	output:
	path "test/index"

	shell:
	'''
	mkdir -p test/index
	zcat < !{genomeanns[0]} > genome.fa
	zcat < !{genomeanns[1]} > annotations.gtf
	mv genome.fa annotations.gtf test/index
	cd test/index
	STAR \\
		--runThreadN 12 \\
		--runMode genomeGenerate \\
		--genomeDir ./ \\
		--genomeFastaFiles genome.fa \\
		--sjdbGTFfile annotations.gtf \\
		--genomeSAindexNbases 12
	'''
}