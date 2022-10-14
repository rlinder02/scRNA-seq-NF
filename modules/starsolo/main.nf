params.outdir = 'results'

process STARSOLO {
	tag "STARsolo run"
	label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	tuple val(sample_id), path(reads)
	path index 
	path wlist

	output:
	path "${sample_id}"

	shell:

	'''
	zcat < !{wlist} > barcodes.txt
	STAR \\
	 --genomeDir !{index} \\
	 --readFilesIn !{reads[1]} !{reads[0]} \\
	 --readFilesCommand zcat \\
	 --runDirPerm All_RWX \\
	 --outSAMtype BAM SortedByCoordinate \\
	 --genomeLoad NoSharedMemory \\
	 --soloStrand Forward \\
	 --soloType CB_UMI_Simple \\
	 --soloCBwhitelist barcodes.txt \\
	 --soloBarcodeReadLength 0 \\
	 --soloCBlen 16 \\
	 --soloUMIstart 17 \\
	 --soloUMIlen 12 \\
	 --soloUMIdedup 1MM_CR \\
	 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\
	 --soloUMIfiltering MultiGeneUMI_CR \\
	 --soloCellFilter EmptyDrops_CR \\
	 --clipAdapterType CellRanger4 \\
	 --outFilterScoreMin 30 \\
	 --soloFeatures Gene GeneFull Velocyto \\
	 --soloOutFileNames !{sample_id}/ features.tsv barcodes.tsv matrix.mtx \\
	 --soloMultiMappers EM \\
	 --limitOutSJcollapsed 5000000 \\
	 --outReadsUnmapped Fastx

	if [[ -s Aligned.sortedByCoord.out.bam ]]
	then
		samtools index -@16 Aligned.sortedByCoord.out.bam
	fi
	
	gzip Unmapped.out.mate1 &
	gzip Unmapped.out.mate2 &

	mv *.bam !{sample_id}/
	mv *.bai !{sample_id}/
	mv *.gz !{sample_id}/

	cd !{sample_id}
	for i in Gene/raw Gene/filtered GeneFull/raw GeneFull/filtered Velocyto/raw Velocyto/filtered
	do 
	  cd $i; for j in *; do gzip $j & done
	  cd ../../
	done

	wait
	echo "ALL DONE!"
	'''
}

process STARSOLO_TEST {
	tag "STARsolo on test data"
	label 'big_mem'
	publishDir params.outdir, mode:'copy'

	input:
	path reads
	path index
	path wlist 

	output:
	path "output"

	shell:
	println reads[0]
	println reads[1]

	'''
	mkdir -p test/alignment
	zcat < !{wlist} > barcodes.txt
	STAR \\
	 --genomeDir !{index} \\
	 --readFilesIn !{reads[1]} !{reads[0]} \\
	 --readFilesCommand zcat \\
	 --runDirPerm All_RWX \\
	 --outSAMtype BAM SortedByCoordinate \\
	 --genomeLoad NoSharedMemory \\
	 --soloStrand Forward \\
	 --soloType CB_UMI_Simple \\
	 --soloCBwhitelist barcodes.txt \\
	 --soloBarcodeReadLength 0 \\
	 --soloCBlen 16 \\
	 --soloUMIstart 17 \\
	 --soloUMIlen 12 \\
	 --soloUMIdedup 1MM_CR \\
	 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\
	 --soloUMIfiltering MultiGeneUMI_CR \\
	 --soloCellFilter EmptyDrops_CR \\
	 --clipAdapterType CellRanger4 \\
	 --outFilterScoreMin 30 \\
	 --soloFeatures Gene GeneFull Velocyto \\
	 --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx \\
	 --soloMultiMappers EM \\
	 --outReadsUnmapped Fastx

	if [[ -s Aligned.sortedByCoord.out.bam ]]
	then
		samtools index -@16 Aligned.sortedByCoord.out.bam
	fi
	
	gzip Unmapped.out.mate1 &
	gzip Unmapped.out.mate2 &

	mv *.bam output/
	mv *.bai output/
	mv *.gz output/

	cd output
	for i in Gene/raw Gene/filtered GeneFull/raw GeneFull/filtered Velocyto/raw Velocyto/filtered
	do 
	  cd $i; for j in *; do gzip $j & done
	  cd ../../
	done

	wait
	echo "ALL DONE!"
	'''
}