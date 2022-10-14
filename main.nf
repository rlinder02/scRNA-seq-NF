#!/usr/bin/env nextflow

/*
 * Proof of concept of a scRNAseq secondary analysis pipeline implemented with Nextflow
 * 
 * Author:
 * - Robert Linder <rlinder02@gmail.com>
 */

nextflow.enable.dsl = 2


/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.reads = "$projectDir/assets/*_R{1,2}_001.{fastq,fq}.gz"
params.genome = "$projectDir/assets/GRCh38.primary_assembly.genome.fa.gz"
params.annotations = "$projectDir/assets/gencode.v32.primary_assembly.annotation.gtf.gz"
params.whitelist10k = "$projectDir/assets/3M-february-2018.txt.gz"
params.outdir = "$projectDir/results"
//change params.runtype to "test" at the command line when launching this nf pipeline to run the full pipeline on a single user-defined chromosome and subsample of reads
params.runtype = "full"
params.testchr = "chr10"
params.subsamplefastq = 200000


log.info """\
 scR N A S E Q - N F   P I P E L I N E
 ===================================
 genome: ${params.genome}
 annotations  : ${params.annotations} 
 reads        : ${params.reads}
 outdir       : ${params.outdir}
 """

// import modules; can only invoke each module once unless it is given multiple names
include {FASTQC; FASTQC_TEST} from './modules/fastqc/main'
include {MULTIQC} from './modules/multiqc/main'
include {TEST_GENOME; TEST_FASTQS} from './modules/dry_run/main'
include {INDEX; INDEX_TEST} from './modules/index/main'
include {STARSOLO; STARSOLO_TEST} from './modules/starsolo/main'

/* 
* main script flow
*/

workflow {
	read_pairs_ch = channel.fromFilePairs(params.reads, checkIfExists: true)
	if (params.runtype == "test") {
		TEST_GENOME(params.genome, params.annotations, params.testchr)
		TEST_FASTQS(read_pairs_ch.first(), params.subsamplefastq)
		FASTQC_TEST(TEST_FASTQS.out)
        INDEX_TEST(TEST_GENOME.out)
		STARSOLO_TEST(TEST_FASTQS.out, INDEX_TEST.out, params.whitelist10k)
        MULTIQC(FASTQC_TEST.out)
    } else {
        FASTQC(read_pairs_ch)
        MULTIQC(FASTQC.out)
        INDEX(params.genome, params.annotations) 
        STARSOLO(read_pairs_ch, INDEX.out.collect(), params.whitelist10k)
    }
}


/* 
 * completion handler
 */
workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}