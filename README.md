# scRNAseq-NF

This is an experimental Nextflow pipeline currently under development for analyzing scRNA-seq droplet-based data from 10x Genomics, although support for other types of scRNA-seq experiments will be added in the future. 

[![homepage](https://img.shields.io/badge/nextflow-%E2%89%A522.04.5-brightgreen.svg)](https://nextflow.io/ "Redirect to nextflow homepage")

## Requirements

- Unix-like operating system (Linux, macOS, etc)
- Java 8 or later
- Docker 20.10.17 or later
- Nextflow 22.04.5 or later

## Components

The following software tools are used for this pipeline:

- [STAR 2.7.10a_alpha_220818](https://github.com/alexdobin/STAR/releases/tag/2.7.10a_alpha_220818) 
- [wget](https://www.gnu.org/software/wget/) 1.20.3
- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 0.11.9
- [multiqc](https://multiqc.info/) 1.13
- [samtools](http://www.htslib.org/) 1.15.1
- [bbmap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/) 38.97
- [rsem](https://github.com/deweylab/RSEM) 1.3.3
- [seqtk](https://github.com/lh3/seqtk) 1.3
- [cdbtools](https://github.com/gpertea/cdbfasta) 0.99

## Quickstart

If you don't already have Nextflow, install by using the following command:

```
curl -s https://get.nextflow.io | bash
```

The Dockerfile and conda.yaml are included if updated or additional software tools are desired to be added to the Docker image used throughout this pipeline. To download the Docker image from dockerhub, use this command:

```
docker pull rlinder02/starsolo-sc-rnaseq:v0.1.0
```

Download the pipeline from GitHub (do not try running directly from GitHub):

```
git clone https://github.com/rlinder02/scRNA-seq-NF.git
```

Run the pipeline:

```
nextflow run main.nf
```
The nextflow.config file ensures that the Docker image is used by default for containerised execution, unless overriden at the command line.


## Pipeline Description

This pipeline uses STARsolo to create count matrices starting from raw fastq files generated by a scRNA-seq experiment. At this moment, only 10X Genomics droplet-based data can be used with the pipeline. In the future, support for additional scRNA-seq experimental setups will be added, as will support for using Scanpy to run tertiary analyses of the scRNA-seq data. Currently, this pipeline is configured to be run either locally or on AWS. Input/output files can be stored either locally or in an S3 bucket. Instructions for configuring an AWS account to work with Nextflow can be found [here](https://staphb.org/resources/2020-04-29-nextflow_batch.html) and more generally [here](https://seqera.io/blog/nextflow-and-aws-batch-inside-the-integration-part-2-of-3/).

## Input files

This pipeline requires the following gzipped input files:

- 10x RNAseq paired end reads, `*.fastq.gz`
- Genome fasta, `*.fa.gz`
- Whitelist of 10x barcodes (should not be gzipped)
- Genome annotation file, `*.gtf.gz`

The RNAseq read file names should match the following pattern: _*_R{1,2}_001.{fastq, fq}.gz_

Ensure that chromosomes are labelled in the same manner between the fasta and gtf files (ie., both use _chr1, chr2, etc..._ or _chrom1, chrom2, etc..._; the format should be the same).

## Pipeline parameters 

At the command line, you can specify several options (documented [here](https://www.nextflow.io/docs/latest/)). The most relevant options for this pipeline include:

`-- runtype`

- This specifies whether or not you are running a test of the pipeline (accepts _test_ or _full_, by default this is _full_). A test takes only the first pair of reads, extracts a user-specified subsample of reads (default is 200000, can be changed by using the `--subsamplefastq` parameter), takes and extracts a single user-specified chromosome from the fasta file (default is _chr10_, this can be changed by using the `--testchr` parameter), and runs the pipeline using this reduced dataset (with the full set of annotations). 

`--reads`

- This specifies the fastq reads to use as input.

`--genome`

- This specifies the fasta file to use.

`--annotations`

- This specifies the genome annotation file to use.

`--whitelist10k`

- This specifies the 10X barcode whitelist to be used.

`--outdir`

- The directory into which all output will be copied.

Currently, default STAR parameters are set up for 10X 3'v3, v3.1 chemistry using the *3M-february-2018.txt* barcode list. Support for automatic detection of the 10X version will be added in the future. However, for a complete and detailed list of parameters for running STARsolo on multiple types of 10X chemistries, please see [here](https://github.com/cellgeni/STARsolo).

The indexing and alignment/counting processes are currently labelled as `big_mem`. The number of CPUs and requested memory can be changed in the nextflow.config file and depends on the resources available. It is important to remember that, when running on AWS, the number of vCPUs provisioned for a job will not be the same as the number of CPUs reqested in the config file, as there are normally 2+ vCPUs per CPU (exact specifications for a variety of EC2 instances are listed in table form [here](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/cpu-options-supported-instances-values.html)).

## Pipeline results

The results are copied to a generic `results` folder, with STARsolo genome indexing results stored in the `index` sub-folder and the STARsolo alignment and count matrix files are stored in `sample_id` sub-folders. The `sample_id` is generated from the part of the fastq file-name that comes before the above-mentioned pattern. 

## Credits

The parameters used to set up the STARsolo run were originally specified [here](https://github.com/cellgeni/STARsolo), as well as code for the end of the STARsolo run. The general format of the pipeline was inspired by the work done [here](https://github.com/nextflow-io/rnaseq-nf) and [here](https://github.com/nf-core/rnaseq).


