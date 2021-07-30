# **PLASTER**: Phased Long Allele Sequence Typing with Error Removal

PLASTER is a comprehensive data processing pipeline for allele typing from long amplicon sequencing data generated on the PacBio SMRT platform. Inputs are PacBio subreads in BAM format, as well as sample barcodes and target amplicon details. Outputs are phased BAMs for each sample amplicon and variant calls in VCF format. Additionally the pipeline supports Pharmacogenomic star alelle assignment using the PharmVar database, and gene fusion detection for CYP2D6 and CYP2D7 fusion alleles.

The pipeline is built using [Nextflow](https://nextflow.io/), a workflow tool to run tasks across multiple compute infrastructures in a  portable and efficient manner. Included is a [Docker](https://www.docker.com/) container, making installation trivial and results highly reproducible. 

## Pipeline Overview
<p align="center"><img src="doc/diagram.png"/></p>

## Prerequisites

* Make sure [Nextflow](https://nextflow.io/) and either [Singularity](https://sylabs.io/guides/3.0/user-guide/index.html) or [Docker](https://www.docker.com/) are installed on your system

## Usage

### Pre-processing

* **Running the test dataset**
  ```
  nextflow run bahlolab/PLASTER -profile preproc,test,singularity
  ```
  This command will download the pipeline from GitHub and run the pre-processing stage on a minimal test dataset using singularity to run the software container. Replace "singularity" with "docker" to use docker instead.
* **Running your own dataset**
  ```
  nextflow run bahlolab/PLASTER -profile preproc,singularity -c <my_dataset.config>
  ```
  where `my_dataset.config` is a Nextflow config file specifying the following required parameters:
  ```Nextflow
  params {
    subreads_bam = '/PATH/TO/MY/subreads.bam'
    barcodes_fasta = '/PATH/TO/MY/barcodes.fasta'
    amplicons_json = '/PATH/TO/MY/amplicons.json'
    ref_fasta = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz'
  }
  ```
  See [Pre-processing Parameters](doc/preproc.md) from more details

### Allele-typing

* **Running the test dataset**
  ```
  nextflow run bahlolab/PLASTER -profile typing,test,singularity
  ```
  This command will download the pipeline from GitHub and run the pre-processing stage on a minimal test dataset using singularity to run the software container. Replace "singularity" with "docker" to use docker instead.
* **Running your own dataset**
  ```
  nextflow run bahlolab/PLASTER -profile typing,singularity -c <my_dataset.config>
  ```
  where `my_dataset.config` is a file specifying the following required parameters:
  ```Nextflow
  params {
    manifest = '@PROJECT_DIR@/test/typing/manifest.csv'
    amplicons_json = '@PROJECT_DIR@/test/typing/amplicons.json'
    ref_fasta = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz'
  }
  ```
  See [Allele-typing Parameters](doc/typing.md) from more details


