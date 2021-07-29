# **PLASTER**: Phased Long Allele Sequence Typing with Error Removal

PLASTER is a comprehensive data processing pipeline for allele typing from long amplicon sequencing data generated on the PacBio SMRT platform. Inputs are pacbio subreads in BAM format, as well as sample barcodes and target amplicon details. Outputs are phased BAMs for each sample amplicon and variant calls in VCF format. Additionally the pipeline supports Pharmacogenomic star alelle assignment using the PharmVar database, and gene fusion detection for CYP2D6 and CYP2D7 fusion alleles.

The pipeline is built using [Nextflow](https://nextflow.io/), a workflow tool to run tasks across multiple compute infrastructures in a  portable and efficient manner. Included is a [Docker](https://www.docker.com/) container, making installation trivial and results highly reproducible. 

## Pipeline Overview
<p align="center"><img src="images/diagram.png"/></p>

## Prerequisites

* Make sure [Nextflow](https://nextflow.io/) and either [Singularity](https://sylabs.io/guides/3.0/user-guide/index.html) or [Docker](https://www.docker.com/) are installed on your system

## Usage

#### Pre-processing

* **Running the test dataset**
  ```
  nextflow run bahlolab/PLASTER -profile preproc,test,singularity
  ```
  This command will download the pipeline from GitHub and run the pre-processing stage on a minimal test dataset using singularity to run the software container. Replace "singularity" with "docker" to use docker instead.
* **Running your own dataset**
  ```
  nextflow run bahlolab/PLASTER -profile preproc,singularity -c <my_dataset.config>
  ```
  where `my_dataset.config` is a file specifying the following required parameters:
  ```Nextflow
  params {
    run_id = 'my-dataset'
    subreads_bam = '/PATH/TO/MY/subreads.bam'
    sample_manifest = '/PATH/TO/MY/sample-manifest.tsv'
    amplicons_json = '/PATH/TO/MY/amplicons.json'
    barcodes_fasta = '/PATH/TO/MY/barcodes.fasta'
    ref_fasta = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz'
  }
  ```
  and `sample_manifest` is a TSV file with the following format:
  ```
  sample	barcode	amplicons
  NA07439	BC04	CYP2D6;CYP2D7
  NA10005	BC22	CYP2D6;CYP2D7
  NA17203	BC59	CYP2D6;CYP2D7
  ...  ...  ...
  ```
  and `amplicons_json` is a JSON file with the following format:
  ```JSON
  {
    "CYP2D6": {
      "chrom": "chr22",
      "start": 42125398,
      "end": 42131503,
      "strand": "-",
      "fwd_primer": "TGTGAATATTGTCTTTGTGTGGGTG",
      "rvs_primer": "CAGGACTCAGGTAATCATATGCTCA"
    },
    "CYP2D7": {
      "chrom": "chr22",
      "start": 42137550,
      "end": 42145176,
      "strand": "-",
      "fwd_primer": "TGTGAATATTGTCTTTGTGTGGGTG",
      "rvs_primer": "CAGGACTCAGGTAATCATATGCTCA"
    }
  }
  ```

#### Allele-typing

* **Running the test dataset**
  ```
  nextflow run bahlolab/PLASTER -profile typing,test,singularity
  ```
  This command will download the pipeline from GitHub and run the pre-processing stage on a minimal test dataset using singularity to run the software container. Replace "singularity" with "docker" to use docker instead.


### Details

#### Pre-processing

1. **Consensus**  
Circular consensus sequences are produced using PacBio's CCS tool with the recommend minimum 3 passes and minimum 0.99 accuracy.  
1. **Tag and trim barcodes**  
PacBio's lima tool is used along with the set of input barcode sequences to tag reads with barcode identity and to trim the barcode and adapter sequence. This is done for all CCS reads, as well as subreads that failed the CCS process for QC purposes.  
**3. Alignment**  
CCS reads and subreads are aligned to the reference genome (hg38) using the PacBio's minimap2 wrapper pbmm2.  
1. **Tag, filter and trim targets**  
A custom python script 'bam_annotate_trim_amplicons.py' is used to identify reads from a given target amplicon by alignment position and matching the forward and reverse primers. Reads lacking the correctly oriented forward and reverse primers are filtered out, and remaining reads are trimmed to the primer sites. Bam tags are used to store read metadata.
1. **Collect & report stats**  
Custom Rmarkdown script that generates the 'PacBio AmpSeq Preprocessing Report'. This report shows a breakdown of the numbers reads assigned to each amplicon, reads that produced CCS, reads that were able to be barcoded, and reads that have correct primer pair orientation ('complete' reads).  
1. **Realignment**  
BAM files are realigned after trimming  
1. **Split to sample-target BAMs**  
Monolithic CCS BAM file is split into separate BAM files for each unique sample target amplicon  


#### Allele-typing

1. **Call SNPs**  
GATK HaplotypeCaller is used to generate SNV calls in diploid mode for each sample-amplicon.  
1. **Assign Alleles to Reads**  
The set of alleles carried by individual reads is determined using WhatsHap Haplotag and a custom python script, and output as a table.  
1. **Filter, Cluster and Phase Reads**  
Several filters are applied to reads and variants in an iterative process to remove noise, outliers, chimeras and uninformative sites based variant allele frequencies. Pairwise SNP distance is then calculated between the reads and subsequently used for clustering. Clusters are determined using k-medoids (partitioning around medoids) with k values from 2-4 trialled, with the largest 2 clusters selected in each case (assuming 2 phases). Several heuristic filters are then applied to the clustering results to select the best clustering. If no 2-phase clustering is found the reads are assumed to be single-phase (i.e. Hemi/homozygous). See visualisation below.  
1. **Split Reads into Phases**  
Sample-amplicon BAMs are separated into individual phased BAMs based on the results of step 3
1. **Call short variants**  
SNPs and Indels are jointly-called using GATK HaplotypeCaller and GenotypeGVCFs with ploidy set to 1.  
1. **Call Structural Variants**  
SVs are jointly called using 'pbsv discover' and 'pbsv call'. This tools only produces diploid calls which aren't appropriate for the single-phase input BAMs, so a custom python script 'vcf_dip_to_hap.py' is used to convert heterozygous diploid calls to haploid calls based on a threshold of 0.5 B-allele frequency.  
1. **Merge Calls**  
All calls are merged into a single VCF, with haploid sample phases merged into phased diploid calls.  
1. **Annotate Variants**  
Ensembl VEP is used to provided annotations including reference population allele frequencies and variant effect predictions.  


