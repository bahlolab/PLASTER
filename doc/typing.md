#### Allele-typing Parameters

* **ref_fasta (required)**  
  Path to a FASTA file to be used as the reference genome. Also supports FTP and HTTP(S) if prefixed appropriately (i.e., with `'ftp://'`, `'http://'` or `'https://'`). 
* **manifest (required)**  
  Path to a CSV file, output from the pre-processing stage, with the following format:
  ```
  sample,amplicon,n_reads,bam_file
  NA07439,CYP2D6,218,/PATH/TO/SM-NA07439.AM-CYP2D6.bam
  NA07439,CYP2D6,151,/PATH/TO/SM-NA07439.AM-CYP2D6.bam
  NA07439,CYP2D7,14,/PATH/TO/SM-NA07439.AM-CYP2D7.bam
  ...
  ```
* **amplicons_json (required)**  
  Path to a JSON file describing target amplicons with the following format:
  ```JSON
  {
    "CYP2D6": {
      "chrom": "chr22",
      "start": 42125398,
      "end": 42131503,
      "strand": "-",
      "fusion": "CYP2D7",
      "pharmvar_gene": "CYP2D6",
      "pharmvar_ver": "4.2.6.1",
      "vep_feature":"ENST00000645361"
    },
    "CYP2D7": {
      "chrom": "chr22",
      "start": 42137550,
      "end": 42145176,
      "strand": "-",
    }
  }
  ```
  Note the following fields are optional:
    * `"fusion"` - enables fusion detection by specifying the name of a second amplicon to check for fusions with
    * `"pharmvar_gene"`, `"pharmvar_ver"` - enables PharmVar star allele assignment
    * `"vep_feature"` - used in PharmVar star allele assignment to filter for variants affecting a given Ensembl trasnscript use VEP annotation
* **ploidy (optional, default: 2)**  
  Expected ploidy (i.e., copy number) for all samples.
* **copy_num (optional, default: null)**  
  Path to a CSV file specifying the copy number of any subset of input samples and amplicons. This will override the default ploidy argument for these samples. Should be in the following format:
  ```
  sample,amplicon,copy_num
  NA07439,CYP2D6,3
  NA10005,CYP2D6,2
  NA12244,CYP2D6,2
  ...
  ```
* **max_reads (optional, default: 500)**  
  Maximum number of reads for genotyping. Sample amplicons with greater than this number of reads will be downsampled to this number.
* **min_reads (optional, default: 25)**  
  Minimum number of reads for initial genotyping. Samples amplicons with fewer than this number of reads will be excluded.
* **min_reads_phased (optional, default: 5)**  
  Minimum number of reads for phase genotyping. Samples amplicons phases with fewer than this number of reads will be excluded.
* **qd_1 (optional, default: '2.0')**  
  Minimum GATK QD (QualByDepth) for initial genotyping. This is fairly lenient as initial genotyping  is performed on unfiltered, unphased reads.
* **qd_2 (optional, default: '20.0')**  
  Minimum GATK QD (QualByDepth) for phased genotyping. This is stricter as input read quality should be much higher.
* **vep_assembly (optional, default: null)**  
  `--assembly` parameter passed to Ensembl `VEP`. Only used if `vep_cache_ver` also specified.
* **vep_cache_ver (optional, default: null)**  
  `--cache_version` parameter passed to Ensembl `VEP`. Only used if `vep_assembly` also specified.

  