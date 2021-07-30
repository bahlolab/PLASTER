#### Pre-processing Parameters

* **subreads_bam (required)**  
  Path to a PacBio subreads BAM file with corresponding `.pbi` index file, output from a sequel I/II machine.
* **ref_fasta (required)**  
  Path to a FASTA file to be used as the reference genome. Also supports FTP and HTTP(S) if prefixed appropriately (i.e., with `'ftp://'`, `'http://'` or `'https://'`). 
* **barcodes_fasta (required)**  
  Path to a FASTA file containing barcode sequences and samples names, e.g.:
  ```
  >NA07439
  TGTGTATCAGTACATG
  >NA12244
  ACACGCATGACACACT
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