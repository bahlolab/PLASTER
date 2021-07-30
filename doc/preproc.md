### Pre-processing Parameters

* **`subreads_bam` (required)**  
  Path to a PacBio subreads BAM file with corresponding `.pbi` index file, output from a sequel I/II machine.
* **`ref_fasta` (required)**  
  Path to a FASTA file to be used as the reference genome. Also supports FTP and HTTP(S) if prefixed appropriately (i.e., with `'ftp://'`, `'http://'` or `'https://'`). 
* **`barcodes_fasta` (required)**  
  Path to a FASTA file containing barcode sequences and samples names, e.g.:
  ```
  >NA07439
  TGTGTATCAGTACATG
  >NA12244
  ACACGCATGACACACT
  ...
  ```
* **`amplicons_json` (required)**  
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
* **`run_id` (optional, default: "preproc-run")**  
  Identifier used for naming output files.
* **`ccs_min_passes` (optional, default: 3)**  
  Parameter passed to PacBio ccs as `--min-passes`.
* **`ccs_min_acc` (optional, default: 0.99)**  
  Parameter passed to PacBio ccs as `--min-snr`.
* **`ccs_min_len` (optional, default: 250)**  
  Parameter passed to PacBio ccs as `--min-length`.
* **`ccs_max_len` (optional, default: 25000)**  
  Parameter passed to PacBio ccs as `--max-length`.
* **c`cs_n_parallel` (optional, default: 100)**  
  Number of parallel jobs to run for PacBio `ccs`.
  