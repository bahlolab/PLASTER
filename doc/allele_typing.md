### PLASTER Allele-typing parameters

**Required parameters**

* **sample_manifest** - A TSV file with sample names, barcode sequence names and amplicon names run for each sample, with the following format:
  ```
  sample	barcode	amplicons
  NA07439	BC04	CYP2D6;CYP2D7
  NA10005	BC22	CYP2D6;CYP2D7
  NA17203	BC59	CYP2D6;CYP2D7
  ...  ...  ...
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
  and `barcodes_fasta` is a FASTA file with sequence names matching the "barcode"" column in `sample_manifest`. Note that only symmetric barcoding is supported (same barcode on forward and reverse primers).