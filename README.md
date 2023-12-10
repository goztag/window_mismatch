# window-mismatch

Calculate pairwise mismatch in _n_ SNP windows from VCF files.

## Dependencies
pysam
itertools
pandas

## Input 
Population VCF file with biallelic SNPs only and no missing data.

## Options

```bash

-vcf    <filename>      "Input VCF file"
-win    <INT>           "Window size (#SNPs)"
-tot    <INT>           "Total SNP count"
-out    <filename>      "Output csv file"

```



