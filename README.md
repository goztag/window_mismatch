# window-mismatch

Calculate pairwise mismatch in _n_ SNP windows from VCF files.

## Dependencies
- pysam
- itertools
- pandas

## Input 
Population VCF file with biallelic SNPs only and no missing data.

## Options

```bash

-vcf    <filename>      "Input VCF file"
-win    <INT>           "Window size (#SNPs)"
-tot    <INT>           "Total SNP count"
-out    <filename>      "Output csv file"

```

## Example

```bash
python3 window_mismatch.py -vcf in.vcf -win 1000 -tot 100000 -out out.csv
```

## Output

Table of average pairwise mismatch values per window. 


