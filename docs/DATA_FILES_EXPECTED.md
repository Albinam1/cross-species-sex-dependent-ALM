# Expected external input files

The scripts use explicit relative filenames. Download or generate the corresponding resources before running the workflow, or edit the SETTINGS blocks in the scripts.

Examples include:
- IMPC export for abnormal lean body mass.
- Ensembl BioMart orthology export (`mart.txt`).
- GRCh37 BioMart coordinate export (`mart_GRCh37_coordinates.txt`).
- ALM GWAS VCFs: `ebi-a-GCST90000025.vcf.gz`, `ebi-a-GCST90000026.vcf.gz`, `ebi-a-GCST90000027.vcf.gz`.
- GTEx v8 sex-biased eQTL archive.
- 1000 Genomes Phase 3 EUR PLINK reference files for LD clumping.

The manuscript-level workflow also depends on the output of a regional EPB41L1 analysis step (`05_EPB41L1_REGIONAL_SEX_ANALYSIS.xlsx`). The code that creates this file was not among the scripts supplied when this repository package was assembled and should be added before the repository is made public.
