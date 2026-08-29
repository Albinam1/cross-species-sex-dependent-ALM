# Cross-species analysis of sex-dependent genetic effects on appendicular lean mass

This repository contains the analysis code used for a cross-species study integrating mouse knockout phenotypes, mouse-to-human orthology, sex-stratified human appendicular lean mass (ALM) GWAS data, linkage disequilibrium (LD) analysis, and GTEx v8 regulatory data.

## Study workflow

1. **IMPC phenotype screening** — identify homozygous knockout genes associated with altered lean body mass and classify detectable phenotypes as male-only, female-only, sex-conserved, or sex-opposite.
2. **Mouse-to-human orthology** — map mouse genes to high-confidence one-to-one human orthologs using Ensembl BioMart.
3. **GRCh37 locus preparation** — obtain human gene coordinates and construct gene-centered ±500 kb candidate regions.
4. **Sex-stratified ALM GWAS screen** — extract and harmonize combined, male-specific, and female-specific GWAS variants and test male–female effect heterogeneity.
5. **Regional EPB41L1 analysis** — characterize the prioritized chromosome 20 locus and apply region-wide multiple-testing correction. **The step-05 script was not included in the files used to assemble this repository and should be added before public release.**
6. **LD clumping** — characterize correlation among sex-heterogeneous EPB41L1 variants using PLINK v1.9 and a 1000 Genomes Phase 3 European reference panel.
7. **GTEx v8 cis-eQTL analysis** — retrieve significant EPB41L1 single-tissue cis-eQTLs and compare them with prioritized GWAS variants.
8. **GTEx sex-biased eQTL analysis** — evaluate available sex-biased cis-eQTL evidence for EPB41L1.

## Repository structure

```text
scripts/
  main/
    00_prepare_biomart_table.py
    01_IMPC_sex_specific_analysis.py
    02_mouse_to_human_orthologs_all_genes.py
    03a_prepare_human_genes_for_GRCh37.py
    03b_prepare_GRCh37_coordinates.py
    04_all_genes_sex_stratified_GWAS_screen.py
    06_prepare_EPB41L1_LD_clumping.py
    06b_prepare_significant_LD_clumping.py
    07_EPB41L1_GTEx_significant_eQTL.py
    08_EPB41L1_GTEx_sex_biased_eQTL.py
  exploratory/
    02_mouse_to_human_orthologs.py
    03_extract_sex_stratified_GWAS_loci.py
    04_GWAS_locus_QC_and_index_SNP_analysis.py
```

The `exploratory/` directory contains earlier focused analyses of predefined genes and is retained for provenance; the manuscript-level systematic screen is implemented in `04_all_genes_sex_stratified_GWAS_screen.py`.

## Main statistical test

For male- and female-specific GWAS performed in mutually exclusive samples, sex heterogeneity was evaluated as

```text
Z_sex = (beta_male - beta_female) / sqrt(SE_male^2 + SE_female^2)
P_sex = 2 * Phi(-abs(Z_sex))
```

Multiple testing was controlled using Benjamini–Hochberg false discovery rate and Bonferroni correction.

## External data sources

Raw data are not redistributed in this repository. The analysis uses:

- International Mouse Phenotyping Consortium (IMPC): abnormal lean body mass phenotype (MP:0003959).
- Ensembl BioMart for mouse-to-human orthology and GRCh37 gene coordinates.
- Pei et al. appendicular lean mass GWAS summary statistics: GCST90000025 (combined), GCST90000026 (male), and GCST90000027 (female).
- 1000 Genomes Project Phase 3 European ancestry samples for LD reference calculations.
- GTEx Analysis v8 single-tissue cis-eQTL and sex-biased eQTL resources.

Users should download the relevant source data from the original providers and comply with their terms of use.

## Software

Python dependencies are listed in `requirements.txt`. LD clumping additionally requires **PLINK v1.9**.

## Reproducibility notes

- GWAS and LD coordinates are handled in **GRCh37**.
- GTEx v8 uses **GRCh38**; cross-resource comparisons are therefore based on verified variant identifiers rather than unverified direct coordinate matching.
- The primary human screen selects the strongest combined-sex ALM association within each predefined gene ±500 kb region before testing sex heterogeneity, reducing direct selection on the heterogeneity statistic.
- LD clumping is used to describe correlation among reference-covered variants and is not interpreted as formal fine-mapping.
- Direct GWAS/eQTL overlap is considered regulatory evidence, not proof of causal colocalization.

## Citation

If you use this code, please cite the associated manuscript once its final bibliographic information is available.

## Figure generation

The repository includes reproducible plotting scripts in `scripts/figures/`:

- `figure_2_candidate_gene_sex_heterogeneity.py` — candidate-gene sex-heterogeneity plot from `04_ALL_GENES_SEX_STRATIFIED_GWAS_SCREEN.xlsx`.
- `figure_3_EPB41L1_regional_sex_heterogeneity.py` — regional EPB41L1 sex-heterogeneity plot from `05_EPB41L1_REGIONAL_SEX_ANALYSIS.xlsx`.
- `figure_4_EPB41L1_genetic_regulatory_characterization.py` — combined sex-stratified effect and GTEx regulatory figure from the regional, conventional GTEx, and (when present) sex-biased GTEx outputs.

Each plotting script saves PNG (600 dpi), PDF, and SVG files to a local `figures/` directory. The plots are generated directly from the pipeline output workbooks rather than from hard-coded numerical values.

`Figure 1` is a conceptual study-workflow schematic rather than a statistical plot and is therefore not generated from an analysis table.
