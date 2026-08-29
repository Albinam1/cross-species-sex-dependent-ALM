from pathlib import Path
import pandas as pd


# ============================================================
# SETTINGS
# ============================================================

COORD_FILE = Path(
    "mart_GRCh37_coordinates.txt"
)

GENE_FILE = Path(
    "03_human_genes_for_GRCh37_lookup.xlsx"
)

OUTPUT_FILE = Path(
    "03_HUMAN_GENE_COORDINATES_GRCh37.xlsx"
)


# ============================================================
# LOAD BIOMART AUTOMATICALLY
# ============================================================

print("=" * 80)
print("1/6 Loading GRCh37 BioMart coordinates...")
print("=" * 80)

coord = pd.read_csv(
    COORD_FILE,
    sep=None,
    engine="python",
    dtype=str
)

print("Rows:", len(coord))

print("\nColumns:")
for i, c in enumerate(coord.columns):
    print(i, repr(c))


# ============================================================
# IDENTIFY EXPECTED COLUMNS
# ============================================================

required = [
    "Gene stable ID",
    "Chromosome/scaffold name",
    "Gene start (bp)",
    "Gene end (bp)",
]

missing = [
    x for x in required
    if x not in coord.columns
]

if missing:
    raise ValueError(
        "Missing GRCh37 columns:\n"
        + "\n".join(missing)
    )


# ============================================================
# STANDARDIZE
# ============================================================

rename = {
    "Gene stable ID":
        "Human_Ensembl",

    "Gene name":
        "GRCh37_gene_name",

    "Chromosome/scaffold name":
        "Chromosome",

    "Gene start (bp)":
        "Gene_start_GRCh37",

    "Gene end (bp)":
        "Gene_end_GRCh37",

    "Gene type":
        "Gene_type",

    "Strand":
        "Strand",
}

coord = coord.rename(
    columns={
        k: v
        for k, v in rename.items()
        if k in coord.columns
    }
)


coord["Human_Ensembl"] = (
    coord["Human_Ensembl"]
    .astype(str)
    .str.strip()
    .str.replace(
        r"\.\d+$",
        "",
        regex=True
    )
)


coord["Gene_start_GRCh37"] = (
    pd.to_numeric(
        coord["Gene_start_GRCh37"],
        errors="coerce"
    )
)

coord["Gene_end_GRCh37"] = (
    pd.to_numeric(
        coord["Gene_end_GRCh37"],
        errors="coerce"
    )
)


# ============================================================
# REMOVE DUPLICATES
# ============================================================

coord = (
    coord
    .drop_duplicates(
        subset=[
            "Human_Ensembl",
            "Chromosome",
            "Gene_start_GRCh37",
            "Gene_end_GRCh37",
        ]
    )
)


# ============================================================
# LOAD OUR 441 GENES
# ============================================================

genes = pd.read_excel(
    GENE_FILE,
    engine="openpyxl"
)


genes["Human_Ensembl"] = (
    genes["Human_Ensembl"]
    .astype(str)
    .str.strip()
    .str.replace(
        r"\.\d+$",
        "",
        regex=True
    )
)


# ============================================================
# MERGE
# ============================================================

print("\n" + "=" * 80)
print("2/6 Merging with IMPC-derived human orthologs...")
print("=" * 80)

merged = genes.merge(
    coord,
    on="Human_Ensembl",
    how="left"
)


# ============================================================
# CHROMOSOME QC
# ============================================================

print("\n" + "=" * 80)
print("3/6 Chromosome QC...")
print("=" * 80)

# Main autosomes + X/Y
canonical_chr = (
    [str(x) for x in range(1, 23)]
    + ["X", "Y"]
)

merged["Canonical_chromosome"] = (
    merged["Chromosome"]
    .astype(str)
    .isin(canonical_chr)
)


# ============================================================
# VALID COORDINATES
# ============================================================

merged["Valid_GRCh37_coordinates"] = (
    merged["Chromosome"].notna()
    &
    merged["Gene_start_GRCh37"].notna()
    &
    merged["Gene_end_GRCh37"].notna()
)


valid = merged[
    merged["Valid_GRCh37_coordinates"]
    &
    merged["Canonical_chromosome"]
].copy()


missing_coord = merged[
    ~merged["Valid_GRCh37_coordinates"]
].copy()


noncanonical = merged[
    merged["Valid_GRCh37_coordinates"]
    &
    ~merged["Canonical_chromosome"]
].copy()


# ============================================================
# LOCUS WINDOW
# ============================================================

WINDOW_BP = 500_000

valid["Locus_start_GRCh37"] = (
    valid["Gene_start_GRCh37"]
    - WINDOW_BP
).clip(lower=1)

valid["Locus_end_GRCh37"] = (
    valid["Gene_end_GRCh37"]
    + WINDOW_BP
)


# ============================================================
# SUMMARY
# ============================================================

summary = pd.DataFrame({
    "Metric": [
        "Human orthologs expected",
        "Genes with GRCh37 coordinates",
        "Genes with canonical-chromosome coordinates",
        "Genes missing GRCh37 coordinates",
        "Genes on non-canonical contigs",
    ],

    "N": [
        genes["Human_Ensembl"].nunique(),

        merged.loc[
            merged["Valid_GRCh37_coordinates"],
            "Human_Ensembl"
        ].nunique(),

        valid["Human_Ensembl"].nunique(),

        missing_coord["Human_Ensembl"].nunique(),

        noncanonical["Human_Ensembl"].nunique(),
    ]
})


print("\n" + "=" * 80)
print("FINAL SUMMARY")
print("=" * 80)

print(
    summary.to_string(
        index=False
    )
)


# ============================================================
# CHECK FOCAL GENES
# ============================================================

focal = valid[
    valid["Human_gene"].isin(
        [
            "CLPP",
            "MTA1",
            "RSPO1",
            "SNAP47",
            "LINGO2"
        ]
    )
].copy()


print("\nFOCAL GENES")
print("-" * 80)

print(
    focal[
        [
            "Human_gene",
            "Human_Ensembl",
            "Chromosome",
            "Gene_start_GRCh37",
            "Gene_end_GRCh37",
            "Locus_start_GRCh37",
            "Locus_end_GRCh37",
        ]
    ]
    .to_string(index=False)
)


# ============================================================
# EXPORT
# ============================================================

print("\n" + "=" * 80)
print("6/6 Saving workbook...")
print("=" * 80)

with pd.ExcelWriter(
    OUTPUT_FILE,
    engine="openpyxl"
) as writer:

    summary.to_excel(
        writer,
        sheet_name="S0_Summary",
        index=False
    )

    valid.to_excel(
        writer,
        sheet_name="S1_GWAS_loci_GRCh37",
        index=False
    )

    merged.to_excel(
        writer,
        sheet_name="S2_All_mapping",
        index=False
    )

    missing_coord.to_excel(
        writer,
        sheet_name="S3_Missing_coordinates",
        index=False
    )

    noncanonical.to_excel(
        writer,
        sheet_name="S4_Noncanonical",
        index=False
    )

    focal.to_excel(
        writer,
        sheet_name="S5_Focal_genes",
        index=False
    )


print(
    "\nSaved:",
    OUTPUT_FILE.resolve()
)

print("\nDONE.")