from pathlib import Path
import pandas as pd


# ============================================================
# SETTINGS
# ============================================================

INPUT_FILE = Path(
    "02_MOUSE_HUMAN_ORTHOLOGS_ALL_GENES.xlsx"
)

INPUT_SHEET = "S3_GWAS_input"

OUTPUT_TXT = Path(
    "03_human_ensembl_ids_for_GRCh37.txt"
)

OUTPUT_XLSX = Path(
    "03_human_genes_for_GRCh37_lookup.xlsx"
)


# ============================================================
# LOAD
# ============================================================

if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"File not found:\n{INPUT_FILE.resolve()}"
    )

df = pd.read_excel(
    INPUT_FILE,
    sheet_name=INPUT_SHEET,
    engine="openpyxl"
)

print("=" * 80)
print("INPUT")
print("=" * 80)

print("Rows:", len(df))
print("Unique mouse genes:", df["Gene"].nunique())
print("Unique human genes:", df["Human_gene"].nunique())
print("Unique Human Ensembl IDs:", df["Human_Ensembl"].nunique())


# ============================================================
# CLEAN ENSEMBL IDS
# ============================================================

genes = df[
    [
        "Gene",
        "Human_gene",
        "Human_Ensembl",
        "Sex_pattern",
        "Male_KO_direction",
        "Female_KO_direction",
    ]
].copy()


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


genes = (
    genes
    .dropna(
        subset=["Human_Ensembl"]
    )
    .drop_duplicates(
        subset=["Human_Ensembl"]
    )
    .sort_values("Human_Ensembl")
)


# ============================================================
# EXPORT LIST FOR BIOMART
# ============================================================

with open(
    OUTPUT_TXT,
    "w",
    encoding="utf-8"
) as f:

    for ens in genes["Human_Ensembl"]:
        f.write(str(ens) + "\n")


genes.to_excel(
    OUTPUT_XLSX,
    index=False
)


# ============================================================
# QC
# ============================================================

print("\n" + "=" * 80)
print("FINAL SUMMARY")
print("=" * 80)

print(
    "Human Ensembl IDs exported:",
    genes["Human_Ensembl"].nunique()
)

print(
    "Sex patterns:"
)

print(
    genes["Sex_pattern"]
    .value_counts()
)

print("\nSaved:")
print(OUTPUT_TXT.resolve())
print(OUTPUT_XLSX.resolve())

print("\nDONE.")