from pathlib import Path
import pandas as pd


# ============================================================
# SETTINGS
# ============================================================

INPUT_FILE = Path("mart.txt")

OUTPUT_XLSX = Path(
    "mouse_human_orthologs_all_genes.xlsx"
)

OUTPUT_CSV = Path(
    "mouse_human_orthologs_all_genes.csv"
)


# ============================================================
# LOAD BIOMART TXT
# ============================================================

print("=" * 80)
print("1/6 Loading BioMart mart.txt...")
print("=" * 80)

if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"File not found:\n{INPUT_FILE.resolve()}"
    )


df = pd.read_csv(
    INPUT_FILE,
    sep=None,
    engine="python",
    dtype=str
)

print(f"Raw rows: {len(df):,}")

print("\nColumns:")
for i, col in enumerate(df.columns):
    print(i, repr(col))


# ============================================================
# EXPECTED COLUMNS
# ============================================================

required = [
    "Gene stable ID",
    "Gene name",
    "Human gene stable ID",
    "Human gene name",
    "Human homology type",
    "Human orthology confidence [0 low, 1 high]",
]

missing = [
    col
    for col in required
    if col not in df.columns
]

if missing:
    raise ValueError(
        "\nRequired columns are missing:\n"
        + "\n".join(missing)
        + "\n\nYour file contains:\n"
        + "\n".join(df.columns.astype(str))
    )


# ============================================================
# SELECT ONLY GENE-LEVEL INFORMATION
# ============================================================

print("\n" + "=" * 80)
print("2/6 Selecting gene-level ortholog fields...")
print("=" * 80)

columns_to_keep = [
    "Gene stable ID",
    "Gene name",
    "Human gene stable ID",
    "Human gene name",
    "Human homology type",
    "Human orthology confidence [0 low, 1 high]",
]

# Identity columns are optional
optional_columns = [
    "%id. target Human gene identical to query gene",
    "%id. query gene identical to target Human gene",
]

for col in optional_columns:
    if col in df.columns:
        columns_to_keep.append(col)


clean = df[
    columns_to_keep
].copy()


# ============================================================
# CLEAN STRINGS
# ============================================================

print("\n" + "=" * 80)
print("3/6 Cleaning values...")
print("=" * 80)

for col in clean.columns:
    clean[col] = (
        clean[col]
        .astype(str)
        .str.strip()
    )

    clean[col] = clean[col].replace(
        {
            "": pd.NA,
            "nan": pd.NA,
            "NaN": pd.NA,
            "None": pd.NA,
            "-": pd.NA,
        }
    )


# ============================================================
# REMOVE TRANSCRIPT-DRIVEN DUPLICATES
# ============================================================

print("\n" + "=" * 80)
print("4/6 Collapsing duplicate transcript rows...")
print("=" * 80)

before = len(clean)

clean = clean.drop_duplicates()

print(
    f"Duplicate rows removed: "
    f"{before - len(clean):,}"
)

print(
    f"Gene-level rows remaining: "
    f"{len(clean):,}"
)


# ============================================================
# RENAME COLUMNS FOR ANALYSIS
# ============================================================

rename_map = {
    "Gene stable ID":
        "Mouse_Ensembl",

    "Gene name":
        "Mouse_gene",

    "Human gene stable ID":
        "Human_Ensembl",

    "Human gene name":
        "Human_gene",

    "Human homology type":
        "Orthology_type",

    "Human orthology confidence [0 low, 1 high]":
        "Orthology_confidence",

    "%id. target Human gene identical to query gene":
        "Human_identity_to_mouse_pct",

    "%id. query gene identical to target Human gene":
        "Mouse_identity_to_human_pct",
}

clean = clean.rename(
    columns=rename_map
)


# Convert numeric columns
numeric_cols = [
    "Orthology_confidence",
    "Human_identity_to_mouse_pct",
    "Mouse_identity_to_human_pct",
]

for col in numeric_cols:
    if col in clean.columns:
        clean[col] = pd.to_numeric(
            clean[col],
            errors="coerce"
        )


# ============================================================
# QC CLASSIFICATION
# ============================================================

print("\n" + "=" * 80)
print("5/6 Classifying ortholog mappings...")
print("=" * 80)


def classify_orthology(x):

    if pd.isna(x):
        return "unmapped"

    x = str(x).lower()

    if "one2one" in x:
        return "one_to_one"

    if "one2many" in x:
        return "one_to_many"

    if "many2many" in x:
        return "many_to_many"

    return "other"


clean["Orthology_class"] = (
    clean["Orthology_type"]
    .apply(classify_orthology)
)


# Number of unique human orthologs for each mouse gene
ortholog_counts = (
    clean[
        clean["Human_gene"].notna()
    ]
    .groupby("Mouse_gene")[
        "Human_gene"
    ]
    .nunique()
    .rename("N_human_orthologs")
)


clean = clean.merge(
    ortholog_counts,
    on="Mouse_gene",
    how="left"
)

clean["N_human_orthologs"] = (
    clean["N_human_orthologs"]
    .fillna(0)
    .astype(int)
)


# Primary high-confidence one-to-one set
clean["Primary_one2one"] = (
    clean["Human_gene"].notna()
    &
    (clean["Orthology_class"] == "one_to_one")
    &
    (clean["Orthology_confidence"] == 1)
    &
    (clean["N_human_orthologs"] == 1)
)


# ============================================================
# VALIDATION OF KEY GENES
# ============================================================

key_genes = [
    "Clpp",
    "Mta1",
    "Rspo1",
    "Snap47",
    "Lingo2",
]

print("\nKey-gene validation:")
print("-" * 80)

key_check = clean[
    clean["Mouse_gene"].isin(key_genes)
][
    [
        "Mouse_gene",
        "Mouse_Ensembl",
        "Human_gene",
        "Human_Ensembl",
        "Orthology_type",
        "Orthology_confidence",
        "Primary_one2one",
    ]
].drop_duplicates()


print(
    key_check
    .sort_values("Mouse_gene")
    .to_string(index=False)
)


missing_key = (
    set(key_genes)
    -
    set(key_check["Mouse_gene"])
)

if missing_key:
    print(
        "\nWARNING: missing focal genes:",
        sorted(missing_key)
    )


# ============================================================
# SUMMARY
# ============================================================

summary = pd.DataFrame(
    {
        "Metric": [
            "Raw BioMart rows",
            "Gene-level unique mappings",
            "Unique mouse genes",
            "Mouse genes with human ortholog",
            "High-confidence one-to-one mouse genes",
            "One-to-many mouse genes",
            "Many-to-many mouse genes",
        ],

        "N": [
            len(df),

            len(clean),

            clean["Mouse_gene"].nunique(),

            clean.loc[
                clean["Human_gene"].notna(),
                "Mouse_gene"
            ].nunique(),

            clean.loc[
                clean["Primary_one2one"],
                "Mouse_gene"
            ].nunique(),

            clean.loc[
                clean["Orthology_class"]
                == "one_to_many",
                "Mouse_gene"
            ].nunique(),

            clean.loc[
                clean["Orthology_class"]
                == "many_to_many",
                "Mouse_gene"
            ].nunique(),
        ]
    }
)


print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)

print(
    summary.to_string(
        index=False
    )
)


# ============================================================
# EXPORT
# ============================================================

print("\n" + "=" * 80)
print("6/6 Writing output files...")
print("=" * 80)

primary = clean[
    clean["Primary_one2one"]
].copy()

ambiguous = clean[
    clean["Human_gene"].notna()
    &
    ~clean["Primary_one2one"]
].copy()

unmapped = clean[
    clean["Human_gene"].isna()
].copy()


with pd.ExcelWriter(
    OUTPUT_XLSX,
    engine="openpyxl"
) as writer:

    summary.to_excel(
        writer,
        sheet_name="S0_Summary",
        index=False
    )

    clean.to_excel(
        writer,
        sheet_name="S1_All_gene_mappings",
        index=False
    )

    primary.to_excel(
        writer,
        sheet_name="S2_Primary_one2one",
        index=False
    )

    ambiguous.to_excel(
        writer,
        sheet_name="S3_Ambiguous",
        index=False
    )

    unmapped.to_excel(
        writer,
        sheet_name="S4_Unmapped",
        index=False
    )

    key_check.to_excel(
        writer,
        sheet_name="S5_Key_genes_check",
        index=False
    )


clean.to_csv(
    OUTPUT_CSV,
    index=False
)


print(
    "\nSaved Excel:",
    OUTPUT_XLSX.resolve()
)

print(
    "Saved CSV:",
    OUTPUT_CSV.resolve()
)

print("\nDONE.")