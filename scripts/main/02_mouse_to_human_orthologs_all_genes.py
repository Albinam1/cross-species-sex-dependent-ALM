from pathlib import Path
import pandas as pd


# ============================================================
# SETTINGS
# ============================================================

IMPC_FILE = Path("01_IMPC_sex_specific_analysis.xlsx")

BIOMART_FILE = Path(
    "mouse_human_orthologs_all_genes.xlsx"
)

OUTPUT_FILE = Path(
    "02_MOUSE_HUMAN_ORTHOLOGS_ALL_GENES.xlsx"
)

IMPC_SHEET = "S4_Sex_comparison"
BIOMART_SHEET = "S1_All_gene_mappings"


# ============================================================
# 1. CHECK FILES
# ============================================================

print("=" * 80)
print("1/8 Checking input files...")
print("=" * 80)

for f in [IMPC_FILE, BIOMART_FILE]:
    if not f.exists():
        raise FileNotFoundError(
            f"File not found:\n{f.resolve()}"
        )

print("IMPC:", IMPC_FILE.resolve())
print("BioMart:", BIOMART_FILE.resolve())


# ============================================================
# 2. LOAD IMPC
# ============================================================

print("\n" + "=" * 80)
print("2/8 Loading IMPC sex-specific results...")
print("=" * 80)

impc = pd.read_excel(
    IMPC_FILE,
    sheet_name=IMPC_SHEET,
    engine="openpyxl"
)

print(f"Raw IMPC rows: {len(impc):,}")
print(f"Unique IMPC genes: {impc['Gene'].nunique():,}")

print("\nIMPC columns:")
for col in impc.columns:
    print(" -", repr(col))


# ------------------------------------------------------------
# Clean gene symbols
# ------------------------------------------------------------

impc["Gene"] = (
    impc["Gene"]
    .astype(str)
    .str.strip()
)

impc = impc[
    impc["Gene"].notna()
    &
    (impc["Gene"] != "")
    &
    (impc["Gene"] != "nan")
].copy()


# ============================================================
# 3. LOAD PREPARED BIOMART
# ============================================================

print("\n" + "=" * 80)
print("3/8 Loading prepared Ensembl BioMart mapping...")
print("=" * 80)

bm = pd.read_excel(
    BIOMART_FILE,
    sheet_name=BIOMART_SHEET,
    engine="openpyxl"
)

print(f"Prepared BioMart rows: {len(bm):,}")

print("\nBioMart columns:")
for col in bm.columns:
    print(" -", repr(col))


required_biomart = [
    "Mouse_Ensembl",
    "Mouse_gene",
    "Human_Ensembl",
    "Human_gene",
    "Orthology_type",
    "Orthology_confidence",
    "Orthology_class",
    "N_human_orthologs",
    "Primary_one2one",
]

missing = [
    col
    for col in required_biomart
    if col not in bm.columns
]

if missing:
    raise ValueError(
        "\nRequired BioMart columns are missing:\n"
        + "\n".join(missing)
    )


# ============================================================
# 4. CLEAN BIOMART
# ============================================================

print("\n" + "=" * 80)
print("4/8 Cleaning BioMart mapping...")
print("=" * 80)

for col in [
    "Mouse_gene",
    "Human_gene",
    "Mouse_Ensembl",
    "Human_Ensembl",
]:
    bm[col] = (
        bm[col]
        .astype("string")
        .str.strip()
    )


# ------------------------------------------------------------
# Convert Primary_one2one robustly to Boolean
# ------------------------------------------------------------

if bm["Primary_one2one"].dtype != bool:

    bm["Primary_one2one"] = (
        bm["Primary_one2one"]
        .astype(str)
        .str.strip()
        .str.lower()
        .map({
            "true": True,
            "false": False,
            "1": True,
            "0": False,
        })
        .fillna(False)
    )


# ------------------------------------------------------------
# Remove exact duplicate mappings
# ------------------------------------------------------------

before = len(bm)

bm = bm.drop_duplicates(
    subset=[
        "Mouse_gene",
        "Mouse_Ensembl",
        "Human_gene",
        "Human_Ensembl",
        "Orthology_type",
    ]
).copy()

print(
    "Duplicate mapping rows removed:",
    f"{before - len(bm):,}"
)

print(
    "Unique mouse genes in BioMart:",
    f"{bm['Mouse_gene'].nunique():,}"
)


# ============================================================
# 5. MERGE IMPC WITH BIOMART
# ============================================================

print("\n" + "=" * 80)
print("5/8 Mapping IMPC genes to human orthologs...")
print("=" * 80)

mapped = impc.merge(
    bm,
    left_on="Gene",
    right_on="Mouse_gene",
    how="left",
    validate="one_to_many"
)


n_impc = impc["Gene"].nunique()

n_with_human = mapped.loc[
    mapped["Human_gene"].notna(),
    "Gene"
].nunique()

n_without_human = (
    n_impc - n_with_human
)


print(f"IMPC genes entering mapping: {n_impc}")
print(f"Genes with >=1 human ortholog: {n_with_human}")
print(f"Genes without human ortholog: {n_without_human}")


# ============================================================
# 6. DEFINE PRIMARY ONE-TO-ONE SET
# ============================================================

print("\n" + "=" * 80)
print("6/8 Defining primary high-confidence one-to-one set...")
print("=" * 80)

primary = mapped[
    mapped["Primary_one2one"] == True
].copy()


# One row per mouse-human gene pair
primary = primary.drop_duplicates(
    subset=[
        "Gene",
        "Human_gene"
    ]
).copy()


print(
    "High-confidence one-to-one IMPC genes:",
    primary["Gene"].nunique()
)

print(
    "Unique human orthologs:",
    primary["Human_gene"].nunique()
)


# ============================================================
# SEX-PATTERN SUBSETS
# ============================================================

sex_conserved = primary[
    primary["Sex_pattern"] == "sex_conserved"
].copy()

sex_opposite = primary[
    primary["Sex_pattern"] == "sex_opposite"
].copy()

male_only = primary[
    primary["Sex_pattern"] == "male_only"
].copy()

female_only = primary[
    primary["Sex_pattern"] == "female_only"
].copy()


# ============================================================
# NON-PRIMARY / UNMAPPED
# ============================================================

primary_genes = set(
    primary["Gene"].dropna()
)

ambiguous = mapped[
    mapped["Human_gene"].notna()
    &
    ~mapped["Gene"].isin(primary_genes)
].copy()


mapped_genes = set(
    mapped.loc[
        mapped["Human_gene"].notna(),
        "Gene"
    ]
)

unmapped_gene_names = (
    set(impc["Gene"])
    -
    mapped_genes
)

unmapped = impc[
    impc["Gene"].isin(unmapped_gene_names)
].copy()


# ============================================================
# 7. CHECK FIVE PREDEFINED FOCAL GENES
# ============================================================

print("\n" + "=" * 80)
print("7/8 Checking predefined focal genes...")
print("=" * 80)

focal_genes = [
    "Clpp",
    "Mta1",
    "Rspo1",
    "Snap47",
    "Lingo2",
]


focal_columns = [
    "Gene",
    "Mouse_Ensembl",
    "Human_gene",
    "Human_Ensembl",
    "Sex_pattern",
    "Male_KO_direction",
    "Female_KO_direction",
    "Male_P_IMPC_min",
    "Female_P_IMPC_min",
    "Orthology_type",
    "Orthology_confidence",
    "Human_identity_to_mouse_pct",
    "Mouse_identity_to_human_pct",
    "Primary_one2one",
]

focal_columns = [
    c for c in focal_columns
    if c in primary.columns
]


focal = primary[
    primary["Gene"].isin(focal_genes)
][focal_columns].copy()


focal = focal.sort_values("Gene")


if focal.empty:

    print("WARNING: No focal genes found.")

else:

    print(
        focal.to_string(index=False)
    )


missing_focal = (
    set(focal_genes)
    -
    set(focal["Gene"])
)

if missing_focal:

    print(
        "\nWARNING: focal genes missing from primary mapping:"
    )

    for gene in sorted(missing_focal):
        print(" -", gene)

else:

    print(
        "\nAll 5 focal genes successfully retained."
    )


# ============================================================
# SUMMARY
# ============================================================

summary = pd.DataFrame({

    "Metric": [

        "IMPC genes entering ortholog mapping",

        "Mouse genes with >=1 human ortholog",

        "Mouse genes without mapped human ortholog",

        "High-confidence one-to-one ortholog genes",

        "Unique human genes in primary set",

        "Ambiguous/non-primary mapped genes",

        "Primary sex-conserved genes",

        "Primary sex-opposite genes",

        "Primary male-only genes",

        "Primary female-only genes",

        "Predefined focal genes retained",
    ],

    "N": [

        n_impc,

        n_with_human,

        n_without_human,

        primary["Gene"].nunique(),

        primary["Human_gene"].nunique(),

        ambiguous["Gene"].nunique(),

        sex_conserved["Gene"].nunique(),

        sex_opposite["Gene"].nunique(),

        male_only["Gene"].nunique(),

        female_only["Gene"].nunique(),

        focal["Gene"].nunique(),
    ]
})


print("\n" + "=" * 80)
print("FINAL SUMMARY")
print("=" * 80)

print(
    summary.to_string(index=False)
)


# ============================================================
# CREATE GWAS INPUT
# ============================================================

gwas_columns = [
    "Gene",
    "Mouse_Ensembl",
    "Human_gene",
    "Human_Ensembl",

    "Sex_pattern",

    "Male_KO_direction",
    "Male_P_IMPC_min",

    "Female_KO_direction",
    "Female_P_IMPC_min",

    "Orthology_type",
    "Orthology_confidence",

    "Human_identity_to_mouse_pct",
    "Mouse_identity_to_human_pct",
]


gwas_columns = [
    c
    for c in gwas_columns
    if c in primary.columns
]


gwas_input = primary[
    gwas_columns
].copy()


gwas_input = (
    gwas_input
    .drop_duplicates(
        subset=[
            "Gene",
            "Human_gene"
        ]
    )
    .sort_values(
        [
            "Sex_pattern",
            "Human_gene"
        ]
    )
)


# ============================================================
# 8. WRITE OUTPUT
# ============================================================

print("\n" + "=" * 80)
print("8/8 Writing output workbook...")
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

    mapped.to_excel(
        writer,
        sheet_name="S1_All_mapping",
        index=False
    )

    primary.to_excel(
        writer,
        sheet_name="S2_Primary_one2one",
        index=False
    )

    gwas_input.to_excel(
        writer,
        sheet_name="S3_GWAS_input",
        index=False
    )

    sex_conserved.to_excel(
        writer,
        sheet_name="S4_Sex_conserved",
        index=False
    )

    sex_opposite.to_excel(
        writer,
        sheet_name="S5_Sex_opposite",
        index=False
    )

    male_only.to_excel(
        writer,
        sheet_name="S6_Male_only",
        index=False
    )

    female_only.to_excel(
        writer,
        sheet_name="S7_Female_only",
        index=False
    )

    ambiguous.to_excel(
        writer,
        sheet_name="S8_Ambiguous",
        index=False
    )

    unmapped.to_excel(
        writer,
        sheet_name="S9_Unmapped",
        index=False
    )

    focal.to_excel(
        writer,
        sheet_name="S10_Focal_5_genes",
        index=False
    )


print(
    "\nSaved:",
    OUTPUT_FILE.resolve()
)

print("\nDONE.")