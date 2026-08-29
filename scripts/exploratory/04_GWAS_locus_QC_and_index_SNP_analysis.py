from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


# ============================================================
# SETTINGS
# ============================================================

INPUT_FILE = Path(
    "03_SEX_STRATIFIED_GWAS_LOCI.xlsx"
)

OUTPUT_FILE = Path(
    "04_GWAS_LOCUS_QC_INDEX_SNP.xlsx"
)

INPUT_SHEET = "S3_All_harmonized"

# Common variant threshold.
# This avoids prioritizing very rare variants with unstable effects.
MIN_AF = 0.01

# Require combined GWAS result for index-SNP selection.
REQUIRE_COMBINED = True


# ============================================================
# CHECK INPUT
# ============================================================

if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"Input file not found:\n{INPUT_FILE.resolve()}"
    )


# ============================================================
# LOAD DATA
# ============================================================

print("=" * 80)
print("1/8 Loading sex-stratified GWAS locus data...")
print("=" * 80)

df = pd.read_excel(
    INPUT_FILE,
    sheet_name=INPUT_SHEET,
    engine="openpyxl"
)

print(f"Rows loaded: {len(df):,}")
print(f"Genes/loci: {df['Gene'].nunique()}")

print("\nGenes:")
print(
    sorted(
        df["Gene"]
        .dropna()
        .unique()
        .tolist()
    )
)


# ============================================================
# INITIAL QC
# ============================================================

print("\n" + "=" * 80)
print("2/8 Variant-level QC...")
print("=" * 80)

data = df.copy()


# ------------------------------------------------------------
# Remove duplicate records
# ------------------------------------------------------------

before = len(data)

data = (
    data
    .drop_duplicates(
        subset=[
            "Gene",
            "Variant_key"
        ]
    )
    .copy()
)

print(
    f"Duplicate records removed: "
    f"{before - len(data):,}"
)


# ------------------------------------------------------------
# Biallelic check
# ------------------------------------------------------------

# VCF records with ALT containing comma represent multiallelic sites.
data["Is_biallelic"] = (
    ~data["ALT"]
    .astype(str)
    .str.contains(",", regex=False)
)

print(
    "Multiallelic records:",
    (~data["Is_biallelic"]).sum()
)


# ------------------------------------------------------------
# Valid SE
# ------------------------------------------------------------

data["Valid_SE"] = (
    data["SE_Male"].notna()
    &
    data["SE_Female"].notna()
    &
    (data["SE_Male"] > 0)
    &
    (data["SE_Female"] > 0)
)


# ------------------------------------------------------------
# Valid effect estimates
# ------------------------------------------------------------

data["Valid_effects"] = (
    data["Beta_Male"].notna()
    &
    data["Beta_Female"].notna()
)


# ------------------------------------------------------------
# Combined GWAS availability
# ------------------------------------------------------------

data["Has_combined"] = (
    data["Beta_Combined"].notna()
    &
    data["P_Combined"].notna()
)


# ============================================================
# ALLELE FREQUENCY QC
# ============================================================

print("\n" + "=" * 80)
print("3/8 Allele-frequency QC...")
print("=" * 80)


def maf_from_af(af):
    if pd.isna(af):
        return np.nan

    if af < 0 or af > 1:
        return np.nan

    return min(af, 1 - af)


for dataset in [
    "Male",
    "Female",
    "Combined"
]:
    af_col = f"AF_{dataset}"
    maf_col = f"MAF_{dataset}"

    if af_col in data.columns:
        data[maf_col] = (
            data[af_col]
            .apply(maf_from_af)
        )


# Main MAF criterion:
# require >=1% in both male and female datasets when AF exists.

if (
    "MAF_Male" in data.columns
    and
    "MAF_Female" in data.columns
):

    data["Common_variant"] = (
        (
            data["MAF_Male"].isna()
            |
            (data["MAF_Male"] >= MIN_AF)
        )
        &
        (
            data["MAF_Female"].isna()
            |
            (data["MAF_Female"] >= MIN_AF)
        )
    )

else:

    data["Common_variant"] = True


print(
    f"Variants failing MAF >= {MIN_AF}:",
    (~data["Common_variant"]).sum()
)


# ============================================================
# FINAL QC FLAG
# ============================================================

data["Pass_QC"] = (
    data["Is_biallelic"]
    &
    data["Valid_SE"]
    &
    data["Valid_effects"]
    &
    data["Common_variant"]
)

if REQUIRE_COMBINED:
    data["Pass_QC"] &= data["Has_combined"]


qc_pass = data[
    data["Pass_QC"]
].copy()


print(
    "\nVariants passing QC:",
    f"{len(qc_pass):,}"
)

print("\nPassing QC by locus:")

print(
    qc_pass
    .groupby("Gene")
    .size()
    .to_string()
)


# ============================================================
# RECALCULATE SEX HETEROGENEITY
# ============================================================

print("\n" + "=" * 80)
print("4/8 Calculating male-female heterogeneity...")
print("=" * 80)


qc_pass["Z_sex"] = (
    (
        qc_pass["Beta_Male"]
        -
        qc_pass["Beta_Female"]
    )
    /
    np.sqrt(
        qc_pass["SE_Male"] ** 2
        +
        qc_pass["SE_Female"] ** 2
    )
)


qc_pass["P_sex"] = (
    2
    *
    norm.sf(
        np.abs(
            qc_pass["Z_sex"]
        )
    )
)


def sign_pattern(row):

    bm = row["Beta_Male"]
    bf = row["Beta_Female"]

    if bm > 0 and bf > 0:
        return "same_positive"

    if bm < 0 and bf < 0:
        return "same_negative"

    if bm > 0 and bf < 0:
        return "male_positive_female_negative"

    if bm < 0 and bf > 0:
        return "male_negative_female_positive"

    return "other"


qc_pass["Human_beta_pattern"] = (
    qc_pass.apply(
        sign_pattern,
        axis=1
    )
)


# ============================================================
# SELECT INDEX SNP USING COMBINED GWAS
# ============================================================

print("\n" + "=" * 80)
print("5/8 Selecting combined-GWAS index SNPs...")
print("=" * 80)


def select_lead(group):

    x = group[
        group["P_Combined"].notna()
    ].copy()

    if x.empty:
        return None

    x = x.sort_values(
        [
            "P_Combined",
            "Distance_to_gene"
        ],
        ascending=[
            True,
            True
        ]
    )

    return x.iloc[0]


index_rows = []

for gene, group in qc_pass.groupby(
    "Gene",
    sort=True
):

    lead = select_lead(group)

    if lead is not None:
        index_rows.append(
            lead
        )


index_snps = pd.DataFrame(
    index_rows
).reset_index(drop=True)


# ============================================================
# ADD WITHIN-GENE LEAD SNP
# ============================================================

inside_rows = []

for gene, group in qc_pass.groupby(
    "Gene",
    sort=True
):

    x = group[
        group["Inside_gene"] == True
    ].copy()

    if x.empty:
        continue

    lead = select_lead(x)

    if lead is not None:
        inside_rows.append(
            lead
        )


inside_gene_leads = pd.DataFrame(
    inside_rows
).reset_index(drop=True)


# ============================================================
# CORRECT HETEROGENEITY TESTS FOR 5 PREDEFINED INDEX SNPs
# ============================================================

print("\n" + "=" * 80)
print("6/8 Multiple-testing correction for index SNPs...")
print("=" * 80)


index_snps[
    "P_sex_Bonferroni"
] = np.nan

index_snps[
    "P_sex_FDR"
] = np.nan


valid_p = (
    index_snps["P_sex"]
    .notna()
)


if valid_p.sum() > 0:

    pvals = (
        index_snps.loc[
            valid_p,
            "P_sex"
        ]
        .values
    )

    # Bonferroni across 5 predefined loci
    index_snps.loc[
        valid_p,
        "P_sex_Bonferroni"
    ] = np.minimum(
        pvals
        *
        len(pvals),
        1.0
    )

    index_snps.loc[
        valid_p,
        "P_sex_FDR"
    ] = multipletests(
        pvals,
        method="fdr_bh"
    )[1]


index_snps[
    "Nominal_sex_difference"
] = (
    index_snps["P_sex"]
    < 0.05
)

index_snps[
    "Bonferroni_sex_difference"
] = (
    index_snps[
        "P_sex_Bonferroni"
    ]
    < 0.05
)

index_snps[
    "FDR_sex_difference"
] = (
    index_snps[
        "P_sex_FDR"
    ]
    < 0.05
)


# ============================================================
# DISTANCE CLASSIFICATION
# ============================================================

def location_class(row):

    if row["Inside_gene"]:
        return "within_gene"

    d = row["Distance_to_gene"]

    if pd.isna(d):
        return "unknown"

    if d <= 50_000:
        return "within_50kb"

    if d <= 100_000:
        return "within_100kb"

    if d <= 250_000:
        return "within_250kb"

    return "within_500kb"


index_snps[
    "Index_variant_location"
] = (
    index_snps.apply(
        location_class,
        axis=1
    )
)


# ============================================================
# MOUSE INTERPRETATION
# ============================================================

def expected_role(direction):

    if direction == "increased":
        return "KO_increases_lean_mass"

    if direction == "decreased":
        return "KO_decreases_lean_mass"

    return "unknown"


index_snps[
    "Mouse_male_effect"
] = (
    index_snps[
        "Male_KO_direction"
    ]
    .apply(expected_role)
)

index_snps[
    "Mouse_female_effect"
] = (
    index_snps[
        "Female_KO_direction"
    ]
    .apply(expected_role)
)


# IMPORTANT:
# We DO NOT classify GWAS allele direction as concordant with
# mouse KO yet. A GWAS allele does not automatically indicate
# increased/decreased gene activity. That requires eQTL or other
# functional annotation.


# ============================================================
# PUBLICATION TABLE
# ============================================================

publication_columns = [
    "Gene",
    "Mouse_gene",
    "Mouse_sex_pattern",

    "Male_KO_direction",
    "Female_KO_direction",

    "rsID",
    "Variant_key",

    "REF",
    "ALT",
    "Effect_allele",

    "Inside_gene",
    "Distance_to_gene",
    "Index_variant_location",

    "AF_Male",
    "AF_Female",
    "AF_Combined",

    "Beta_Male",
    "SE_Male",
    "P_Male",

    "Beta_Female",
    "SE_Female",
    "P_Female",

    "Beta_Combined",
    "P_Combined",

    "Human_beta_pattern",

    "Z_sex",
    "P_sex",

    "P_sex_Bonferroni",
    "P_sex_FDR",

    "Nominal_sex_difference",
    "Bonferroni_sex_difference",
    "FDR_sex_difference",
]


publication_columns = [
    c
    for c in publication_columns
    if c in index_snps.columns
]


publication = (
    index_snps[
        publication_columns
    ]
    .copy()
    .sort_values("Gene")
)


# ============================================================
# LOCUS-LEVEL QC SUMMARY
# ============================================================

locus_qc = (
    data
    .groupby("Gene")
    .agg(
        N_total=(
            "Variant_key",
            "nunique"
        ),

        N_pass_QC=(
            "Pass_QC",
            "sum"
        ),

        N_inside_gene=(
            "Inside_gene",
            "sum"
        ),

        N_multiallelic=(
            "Is_biallelic",
            lambda x:
                (~x).sum()
        ),

        N_common_variants=(
            "Common_variant",
            "sum"
        ),
    )
    .reset_index()
)


# ============================================================
# SUMMARY
# ============================================================

summary = pd.DataFrame({
    "Metric": [
        "Candidate loci",
        "Raw locus variants",
        "Variants passing QC",
        "Combined-GWAS index SNPs",
        "Index SNPs with opposite beta signs",
        "Index SNPs with nominal sex difference P<0.05",
        "Index SNPs significant after Bonferroni",
        "Index SNPs significant after FDR",
    ],

    "N": [
        data["Gene"].nunique(),
        data["Variant_key"].nunique(),
        qc_pass["Variant_key"].nunique(),
        index_snps["Gene"].nunique(),

        index_snps[
            "Human_beta_pattern"
        ].isin([
            "male_positive_female_negative",
            "male_negative_female_positive"
        ]).sum(),

        index_snps[
            "Nominal_sex_difference"
        ].sum(),

        index_snps[
            "Bonferroni_sex_difference"
        ].sum(),

        index_snps[
            "FDR_sex_difference"
        ].sum(),
    ]
})


# ============================================================
# PRINT RESULTS
# ============================================================

print("\n" + "=" * 80)
print("FINAL SUMMARY")
print("=" * 80)

print(
    summary.to_string(
        index=False
    )
)


print("\n" + "=" * 80)
print("INDEX SNP PUBLICATION TABLE")
print("=" * 80)

pd.set_option(
    "display.max_columns",
    100
)

print(
    publication.to_string(
        index=False
    )
)


print("\n" + "=" * 80)
print("WITHIN-GENE LEAD VARIANTS")
print("=" * 80)

if inside_gene_leads.empty:

    print(
        "No within-gene lead variants identified."
    )

else:

    cols = [
        "Gene",
        "rsID",
        "Variant_key",
        "Effect_allele",
        "Beta_Male",
        "P_Male",
        "Beta_Female",
        "P_Female",
        "Beta_Combined",
        "P_Combined",
        "P_sex",
    ]

    cols = [
        c
        for c in cols
        if c in inside_gene_leads.columns
    ]

    print(
        inside_gene_leads[
            cols
        ]
        .sort_values("Gene")
        .to_string(index=False)
    )


# ============================================================
# WRITE OUTPUT
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

    locus_qc.to_excel(
        writer,
        sheet_name="S1_Locus_QC",
        index=False
    )

    publication.to_excel(
        writer,
        sheet_name="S2_Index_SNPs",
        index=False
    )

    inside_gene_leads.to_excel(
        writer,
        sheet_name="S3_Within_gene_leads",
        index=False
    )

    qc_pass.to_excel(
        writer,
        sheet_name="S4_All_QC_variants",
        index=False
    )

    data[
        ~data["Pass_QC"]
    ].to_excel(
        writer,
        sheet_name="S5_QC_failed",
        index=False
    )


print()
print(
    "Saved:",
    OUTPUT_FILE.resolve()
)

print("\nDONE.")