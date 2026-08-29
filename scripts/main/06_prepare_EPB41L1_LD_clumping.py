from pathlib import Path
import pandas as pd
import numpy as np

# ============================================================
# SETTINGS
# ============================================================

INPUT_FILE = Path(
    "05_EPB41L1_REGIONAL_SEX_ANALYSIS.xlsx"
)

INPUT_SHEET = "S5_All_region_variants"

OUTPUT_ALL = Path(
    "06_EPB41L1_LD_clumping_input.tsv"
)

OUTPUT_STRONG = Path(
    "06_EPB41L1_strong_sex_variants.tsv"
)

OUTPUT_RSID = Path(
    "06_EPB41L1_rsIDs_for_LD.txt"
)

# Variant-level screening threshold used only for LD clumping.
P_THRESHOLD = 0.05

# Strong regional signal subset
STRONG_P_THRESHOLD = 1e-4


# ============================================================
# LOAD
# ============================================================

df = pd.read_excel(
    INPUT_FILE,
    sheet_name=INPUT_SHEET,
    engine="openpyxl"
)

print("=" * 80)
print("EPB41L1 LD CLUMPING PREPARATION")
print("=" * 80)

print("Input variants:", len(df))


# ============================================================
# QC
# ============================================================

data = df.copy()

# valid rsID only
data["Valid_rsID"] = (
    data["rsID"]
    .astype(str)
    .str.match(r"^rs\d+$")
)

# SNP/indel must have heterogeneity result
data["Valid_for_LD"] = (
    data["Valid_rsID"]
    &
    data["P_sex"].notna()
    &
    data["Beta_Male"].notna()
    &
    data["Beta_Female"].notna()
)


valid = data[
    data["Valid_for_LD"]
].copy()


# ============================================================
# SEX-DIFFERENTIAL SET
# ============================================================

sex_candidates = valid[
    valid["P_sex"] < P_THRESHOLD
].copy()

strong = valid[
    valid["P_sex"] < STRONG_P_THRESHOLD
].copy()


# ============================================================
# CREATE PLINK CLUMPING TABLE
# ============================================================

clump = sex_candidates[
    [
        "rsID",
        "CHROM",
        "POS",
        "P_sex",
        "Beta_Male",
        "Beta_Female",
        "Beta_difference_M_minus_F",
        "P_Combined",
        "Inside_EPB41L1",
        "Distance_to_EPB41L1"
    ]
].copy()


clump = clump.rename(
    columns={
        "rsID": "SNP",
        "P_sex": "P"
    }
)


clump = (
    clump
    .sort_values("P")
    .drop_duplicates("SNP")
)


# ============================================================
# STRONG SIGNAL TABLE
# ============================================================

strong = (
    strong
    .sort_values("P_sex")
    .drop_duplicates("rsID")
)


# ============================================================
# EXPORT
# ============================================================

clump.to_csv(
    OUTPUT_ALL,
    sep="\t",
    index=False
)

strong.to_csv(
    OUTPUT_STRONG,
    sep="\t",
    index=False
)

with open(
    OUTPUT_RSID,
    "w",
    encoding="utf-8"
) as f:

    for rsid in strong["rsID"]:
        f.write(str(rsid) + "\n")


# ============================================================
# SUMMARY
# ============================================================

print("\nSUMMARY")
print("-" * 80)

print("Valid rsID variants:", len(valid))
print("P_sex < 0.05:", len(sex_candidates))
print("P_sex < 1e-4:", len(strong))

print("\nTop 20 variants:")

print(
    strong[
        [
            "rsID",
            "POS",
            "Beta_Male",
            "Beta_Female",
            "P_sex",
            "P_sex_FDR_region",
            "P_sex_Bonferroni_region"
        ]
    ]
    .head(20)
    .to_string(index=False)
)

print("\nSaved:")
print(OUTPUT_ALL.resolve())
print(OUTPUT_STRONG.resolve())
print(OUTPUT_RSID.resolve())

print("\nDONE.")