from pathlib import Path
import pandas as pd


INPUT_FILE = Path(
    "05_EPB41L1_REGIONAL_SEX_ANALYSIS.xlsx"
)

OUTPUT_FILE = Path(
    "06_EPB41L1_Bonferroni_LD_input.tsv"
)

df = pd.read_excel(
    INPUT_FILE,
    sheet_name="S5_All_region_variants",
    engine="openpyxl"
)

# ------------------------------------------------------------
# Regional Bonferroni significance
# ------------------------------------------------------------

sig = df[
    (
        df["P_sex_Bonferroni_region"] < 0.05
    )
    &
    (
        df["rsID"]
        .astype(str)
        .str.match(r"^rs\d+$")
    )
].copy()


sig = sig[
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
        "Distance_to_EPB41L1",
    ]
]


sig = sig.rename(
    columns={
        "rsID": "SNP",
        "P_sex": "P"
    }
)


sig = (
    sig
    .sort_values("P")
    .drop_duplicates("SNP")
)


sig.to_csv(
    OUTPUT_FILE,
    sep="\t",
    index=False
)


print(
    "Bonferroni-significant variants:",
    len(sig)
)

print(
    sig.head(20).to_string(index=False)
)

print(
    "\nSaved:",
    OUTPUT_FILE.resolve()
)