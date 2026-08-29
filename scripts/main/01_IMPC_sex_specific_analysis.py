from pathlib import Path
import numpy as np
import pandas as pd


# ============================================================
# SETTINGS
# ============================================================

INPUT_FILE = Path(
    "abnormal lean body mass-associations (2).xlsx"
)

OUTPUT_FILE = Path(
    "01_IMPC_sex_specific_analysis.xlsx"
)

SHEET_NAME = "Data sheet"


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

print("=" * 75)
print("1/7 Loading IMPC data...")
print("=" * 75)

df = pd.read_excel(
    INPUT_FILE,
    sheet_name=SHEET_NAME,
    engine="openpyxl"
)

print(f"Raw rows: {len(df)}")
print(f"Raw unique genes: {df['Gene'].nunique()}")

print("\nColumns:")
for col in df.columns:
    print(" -", col)


# ============================================================
# CLEAN BASIC FIELDS
# ============================================================

print("\n" + "=" * 75)
print("2/7 Cleaning fields...")
print("=" * 75)

data = df.copy()

# clean strings
for col in [
    "Gene",
    "Allele",
    "Phenotype",
    "Zygosity",
    "Sex",
    "Life stage",
    "Parameter",
    "Phenotyping center"
]:
    data[col] = (
        data[col]
        .astype(str)
        .str.strip()
    )

# normalize key columns
data["Sex_clean"] = (
    data["Sex"]
    .str.lower()
    .str.strip()
)

data["Zygosity_clean"] = (
    data["Zygosity"]
    .str.lower()
    .str.strip()
)

data["Phenotype_clean"] = (
    data["Phenotype"]
    .str.lower()
    .str.strip()
)

data["Parameter_clean"] = (
    data["Parameter"]
    .str.lower()
    .str.strip()
)

# numeric P
data["P_IMPC"] = pd.to_numeric(
    data["Most significant P-Value"],
    errors="coerce"
)


# ============================================================
# FILTER TO HOMOZYGOUS MALE/FEMALE KO
# ============================================================

print("\n" + "=" * 75)
print("3/7 Filtering homozygous male/female records...")
print("=" * 75)

filtered = data[
    (data["Zygosity_clean"] == "homozygote")
    &
    (data["Sex_clean"].isin(["male", "female"]))
].copy()

print(
    "Rows after homozygous + known-sex filtering:",
    len(filtered)
)

print(
    "Male rows:",
    (filtered["Sex_clean"] == "male").sum()
)

print(
    "Female rows:",
    (filtered["Sex_clean"] == "female").sum()
)


# ============================================================
# CLASSIFY PHENOTYPE DIRECTION
# ============================================================

print("\n" + "=" * 75)
print("4/7 Classifying knockout effects...")
print("=" * 75)

def classify_direction(x):

    x = str(x).lower().strip()

    if x == "increased lean body mass":
        return "increased"

    elif x == "decreased lean body mass":
        return "decreased"

    else:
        return "unspecified"


filtered["KO_direction"] = (
    filtered["Phenotype_clean"]
    .apply(classify_direction)
)

print(
    filtered["KO_direction"]
    .value_counts(dropna=False)
)


# ============================================================
# REMOVE UNSPECIFIED DIRECTION FROM DIRECTIONAL ANALYSIS
# BUT KEEP IT IN A SEPARATE SHEET
# ============================================================

unspecified = filtered[
    filtered["KO_direction"] == "unspecified"
].copy()

directional = filtered[
    filtered["KO_direction"].isin(
        ["increased", "decreased"]
    )
].copy()

print(
    "\nDirectional records:",
    len(directional)
)

print(
    "Unspecified records:",
    len(unspecified)
)


# ============================================================
# SUMMARIZE EACH GENE WITHIN EACH SEX
# ============================================================

print("\n" + "=" * 75)
print("5/7 Resolving within-sex records...")
print("=" * 75)


def summarize_gene_sex(group):

    directions = sorted(
        group["KO_direction"]
        .dropna()
        .unique()
        .tolist()
    )

    # gene has contradictory phenotype directions
    # WITHIN the same sex
    if len(directions) > 1:

        status = "within_sex_conflict"
        final_direction = "conflicting"

    elif len(directions) == 1:

        status = "consistent"
        final_direction = directions[0]

    else:

        status = "no_direction"
        final_direction = np.nan

    # best statistical record
    valid_p = group.dropna(
        subset=["P_IMPC"]
    )

    if not valid_p.empty:
        best_idx = valid_p["P_IMPC"].idxmin()
        best = group.loc[best_idx]
    else:
        best = group.iloc[0]

    return pd.Series({
        "KO_direction": final_direction,
        "Within_sex_status": status,
        "P_IMPC_min": group["P_IMPC"].min(),
        "N_records": len(group),
        "Allele_best": best["Allele"],
        "Phenotype_best": best["Phenotype"],
        "Parameter_best": best["Parameter"],
        "Life_stage_best": best["Life stage"],
        "Phenotyping_center_best":
            best["Phenotyping center"]
    })


gene_sex = (
    directional
    .groupby(
        ["Gene", "Sex_clean"],
        sort=True
    )
    .apply(summarize_gene_sex)
    .reset_index()
)

print(
    "\nGene-sex combinations:",
    len(gene_sex)
)

print(
    "\nWithin-sex status:"
)

print(
    gene_sex["Within_sex_status"]
    .value_counts()
)


# ============================================================
# SEPARATE WITHIN-SEX CONFLICTS
# ============================================================

within_conflicts = gene_sex[
    gene_sex["Within_sex_status"]
    == "within_sex_conflict"
].copy()

clean_gene_sex = gene_sex[
    gene_sex["Within_sex_status"]
    == "consistent"
].copy()

print(
    "\nWithin-sex conflicts:",
    len(within_conflicts)
)


# ============================================================
# CREATE MALE/FEMALE WIDE TABLE
# ============================================================

print("\n" + "=" * 75)
print("6/7 Comparing males and females...")
print("=" * 75)


male = (
    clean_gene_sex[
        clean_gene_sex["Sex_clean"] == "male"
    ]
    .drop(columns=["Sex_clean"])
    .copy()
)

female = (
    clean_gene_sex[
        clean_gene_sex["Sex_clean"] == "female"
    ]
    .drop(columns=["Sex_clean"])
    .copy()
)


male = male.rename(
    columns={
        col: f"Male_{col}"
        for col in male.columns
        if col != "Gene"
    }
)

female = female.rename(
    columns={
        col: f"Female_{col}"
        for col in female.columns
        if col != "Gene"
    }
)


comparison = pd.merge(
    male,
    female,
    on="Gene",
    how="outer"
)


# ============================================================
# SEX CLASSIFICATION
# ============================================================

def classify_sex_pattern(row):

    male_dir = row.get(
        "Male_KO_direction",
        np.nan
    )

    female_dir = row.get(
        "Female_KO_direction",
        np.nan
    )

    male_present = pd.notna(male_dir)
    female_present = pd.notna(female_dir)

    # both sexes
    if male_present and female_present:

        if male_dir == female_dir:

            return "sex_conserved"

        else:

            return "sex_opposite"

    # male only
    elif male_present and not female_present:

        return "male_only"

    # female only
    elif female_present and not male_present:

        return "female_only"

    else:

        return "unclassified"


comparison["Sex_pattern"] = (
    comparison.apply(
        classify_sex_pattern,
        axis=1
    )
)


# ============================================================
# BIOLOGICAL INTERPRETATION
# ============================================================

def biological_role(direction):

    if direction == "increased":
        return (
            "potential_negative_regulator"
        )

    elif direction == "decreased":
        return (
            "potential_positive_regulator"
        )

    return np.nan


comparison["Male_expected_role"] = (
    comparison["Male_KO_direction"]
    .apply(biological_role)
)

comparison["Female_expected_role"] = (
    comparison["Female_KO_direction"]
    .apply(biological_role)
)


# ============================================================
# CREATE GROUP TABLES
# ============================================================

sex_conserved = comparison[
    comparison["Sex_pattern"]
    == "sex_conserved"
].copy()

sex_opposite = comparison[
    comparison["Sex_pattern"]
    == "sex_opposite"
].copy()

male_only = comparison[
    comparison["Sex_pattern"]
    == "male_only"
].copy()

female_only = comparison[
    comparison["Sex_pattern"]
    == "female_only"
].copy()


# ============================================================
# SUMMARY
# ============================================================

summary = pd.DataFrame({

    "Metric": [

        "Raw IMPC records",

        "Homozygous male/female records",

        "Directional homozygous records",

        "Unique genes after within-sex conflict removal",

        "Male genes",

        "Female genes",

        "Sex-conserved genes",

        "Sex-opposite genes",

        "Male-only genes",

        "Female-only genes",

        "Within-sex conflicting gene-sex records",

        "Unspecified phenotype records"
    ],

    "N": [

        len(df),

        len(filtered),

        len(directional),

        comparison["Gene"].nunique(),

        male["Gene"].nunique(),

        female["Gene"].nunique(),

        sex_conserved["Gene"].nunique(),

        sex_opposite["Gene"].nunique(),

        male_only["Gene"].nunique(),

        female_only["Gene"].nunique(),

        len(within_conflicts),

        len(unspecified)
    ]
})


print("\nSUMMARY")
print("-" * 75)

print(
    summary.to_string(
        index=False
    )
)


# ============================================================
# ADD MORE DETAILED COUNTS
# ============================================================

print("\nSex-conserved directions:")

if not sex_conserved.empty:

    conserved_direction_counts = (
        sex_conserved[
            "Male_KO_direction"
        ]
        .value_counts()
    )

    print(
        conserved_direction_counts
    )


print("\nSex-opposite combinations:")

if not sex_opposite.empty:

    opposite_counts = (
        sex_opposite
        .groupby(
            [
                "Male_KO_direction",
                "Female_KO_direction"
            ]
        )
        .size()
    )

    print(opposite_counts)


# ============================================================
# WRITE EXCEL
# ============================================================

print("\n" + "=" * 75)
print("7/7 Writing output workbook...")
print("=" * 75)


with pd.ExcelWriter(
    OUTPUT_FILE,
    engine="openpyxl"
) as writer:

    summary.to_excel(
        writer,
        sheet_name="S0_Summary",
        index=False
    )

    filtered.to_excel(
        writer,
        sheet_name="S1_Filtered_records",
        index=False
    )

    directional.to_excel(
        writer,
        sheet_name="S2_Directional_records",
        index=False
    )

    gene_sex.to_excel(
        writer,
        sheet_name="S3_Gene_sex_summary",
        index=False
    )

    comparison.to_excel(
        writer,
        sheet_name="S4_Sex_comparison",
        index=False
    )

    sex_conserved.to_excel(
        writer,
        sheet_name="S5_Sex_conserved",
        index=False
    )

    sex_opposite.to_excel(
        writer,
        sheet_name="S6_Sex_opposite",
        index=False
    )

    male_only.to_excel(
        writer,
        sheet_name="S7_Male_only",
        index=False
    )

    female_only.to_excel(
        writer,
        sheet_name="S8_Female_only",
        index=False
    )

    within_conflicts.to_excel(
        writer,
        sheet_name="S9_WithinSex_conflicts",
        index=False
    )

    unspecified.to_excel(
        writer,
        sheet_name="S10_Unspecified",
        index=False
    )


print()
print(f"Saved: {OUTPUT_FILE.resolve()}")


# ============================================================
# PRINT IMPORTANT GENES
# ============================================================

print("\n" + "=" * 75)
print("SEX-OPPOSITE GENES")
print("=" * 75)

if sex_opposite.empty:

    print(
        "No sex-opposite genes found."
    )

else:

    print(
        sex_opposite[
            [
                "Gene",
                "Male_KO_direction",
                "Female_KO_direction",
                "Male_P_IMPC_min",
                "Female_P_IMPC_min"
            ]
        ]
        .sort_values("Gene")
        .to_string(index=False)
    )


print("\n" + "=" * 75)
print("SEX-CONSERVED GENES")
print("=" * 75)

print(
    sex_conserved[
        [
            "Gene",
            "Male_KO_direction",
            "Female_KO_direction",
            "Male_P_IMPC_min",
            "Female_P_IMPC_min"
        ]
    ]
    .sort_values("Gene")
    .head(50)
    .to_string(index=False)
)


print("\nDONE.")