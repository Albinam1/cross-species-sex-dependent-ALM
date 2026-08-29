from pathlib import Path
import gzip
import numpy as np
import pandas as pd

from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


# ============================================================
# SETTINGS
# ============================================================

MALE_FILE = Path(
    "ebi-a-GCST90000026.vcf.gz"
)

FEMALE_FILE = Path(
    "ebi-a-GCST90000027.vcf.gz"
)

COMBINED_FILE = Path(
    "ebi-a-GCST90000025.vcf.gz"
)


OUTPUT_XLSX = Path(
    "05_EPB41L1_REGIONAL_SEX_ANALYSIS.xlsx"
)

OUTPUT_TSV = Path(
    "05_EPB41L1_all_variants.tsv.gz"
)


# ============================================================
# GRCh37 REGION
# ============================================================

CHROM = "20"

GENE_START = 34679426
GENE_END   = 34820721

REGION_START = 34179426
REGION_END   = 35320721

GENE = "EPB41L1"


# ============================================================
# QC
# ============================================================

MIN_MAF = 0.01


# ============================================================
# FUNCTIONS
# ============================================================

def safe_float(x):

    if x is None:
        return np.nan

    if x in {"", ".", "NA", "nan"}:
        return np.nan

    try:
        return float(x)

    except Exception:
        return np.nan


def maf_from_af(af):

    if pd.isna(af):
        return np.nan

    if af < 0 or af > 1:
        return np.nan

    return min(
        af,
        1 - af
    )


def lp_to_p(lp):

    if pd.isna(lp):
        return np.nan

    # ordinary floating point representation
    if lp <= 300:
        return 10 ** (-lp)

    return 0.0


def parse_format(
    format_string,
    sample_string
):

    keys = format_string.split(":")
    vals = sample_string.split(":")

    return dict(
        zip(keys, vals)
    )


def clean_chrom(x):

    return (
        str(x)
        .replace("chr", "")
        .strip()
    )


def variant_key(
    chrom,
    pos,
    ref,
    alt
):

    return (
        f"{chrom}:{pos}:{ref}:{alt}"
    )


# ============================================================
# EXTRACT COMPLETE REGION
# ============================================================

def extract_region(
    filepath,
    dataset
):

    print(
        f"\nReading {dataset}: "
        f"{filepath.name}"
    )

    records = []

    scanned = 0
    region_records = 0


    with gzip.open(
        filepath,
        "rt",
        encoding="utf-8",
        errors="replace"
    ) as handle:

        for line in handle:

            if line.startswith("#"):
                continue


            scanned += 1


            parts = (
                line.rstrip("\n")
                .split("\t")
            )


            if len(parts) < 10:
                continue


            chrom = clean_chrom(
                parts[0]
            )


            if chrom != CHROM:
                continue


            try:

                pos = int(
                    parts[1]
                )

            except Exception:

                continue


            if (
                pos < REGION_START
                or
                pos > REGION_END
            ):

                continue


            ref = parts[3]
            alt = parts[4]


            # Keep only biallelic variants
            if "," in alt:
                continue

            if (
                "<" in alt
                or
                ">" in alt
            ):
                continue


            fmt = parse_format(
                parts[8],
                parts[9]
            )


            beta = safe_float(
                fmt.get("ES")
            )

            se = safe_float(
                fmt.get("SE")
            )

            lp = safe_float(
                fmt.get("LP")
            )

            af = safe_float(
                fmt.get("AF")
            )

            maf = maf_from_af(
                af
            )


            key = variant_key(
                chrom,
                pos,
                ref,
                alt
            )


            records.append(
                {
                    "Variant_key": key,

                    "CHROM": chrom,
                    "POS": pos,

                    "rsID": parts[2],

                    "REF": ref,
                    "ALT": alt,

                    f"Beta_{dataset}":
                        beta,

                    f"SE_{dataset}":
                        se,

                    f"LP_{dataset}":
                        lp,

                    f"P_{dataset}":
                        lp_to_p(lp),

                    f"AF_{dataset}":
                        af,

                    f"MAF_{dataset}":
                        maf,
                }
            )


            region_records += 1


    print(
        f"  variants extracted: "
        f"{region_records:,}"
    )


    return pd.DataFrame(
        records
    )


# ============================================================
# LOAD ALL THREE DATASETS
# ============================================================

print("=" * 90)
print("EPB41L1 REGIONAL SEX-STRATIFIED ANALYSIS")
print("=" * 90)


combined = extract_region(
    COMBINED_FILE,
    "Combined"
)

male = extract_region(
    MALE_FILE,
    "Male"
)

female = extract_region(
    FEMALE_FILE,
    "Female"
)


# ============================================================
# MERGE SAME VARIANTS
# ============================================================

print("\nMerging datasets...")


variants = (
    combined
    .merge(
        male,
        on=[
            "Variant_key",
            "CHROM",
            "POS",
            "rsID",
            "REF",
            "ALT",
        ],
        how="inner"
    )
    .merge(
        female,
        on=[
            "Variant_key",
            "CHROM",
            "POS",
            "rsID",
            "REF",
            "ALT",
        ],
        how="inner"
    )
)


print(
    "Variants present in all three GWAS:",
    len(variants)
)


# ============================================================
# QC
# ============================================================

variants["Pass_QC"] = (

    variants["Beta_Male"].notna()
    &
    variants["Beta_Female"].notna()
    &
    variants["Beta_Combined"].notna()

    &
    variants["SE_Male"].notna()
    &
    variants["SE_Female"].notna()

    &
    (variants["SE_Male"] > 0)
    &
    (variants["SE_Female"] > 0)

    &
    variants["MAF_Male"].notna()
    &
    variants["MAF_Female"].notna()
    &
    variants["MAF_Combined"].notna()

    &
    (variants["MAF_Male"] >= MIN_MAF)
    &
    (variants["MAF_Female"] >= MIN_MAF)
    &
    (variants["MAF_Combined"] >= MIN_MAF)
)


qc = variants[
    variants["Pass_QC"]
].copy()


print(
    "Variants passing QC:",
    len(qc)
)


# ============================================================
# POSITION RELATIVE TO EPB41L1
# ============================================================

qc["Inside_EPB41L1"] = (

    (qc["POS"] >= GENE_START)
    &
    (qc["POS"] <= GENE_END)
)


def distance_to_gene(pos):

    if GENE_START <= pos <= GENE_END:
        return 0

    if pos < GENE_START:
        return GENE_START - pos

    return pos - GENE_END


qc["Distance_to_EPB41L1"] = (
    qc["POS"]
    .apply(distance_to_gene)
)


# ============================================================
# SEX HETEROGENEITY
# ============================================================

qc["Beta_difference_M_minus_F"] = (

    qc["Beta_Male"]
    -
    qc["Beta_Female"]
)


qc["SE_difference"] = np.sqrt(

    qc["SE_Male"] ** 2
    +
    qc["SE_Female"] ** 2
)


qc["Z_sex"] = (

    qc["Beta_difference_M_minus_F"]
    /
    qc["SE_difference"]
)


qc["P_sex"] = (

    2
    *
    norm.sf(
        np.abs(
            qc["Z_sex"]
        )
    )
)


# ============================================================
# MULTIPLE TESTING WITHIN REGION
# ============================================================

qc["P_sex_FDR_region"] = (
    multipletests(
        qc["P_sex"],
        method="fdr_bh"
    )[1]
)


qc["P_sex_Bonferroni_region"] = np.minimum(

    qc["P_sex"]
    *
    len(qc),

    1
)


# ============================================================
# EFFECT PATTERN
# ============================================================

def pattern(row):

    bm = row["Beta_Male"]
    bf = row["Beta_Female"]

    if bm > 0 and bf > 0:

        if bm > bf:
            return "positive_stronger_male"

        if bf > bm:
            return "positive_stronger_female"

        return "same_positive"


    if bm < 0 and bf < 0:

        if abs(bm) > abs(bf):
            return "negative_stronger_male"

        if abs(bf) > abs(bm):
            return "negative_stronger_female"

        return "same_negative"


    if bm > 0 and bf < 0:
        return "opposite_male_positive"


    if bm < 0 and bf > 0:
        return "opposite_female_positive"


    return "other"


qc["Sex_effect_pattern"] = (
    qc.apply(
        pattern,
        axis=1
    )
)


# ============================================================
# FLAG IMPORTANT VARIANTS
# ============================================================

qc["Is_rs532201406"] = (
    qc["rsID"] == "rs532201406"
)


qc["Is_rs1006296"] = (
    qc["rsID"] == "rs1006296"
)


# ============================================================
# SORT
# ============================================================

qc = qc.sort_values(
    "P_sex"
).reset_index(
    drop=True
)


qc["Sex_heterogeneity_rank"] = (
    np.arange(
        1,
        len(qc) + 1
    )
)


# ============================================================
# SUMMARY
# ============================================================

summary = pd.DataFrame(
    {
        "Metric": [

            "Region",

            "Variants present in all 3 GWAS",

            "Variants passing QC",

            "Variants inside EPB41L1",

            "Variants with P_sex < 0.05",

            "Variants with P_sex < 0.01",

            "Variants with P_sex < 0.001",

            "Variants with P_sex < 1e-4",

            "Variants with FDR_region < 0.05",

            "Variants with Bonferroni_region < 0.05",

            "Variants with beta_male > beta_female",

            "Same-positive variants stronger in males",
        ],

        "Value": [

            f"chr{CHROM}:{REGION_START}-{REGION_END}",

            len(variants),

            len(qc),

            int(
                qc[
                    "Inside_EPB41L1"
                ].sum()
            ),

            int(
                (
                    qc["P_sex"] < 0.05
                ).sum()
            ),

            int(
                (
                    qc["P_sex"] < 0.01
                ).sum()
            ),

            int(
                (
                    qc["P_sex"] < 0.001
                ).sum()
            ),

            int(
                (
                    qc["P_sex"] < 1e-4
                ).sum()
            ),

            int(
                (
                    qc[
                        "P_sex_FDR_region"
                    ] < 0.05
                ).sum()
            ),

            int(
                (
                    qc[
                        "P_sex_Bonferroni_region"
                    ] < 0.05
                ).sum()
            ),

            int(
                (
                    qc[
                        "Beta_difference_M_minus_F"
                    ] > 0
                ).sum()
            ),

            int(
                (
                    qc[
                        "Sex_effect_pattern"
                    ]
                    ==
                    "positive_stronger_male"
                ).sum()
            ),
        ]
    }
)


# ============================================================
# TOP SEX-DIFFERENTIAL VARIANTS
# ============================================================

top_sex = (
    qc
    .sort_values(
        "P_sex"
    )
    .head(100)
    .copy()
)


# ============================================================
# TOP COMBINED GWAS VARIANTS
# ============================================================

top_combined = (
    qc
    .sort_values(
        "LP_Combined",
        ascending=False
    )
    .head(100)
    .copy()
)


# ============================================================
# INTRAGENIC EPB41L1 VARIANTS
# ============================================================

inside_gene = (
    qc[
        qc[
            "Inside_EPB41L1"
        ]
    ]
    .sort_values(
        "P_sex"
    )
    .copy()
)


# ============================================================
# KEY VARIANTS
# ============================================================

key_variants = (
    qc[
        qc["rsID"].isin(
            [
                "rs532201406",
                "rs1006296"
            ]
        )
    ]
    .copy()
)


# ============================================================
# PRINT
# ============================================================

print("\n" + "=" * 90)
print("FINAL SUMMARY")
print("=" * 90)

print(
    summary.to_string(
        index=False
    )
)


print("\n" + "=" * 90)
print("KEY VARIANTS")
print("=" * 90)


cols = [

    "rsID",
    "Variant_key",
    "POS",

    "Inside_EPB41L1",
    "Distance_to_EPB41L1",

    "Beta_Male",
    "SE_Male",
    "P_Male",

    "Beta_Female",
    "SE_Female",
    "P_Female",

    "Beta_Combined",
    "P_Combined",

    "Beta_difference_M_minus_F",

    "Z_sex",
    "P_sex",

    "P_sex_FDR_region",
    "P_sex_Bonferroni_region",

    "Sex_effect_pattern",
]


print(
    key_variants[
        cols
    ]
    .to_string(
        index=False
    )
)


print("\n" + "=" * 90)
print("TOP 20 SEX-HETEROGENEITY VARIANTS")
print("=" * 90)

print(
    top_sex[
        cols
    ]
    .head(20)
    .to_string(
        index=False
    )
)


# ============================================================
# EXPORT
# ============================================================

qc.to_csv(
    OUTPUT_TSV,
    sep="\t",
    index=False,
    compression="gzip"
)


with pd.ExcelWriter(
    OUTPUT_XLSX,
    engine="openpyxl"
) as writer:

    summary.to_excel(
        writer,
        sheet_name="S0_Summary",
        index=False
    )

    key_variants.to_excel(
        writer,
        sheet_name="S1_Key_variants",
        index=False
    )

    top_sex.to_excel(
        writer,
        sheet_name="S2_Top_sex_heterogeneity",
        index=False
    )

    top_combined.to_excel(
        writer,
        sheet_name="S3_Top_combined_GWAS",
        index=False
    )

    inside_gene.to_excel(
        writer,
        sheet_name="S4_Inside_EPB41L1",
        index=False
    )

    qc.to_excel(
        writer,
        sheet_name="S5_All_region_variants",
        index=False
    )


print("\nSaved:")
print(OUTPUT_XLSX.resolve())
print(OUTPUT_TSV.resolve())

print("\nDONE.")