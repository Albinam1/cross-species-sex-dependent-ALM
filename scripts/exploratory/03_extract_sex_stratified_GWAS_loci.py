from pathlib import Path
import gzip
import math

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

OUTPUT_FILE = Path(
    "03_SEX_STRATIFIED_GWAS_LOCI.xlsx"
)


# ============================================================
# GENOME BUILD
# ============================================================
#
# Pei2020 / OpenGWAS files:
# GRCh37 / hg19
#
# Coordinates below are HUMAN gene coordinates in GRCh37.
#
# We analyze gene boundaries +/- 500 kb.
#
# This is deliberately a LOCUS analysis, not "nearest gene"
# annotation.
# ============================================================

WINDOW_BP = 500_000


GENES = {

    "CLPP": {
        "chr": "19",
        "start": 6_361_463,
        "end": 6_368_919,
        "mouse_gene": "Clpp",
        "mouse_pattern": "sex_opposite",
        "male_KO": "increased",
        "female_KO": "decreased",
    },

    "MTA1": {
        "chr": "14",
        "start": 105_886_159,
        "end": 105_937_066,
        "mouse_gene": "Mta1",
        "mouse_pattern": "sex_opposite",
        "male_KO": "decreased",
        "female_KO": "increased",
    },

    "RSPO1": {
        "chr": "1",
        "start": 38_076_951,
        "end": 38_100_595,
        "mouse_gene": "Rspo1",
        "mouse_pattern": "sex_opposite",
        "male_KO": "increased",
        "female_KO": "decreased",
    },

    "SNAP47": {
        "chr": "1",
        "start": 227_915_869,
        "end": 227_968_927,
        "mouse_gene": "Snap47",
        "mouse_pattern": "sex_opposite",
        "male_KO": "increased",
        "female_KO": "decreased",
    },

    "LINGO2": {
        "chr": "9",
        "start": 27_948_076,
        "end": 28_670_283,
        "mouse_gene": "Lingo2",
        "mouse_pattern": "sex_conserved",
        "male_KO": "increased",
        "female_KO": "increased",
    },
}


# ============================================================
# CHECK FILES
# ============================================================

print("=" * 78)
print("SEX-STRATIFIED APPENDICULAR LEAN-MASS GWAS ANALYSIS")
print("=" * 78)

for f in [
    MALE_FILE,
    FEMALE_FILE,
    COMBINED_FILE
]:

    if not f.exists():

        raise FileNotFoundError(
            f"\nFile not found:\n{f.resolve()}"
        )

    print(
        f"{f.name}: "
        f"{f.stat().st_size / 1024**3:.2f} GB"
    )


# ============================================================
# PREPARE LOCI
# ============================================================

loci = []

for gene, x in GENES.items():

    locus_start = max(
        1,
        x["start"] - WINDOW_BP
    )

    locus_end = (
        x["end"] + WINDOW_BP
    )

    loci.append({
        "Gene": gene,
        "Mouse_gene": x["mouse_gene"],
        "Chromosome": x["chr"],
        "Gene_start_GRCh37": x["start"],
        "Gene_end_GRCh37": x["end"],
        "Locus_start": locus_start,
        "Locus_end": locus_end,
        "Window_bp": WINDOW_BP,
        "Mouse_sex_pattern":
            x["mouse_pattern"],
        "Male_KO_direction":
            x["male_KO"],
        "Female_KO_direction":
            x["female_KO"],
    })


loci_df = pd.DataFrame(loci)


print("\nCandidate loci:")

print(
    loci_df[
        [
            "Gene",
            "Chromosome",
            "Locus_start",
            "Locus_end",
            "Mouse_sex_pattern"
        ]
    ].to_string(index=False)
)


# ============================================================
# FAST LOOKUP BY CHROMOSOME
# ============================================================

loci_by_chr = {}

for row in loci:

    chrom = row["Chromosome"]

    loci_by_chr.setdefault(
        chrom,
        []
    ).append(row)


# ============================================================
# PARSE FORMAT FIELD
# ============================================================

def parse_format(
    format_string,
    sample_string
):

    keys = format_string.split(":")
    values = sample_string.split(":")

    return dict(
        zip(keys, values)
    )


def safe_float(x):

    if x is None:
        return np.nan

    if x in {
        "",
        ".",
        "NA",
        "nan"
    }:
        return np.nan

    try:
        return float(x)

    except Exception:
        return np.nan


# ============================================================
# LP -> P VALUE
# ============================================================

def lp_to_p(lp):

    """
    OpenGWAS GWAS-VCF commonly stores
    LP = -log10(P).
    """

    if pd.isna(lp):
        return np.nan

    try:

        if lp > 323:
            # smaller than normal Python float range
            return 0.0

        return 10 ** (-lp)

    except Exception:
        return np.nan


# ============================================================
# EXTRACT ONE VCF
# ============================================================

def extract_vcf(
    filepath,
    dataset_name
):

    print("\n" + "=" * 78)
    print(f"Reading {dataset_name}: {filepath.name}")
    print("=" * 78)

    records = []

    n_total = 0
    n_extracted = 0

    sample_column_name = None

    with gzip.open(
        filepath,
        "rt",
        encoding="utf-8",
        errors="replace"
    ) as handle:

        for line in handle:

            # ------------------------------------------------
            # metadata
            # ------------------------------------------------

            if line.startswith("##"):

                continue


            # ------------------------------------------------
            # VCF header
            # ------------------------------------------------

            if line.startswith("#CHROM"):

                header = (
                    line.rstrip("\n")
                    .split("\t")
                )

                if len(header) < 10:

                    raise ValueError(
                        f"{filepath.name}: "
                        "VCF does not contain a sample column."
                    )

                sample_column_name = header[9]

                print(
                    "VCF sample column:",
                    sample_column_name
                )

                continue


            if line.startswith("#"):
                continue


            n_total += 1


            # ------------------------------------------------
            # progress
            # ------------------------------------------------

            if n_total % 2_000_000 == 0:

                print(
                    f"  scanned "
                    f"{n_total:,} variants..."
                )


            parts = (
                line.rstrip("\n")
                .split("\t")
            )

            if len(parts) < 10:
                continue


            chrom = (
                parts[0]
                .replace("chr", "")
            )

            # skip chromosomes with no loci
            if chrom not in loci_by_chr:
                continue


            try:
                pos = int(parts[1])

            except Exception:
                continue


            # ------------------------------------------------
            # check locus
            # ------------------------------------------------

            matching_loci = []

            for locus in loci_by_chr[chrom]:

                if (
                    locus["Locus_start"]
                    <= pos
                    <= locus["Locus_end"]
                ):

                    matching_loci.append(
                        locus
                    )


            if not matching_loci:
                continue


            # ------------------------------------------------
            # variant fields
            # ------------------------------------------------

            rsid = parts[2]

            ref = parts[3]
            alt = parts[4]

            format_string = parts[8]
            sample_string = parts[9]


            fmt = parse_format(
                format_string,
                sample_string
            )


            # GWAS-VCF fields
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

            internal_id = fmt.get(
                "ID",
                None
            )


            p_value = lp_to_p(lp)


            # ------------------------------------------------
            # one variant can theoretically fall in overlapping
            # candidate windows
            # ------------------------------------------------

            for locus in matching_loci:

                records.append({
                    "Dataset":
                        dataset_name,

                    "Gene":
                        locus["Gene"],

                    "Mouse_gene":
                        locus["Mouse_gene"],

                    "Mouse_sex_pattern":
                        locus[
                            "Mouse_sex_pattern"
                        ],

                    "Male_KO_direction":
                        locus[
                            "Male_KO_direction"
                        ],

                    "Female_KO_direction":
                        locus[
                            "Female_KO_direction"
                        ],

                    "CHROM":
                        chrom,

                    "POS":
                        pos,

                    "rsID":
                        rsid,

                    "REF":
                        ref,

                    "ALT":
                        alt,

                    # GWAS-VCF ES is interpreted
                    # relative to ALT
                    "Effect_allele":
                        alt,

                    "Other_allele":
                        ref,

                    "BETA_ALT":
                        beta,

                    "SE":
                        se,

                    "LP":
                        lp,

                    "P":
                        p_value,

                    "AF_ALT":
                        af,

                    "VCF_internal_ID":
                        internal_id,

                    "Distance_to_gene":
                        (
                            locus[
                                "Gene_start_GRCh37"
                            ] - pos
                            if pos <
                            locus[
                                "Gene_start_GRCh37"
                            ]

                            else

                            pos -
                            locus[
                                "Gene_end_GRCh37"
                            ]
                            if pos >
                            locus[
                                "Gene_end_GRCh37"
                            ]

                            else 0
                        ),

                    "Inside_gene":
                        (
                            locus[
                                "Gene_start_GRCh37"
                            ]
                            <= pos
                            <= locus[
                                "Gene_end_GRCh37"
                            ]
                        )
                })

                n_extracted += 1


    out = pd.DataFrame(records)

    print(
        f"\nScanned variants: "
        f"{n_total:,}"
    )

    print(
        f"Extracted locus records: "
        f"{n_extracted:,}"
    )


    if not out.empty:

        print("\nVariants per locus:")

        print(
            out.groupby("Gene")
            .size()
            .to_string()
        )


    return out


# ============================================================
# EXTRACT ALL THREE GWAS
# ============================================================

male = extract_vcf(
    MALE_FILE,
    "Male"
)

female = extract_vcf(
    FEMALE_FILE,
    "Female"
)

combined = extract_vcf(
    COMBINED_FILE,
    "Combined"
)


# ============================================================
# QC
# ============================================================

print("\n" + "=" * 78)
print("VCF QC")
print("=" * 78)


def qc_table(df, label):

    return {
        "Dataset": label,

        "N_records":
            len(df),

        "N_unique_variants":
            df[
                [
                    "CHROM",
                    "POS",
                    "REF",
                    "ALT"
                ]
            ]
            .drop_duplicates()
            .shape[0],

        "Missing_beta":
            df["BETA_ALT"]
            .isna()
            .sum(),

        "Missing_SE":
            df["SE"]
            .isna()
            .sum(),

        "Missing_P":
            df["P"]
            .isna()
            .sum(),
    }


qc = pd.DataFrame([
    qc_table(
        male,
        "Male"
    ),

    qc_table(
        female,
        "Female"
    ),

    qc_table(
        combined,
        "Combined"
    )
])

print(qc.to_string(index=False))


# ============================================================
# CREATE UNIQUE VARIANT KEY
# ============================================================

def add_variant_key(df):

    df = df.copy()

    df["Variant_key"] = (
        df["CHROM"].astype(str)
        + ":"
        + df["POS"].astype(str)
        + ":"
        + df["REF"].astype(str)
        + ":"
        + df["ALT"].astype(str)
    )

    return df


male = add_variant_key(male)
female = add_variant_key(female)
combined = add_variant_key(combined)


# ============================================================
# PREP FOR MERGE
# ============================================================

base_cols = [
    "Gene",
    "Mouse_gene",
    "Mouse_sex_pattern",
    "Male_KO_direction",
    "Female_KO_direction",
    "CHROM",
    "POS",
    "rsID",
    "REF",
    "ALT",
    "Effect_allele",
    "Other_allele",
    "Inside_gene",
    "Distance_to_gene",
    "Variant_key"
]


def rename_dataset(
    df,
    suffix
):

    keep = (
        base_cols
        + [
            "BETA_ALT",
            "SE",
            "P",
            "LP",
            "AF_ALT"
        ]
    )

    x = df[
        [c for c in keep if c in df.columns]
    ].copy()


    rename = {
        "BETA_ALT":
            f"Beta_{suffix}",

        "SE":
            f"SE_{suffix}",

        "P":
            f"P_{suffix}",

        "LP":
            f"LP_{suffix}",

        "AF_ALT":
            f"AF_{suffix}",
    }

    return x.rename(
        columns=rename
    )


male_m = rename_dataset(
    male,
    "Male"
)

female_m = rename_dataset(
    female,
    "Female"
)

combined_m = rename_dataset(
    combined,
    "Combined"
)


# ============================================================
# MERGE MALE + FEMALE
# ============================================================

merge_keys = [
    "Gene",
    "Mouse_gene",
    "Mouse_sex_pattern",
    "Male_KO_direction",
    "Female_KO_direction",
    "CHROM",
    "POS",
    "rsID",
    "REF",
    "ALT",
    "Effect_allele",
    "Other_allele",
    "Inside_gene",
    "Distance_to_gene",
    "Variant_key"
]


sex = male_m.merge(
    female_m,
    on=merge_keys,
    how="outer"
)


# ============================================================
# ADD COMBINED
# ============================================================

sex = sex.merge(
    combined_m,
    on=merge_keys,
    how="outer"
)


# ============================================================
# FORMAL SEX HETEROGENEITY TEST
# ============================================================
#
# Male and female samples are non-overlapping,
# therefore:
#
# Z = (beta_male - beta_female)
#     / sqrt(SE_male^2 + SE_female^2)
#
# ============================================================

valid = (
    sex["Beta_Male"].notna()
    &
    sex["Beta_Female"].notna()
    &
    sex["SE_Male"].notna()
    &
    sex["SE_Female"].notna()
    &
    (sex["SE_Male"] > 0)
    &
    (sex["SE_Female"] > 0)
)


sex["Z_sex_difference"] = np.nan
sex["P_sex_difference"] = np.nan


sex.loc[
    valid,
    "Z_sex_difference"
] = (

    (
        sex.loc[
            valid,
            "Beta_Male"
        ]

        -

        sex.loc[
            valid,
            "Beta_Female"
        ]
    )

    /

    np.sqrt(
        sex.loc[
            valid,
            "SE_Male"
        ] ** 2

        +

        sex.loc[
            valid,
            "SE_Female"
        ] ** 2
    )
)


sex.loc[
    valid,
    "P_sex_difference"
] = (

    2
    *
    norm.sf(
        np.abs(
            sex.loc[
                valid,
                "Z_sex_difference"
            ]
        )
    )
)


# ============================================================
# MULTIPLE TESTING FOR SEX DIFFERENCE
# ============================================================

sex[
    "FDR_sex_difference"
] = np.nan


valid_p = (
    sex["P_sex_difference"]
    .notna()
)


if valid_p.sum() > 0:

    sex.loc[
        valid_p,
        "FDR_sex_difference"
    ] = multipletests(
        sex.loc[
            valid_p,
            "P_sex_difference"
        ],
        method="fdr_bh"
    )[1]


# ============================================================
# EFFECT DIRECTION CLASSIFICATION
# ============================================================

def effect_direction(beta):

    if pd.isna(beta):
        return "missing"

    if beta > 0:
        return "increased_lean_mass"

    if beta < 0:
        return "decreased_lean_mass"

    return "zero"


sex["Male_GWAS_direction"] = (
    sex["Beta_Male"]
    .apply(effect_direction)
)

sex["Female_GWAS_direction"] = (
    sex["Beta_Female"]
    .apply(effect_direction)
)


# ============================================================
# SEX EFFECT PATTERN
# ============================================================

def classify_human_pattern(row):

    bm = row["Beta_Male"]
    bf = row["Beta_Female"]

    if pd.isna(bm) or pd.isna(bf):
        return "incomplete"

    if bm > 0 and bf > 0:
        return "same_positive"

    if bm < 0 and bf < 0:
        return "same_negative"

    if bm > 0 and bf < 0:
        return "opposite_male_positive"

    if bm < 0 and bf > 0:
        return "opposite_female_positive"

    return "other"


sex["Human_effect_pattern"] = (
    sex.apply(
        classify_human_pattern,
        axis=1
    )
)


# ============================================================
# FLAG POTENTIALLY INTERESTING VARIANTS
# ============================================================

sex[
    "Nominal_sex_heterogeneity"
] = (
    sex["P_sex_difference"]
    < 0.05
)

sex[
    "FDR_sex_heterogeneity"
] = (
    sex["FDR_sex_difference"]
    < 0.05
)

sex[
    "Opposite_GWAS_direction"
] = (
    sex[
        "Human_effect_pattern"
    ].isin([
        "opposite_male_positive",
        "opposite_female_positive"
    ])
)


# ============================================================
# LEAD SNP PER GENE / SEX
# ============================================================

def get_lead(
    df,
    p_col,
    beta_col,
    label
):

    x = df[
        df[p_col].notna()
    ].copy()

    if x.empty:

        return pd.DataFrame()


    x = (
        x.sort_values(
            [
                "Gene",
                p_col
            ]
        )
        .groupby(
            "Gene",
            as_index=False
        )
        .first()
    )


    columns = [
        "Gene",
        "rsID",
        "Variant_key",
        "Effect_allele",
        "Other_allele",
        beta_col,
        p_col
    ]


    x = x[columns]


    x = x.rename(
        columns={
            "rsID":
                f"Lead_rsID_{label}",

            "Variant_key":
                f"Lead_variant_{label}",

            "Effect_allele":
                f"Effect_allele_{label}",

            "Other_allele":
                f"Other_allele_{label}",

            beta_col:
                f"Lead_beta_{label}",

            p_col:
                f"Lead_P_{label}"
        }
    )


    return x


lead_male = get_lead(
    sex,
    "P_Male",
    "Beta_Male",
    "Male"
)

lead_female = get_lead(
    sex,
    "P_Female",
    "Beta_Female",
    "Female"
)

lead_combined = get_lead(
    sex,
    "P_Combined",
    "Beta_Combined",
    "Combined"
)


gene_summary = (
    loci_df[
        [
            "Gene",
            "Mouse_gene",
            "Mouse_sex_pattern",
            "Male_KO_direction",
            "Female_KO_direction"
        ]
    ]
    .merge(
        lead_male,
        on="Gene",
        how="left"
    )
    .merge(
        lead_female,
        on="Gene",
        how="left"
    )
    .merge(
        lead_combined,
        on="Gene",
        how="left"
    )
)


# ============================================================
# LOCUS COUNTS
# ============================================================

counts = (
    sex.groupby("Gene")
    .agg(
        N_variants=(
            "Variant_key",
            "nunique"
        ),

        N_male=(
            "Beta_Male",
            lambda x:
                x.notna().sum()
        ),

        N_female=(
            "Beta_Female",
            lambda x:
                x.notna().sum()
        ),

        N_combined=(
            "Beta_Combined",
            lambda x:
                x.notna().sum()
        ),

        N_opposite_direction=(
            "Opposite_GWAS_direction",
            "sum"
        ),

        N_nominal_sex_diff=(
            "Nominal_sex_heterogeneity",
            "sum"
        ),

        N_FDR_sex_diff=(
            "FDR_sex_heterogeneity",
            "sum"
        )
    )
    .reset_index()
)


gene_summary = gene_summary.merge(
    counts,
    on="Gene",
    how="left"
)


# ============================================================
# MOST SEX-DIFFERENT VARIANT PER GENE
# ============================================================

het = (
    sex[
        sex["P_sex_difference"]
        .notna()
    ]
    .sort_values(
        [
            "Gene",
            "P_sex_difference"
        ]
    )
    .groupby(
        "Gene",
        as_index=False
    )
    .first()
)


if not het.empty:

    het = het[
        [
            "Gene",
            "rsID",
            "Variant_key",
            "Beta_Male",
            "SE_Male",
            "P_Male",
            "Beta_Female",
            "SE_Female",
            "P_Female",
            "Z_sex_difference",
            "P_sex_difference",
            "FDR_sex_difference",
            "Human_effect_pattern"
        ]
    ]

    het = het.rename(
        columns={
            "rsID":
                "Top_heterogeneity_rsID",

            "Variant_key":
                "Top_heterogeneity_variant"
        }
    )

    gene_summary = (
        gene_summary.merge(
            het,
            on="Gene",
            how="left"
        )
    )


# ============================================================
# SUMMARY
# ============================================================

summary = pd.DataFrame({
    "Metric": [

        "Candidate genes",

        "Male GWAS locus records",

        "Female GWAS locus records",

        "Combined GWAS locus records",

        "Merged unique locus variants",

        "Variants testable for male-female heterogeneity",

        "Variants with opposite beta signs",

        "Variants with nominal sex heterogeneity P<0.05",

        "Variants with FDR sex heterogeneity <0.05"
    ],

    "N": [

        len(GENES),

        len(male),

        len(female),

        len(combined),

        sex["Variant_key"].nunique(),

        valid.sum(),

        sex[
            "Opposite_GWAS_direction"
        ].sum(),

        sex[
            "Nominal_sex_heterogeneity"
        ].sum(),

        sex[
            "FDR_sex_heterogeneity"
        ].sum()
    ]
})


print("\n" + "=" * 78)
print("FINAL SUMMARY")
print("=" * 78)

print(
    summary.to_string(
        index=False
    )
)


print("\nGENE-LEVEL SUMMARY")
print("-" * 78)

print(
    gene_summary.to_string(
        index=False
    )
)


# ============================================================
# HIGH-INTEREST VARIANTS
# ============================================================

priority_variants = sex[
    (
        sex[
            "Opposite_GWAS_direction"
        ]
    )
    |
    (
        sex[
            "P_sex_difference"
        ]
        < 0.05
    )
].copy()


priority_variants = (
    priority_variants
    .sort_values(
        [
            "Gene",
            "P_sex_difference"
        ],
        na_position="last"
    )
)


# ============================================================
# SAVE EXCEL
# ============================================================

print("\n" + "=" * 78)
print("Writing output workbook...")
print("=" * 78)


with pd.ExcelWriter(
    OUTPUT_FILE,
    engine="openpyxl"
) as writer:

    summary.to_excel(
        writer,
        sheet_name="S0_Summary",
        index=False
    )

    loci_df.to_excel(
        writer,
        sheet_name="S1_Loci_GRCh37",
        index=False
    )

    gene_summary.to_excel(
        writer,
        sheet_name="S2_Gene_summary",
        index=False
    )

    sex.to_excel(
        writer,
        sheet_name="S3_All_harmonized",
        index=False
    )

    priority_variants.to_excel(
        writer,
        sheet_name="S4_Priority_variants",
        index=False
    )

    male.to_excel(
        writer,
        sheet_name="S5_Male_raw_loci",
        index=False
    )

    female.to_excel(
        writer,
        sheet_name="S6_Female_raw_loci",
        index=False
    )

    combined.to_excel(
        writer,
        sheet_name="S7_Combined_raw_loci",
        index=False
    )

    qc.to_excel(
        writer,
        sheet_name="S8_QC",
        index=False
    )


print(
    "\nSaved:",
    OUTPUT_FILE.resolve()
)

print("\nDONE.")