from pathlib import Path
import tarfile
import io
import re
import pandas as pd
import numpy as np


# ============================================================
# SETTINGS
# ============================================================

ARCHIVE_FILE = Path(
    "GTEx_Analysis_v8_sbeQTLs.tar.gz"
)

OUTPUT_FILE = Path(
    "08_EPB41L1_GTEx_sex_biased_eQTL.xlsx"
)

GENE_SYMBOL = "EPB41L1"

GENE_ENSEMBL = "ENSG00000088367"

KEY_GWAS_SNPS = {
    "rs6060632",
    "rs1006296",
    "rs113092336",
    "rs532201406",
}


# ============================================================
# CHECK INPUT
# ============================================================

if not ARCHIVE_FILE.exists():
    raise FileNotFoundError(
        f"Archive not found:\n{ARCHIVE_FILE.resolve()}"
    )


print("=" * 90)
print("GTEx v8 SEX-BIASED eQTL ANALYSIS: EPB41L1")
print("=" * 90)

print(
    "Archive:",
    ARCHIVE_FILE.resolve()
)


# ============================================================
# OPEN ARCHIVE
# ============================================================

tar = tarfile.open(
    ARCHIVE_FILE,
    mode="r:gz"
)

members = [
    m
    for m in tar.getmembers()
    if m.isfile()
]


print(
    "\nFiles inside archive:",
    len(members)
)


for m in members:
    print(" -", m.name)


# ============================================================
# READ README FILES
# ============================================================

print("\n" + "=" * 90)
print("README CONTENT")
print("=" * 90)

readme_texts = []

for member in members:

    name_lower = (
        member.name
        .lower()
    )

    if "readme" in name_lower:

        f = tar.extractfile(
            member
        )

        if f is None:
            continue

        raw = f.read()

        try:
            text = raw.decode(
                "utf-8"
            )

        except Exception:
            text = raw.decode(
                "latin1",
                errors="replace"
            )

        readme_texts.append(
            {
                "file": member.name,
                "text": text
            }
        )

        print(
            f"\n--- {member.name} ---\n"
        )

        print(
            text[:12000]
        )


# ============================================================
# IDENTIFY DATA FILES
# ============================================================

data_members = []

for member in members:

    name = member.name.lower()

    if "readme" in name:
        continue

    if (
        name.endswith(".txt")
        or
        name.endswith(".tsv")
        or
        name.endswith(".csv")
        or
        name.endswith(".txt.gz")
        or
        name.endswith(".tsv.gz")
        or
        name.endswith(".csv.gz")
    ):

        data_members.append(
            member
        )


print("\n" + "=" * 90)
print("CANDIDATE DATA FILES")
print("=" * 90)

for member in data_members:
    print(
        member.name
    )


# ============================================================
# HELPERS
# ============================================================

def normalize_gene_id(x):

    if pd.isna(x):
        return ""

    return (
        str(x)
        .strip()
        .split(".")[0]
    )


def detect_separator_from_name(
    name
):

    name = name.lower()

    if ".csv" in name:
        return ","

    return "\t"


def read_member_table(
    tar_obj,
    member
):

    f = tar_obj.extractfile(
        member
    )

    if f is None:
        return None

    raw = f.read()


    # --------------------------------------------------------
    # decompression for gz file stored INSIDE tar
    # --------------------------------------------------------

    if member.name.lower().endswith(
        ".gz"
    ):

        import gzip

        raw = gzip.decompress(
            raw
        )


    sep = detect_separator_from_name(
        member.name
    )


    try:

        df = pd.read_csv(
            io.BytesIO(raw),
            sep=sep,
            dtype=str,
            low_memory=False
        )

    except Exception:

        # fallback auto detection
        df = pd.read_csv(
            io.BytesIO(raw),
            sep=None,
            engine="python",
            dtype=str
        )


    return df


def find_column(
    columns,
    candidate_patterns
):

    cols = list(
        columns
    )


    # exact case-insensitive match first
    lower_map = {
        str(c).lower(): c
        for c in cols
    }


    for candidate in candidate_patterns:

        if candidate.lower() in lower_map:

            return lower_map[
                candidate.lower()
            ]


    # substring match
    for col in cols:

        low = str(col).lower()

        for candidate in candidate_patterns:

            if (
                candidate.lower()
                in low
            ):

                return col


    return None


# ============================================================
# SCAN FILES FOR EPB41L1
# ============================================================

print("\n" + "=" * 90)
print("SCANNING sb-eQTL FILES FOR EPB41L1")
print("=" * 90)


hits = []

file_summary = []


for i, member in enumerate(
    data_members,
    start=1
):

    print(
        f"\n[{i}/{len(data_members)}] "
        f"{member.name}"
    )


    try:

        df = read_member_table(
            tar,
            member
        )

    except Exception as e:

        print(
            "  Could not read:",
            repr(e)
        )

        continue


    if (
        df is None
        or
        df.empty
    ):

        print(
            "  Empty table."
        )

        continue


    print(
        "  Rows:",
        len(df)
    )

    print(
        "  Columns:",
        list(df.columns)
    )


    # --------------------------------------------------------
    # find gene column
    # --------------------------------------------------------

    gene_col = find_column(
        df.columns,
        [
            "gene_id",
            "gene",
            "gencode_id",
            "phenotype_id",
            "gene_symbol",
            "symbol"
        ]
    )


    variant_col = find_column(
        df.columns,
        [
            "variant_id",
            "variant",
            "snp",
            "snp_id",
            "rsid",
            "rs_id"
        ]
    )


    p_col = find_column(
        df.columns,
        [
            "pval",
            "p_value",
            "pvalue",
            "p_nominal",
            "interaction_p",
            "sex_p"
        ]
    )


    slope_col = find_column(
        df.columns,
        [
            "slope",
            "beta",
            "effect",
            "interaction_effect",
            "sex_effect"
        ]
    )


    fdr_col = find_column(
        df.columns,
        [
            "qval",
            "q_value",
            "fdr",
            "adj_p",
            "padj"
        ]
    )


    # --------------------------------------------------------
    # detect tissue from filename
    # --------------------------------------------------------

    tissue = (
        Path(
            member.name
        )
        .name
    )


    tissue = re.sub(
        r"\.(txt|tsv|csv)(\.gz)?$",
        "",
        tissue,
        flags=re.IGNORECASE
    )


    # --------------------------------------------------------
    # locate EPB41L1
    # --------------------------------------------------------

    mask = pd.Series(
        False,
        index=df.index
    )


    if gene_col is not None:

        gene_values = (
            df[gene_col]
            .astype(str)
            .str.strip()
        )


        normalized = (
            gene_values
            .apply(
                normalize_gene_id
            )
        )


        mask = (
            gene_values
            .str.upper()
            .eq(
                GENE_SYMBOL
            )
            |
            normalized.eq(
                GENE_ENSEMBL
            )
        )


    # fallback: scan all columns for symbol/id
    if not mask.any():

        for col in df.columns:

            vals = (
                df[col]
                .astype(str)
                .str.strip()
            )


            normalized = (
                vals
                .apply(
                    normalize_gene_id
                )
            )


            col_mask = (
                vals
                .str.upper()
                .eq(
                    GENE_SYMBOL
                )
                |
                normalized.eq(
                    GENE_ENSEMBL
                )
            )


            if col_mask.any():

                mask = (
                    mask
                    |
                    col_mask
                )


    n_hits = int(
        mask.sum()
    )


    file_summary.append(
        {
            "File": member.name,
            "Tissue_guess": tissue,
            "Rows": len(df),
            "Gene_column": gene_col,
            "Variant_column": variant_col,
            "P_column": p_col,
            "Effect_column": slope_col,
            "FDR_column": fdr_col,
            "EPB41L1_hits": n_hits,
        }
    )


    if n_hits == 0:

        print(
            "  EPB41L1: no records"
        )

        continue


    print(
        f"  EPB41L1 records: "
        f"{n_hits}"
    )


    sub = df.loc[
        mask
    ].copy()


    sub[
        "Source_file"
    ] = member.name

    sub[
        "Tissue_guess"
    ] = tissue


    # standard fields
    sub[
        "_gene_col"
    ] = gene_col

    sub[
        "_variant_col"
    ] = variant_col

    sub[
        "_p_col"
    ] = p_col

    sub[
        "_effect_col"
    ] = slope_col

    sub[
        "_fdr_col"
    ] = fdr_col


    hits.append(
        sub
    )


# ============================================================
# COMBINE EPB41L1 RESULTS
# ============================================================

if hits:

    all_hits = pd.concat(
        hits,
        ignore_index=True,
        sort=False
    )

else:

    all_hits = pd.DataFrame()


file_summary_df = pd.DataFrame(
    file_summary
)


print("\n" + "=" * 90)
print("EPB41L1 SUMMARY")
print("=" * 90)


print(
    "Files containing EPB41L1:",
    int(
        (
            file_summary_df[
                "EPB41L1_hits"
            ] > 0
        ).sum()
    )
)


print(
    "Total EPB41L1 sb-eQTL records:",
    len(all_hits)
)


# ============================================================
# STANDARDIZE EPB41L1 TABLE
# ============================================================

standard_rows = []


if not all_hits.empty:

    for _, row in all_hits.iterrows():

        variant_col = row.get(
            "_variant_col"
        )

        p_col = row.get(
            "_p_col"
        )

        effect_col = row.get(
            "_effect_col"
        )

        fdr_col = row.get(
            "_fdr_col"
        )


        variant = (
            row.get(
                variant_col
            )
            if (
                variant_col
                and
                variant_col in row.index
            )
            else np.nan
        )


        pval = (
            row.get(
                p_col
            )
            if (
                p_col
                and
                p_col in row.index
            )
            else np.nan
        )


        effect = (
            row.get(
                effect_col
            )
            if (
                effect_col
                and
                effect_col in row.index
            )
            else np.nan
        )


        fdr = (
            row.get(
                fdr_col
            )
            if (
                fdr_col
                and
                fdr_col in row.index
            )
            else np.nan
        )


        standard_rows.append(
            {
                "Gene":
                    GENE_SYMBOL,

                "Gene_Ensembl":
                    GENE_ENSEMBL,

                "Tissue":
                    row.get(
                        "Tissue_guess"
                    ),

                "Source_file":
                    row.get(
                        "Source_file"
                    ),

                "Variant":
                    variant,

                "Effect":
                    pd.to_numeric(
                        effect,
                        errors="coerce"
                    ),

                "P_value":
                    pd.to_numeric(
                        pval,
                        errors="coerce"
                    ),

                "FDR_or_qvalue":
                    pd.to_numeric(
                        fdr,
                        errors="coerce"
                    ),
            }
        )


standard = pd.DataFrame(
    standard_rows
)


# ============================================================
# TRY TO EXTRACT rsID
# ============================================================

def extract_rsid(x):

    if pd.isna(x):
        return np.nan


    text = str(x)


    m = re.search(
        r"(rs\d+)",
        text,
        flags=re.IGNORECASE
    )


    if m:
        return m.group(1)


    return np.nan


if not standard.empty:

    standard[
        "rsID"
    ] = (
        standard[
            "Variant"
        ]
        .apply(
            extract_rsid
        )
    )


    standard[
        "Key_GWAS_variant"
    ] = (
        standard[
            "rsID"
        ]
        .isin(
            KEY_GWAS_SNPS
        )
    )


# ============================================================
# SIGNIFICANCE CLASSIFICATION
# ============================================================

if not standard.empty:

    standard[
        "Significant_FDR_0.05"
    ] = (
        standard[
            "FDR_or_qvalue"
        ].notna()
        &
        (
            standard[
                "FDR_or_qvalue"
            ] < 0.05
        )
    )


    standard[
        "Nominal_P_lt_0.05"
    ] = (
        standard[
            "P_value"
        ].notna()
        &
        (
            standard[
                "P_value"
            ] < 0.05
        )
    )


# ============================================================
# SKELETAL MUSCLE
# ============================================================

if not standard.empty:

    muscle = standard[
        standard[
            "Tissue"
        ]
        .astype(str)
        .str.lower()
        .str.contains(
            "muscle"
        )
    ].copy()

else:

    muscle = pd.DataFrame()


# ============================================================
# DIRECT KEY SNP OVERLAP
# ============================================================

if not standard.empty:

    direct_overlap = standard[
        standard[
            "Key_GWAS_variant"
        ]
    ].copy()

else:

    direct_overlap = pd.DataFrame()


# ============================================================
# TISSUE SUMMARY
# ============================================================

if not standard.empty:

    tissue_summary = (
        standard
        .groupby(
            "Tissue",
            dropna=False
        )
        .agg(
            N_EPB41L1_sb_eQTL=(
                "Variant",
                "count"
            ),

            Min_P=(
                "P_value",
                "min"
            ),

            Min_FDR=(
                "FDR_or_qvalue",
                "min"
            ),

            N_FDR_lt_005=(
                "Significant_FDR_0.05",
                "sum"
            ),

            N_nominal_P_lt_005=(
                "Nominal_P_lt_0.05",
                "sum"
            ),

            N_key_GWAS_overlap=(
                "Key_GWAS_variant",
                "sum"
            ),
        )
        .reset_index()
        .sort_values(
            [
                "Min_FDR",
                "Min_P"
            ],
            na_position="last"
        )
    )

else:

    tissue_summary = pd.DataFrame()


# ============================================================
# GLOBAL SUMMARY
# ============================================================

summary_rows = [
    {
        "Metric":
            "EPB41L1 present in sb-eQTL archive",

        "Value":
            "YES"
            if len(all_hits) > 0
            else "NO"
    },

    {
        "Metric":
            "Files/tissues containing EPB41L1",

        "Value":
            int(
                (
                    file_summary_df[
                        "EPB41L1_hits"
                    ] > 0
                ).sum()
            )
    },

    {
        "Metric":
            "Total EPB41L1 sb-eQTL records",

        "Value":
            len(standard)
    },

    {
        "Metric":
            "EPB41L1 sb-eQTL records with FDR/q < 0.05",

        "Value":
            int(
                standard[
                    "Significant_FDR_0.05"
                ].sum()
            )
            if not standard.empty
            else 0
    },

    {
        "Metric":
            "EPB41L1 skeletal-muscle sb-eQTL records",

        "Value":
            len(muscle)
    },

    {
        "Metric":
            "Direct overlaps with key GWAS variants",

        "Value":
            len(
                direct_overlap
            )
    },
]


summary = pd.DataFrame(
    summary_rows
)


# ============================================================
# PRINT RESULTS
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
print("TISSUE SUMMARY")
print("=" * 90)

if tissue_summary.empty:

    print(
        "No EPB41L1 sex-biased eQTL records found."
    )

else:

    print(
        tissue_summary
        .head(30)
        .to_string(
            index=False
        )
    )


print("\n" + "=" * 90)
print("SKELETAL MUSCLE")
print("=" * 90)

if muscle.empty:

    print(
        "No EPB41L1 sb-eQTL records detected "
        "for skeletal muscle."
    )

else:

    print(
        muscle
        .sort_values(
            [
                "FDR_or_qvalue",
                "P_value"
            ],
            na_position="last"
        )
        .head(30)
        .to_string(
            index=False
        )
    )


print("\n" + "=" * 90)
print("DIRECT KEY GWAS VARIANT OVERLAP")
print("=" * 90)

if direct_overlap.empty:

    print(
        "No direct rsID overlap detected "
        "for rs6060632 / rs1006296 / "
        "rs113092336 / rs532201406."
    )

else:

    print(
        direct_overlap
        .sort_values(
            [
                "FDR_or_qvalue",
                "P_value"
            ],
            na_position="last"
        )
        .to_string(
            index=False
        )
    )


# ============================================================
# SAVE EXCEL
# ============================================================

with pd.ExcelWriter(
    OUTPUT_FILE,
    engine="openpyxl"
) as writer:

    summary.to_excel(
        writer,
        sheet_name="S0_Summary",
        index=False
    )

    file_summary_df.to_excel(
        writer,
        sheet_name="S1_File_scan",
        index=False
    )

    standard.to_excel(
        writer,
        sheet_name="S2_EPB41L1_sbeQTL",
        index=False
    )

    tissue_summary.to_excel(
        writer,
        sheet_name="S3_Tissue_summary",
        index=False
    )

    muscle.to_excel(
        writer,
        sheet_name="S4_Skeletal_muscle",
        index=False
    )

    direct_overlap.to_excel(
        writer,
        sheet_name="S5_Key_GWAS_overlap",
        index=False
    )


# README text into workbook
readme_df = pd.DataFrame(
    readme_texts
)

readme_df.to_excel(
    OUTPUT_FILE,
    sheet_name="README",
    index=False
)


print("\nSaved:")
print(
    OUTPUT_FILE.resolve()
)

print("\nDONE.")