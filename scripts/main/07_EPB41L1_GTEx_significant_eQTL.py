from pathlib import Path
import requests
import pandas as pd
import numpy as np
import time


# ============================================================
# SETTINGS
# ============================================================

GENE_SYMBOL = "EPB41L1"

# Stable Ensembl ID
ENSEMBL_ID = "ENSG00000088367"

DATASET = "gtex_v8"

OUTPUT_XLSX = Path(
    "07_EPB41L1_GTEx_significant_eQTL.xlsx"
)

OUTPUT_TSV = Path(
    "07_EPB41L1_GTEx_significant_eQTL.tsv"
)


# ============================================================
# GTEx API
# ============================================================

BASE_URL = (
    "https://gtexportal.org/api/v2/"
    "association/singleTissueEqtl"
)


# ============================================================
# GET VERSIONED GENCODE ID
# ============================================================
#
# GTEx API often requires the versioned GENCODE identifier.
# We first query the gene endpoint.
# ============================================================

GENE_URL = (
    "https://gtexportal.org/api/v2/"
    "reference/gene"
)


print("=" * 80)
print("EPB41L1 GTEx v8 eQTL ANALYSIS")
print("=" * 80)

print("\nFinding GTEx GENCODE identifier...")


params = {
    "geneId": ENSEMBL_ID,
    "datasetId": DATASET,
}


r = requests.get(
    GENE_URL,
    params=params,
    timeout=60
)

r.raise_for_status()

gene_json = r.json()


# ============================================================
# FIND VERSIONED ID
# ============================================================

gene_records = gene_json.get(
    "data",
    []
)


if not gene_records:

    raise RuntimeError(
        "EPB41L1 was not found in the GTEx gene reference."
    )


print(
    f"Gene reference records returned: "
    f"{len(gene_records)}"
)


gencode_id = None

for record in gene_records:

    symbol = str(
        record.get(
            "geneSymbol",
            ""
        )
    ).upper()

    gid = record.get(
        "gencodeId"
    )

    if (
        symbol == GENE_SYMBOL
        and gid is not None
    ):

        gencode_id = gid
        break


if gencode_id is None:

    # fallback
    gencode_id = (
        gene_records[0]
        .get("gencodeId")
    )


if gencode_id is None:

    raise RuntimeError(
        "Could not determine versioned GTEx GENCODE ID."
    )


print(
    "GTEx GENCODE ID:",
    gencode_id
)


# ============================================================
# QUERY SIGNIFICANT eQTLs
# ============================================================

print(
    "\nDownloading significant single-tissue eQTLs..."
)


all_records = []

page = 0
items_per_page = 250


while True:

    params = {
        "gencodeId": gencode_id,
        "datasetId": DATASET,
        "page": page,
        "itemsPerPage": items_per_page,
    }


    response = requests.get(
        BASE_URL,
        params=params,
        timeout=120
    )

    response.raise_for_status()

    result = response.json()


    records = result.get(
        "data",
        []
    )


    if not records:
        break


    all_records.extend(
        records
    )


    paging = result.get(
        "paging_info",
        {}
    )


    total_pages = paging.get(
        "numberOfPages"
    )


    print(
        f"  page {page + 1}: "
        f"{len(records)} records"
    )


    page += 1


    if (
        total_pages is not None
        and page >= total_pages
    ):
        break


    if len(records) < items_per_page:
        break


    time.sleep(0.1)


print(
    "\nTotal significant GTEx eQTL records:",
    len(all_records)
)


if not all_records:

    raise RuntimeError(
        "No significant GTEx v8 eQTLs returned for EPB41L1."
    )


# ============================================================
# DATAFRAME
# ============================================================

df = pd.json_normalize(
    all_records
)


print("\nReturned columns:")

for col in df.columns:
    print(" -", col)


# ============================================================
# STANDARDIZE IMPORTANT FIELDS
# ============================================================

def find_column(
    dataframe,
    candidates
):

    lower_map = {
        c.lower(): c
        for c in dataframe.columns
    }

    for x in candidates:

        if x.lower() in lower_map:
            return lower_map[
                x.lower()
            ]

    return None


col_variant = find_column(
    df,
    [
        "variantId",
        "variant_id"
    ]
)

col_rsid = find_column(
    df,
    [
        "snpId",
        "snp_id",
        "rsId"
    ]
)

col_tissue = find_column(
    df,
    [
        "tissueSiteDetailId",
        "tissueSiteDetail"
    ]
)

col_p = find_column(
    df,
    [
        "pValue",
        "pvalue",
        "pval_nominal"
    ]
)

col_nes = find_column(
    df,
    [
        "nes",
        "normalizedEffectSize"
    ]
)

col_gene = find_column(
    df,
    [
        "geneSymbol",
        "gene_symbol"
    ]
)


# ============================================================
# CREATE CLEAN TABLE
# ============================================================

clean = pd.DataFrame()


clean["Gene"] = GENE_SYMBOL
clean["GTEx_GENCODE_ID"] = gencode_id


if col_variant:
    clean["GTEx_variant_ID"] = df[
        col_variant
    ]
else:
    clean["GTEx_variant_ID"] = np.nan


if col_rsid:
    clean["rsID"] = df[
        col_rsid
    ]
else:
    clean["rsID"] = np.nan


if col_tissue:
    clean["Tissue"] = df[
        col_tissue
    ]
else:
    clean["Tissue"] = np.nan


if col_nes:
    clean["NES"] = pd.to_numeric(
        df[col_nes],
        errors="coerce"
    )
else:
    clean["NES"] = np.nan


if col_p:
    clean["P_GTEx"] = pd.to_numeric(
        df[col_p],
        errors="coerce"
    )
else:
    clean["P_GTEx"] = np.nan


# Preserve all API information too
raw = df.copy()


# ============================================================
# EFFECT DIRECTION
# ============================================================

clean[
    "Expression_effect"
] = np.where(
    clean["NES"] > 0,
    "increased_expression",
    np.where(
        clean["NES"] < 0,
        "decreased_expression",
        "unknown"
    )
)


# ============================================================
# KEY GWAS VARIANTS
# ============================================================

KEY_GWAS_SNPS = {
    "rs113092336",
    "rs6060632",
    "rs532201406",
    "rs1006296",
}


clean[
    "Key_GWAS_variant"
] = (
    clean["rsID"]
    .astype(str)
    .isin(
        KEY_GWAS_SNPS
    )
)


# ============================================================
# SUMMARY BY TISSUE
# ============================================================

tissue_summary = (
    clean
    .groupby(
        "Tissue",
        dropna=False
    )
    .agg(
        N_significant_eQTLs=(
            "GTEx_variant_ID",
            "count"
        ),

        Min_P_GTEx=(
            "P_GTEx",
            "min"
        ),

        Median_NES=(
            "NES",
            "median"
        ),

        N_positive_NES=(
            "NES",
            lambda x:
                (x > 0).sum()
        ),

        N_negative_NES=(
            "NES",
            lambda x:
                (x < 0).sum()
        ),
    )
    .reset_index()
    .sort_values(
        "Min_P_GTEx"
    )
)


# ============================================================
# TOP eQTL PER TISSUE
# ============================================================

top_per_tissue = (
    clean[
        clean["P_GTEx"]
        .notna()
    ]
    .sort_values(
        [
            "Tissue",
            "P_GTEx"
        ]
    )
    .groupby(
        "Tissue",
        as_index=False
    )
    .first()
)


# ============================================================
# KEY GWAS VARIANT OVERLAP
# ============================================================

key_overlap = clean[
    clean[
        "Key_GWAS_variant"
    ]
].copy()


# ============================================================
# GLOBAL SUMMARY
# ============================================================

summary = pd.DataFrame(
    {
        "Metric": [
            "GTEx GENCODE ID",
            "Significant eQTL records",
            "Tissues with significant EPB41L1 eQTLs",
            "Unique significant eVariants",
            "Positive NES associations",
            "Negative NES associations",
            "Direct overlaps with key GWAS rsIDs",
        ],

        "Value": [
            gencode_id,

            len(clean),

            clean[
                "Tissue"
            ].nunique(),

            clean[
                "GTEx_variant_ID"
            ].nunique(),

            int(
                (
                    clean["NES"] > 0
                ).sum()
            ),

            int(
                (
                    clean["NES"] < 0
                ).sum()
            ),

            len(
                key_overlap
            ),
        ],
    }
)


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
print("TOP GTEx TISSUES")
print("=" * 80)

print(
    tissue_summary
    .head(20)
    .to_string(
        index=False
    )
)


print("\n" + "=" * 80)
print("DIRECT OVERLAP WITH KEY GWAS VARIANTS")
print("=" * 80)


if key_overlap.empty:

    print(
        "No direct rsID overlap found."
    )

else:

    print(
        key_overlap
        .sort_values(
            "P_GTEx"
        )
        .to_string(
            index=False
        )
    )


# ============================================================
# EXPORT
# ============================================================

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
        sheet_name="S1_All_significant_eQTL",
        index=False
    )

    tissue_summary.to_excel(
        writer,
        sheet_name="S2_Tissue_summary",
        index=False
    )

    top_per_tissue.to_excel(
        writer,
        sheet_name="S3_Top_per_tissue",
        index=False
    )

    key_overlap.to_excel(
        writer,
        sheet_name="S4_Key_GWAS_overlap",
        index=False
    )

    raw.to_excel(
        writer,
        sheet_name="S5_Raw_API",
        index=False
    )


clean.to_csv(
    OUTPUT_TSV,
    sep="\t",
    index=False
)


print("\nSaved:")
print(OUTPUT_XLSX.resolve())
print(OUTPUT_TSV.resolve())

print("\nDONE.")