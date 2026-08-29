from pathlib import Path
import gzip
import heapq
from collections import defaultdict

import numpy as np
import pandas as pd

from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


# ============================================================
# SETTINGS
# ============================================================

COORD_FILE = Path(
    "03_HUMAN_GENE_COORDINATES_GRCh37.xlsx"
)

COORD_SHEET = "S1_GWAS_loci_GRCh37"

MALE_FILE = Path(
    "ebi-a-GCST90000026.vcf.gz"
)

FEMALE_FILE = Path(
    "ebi-a-GCST90000027.vcf.gz"
)

COMBINED_FILE = Path(
    "ebi-a-GCST90000025.vcf.gz"
)


OUTPUT_EXCEL = Path(
    "04_ALL_GENES_SEX_STRATIFIED_GWAS_SCREEN.xlsx"
)

OUTPUT_VARIANTS = Path(
    "04_candidate_variants.tsv.gz"
)


# ============================================================
# ANALYSIS PARAMETERS
# ============================================================

# Already defined in 03b:
WINDOW_BP = 500_000

# Common-variant analysis
MIN_MAF = 0.01

# Keep the top N combined-GWAS variants per gene.
# Then choose the strongest one that also passes male/female QC.
TOP_K = 50

# Conventional GWAS threshold
GWAS_P_THRESHOLD = 5e-8

FOCAL_GENES = {
    "CLPP",
    "MTA1",
    "RSPO1",
    "SNAP47",
    "LINGO2",
}


# ============================================================
# BASIC FUNCTIONS
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

    return min(af, 1 - af)


def lp_to_p(lp):

    if pd.isna(lp):
        return np.nan

    if lp > 320:
        return 0.0

    return 10 ** (-lp)


def parse_format(format_string, sample_string):

    keys = format_string.split(":")
    values = sample_string.split(":")

    return dict(zip(keys, values))


def variant_key(chrom, pos, ref, alt):

    return (
        f"{chrom}:{pos}:{ref}:{alt}"
    )


def clean_chrom(x):

    return (
        str(x)
        .replace("chr", "")
        .strip()
    )


def is_valid_biallelic(ref, alt):

    if pd.isna(ref) or pd.isna(alt):
        return False

    if "," in str(alt):
        return False

    # symbolic structural variants
    if "<" in str(alt) or ">" in str(alt):
        return False

    return True


# ============================================================
# CHECK INPUT FILES
# ============================================================

print("=" * 90)
print("SYSTEMATIC SEX-STRATIFIED GWAS SCREEN")
print("=" * 90)

for f in [
    COORD_FILE,
    MALE_FILE,
    FEMALE_FILE,
    COMBINED_FILE
]:

    if not f.exists():

        raise FileNotFoundError(
            f"File not found:\n{f.resolve()}"
        )

    if f.suffix == ".gz":

        print(
            f"{f.name}: "
            f"{f.stat().st_size / 1024**3:.2f} GB"
        )

    else:

        print(f"{f.name}: found")


# ============================================================
# 1. LOAD 439 HUMAN GENE LOCI
# ============================================================

print("\n" + "=" * 90)
print("1/10 Loading GRCh37 candidate loci...")
print("=" * 90)

genes = pd.read_excel(
    COORD_FILE,
    sheet_name=COORD_SHEET,
    engine="openpyxl"
)


required_cols = [
    "Gene",
    "Human_gene",
    "Human_Ensembl",
    "Sex_pattern",
    "Male_KO_direction",
    "Female_KO_direction",
    "Chromosome",
    "Gene_start_GRCh37",
    "Gene_end_GRCh37",
    "Locus_start_GRCh37",
    "Locus_end_GRCh37",
]

missing = [
    c
    for c in required_cols
    if c not in genes.columns
]

if missing:

    raise ValueError(
        "Missing coordinate columns:\n"
        + "\n".join(missing)
    )


# Numeric
for c in [
    "Gene_start_GRCh37",
    "Gene_end_GRCh37",
    "Locus_start_GRCh37",
    "Locus_end_GRCh37",
]:

    genes[c] = pd.to_numeric(
        genes[c],
        errors="coerce"
    )


genes["Chromosome"] = (
    genes["Chromosome"]
    .apply(clean_chrom)
)


genes = genes.dropna(
    subset=[
        "Chromosome",
        "Gene_start_GRCh37",
        "Gene_end_GRCh37",
        "Locus_start_GRCh37",
        "Locus_end_GRCh37",
    ]
).copy()


for c in [
    "Gene_start_GRCh37",
    "Gene_end_GRCh37",
    "Locus_start_GRCh37",
    "Locus_end_GRCh37",
]:

    genes[c] = genes[c].astype(int)


genes = (
    genes
    .drop_duplicates(
        subset=[
            "Human_Ensembl"
        ]
    )
    .reset_index(drop=True)
)


genes["Gene_ID"] = [
    f"G{i:04d}"
    for i in range(1, len(genes) + 1)
]


print(
    "Candidate human genes:",
    len(genes)
)

print(
    "\nMouse sex patterns:"
)

print(
    genes["Sex_pattern"]
    .value_counts()
)


# ============================================================
# 2. MERGE OVERLAPPING ±500 kb WINDOWS
# ============================================================

print("\n" + "=" * 90)
print("2/10 Merging overlapping genomic intervals...")
print("=" * 90)


merged_intervals = []


for chrom, chr_df in genes.groupby(
    "Chromosome"
):

    chr_df = chr_df.sort_values(
        "Locus_start_GRCh37"
    )

    current_start = None
    current_end = None
    current_genes = []

    for _, row in chr_df.iterrows():

        start = int(
            row["Locus_start_GRCh37"]
        )

        end = int(
            row["Locus_end_GRCh37"]
        )

        gene_id = row["Gene_ID"]

        if current_start is None:

            current_start = start
            current_end = end
            current_genes = [gene_id]

        elif start <= current_end:

            current_end = max(
                current_end,
                end
            )

            current_genes.append(
                gene_id
            )

        else:

            merged_intervals.append(
                {
                    "Chromosome": chrom,
                    "Merged_start": current_start,
                    "Merged_end": current_end,
                    "Gene_IDs": current_genes.copy(),
                }
            )

            current_start = start
            current_end = end
            current_genes = [gene_id]


    if current_start is not None:

        merged_intervals.append(
            {
                "Chromosome": chrom,
                "Merged_start": current_start,
                "Merged_end": current_end,
                "Gene_IDs": current_genes.copy(),
            }
        )


merged = pd.DataFrame(
    merged_intervals
)


merged["Merged_Locus_ID"] = [
    f"L{i:04d}"
    for i in range(1, len(merged) + 1)
]


gene_id_to_name = dict(
    zip(
        genes["Gene_ID"],
        genes["Human_gene"]
    )
)


merged["Genes"] = (
    merged["Gene_IDs"]
    .apply(
        lambda xs:
        ";".join(
            gene_id_to_name[x]
            for x in xs
        )
    )
)


merged["N_genes"] = (
    merged["Gene_IDs"]
    .apply(len)
)


print(
    "439 gene windows collapsed into",
    len(merged),
    "non-overlapping genomic regions."
)

print(
    "Merged regions containing >1 candidate gene:",
    (merged["N_genes"] > 1).sum()
)


# ============================================================
# LOOKUP STRUCTURES
# ============================================================

genes_by_chr = {}

for chrom, x in genes.groupby(
    "Chromosome"
):

    records = (
        x.sort_values(
            "Locus_start_GRCh37"
        )
        .to_dict("records")
    )

    genes_by_chr[chrom] = records


merged_by_chr = {}

for chrom, x in merged.groupby(
    "Chromosome"
):

    records = (
        x.sort_values(
            "Merged_start"
        )
        .to_dict("records")
    )

    merged_by_chr[chrom] = records


# ============================================================
# TOP-K HEAPS
# ============================================================

# Stores best combined variants per gene
gene_heaps = defaultdict(list)

# Stores best variants physically inside gene boundaries
inside_gene_heaps = defaultdict(list)

# Stores best variants per merged genomic region
locus_heaps = defaultdict(list)


counter = 0


def add_to_heap(
    heap_dict,
    key,
    record,
    score
):

    global counter

    counter += 1

    heap = heap_dict[key]

    item = (
        score,
        counter,
        record
    )

    if len(heap) < TOP_K:

        heapq.heappush(
            heap,
            item
        )

    elif score > heap[0][0]:

        heapq.heapreplace(
            heap,
            item
        )


# ============================================================
# 3. SCAN COMBINED GWAS ONCE
# ============================================================

print("\n" + "=" * 90)
print("3/10 Scanning combined GWAS...")
print("=" * 90)

n_total = 0
n_union = 0
n_qc = 0


current_chrom = None

gene_idx = 0
active_genes = []

locus_idx = 0
current_loci = []


with gzip.open(
    COMBINED_FILE,
    "rt",
    encoding="utf-8",
    errors="replace"
) as handle:

    for line in handle:

        if line.startswith("#"):
            continue

        n_total += 1

        if n_total % 2_000_000 == 0:

            print(
                f"  scanned {n_total:,} variants"
            )


        parts = (
            line.rstrip("\n")
            .split("\t")
        )

        if len(parts) < 10:
            continue


        chrom = clean_chrom(
            parts[0]
        )


        if chrom not in merged_by_chr:
            continue


        try:
            pos = int(parts[1])
        except Exception:
            continue


        # ----------------------------------------------------
        # Reset sweep structures when chromosome changes
        # ----------------------------------------------------

        if chrom != current_chrom:

            current_chrom = chrom

            gene_idx = 0
            active_genes = []

            locus_idx = 0

            chr_genes = (
                genes_by_chr.get(
                    chrom,
                    []
                )
            )

            chr_loci = (
                merged_by_chr.get(
                    chrom,
                    []
                )
            )


        # ----------------------------------------------------
        # Locate merged candidate region
        # ----------------------------------------------------

        while (
            locus_idx < len(chr_loci)
            and
            chr_loci[locus_idx]["Merged_end"] < pos
        ):

            locus_idx += 1


        if locus_idx >= len(chr_loci):
            continue


        locus = chr_loci[locus_idx]


        if not (
            locus["Merged_start"]
            <= pos
            <= locus["Merged_end"]
        ):

            continue


        n_union += 1


        # ----------------------------------------------------
        # Variant fields
        # ----------------------------------------------------

        rsid = parts[2]
        ref = parts[3]
        alt = parts[4]

        if not is_valid_biallelic(
            ref,
            alt
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


        if (
            pd.isna(beta)
            or
            pd.isna(se)
            or
            se <= 0
            or
            pd.isna(lp)
            or
            pd.isna(maf)
            or
            maf < MIN_MAF
        ):

            continue


        n_qc += 1


        key = variant_key(
            chrom,
            pos,
            ref,
            alt
        )


        base_record = {

            "Variant_key": key,

            "CHROM": chrom,
            "POS": pos,

            "rsID": rsid,

            "REF": ref,
            "ALT": alt,

            "Effect_allele": alt,
            "Other_allele": ref,

            "Beta_Combined": beta,
            "SE_Combined": se,

            "LP_Combined": lp,
            "P_Combined": lp_to_p(lp),

            "AF_Combined": af,
            "MAF_Combined": maf,

            "Merged_Locus_ID":
                locus["Merged_Locus_ID"],
        }


        # ----------------------------------------------------
        # Update merged-locus top variants
        # ----------------------------------------------------

        add_to_heap(
            locus_heaps,
            locus["Merged_Locus_ID"],
            base_record.copy(),
            lp
        )


        # ----------------------------------------------------
        # Maintain active candidate genes
        # ----------------------------------------------------

        while (
            gene_idx < len(chr_genes)
            and
            chr_genes[gene_idx][
                "Locus_start_GRCh37"
            ] <= pos
        ):

            active_genes.append(
                chr_genes[gene_idx]
            )

            gene_idx += 1


        active_genes = [
            g
            for g in active_genes
            if g[
                "Locus_end_GRCh37"
            ] >= pos
        ]


        # ----------------------------------------------------
        # Assign variant to all overlapping candidate genes
        # ----------------------------------------------------

        for g in active_genes:

            if not (
                g["Locus_start_GRCh37"]
                <= pos
                <= g["Locus_end_GRCh37"]
            ):
                continue


            record = (
                base_record.copy()
            )


            record["Gene_ID"] = (
                g["Gene_ID"]
            )

            record["Gene"] = (
                g["Human_gene"]
            )

            record["Human_Ensembl"] = (
                g["Human_Ensembl"]
            )


            inside_gene = (
                g["Gene_start_GRCh37"]
                <= pos
                <= g["Gene_end_GRCh37"]
            )


            record["Inside_gene"] = (
                inside_gene
            )


            if pos < g["Gene_start_GRCh37"]:

                distance = (
                    g["Gene_start_GRCh37"]
                    - pos
                )

            elif pos > g["Gene_end_GRCh37"]:

                distance = (
                    pos
                    - g["Gene_end_GRCh37"]
                )

            else:

                distance = 0


            record["Distance_to_gene"] = (
                distance
            )


            add_to_heap(
                gene_heaps,
                g["Gene_ID"],
                record.copy(),
                lp
            )


            if inside_gene:

                add_to_heap(
                    inside_gene_heaps,
                    g["Gene_ID"],
                    record.copy(),
                    lp
                )


print(
    "\nCombined GWAS variants scanned:",
    f"{n_total:,}"
)

print(
    "Variants inside union of candidate regions:",
    f"{n_union:,}"
)

print(
    "Variants passing combined QC:",
    f"{n_qc:,}"
)


# ============================================================
# CONVERT HEAPS TO CANDIDATE TABLE
# ============================================================

def heap_to_records(
    heap_dict,
    source
):

    records = []

    for group_id, heap in heap_dict.items():

        sorted_heap = sorted(
            heap,
            key=lambda x: x[0],
            reverse=True
        )

        for rank, (_, _, record) in enumerate(
            sorted_heap,
            start=1
        ):

            x = record.copy()

            x["Candidate_source"] = (
                source
            )

            x["Combined_rank"] = rank

            records.append(x)

    return records


candidate_records = []

candidate_records += heap_to_records(
    gene_heaps,
    "gene_window"
)

candidate_records += heap_to_records(
    inside_gene_heaps,
    "inside_gene"
)

candidate_records += heap_to_records(
    locus_heaps,
    "merged_locus"
)


candidate_df = pd.DataFrame(
    candidate_records
)


candidate_keys = set(
    candidate_df[
        "Variant_key"
    ].dropna()
)


print(
    "\nUnique variants retained for male/female lookup:",
    f"{len(candidate_keys):,}"
)


# ============================================================
# 4. EXTRACT ONLY CANDIDATE VARIANTS FROM SEX-SPECIFIC GWAS
# ============================================================

def extract_selected_variants(
    filepath,
    dataset_name,
    keys
):

    print("\n" + "=" * 90)
    print(
        f"4/10 Extracting selected variants: {dataset_name}"
    )
    print("=" * 90)

    out = {}

    n_total = 0
    n_found = 0


    with gzip.open(
        filepath,
        "rt",
        encoding="utf-8",
        errors="replace"
    ) as handle:

        for line in handle:

            if line.startswith("#"):
                continue


            n_total += 1


            if n_total % 3_000_000 == 0:

                print(
                    f"  scanned {n_total:,}"
                )


            parts = (
                line.rstrip("\n")
                .split("\t")
            )


            if len(parts) < 10:
                continue


            chrom = clean_chrom(
                parts[0]
            )

            try:
                pos = int(parts[1])
            except Exception:
                continue


            ref = parts[3]
            alt = parts[4]


            key = variant_key(
                chrom,
                pos,
                ref,
                alt
            )


            if key not in keys:
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


            out[key] = {

                f"Beta_{dataset_name}":
                    beta,

                f"SE_{dataset_name}":
                    se,

                f"LP_{dataset_name}":
                    lp,

                f"P_{dataset_name}":
                    lp_to_p(lp),

                f"AF_{dataset_name}":
                    af,

                f"MAF_{dataset_name}":
                    maf_from_af(af),
            }


            n_found += 1


    print(
        f"Found {n_found:,} / "
        f"{len(keys):,} candidate variants"
    )

    return out


male_lookup = extract_selected_variants(
    MALE_FILE,
    "Male",
    candidate_keys
)


female_lookup = extract_selected_variants(
    FEMALE_FILE,
    "Female",
    candidate_keys
)


# ============================================================
# 5. ATTACH SEX-SPECIFIC STATISTICS
# ============================================================

print("\n" + "=" * 90)
print("5/10 Harmonizing candidate variant statistics...")
print("=" * 90)


unique_variants = (
    candidate_df[
        [
            "Variant_key",
            "CHROM",
            "POS",
            "rsID",
            "REF",
            "ALT",
            "Effect_allele",
            "Other_allele",
            "Beta_Combined",
            "SE_Combined",
            "LP_Combined",
            "P_Combined",
            "AF_Combined",
            "MAF_Combined",
        ]
    ]
    .drop_duplicates(
        subset=["Variant_key"]
    )
    .copy()
)


male_df = pd.DataFrame.from_dict(
    male_lookup,
    orient="index"
)

male_df["Variant_key"] = (
    male_df.index
)

male_df = male_df.reset_index(
    drop=True
)


female_df = pd.DataFrame.from_dict(
    female_lookup,
    orient="index"
)

female_df["Variant_key"] = (
    female_df.index
)

female_df = female_df.reset_index(
    drop=True
)


variants = (
    unique_variants
    .merge(
        male_df,
        on="Variant_key",
        how="left"
    )
    .merge(
        female_df,
        on="Variant_key",
        how="left"
    )
)


# ============================================================
# COMPLETE VARIANT QC
# ============================================================

variants["Pass_sex_QC"] = (

    variants["Beta_Male"].notna()
    &
    variants["Beta_Female"].notna()

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
    (variants["MAF_Male"] >= MIN_MAF)
    &
    (variants["MAF_Female"] >= MIN_MAF)
)


# ============================================================
# SEX HETEROGENEITY
# ============================================================

print("\nCalculating sex heterogeneity...")

variants["Z_sex"] = np.nan
variants["P_sex"] = np.nan


# Select ROWS that passed QC
valid_mask = (
    variants["Pass_sex_QC"] == True
)

valid = variants.loc[
    valid_mask
].copy()


print(
    "Variants passing male/female QC:",
    f"{len(valid):,}"
)


if len(valid) > 0:

    # --------------------------------------------------------
    # Formal male-female heterogeneity test
    #
    # Male and female samples are independent:
    #
    # Z = (beta_male - beta_female)
    #     / sqrt(SE_male^2 + SE_female^2)
    # --------------------------------------------------------

    valid["Z_sex"] = (
        (
            valid["Beta_Male"]
            -
            valid["Beta_Female"]
        )
        /
        np.sqrt(
            valid["SE_Male"] ** 2
            +
            valid["SE_Female"] ** 2
        )
    )


    valid["P_sex"] = (
        2
        *
        norm.sf(
            np.abs(
                valid["Z_sex"]
            )
        )
    )


    # Write results back into full variant table
    variants.loc[
        valid.index,
        "Z_sex"
    ] = valid["Z_sex"]


    variants.loc[
        valid.index,
        "P_sex"
    ] = valid["P_sex"]


print(
    "Variants with heterogeneity test:",
    variants["P_sex"].notna().sum()
)


# ============================================================
# MALE/FEMALE EFFECT PATTERN
# ============================================================

def beta_pattern(row):

    bm = row["Beta_Male"]
    bf = row["Beta_Female"]

    if pd.isna(bm) or pd.isna(bf):
        return "missing"

    if bm > 0 and bf > 0:
        return "same_positive"

    if bm < 0 and bf < 0:
        return "same_negative"

    if bm > 0 and bf < 0:
        return "male_positive_female_negative"

    if bm < 0 and bf > 0:
        return "male_negative_female_positive"

    return "other"


variants[
    "Human_beta_pattern"
] = variants.apply(
    beta_pattern,
    axis=1
)

# ============================================================
# 6. SELECT ONE INDEX SNP PER GENE
# ============================================================

print("\n" + "=" * 90)
print("6/10 Selecting final combined-GWAS index SNP per gene...")
print("=" * 90)


candidate_df = candidate_df.merge(
    variants[
        [
            "Variant_key",
            "Pass_sex_QC",
            "Beta_Male",
            "SE_Male",
            "P_Male",
            "AF_Male",
            "MAF_Male",
            "Beta_Female",
            "SE_Female",
            "P_Female",
            "AF_Female",
            "MAF_Female",
            "Z_sex",
            "P_sex",
            "Human_beta_pattern",
        ]
    ],
    on="Variant_key",
    how="left"
)


gene_candidates = candidate_df[
    candidate_df[
        "Candidate_source"
    ] == "gene_window"
].copy()


gene_candidates = gene_candidates[
    gene_candidates[
        "Pass_sex_QC"
    ] == True
].copy()


gene_candidates = gene_candidates.sort_values(
    [
        "Gene_ID",
        "Combined_rank"
    ]
)


gene_index = (
    gene_candidates
    .groupby(
        "Gene_ID",
        as_index=False
    )
    .first()
)


# ============================================================
# ADD ORIGINAL GENE INFORMATION
# ============================================================

gene_info_cols = [
    "Gene_ID",
    "Gene",
    "Human_gene",
    "Human_Ensembl",
    "Sex_pattern",
    "Male_KO_direction",
    "Female_KO_direction",
    "Chromosome",
    "Gene_start_GRCh37",
    "Gene_end_GRCh37",
    "Locus_start_GRCh37",
    "Locus_end_GRCh37",
]


gene_info = genes[
    gene_info_cols
].copy()


gene_index = gene_info.merge(
    gene_index.drop(
        columns=[
            "Gene",
            "Human_Ensembl"
        ],
        errors="ignore"
    ),
    on="Gene_ID",
    how="left"
)


# ============================================================
# MULTIPLE TESTING ACROSS SYSTEMATIC GENE SCREEN
# ============================================================

gene_index[
    "P_sex_Bonferroni_all_genes"
] = np.nan

gene_index[
    "P_sex_FDR_all_genes"
] = np.nan


testable = gene_index[
    "P_sex"
].notna()


n_testable_genes = (
    testable.sum()
)


if n_testable_genes > 0:

    p = gene_index.loc[
        testable,
        "P_sex"
    ].values


    gene_index.loc[
        testable,
        "P_sex_Bonferroni_all_genes"
    ] = np.minimum(
        p * n_testable_genes,
        1
    )


    gene_index.loc[
        testable,
        "P_sex_FDR_all_genes"
    ] = multipletests(
        p,
        method="fdr_bh"
    )[1]


gene_index[
    "Sex_difference_Bonferroni"
] = (
    gene_index[
        "P_sex_Bonferroni_all_genes"
    ] < 0.05
)

gene_index[
    "Sex_difference_FDR"
] = (
    gene_index[
        "P_sex_FDR_all_genes"
    ] < 0.05
)


# GWAS significance
gene_index[
    "Combined_GWAS_significant"
] = (
    gene_index[
        "P_Combined"
    ] < GWAS_P_THRESHOLD
)

gene_index[
    "Male_GWAS_significant"
] = (
    gene_index[
        "P_Male"
    ] < GWAS_P_THRESHOLD
)

gene_index[
    "Female_GWAS_significant"
] = (
    gene_index[
        "P_Female"
    ] < GWAS_P_THRESHOLD
)


# ============================================================
# CHECK SHARED INDEX VARIANTS
# ============================================================

shared_counts = (
    gene_index[
        "Variant_key"
    ]
    .value_counts()
)


gene_index[
    "N_candidate_genes_sharing_index_variant"
] = (
    gene_index[
        "Variant_key"
    ]
    .map(shared_counts)
)


# ============================================================
# 7. WITHIN-GENE LEAD VARIANTS
# ============================================================

inside_candidates = candidate_df[
    (
        candidate_df[
            "Candidate_source"
        ] == "inside_gene"
    )
    &
    (
        candidate_df[
            "Pass_sex_QC"
        ] == True
    )
].copy()


inside_candidates = (
    inside_candidates
    .sort_values(
        [
            "Gene_ID",
            "Combined_rank"
        ]
    )
)


within_gene_leads = (
    inside_candidates
    .groupby(
        "Gene_ID",
        as_index=False
    )
    .first()
)


within_gene_leads = (
    gene_info.merge(
        within_gene_leads.drop(
            columns=[
                "Gene",
                "Human_Ensembl"
            ],
            errors="ignore"
        ),
        on="Gene_ID",
        how="left"
    )
)


# ============================================================
# 8. MERGED-LOCUS INDEX SNPs
# ============================================================

locus_candidates = candidate_df[
    (
        candidate_df[
            "Candidate_source"
        ] == "merged_locus"
    )
    &
    (
        candidate_df[
            "Pass_sex_QC"
        ] == True
    )
].copy()


locus_candidates = locus_candidates.sort_values(
    [
        "Merged_Locus_ID",
        "Combined_rank"
    ]
)


locus_index = (
    locus_candidates
    .groupby(
        "Merged_Locus_ID",
        as_index=False
    )
    .first()
)


locus_info = merged[
    [
        "Merged_Locus_ID",
        "Chromosome",
        "Merged_start",
        "Merged_end",
        "Genes",
        "N_genes",
    ]
].copy()


locus_index = locus_info.merge(
    locus_index.drop(
        columns=[
            "CHROM"
        ],
        errors="ignore"
    ),
    on="Merged_Locus_ID",
    how="left"
)


# Multiple testing across independent predefined
# merged candidate regions
locus_index[
    "P_sex_FDR_loci"
] = np.nan

locus_index[
    "P_sex_Bonferroni_loci"
] = np.nan


valid_locus = (
    locus_index["P_sex"].notna()
)


n_testable_loci = (
    valid_locus.sum()
)


if n_testable_loci > 0:

    p = locus_index.loc[
        valid_locus,
        "P_sex"
    ].values


    locus_index.loc[
        valid_locus,
        "P_sex_Bonferroni_loci"
    ] = np.minimum(
        p * n_testable_loci,
        1
    )


    locus_index.loc[
        valid_locus,
        "P_sex_FDR_loci"
    ] = multipletests(
        p,
        method="fdr_bh"
    )[1]


# ============================================================
# 9. PREDEFINED FOCAL FIVE
# ============================================================

print("\n" + "=" * 90)
print("9/10 Focused analysis of five predefined genes...")
print("=" * 90)


focal = gene_index[
    gene_index[
        "Human_gene"
    ].isin(
        FOCAL_GENES
    )
].copy()


focal[
    "P_sex_Bonferroni_focal5"
] = np.nan

focal[
    "P_sex_FDR_focal5"
] = np.nan


valid_focal = (
    focal["P_sex"].notna()
)


if valid_focal.sum() > 0:

    p = focal.loc[
        valid_focal,
        "P_sex"
    ].values


    focal.loc[
        valid_focal,
        "P_sex_Bonferroni_focal5"
    ] = np.minimum(
        p * valid_focal.sum(),
        1
    )


    focal.loc[
        valid_focal,
        "P_sex_FDR_focal5"
    ] = multipletests(
        p,
        method="fdr_bh"
    )[1]


# ============================================================
# TOP SEX-HETEROGENEITY GENES
# ============================================================

top_sex = (
    gene_index[
        gene_index[
            "P_sex"
        ].notna()
    ]
    .sort_values(
        "P_sex"
    )
    .head(50)
    .copy()
)


# ============================================================
# SUMMARY
# ============================================================

summary = pd.DataFrame(
    {
        "Metric": [

            "Human candidate genes",

            "Merged genomic candidate regions",

            "Genes with testable index SNP",

            "Merged regions with testable index SNP",

            "Gene index SNPs genome-wide significant in combined GWAS",

            "Gene index SNPs genome-wide significant in male GWAS",

            "Gene index SNPs genome-wide significant in female GWAS",

            "Gene index SNPs with opposite male/female beta signs",

            "Gene index SNPs with nominal sex heterogeneity P<0.05",

            "Gene index SNPs significant after FDR across genes",

            "Gene index SNPs significant after Bonferroni across genes",

            "Merged-locus index SNPs significant after FDR",

            "Focal genes retained",
        ],

        "N": [

            len(genes),

            len(merged),

            gene_index[
                "P_sex"
            ].notna().sum(),

            locus_index[
                "P_sex"
            ].notna().sum(),

            gene_index[
                "Combined_GWAS_significant"
            ].sum(),

            gene_index[
                "Male_GWAS_significant"
            ].sum(),

            gene_index[
                "Female_GWAS_significant"
            ].sum(),

            gene_index[
                "Human_beta_pattern"
            ].isin(
                [
                    "male_positive_female_negative",
                    "male_negative_female_positive"
                ]
            ).sum(),

            (
                gene_index[
                    "P_sex"
                ] < 0.05
            ).sum(),

            (
                gene_index[
                    "P_sex_FDR_all_genes"
                ] < 0.05
            ).sum(),

            (
                gene_index[
                    "P_sex_Bonferroni_all_genes"
                ] < 0.05
            ).sum(),

            (
                locus_index[
                    "P_sex_FDR_loci"
                ] < 0.05
            ).sum(),

            len(focal),
        ]
    }
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
print("TOP 20 SEX-HETEROGENEITY CANDIDATE GENES")
print("=" * 90)

display_cols = [
    "Gene",
    "Human_gene",
    "Sex_pattern",
    "rsID",
    "Variant_key",
    "Beta_Male",
    "P_Male",
    "Beta_Female",
    "P_Female",
    "Beta_Combined",
    "P_Combined",
    "P_sex",
    "P_sex_FDR_all_genes",
    "P_sex_Bonferroni_all_genes",
    "Human_beta_pattern",
]


print(
    top_sex[
        display_cols
    ]
    .head(20)
    .to_string(
        index=False
    )
)


print("\n" + "=" * 90)
print("FOCAL FIVE GENES")
print("=" * 90)

focal_cols = [
    "Gene",
    "Human_gene",
    "Sex_pattern",
    "Male_KO_direction",
    "Female_KO_direction",
    "rsID",
    "Variant_key",
    "Beta_Male",
    "P_Male",
    "Beta_Female",
    "P_Female",
    "Beta_Combined",
    "P_Combined",
    "P_sex",
    "P_sex_Bonferroni_focal5",
    "P_sex_FDR_focal5",
]


print(
    focal[
        focal_cols
    ]
    .sort_values(
        "Human_gene"
    )
    .to_string(
        index=False
    )
)


# ============================================================
# 10. EXPORT
# ============================================================

print("\n" + "=" * 90)
print("10/10 Writing output files...")
print("=" * 90)


with pd.ExcelWriter(
    OUTPUT_EXCEL,
    engine="openpyxl"
) as writer:

    summary.to_excel(
        writer,
        sheet_name="S0_Summary",
        index=False
    )

    genes.to_excel(
        writer,
        sheet_name="S1_Gene_input",
        index=False
    )

    merged.drop(
        columns=["Gene_IDs"]
    ).to_excel(
        writer,
        sheet_name="S2_Merged_regions",
        index=False
    )

    gene_index.to_excel(
        writer,
        sheet_name="S3_Gene_index_SNPs",
        index=False
    )

    top_sex.to_excel(
        writer,
        sheet_name="S4_Top_sex_heterogeneity",
        index=False
    )

    within_gene_leads.to_excel(
        writer,
        sheet_name="S5_Within_gene_leads",
        index=False
    )

    locus_index.to_excel(
        writer,
        sheet_name="S6_Merged_locus_index",
        index=False
    )

    focal.to_excel(
        writer,
        sheet_name="S7_Focal_5_genes",
        index=False
    )


print(
    "\nSaved Excel:",
    OUTPUT_EXCEL.resolve()
)

print(
    "Saved candidate variants:",
    OUTPUT_VARIANTS.resolve()
)

print("\nDONE.")