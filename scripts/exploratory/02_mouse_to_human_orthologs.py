from pathlib import Path
import pandas as pd


# ============================================================
# FILES
# ============================================================

IMPC_FILE = Path(
    "01_IMPC_sex_specific_analysis.xlsx"
)

ORTHOLOG_FILE = Path(
    "ortholog_mapping_5_key_genes.xlsx"
)

OUTPUT_FILE = Path(
    "02_KEY_GENES_MOUSE_HUMAN_SEX_ANALYSIS.xlsx"
)


# ============================================================
# LOAD IMPC SEX COMPARISON
# ============================================================

impc = pd.read_excel(
    IMPC_FILE,
    sheet_name="S4_Sex_comparison",
    engine="openpyxl"
)

print("IMPC genes:", impc["Gene"].nunique())


# ============================================================
# LOAD VERIFIED ORTHOLOGS
# ============================================================

orth = pd.read_excel(
    ORTHOLOG_FILE,
    sheet_name="Gene_level_mapping",
    engine="openpyxl"
)

print("Verified orthologs:", len(orth))


# ============================================================
# MERGE
# ============================================================

result = impc.merge(
    orth,
    left_on="Gene",
    right_on="Mouse_gene",
    how="inner"
)


# ============================================================
# CHECK EXPECTED GENES
# ============================================================

expected = {
    "Clpp",
    "Mta1",
    "Rspo1",
    "Snap47",
    "Lingo2"
}

found = set(result["Gene"])

missing = expected - found

if missing:
    raise ValueError(
        f"Missing key genes after merge: {sorted(missing)}"
    )


# ============================================================
# ADD BIOLOGICAL INTERPRETATION
# ============================================================

def interpret(direction):

    if direction == "increased":
        return "potential_negative_regulator"

    elif direction == "decreased":
        return "potential_positive_regulator"

    return None


result["Male_expected_role"] = (
    result["Male_KO_direction"]
    .apply(interpret)
)

result["Female_expected_role"] = (
    result["Female_KO_direction"]
    .apply(interpret)
)


# ============================================================
# CLASSIFY DIRECTIONAL DIFFERENCE
# ============================================================

def sex_interpretation(row):

    pattern = row["Sex_pattern"]

    if pattern == "sex_conserved":

        return (
            "Same KO direction in males and females"
        )

    if pattern == "sex_opposite":

        return (
            "Opposite KO direction between sexes"
        )

    return pattern


result["Interpretation"] = (
    result.apply(
        sex_interpretation,
        axis=1
    )
)


# ============================================================
# SELECT PUBLICATION-RELEVANT COLUMNS
# ============================================================

publication = result[
    [
        "Gene",
        "Human_gene",
        "Mouse_Ensembl",
        "Human_Ensembl",

        "Sex_pattern",

        "Male_KO_direction",
        "Male_P_IMPC_min",
        "Male_expected_role",

        "Female_KO_direction",
        "Female_P_IMPC_min",
        "Female_expected_role",

        "Orthology_type",
        "Orthology_confidence",

        "Human_identity_to_mouse_pct",
        "Mouse_identity_to_human_pct",

        "Interpretation"
    ]
].copy()


# ============================================================
# SORT
# ============================================================

order = {
    "sex_opposite": 0,
    "sex_conserved": 1,
    "male_only": 2,
    "female_only": 3
}

publication["sort_order"] = (
    publication["Sex_pattern"]
    .map(order)
)

publication = (
    publication
    .sort_values(
        ["sort_order", "Gene"]
    )
    .drop(columns="sort_order")
)


# ============================================================
# SPLIT GROUPS
# ============================================================

sex_opposite = publication[
    publication["Sex_pattern"]
    == "sex_opposite"
].copy()

sex_conserved = publication[
    publication["Sex_pattern"]
    == "sex_conserved"
].copy()


# ============================================================
# SUMMARY
# ============================================================

summary = pd.DataFrame(
    {
        "Metric": [
            "Key genes tested",
            "Sex-opposite genes",
            "Sex-conserved genes",
            "High-confidence one-to-one orthologs"
        ],

        "N": [
            publication["Gene"].nunique(),
            sex_opposite["Gene"].nunique(),
            sex_conserved["Gene"].nunique(),
            (
                (
                    publication["Orthology_type"]
                    == "ortholog_one2one"
                )
                &
                (
                    publication["Orthology_confidence"]
                    == 1
                )
            ).sum()
        ]
    }
)


# ============================================================
# PRINT
# ============================================================

print("\nSUMMARY")
print(summary.to_string(index=False))

print("\nSEX-OPPOSITE GENES")
print(
    sex_opposite[
        [
            "Gene",
            "Human_gene",
            "Male_KO_direction",
            "Female_KO_direction",
            "Male_P_IMPC_min",
            "Female_P_IMPC_min"
        ]
    ].to_string(index=False)
)

print("\nSEX-CONSERVED GENES")
print(
    sex_conserved[
        [
            "Gene",
            "Human_gene",
            "Male_KO_direction",
            "Female_KO_direction",
            "Male_P_IMPC_min",
            "Female_P_IMPC_min"
        ]
    ].to_string(index=False)
)


# ============================================================
# EXPORT
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

    publication.to_excel(
        writer,
        sheet_name="S1_Key_genes",
        index=False
    )

    sex_opposite.to_excel(
        writer,
        sheet_name="S2_Sex_opposite",
        index=False
    )

    sex_conserved.to_excel(
        writer,
        sheet_name="S3_Sex_conserved",
        index=False
    )


print(
    "\nSaved:",
    OUTPUT_FILE.resolve()
)