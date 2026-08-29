from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

INPUT = Path("04_ALL_GENES_SEX_STRATIFIED_GWAS_SCREEN.xlsx")
SHEET = "S3_Gene_index_SNPs"
OUTDIR = Path("figures")
OUTDIR.mkdir(exist_ok=True)

if not INPUT.exists():
    raise FileNotFoundError(f"Missing input: {INPUT.resolve()}")

df = pd.read_excel(INPUT, sheet_name=SHEET, engine="openpyxl")
df = df[df["P_sex"].notna() & (df["P_sex"] > 0)].copy()

# Genomic ordering (autosomes, then X/Y if ever present).
def chr_order(x):
    s = str(x).replace("chr", "")
    if s == "X": return 23
    if s == "Y": return 24
    try: return int(float(s))
    except Exception: return 99

df["_chr_order"] = df["Chromosome"].map(chr_order)
df = df.sort_values(["_chr_order", "Gene_start_GRCh37", "Human_gene"]).reset_index(drop=True)
df["plot_index"] = np.arange(len(df))
df["minus_log10_psex"] = -np.log10(df["P_sex"].clip(lower=np.nextafter(0, 1)))

n_testable = len(df)
nominal_y = -np.log10(0.05)
bonf_y = -np.log10(0.05 / n_testable)

fig, ax = plt.subplots(figsize=(9.0, 4.6))
ax.scatter(df["plot_index"], df["minus_log10_psex"], s=18, alpha=0.72)
ax.axhline(nominal_y, linestyle="--", linewidth=1.0, label="Nominal P = 0.05")
ax.axhline(bonf_y, linestyle="--", linewidth=1.2, label=f"Bonferroni threshold (n = {n_testable})")

hit = df[df["Human_gene"] == "EPB41L1"]
if not hit.empty:
    r = hit.iloc[0]
    ax.scatter([r["plot_index"]], [r["minus_log10_psex"]], s=62, zorder=5)
    ax.annotate(
        "EPB41L1", (r["plot_index"], r["minus_log10_psex"]),
        xytext=(8, 7), textcoords="offset points", fontweight="bold"
    )

ax.set_xlabel("Candidate genes ordered by genomic location")
ax.set_ylabel(r"$-\log_{10}(P_{sex})$")
ax.set_title("Sex heterogeneity across mouse-derived human candidate genes")
ax.legend(frameon=False, loc="upper left")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()

for ext in ("png", "pdf", "svg"):
    kwargs = {"dpi": 600} if ext == "png" else {}
    fig.savefig(OUTDIR / f"Figure_2_candidate_gene_sex_heterogeneity.{ext}", bbox_inches="tight", **kwargs)

print(f"Saved Figure 2 to {OUTDIR.resolve()}")
