from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

INPUT = Path("05_EPB41L1_REGIONAL_SEX_ANALYSIS.xlsx")
SHEET = "S5_All_region_variants"
OUTDIR = Path("figures")
OUTDIR.mkdir(exist_ok=True)

GENE_START = 34_679_426
GENE_END = 34_820_721
KEY_VARIANTS = ["rs113092336", "rs532201406", "rs1006296", "rs6060632"]

if not INPUT.exists():
    raise FileNotFoundError(f"Missing input: {INPUT.resolve()}")

df = pd.read_excel(INPUT, sheet_name=SHEET, engine="openpyxl")
df = df[df["P_sex"].notna() & (df["P_sex"] > 0)].copy()
df["x_mb"] = df["POS"] / 1e6
df["minus_log10_psex"] = -np.log10(df["P_sex"].clip(lower=np.nextafter(0, 1)))

n_region = len(df)
bonf_y = -np.log10(0.05 / n_region)

fig, ax = plt.subplots(figsize=(9.0, 4.8))
ax.scatter(df["x_mb"], df["minus_log10_psex"], s=15, alpha=0.62)
ax.axhspan(0, 0, alpha=0)  # keeps export deterministic across backends
ax.axvspan(GENE_START / 1e6, GENE_END / 1e6, alpha=0.12)
ax.axhline(bonf_y, linestyle="--", linewidth=1.2,
           label=f"Regional Bonferroni threshold (n = {n_region})")

for rsid in KEY_VARIANTS:
    h = df[df["rsID"] == rsid]
    if h.empty:
        continue
    r = h.iloc[0]
    ax.scatter([r["x_mb"]], [r["minus_log10_psex"]], s=58, zorder=5)
    offset = {
        "rs113092336": (5, 8),
        "rs532201406": (5, -18),
        "rs1006296": (5, 6),
        "rs6060632": (5, 8),
    }[rsid]
    ax.annotate(rsid, (r["x_mb"], r["minus_log10_psex"]),
                xytext=offset, textcoords="offset points")

ax.text((GENE_START + GENE_END) / 2 / 1e6, 0.12, "EPB41L1",
        ha="center", va="bottom", fontstyle="italic")
ax.set_xlabel("Chromosome 20 position (Mb, GRCh37)")
ax.set_ylabel(r"$-\log_{10}(P_{sex})$")
ax.set_title("Regional sex heterogeneity at the EPB41L1 locus")
ax.legend(frameon=False, loc="upper right")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()

for ext in ("png", "pdf", "svg"):
    kwargs = {"dpi": 600} if ext == "png" else {}
    fig.savefig(OUTDIR / f"Figure_3_EPB41L1_regional_sex_heterogeneity.{ext}", bbox_inches="tight", **kwargs)

print(f"Saved Figure 3 to {OUTDIR.resolve()}")
