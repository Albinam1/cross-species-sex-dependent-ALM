from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REGIONAL_FILE = Path("05_EPB41L1_REGIONAL_SEX_ANALYSIS.xlsx")
REGIONAL_SHEET = "S5_All_region_variants"
GTEX_FILE = Path("07_EPB41L1_GTEx_significant_eQTL.xlsx")
GTEX_SHEET = "S2_Tissue_summary"
SBEQTL_FILE = Path("08_EPB41L1_GTEx_sex_biased_eQTL.xlsx")
OUTDIR = Path("figures")
OUTDIR.mkdir(exist_ok=True)

SELECTED = ["rs113092336", "rs6060632", "rs532201406", "rs1006296"]
LABELS = {
    "Esophagus_Mucosa": "Esophagus mucosa",
    "Artery_Tibial": "Tibial artery",
    "Artery_Aorta": "Aorta",
    "Nerve_Tibial": "Tibial nerve",
    "Cells_Cultured_fibroblasts": "Cultured fibroblasts",
    "Muscle_Skeletal": "Skeletal muscle",
    "Testis": "Testis",
    "Brain_Spinal_cord_cervical_c-1": "Spinal cord (cervical C1)",
    "Adipose_Subcutaneous": "Subcutaneous adipose",
    "Esophagus_Gastroesophageal_Junction": "Gastroesophageal junction",
    "Liver": "Liver",
    "Skin_Sun_Exposed_Lower_leg": "Sun-exposed skin (lower leg)",
    "Minor_Salivary_Gland": "Minor salivary gland",
    "Heart_Atrial_Appendage": "Heart atrial appendage",
}

for f in (REGIONAL_FILE, GTEX_FILE):
    if not f.exists():
        raise FileNotFoundError(f"Missing input: {f.resolve()}")

reg = pd.read_excel(REGIONAL_FILE, sheet_name=REGIONAL_SHEET, engine="openpyxl")
reg = reg[reg["rsID"].isin(SELECTED)].drop_duplicates("rsID").set_index("rsID")
missing = [x for x in SELECTED if x not in reg.index]
if missing:
    raise ValueError(f"Selected variants missing from regional analysis: {missing}")
reg = reg.loc[SELECTED].reset_index()

gtex = pd.read_excel(GTEX_FILE, sheet_name=GTEX_SHEET, engine="openpyxl")
gtex = gtex.sort_values("N_significant_eQTLs", ascending=False).reset_index(drop=True)
gtex["Label"] = gtex["Tissue"].map(LABELS).fillna(gtex["Tissue"].astype(str).str.replace("_", " "))

# Exact totals from the GTEx output workbook when available.
summary = pd.read_excel(GTEX_FILE, sheet_name="S0_Summary", engine="openpyxl")
summary_map = dict(zip(summary["Metric"].astype(str), summary["Value"]))
records_total = int(summary_map.get("Significant eQTL records", gtex["N_significant_eQTLs"].sum()))
tissues_total = int(summary_map.get("Tissues with significant EPB41L1 eQTLs", len(gtex)))
unique_evariants = int(summary_map.get("Unique significant eVariants", 0))

sb_records = sb_sig = sb_muscle = None
if SBEQTL_FILE.exists():
    sb = pd.read_excel(SBEQTL_FILE, sheet_name="S0_Summary", engine="openpyxl")
    sb_map = dict(zip(sb["Metric"].astype(str), sb["Value"]))
    sb_records = int(sb_map.get("Total EPB41L1 sb-eQTL records", 0))
    sb_sig = int(sb_map.get("EPB41L1 sb-eQTL records with FDR/q < 0.05", 0))
    sb_muscle = int(sb_map.get("EPB41L1 skeletal-muscle sb-eQTL records", 0))

# Colorblind-safe sex colors; neutral bars for tissue counts.
MALE = "#0072B2"
FEMALE = "#D55E00"
BAR = "#B8B8B8"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 10.5,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8.5,
})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 6.1), gridspec_kw={"width_ratios": [1.0, 1.18]})
fig.subplots_adjust(wspace=0.34, bottom=0.16, top=0.90)

# Panel A: forest plot
y = np.arange(len(reg))
offset = 0.13
male_ci = 1.96 * reg["SE_Male"].to_numpy(float)
female_ci = 1.96 * reg["SE_Female"].to_numpy(float)

ax1.errorbar(reg["Beta_Male"], y - offset, xerr=male_ci,
             fmt="o", capsize=3, markersize=6, linewidth=1.2,
             color=MALE, ecolor=MALE, label="Male")
ax1.errorbar(reg["Beta_Female"], y + offset, xerr=female_ci,
             fmt="s", capsize=3, markersize=5.7, linewidth=1.2,
             color=FEMALE, ecolor=FEMALE, label="Female")
ax1.axvline(0, color="0.35", linestyle="--", linewidth=1.0)
ax1.set_yticks(y)
ax1.set_yticklabels(reg["rsID"])
ax1.invert_yaxis()
ax1.set_xlabel("Effect on appendicular lean mass, β (95% CI)")
ax1.set_title("A  Sex-stratified ALM effects at selected EPB41L1-region variants", loc="left", fontweight="bold")
ax1.legend(frameon=False, loc="lower right")
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

# Panel B: tissue eQTL counts. Counts are not effect magnitudes.
y2 = np.arange(len(gtex))
bars = ax2.barh(y2, gtex["N_significant_eQTLs"], color=BAR, edgecolor="0.25", linewidth=0.45)
ax2.set_yticks(y2)
ax2.set_yticklabels(gtex["Label"])
ax2.invert_yaxis()
ax2.set_xlabel("Number of significant EPB41L1 cis-eQTL records")
ax2.set_title("B  Significant EPB41L1 cis-eQTL records across GTEx v8 tissues", loc="left", fontweight="bold")
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

# Highlight skeletal muscle using hatch rather than a stronger color scale.
muscle_rows = np.where(gtex["Tissue"].eq("Muscle_Skeletal"))[0]
if len(muscle_rows):
    i = int(muscle_rows[0])
    bars[i].set_facecolor("white")
    bars[i].set_edgecolor("black")
    bars[i].set_linewidth(1.0)
    bars[i].set_hatch("///")

mx = float(gtex["N_significant_eQTLs"].max())
for b, v in zip(bars, gtex["N_significant_eQTLs"]):
    ax2.text(float(v) + mx * 0.012, b.get_y() + b.get_height()/2, str(int(v)), va="center", fontsize=7.5)
ax2.set_xlim(0, mx * 1.18)

summary_lines = [
    "GTEx v8 summary",
    f"Tissues with significant cis-eQTLs: {tissues_total}",
    f"Unique significant eVariants: {unique_evariants:,}",
    f"Significant cis-eQTL records: {records_total:,}",
]
if len(muscle_rows):
    summary_lines.append(f"Skeletal muscle records: {int(gtex.loc[muscle_rows[0], 'N_significant_eQTLs'])}")
if sb_records is not None:
    summary_lines += [
        "",
        "Sex-biased cis-eQTL analysis",
        f"EPB41L1 records: {sb_records}",
        f"Significant after correction: {sb_sig}",
        f"Skeletal-muscle records: {sb_muscle}",
    ]
ax2.text(0.98, 0.03, "\n".join(summary_lines), transform=ax2.transAxes,
         ha="right", va="bottom", fontsize=7.2,
         bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="0.5", linewidth=0.8))

fig.text(0.57, 0.045,
         "Bar length represents the number of significant cis-eQTL records and should not be interpreted as tissue-specific effect magnitude.",
         ha="center", va="bottom", fontsize=7.5)

for ext in ("png", "pdf", "svg"):
    kwargs = {"dpi": 600} if ext == "png" else {}
    fig.savefig(OUTDIR / f"Figure_4_EPB41L1_genetic_regulatory_characterization.{ext}", bbox_inches="tight", **kwargs)

print(f"Saved Figure 4 to {OUTDIR.resolve()}")
