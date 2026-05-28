import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

summary_file = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod/TandemMod_prediction_summary/TandemMod_prediction_sample_modification_summary.tsv"

outdir = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod/TandemMod_prediction_summary/plots"
os.makedirs(outdir, exist_ok=True)

df = pd.read_csv(summary_file, sep="\t")

# Make sure sample order is correct
sample_order = ["A1", "A2", "A3", "H1", "H2", "H3"]
mod_order = ["m6A", "m1A", "A_I", "m5C", "hm5C", "m7G", "G_I", "psU"]

df["sample"] = pd.Categorical(df["sample"], categories=sample_order, ordered=True)
df["modification"] = pd.Categorical(df["modification"], categories=mod_order, ordered=True)

# -----------------------------
# Plot 1: number of prediction rows
# -----------------------------
plt.figure(figsize=(10, 5))
sns.barplot(
    data=df,
    x="modification",
    y="n_rows",
    hue="sample",
    order=mod_order,
    hue_order=sample_order
)
plt.ylabel("Number of read-level predictions")
plt.xlabel("Modification type")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(f"{outdir}/QC_n_read_level_predictions.pdf")
plt.savefig(f"{outdir}/QC_n_read_level_predictions.png", dpi=300)
plt.close()

# -----------------------------
# Plot 2: unique transcriptomic sites
# -----------------------------
plt.figure(figsize=(10, 5))
sns.barplot(
    data=df,
    x="modification",
    y="n_unique_transcript_sites",
    hue="sample",
    order=mod_order,
    hue_order=sample_order
)
plt.ylabel("Number of unique transcriptomic sites")
plt.xlabel("Modification type")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(f"{outdir}/QC_n_unique_sites.pdf")
plt.savefig(f"{outdir}/QC_n_unique_sites.png", dpi=300)
plt.close()

# -----------------------------
# Plot 3: median probability
# -----------------------------
plt.figure(figsize=(8, 5))
sns.boxplot(
    data=df,
    x="modification",
    y="prob_median",
    order=mod_order,
    showfliers=False
)
sns.stripplot(
    data=df,
    x="modification",
    y="prob_median",
    order=mod_order,
    color="black",
    size=4,
    jitter=True
)
plt.ylabel("Median prediction probability")
plt.xlabel("Modification type")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(f"{outdir}/QC_median_probability.pdf")
plt.savefig(f"{outdir}/QC_median_probability.png", dpi=300)
plt.close()

# -----------------------------
# Plot 4: high-confidence modified fraction
# -----------------------------
metric = "modified_fraction_total_prob_ge_0p9"

plt.figure(figsize=(8, 5))
sns.boxplot(
    data=df,
    x="modification",
    y=metric,
    hue="condition",
    order=mod_order,
    showfliers=False
)
sns.stripplot(
    data=df,
    x="modification",
    y=metric,
    hue="condition",
    order=mod_order,
    dodge=True,
    color="black",
    size=3,
    jitter=True
)
plt.ylabel("Fraction of high-confidence modified reads\n(probability ≥ 0.9)")
plt.xlabel("Modification type")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(f"{outdir}/QC_modified_fraction_prob_ge_0p9.pdf")
plt.savefig(f"{outdir}/QC_modified_fraction_prob_ge_0p9.png", dpi=300)
plt.close()

# -----------------------------
# Plot 5: heatmap of modified fraction
# -----------------------------
heat = df.pivot_table(
    index="modification",
    columns="sample",
    values=metric
)

plt.figure(figsize=(8, 5))
sns.heatmap(
    heat.loc[mod_order, sample_order],
    cmap="viridis",
    annot=True,
    fmt=".3f"
)
plt.xlabel("Sample")
plt.ylabel("Modification type")
plt.tight_layout()
plt.savefig(f"{outdir}/QC_modified_fraction_heatmap_prob_ge_0p9.pdf")
plt.savefig(f"{outdir}/QC_modified_fraction_heatmap_prob_ge_0p9.png", dpi=300)
plt.close()
