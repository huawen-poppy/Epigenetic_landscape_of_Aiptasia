import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from matplotlib.lines import Line2D


# =========================
# Input / output
# =========================

summary_file = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod/TandemMod_prediction_summary/TandemMod_prediction_sample_modification_summary.tsv"

outdir = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod/TandemMod_prediction_summary/plots"
os.makedirs(outdir, exist_ok=True)

df = pd.read_csv(summary_file, sep="\t")


# =========================
# Settings
# =========================

sample_order = ["A1", "A2", "A3",  "H1", "H2", "H3"]
mod_order = ["m6A", "m1A", "A_I", "m5C", "hm5C", "m7G", "G_I", "psU"]

condition_map = {
    "A1": "Aposymbiotic",
    "A2": "Aposymbiotic",
    "A3": "Aposymbiotic",
    "H1": "Symbiotic",
    "H2": "Symbiotic",
    "H3": "Symbiotic",
}

condition_palette = {
    "Aposymbiotic": "#D8B88A",
    "Symbiotic": "#F28E2B",
}

df["condition_plot"] = df["sample"].map(condition_map)

df["sample"] = pd.Categorical(df["sample"], categories=sample_order, ordered=True)
df["modification"] = pd.Categorical(df["modification"], categories=mod_order, ordered=True)
df["condition_plot"] = pd.Categorical(
    df["condition_plot"],
    categories=["Aposymbiotic", "Symbiotic"],
    ordered=True
)

# Your requested color logic, if you need it directly
sample_conditions = [condition_map[s] for s in sample_order]
colors = [
    "#D8B88A" if cond == "Aposymbiotic" else "#F28E2B"
    for cond in sample_conditions
]


# =========================
# Helper functions
# =========================

def p_to_star(p):
    if pd.isna(p):
        return "n.s."
    elif p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "n.s."


def add_stat_annotations(ax, data, x_col, y_col, group_col, x_order):
    """
    Add apo vs symbiotic statistical comparison for each modification.
    Uses Welch's t-test, then BH-FDR correction across modifications.
    """

    pvals = []
    records = []

    for mod in x_order:
        sub = data[data[x_col] == mod]

        apo = sub[sub[group_col] == "Aposymbiotic"][y_col].dropna()
        sym = sub[sub[group_col] == "Symbiotic"][y_col].dropna()

        if len(apo) >= 2 and len(sym) >= 2:
            stat, p = ttest_ind(apo, sym, equal_var=False)
        else:
            p = np.nan

        pvals.append(p)
        records.append({"modification": mod, "pvalue": p})

    valid = ~pd.isna(pvals)

    fdrs = np.full(len(pvals), np.nan)

    if np.sum(valid) > 0:
        _, corrected, _, _ = multipletests(
            np.array(pvals)[valid],
            method="fdr_bh"
        )
        fdrs[valid] = corrected

    # y-axis range
    ymin, ymax = ax.get_ylim()
    yrange = ymax - ymin

    # Add annotation for each modification
    for i, mod in enumerate(x_order):
        sub = data[data[x_col] == mod]
        y_max = sub[y_col].max()

        if pd.isna(y_max):
            continue

        y = y_max + 0.06 * yrange
        h = 0.025 * yrange

        # positions for hue groups inside each x category
        x1 = i - 0.20
        x2 = i + 0.20

        ax.plot(
            [x1, x1, x2, x2],
            [y, y + h, y + h, y],
            lw=1,
            color="black"
        )

        label = p_to_star(fdrs[i])

        ax.text(
            i,
            y + h + 0.01 * yrange,
            label,
            ha="center",
            va="bottom",
            fontsize=11,
            color="black"
        )

        records[i]["FDR"] = fdrs[i]
        records[i]["significance"] = label

    # Expand y-axis so labels are visible
    ax.set_ylim(ymin, ymax + 0.18 * yrange)

    return pd.DataFrame(records)


def move_legend_outside(ax, title="Condition"):
    handles, labels = ax.get_legend_handles_labels()

    # Remove duplicate labels from boxplot + stripplot
    unique = {}
    for h, l in zip(handles, labels):
        if l not in unique:
            unique[l] = h

    ax.legend(
        unique.values(),
        unique.keys(),
        title=title,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
        borderaxespad=0
    )


# =========================
# Plot 1: number of read-level predictions
# =========================

plt.figure(figsize=(10, 5))
ax = sns.barplot(
    data=df,
    x="modification",
    y="n_rows",
    hue="condition_plot",
    order=mod_order,
    palette=condition_palette,
    errorbar="sd"
)

sns.stripplot(
    data=df,
    x="modification",
    y="n_rows",
    hue="condition_plot",
    order=mod_order,
    dodge=True,
    color="black",
    size=3,
    alpha=0.7
)

plt.ylabel("Number of read-level predictions")
plt.xlabel("Modification type")
plt.xticks(rotation=45, ha="right")
move_legend_outside(ax, title="Condition")
plt.tight_layout()
plt.savefig(f"{outdir}/QC_n_read_level_predictions.pdf", bbox_inches="tight")
plt.savefig(f"{outdir}/QC_n_read_level_predictions.png", dpi=300, bbox_inches="tight")
plt.close()


# =========================
# Plot 2: unique transcriptomic sites
# =========================

plt.figure(figsize=(10, 5))
ax = sns.barplot(
    data=df,
    x="modification",
    y="n_unique_transcript_sites",
    hue="condition_plot",
    order=mod_order,
    palette=condition_palette,
    errorbar="sd"
)

sns.stripplot(
    data=df,
    x="modification",
    y="n_unique_transcript_sites",
    hue="condition_plot",
    order=mod_order,
    dodge=True,
    color="black",
    size=3,
    alpha=0.7
)

plt.ylabel("Number of unique transcriptomic sites")
plt.xlabel("Modification type")
plt.xticks(rotation=45, ha="right")
move_legend_outside(ax, title="Condition")
plt.tight_layout()
plt.savefig(f"{outdir}/QC_n_unique_sites.pdf", bbox_inches="tight")
plt.savefig(f"{outdir}/QC_n_unique_sites.png", dpi=300, bbox_inches="tight")
plt.close()


# =========================
# Plot 3: median probability as line plot
# Different modifications in different colors
# =========================

plt.figure(figsize=(10, 5))

mod_palette = dict(zip(mod_order, sns.color_palette("tab10", n_colors=len(mod_order))))

ax = sns.lineplot(
    data=df,
    x="sample",
    y="prob_median",
    hue="modification",
    hue_order=mod_order,
    palette=mod_palette,
    marker="o",
    linewidth=2
)

plt.ylabel("Median prediction probability")
plt.xlabel("Sample")
plt.xticks(rotation=45, ha="right")
plt.ylim(0, 1.02)

ax.legend(
    title="Modification",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    frameon=False,
    borderaxespad=0
)

plt.tight_layout()
plt.savefig(f"{outdir}/QC_median_probability_lineplot.pdf", bbox_inches="tight")
plt.savefig(f"{outdir}/QC_median_probability_lineplot.png", dpi=300, bbox_inches="tight")
plt.close()


# =========================
# Plot 4: high-confidence modified fraction
# With statistical comparison between conditions
# =========================

metric = "modified_fraction_total_prob_ge_0p9"

plt.figure(figsize=(10, 5))

ax = sns.boxplot(
    data=df,
    x="modification",
    y=metric,
    hue="condition_plot",
    order=mod_order,
    palette=condition_palette,
    showfliers=False,
    width=0.65
)

sns.stripplot(
    data=df,
    x="modification",
    y=metric,
    hue="condition_plot",
    order=mod_order,
    dodge=True,
    color="black",
    size=4,
    alpha=0.8
)

plt.ylabel("Fraction of high-confidence modified reads\n(probability ≥ 0.9)")
plt.xlabel("Modification type")
plt.xticks(rotation=45, ha="right")

stat_df = add_stat_annotations(
    ax=ax,
    data=df,
    x_col="modification",
    y_col=metric,
    group_col="condition_plot",
    x_order=mod_order
)

move_legend_outside(ax, title="Condition")

plt.tight_layout()
plt.savefig(f"{outdir}/QC_modification_fraction_prob_ge_0p9.pdf", bbox_inches="tight")
plt.savefig(f"{outdir}/QC_modification_fraction_prob_ge_0p9.png", dpi=300, bbox_inches="tight")
plt.close()

# Save statistical results
stat_df.to_csv(
    f"{outdir}/QC_modification_fraction_prob_ge_0p9_condition_comparison.tsv",
    sep="\t",
    index=False
)


# =========================
# Plot 5: heatmap of modified fraction
# =========================

heat = df.pivot_table(
    index="modification",
    columns="sample",
    values=metric
)

plt.figure(figsize=(8, 5))
ax = sns.heatmap(
    heat.loc[mod_order, sample_order],
    cmap="viridis",
    annot=True,
    fmt=".3f"
)

plt.xlabel("Sample")
plt.ylabel("Modification type")
plt.tight_layout()
plt.savefig(f"{outdir}/QC_modified_fraction_heatmap_prob_ge_0p9.pdf", bbox_inches="tight")
plt.savefig(f"{outdir}/QC_modified_fraction_heatmap_prob_ge_0p9.png", dpi=300, bbox_inches="tight")
plt.close()
