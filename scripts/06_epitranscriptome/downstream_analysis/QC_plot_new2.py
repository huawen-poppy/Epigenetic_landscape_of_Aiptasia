import os
import glob
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


# =========================
# Input / output
# =========================

summary_file = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod/TandemMod_prediction_summary/TandemMod_prediction_sample_modification_summary.tsv"

pred_dir = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod/TandemMod_prediction"

outdir = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod/TandemMod_prediction_summary/plots_updated"
os.makedirs(outdir, exist_ok=True)


# =========================
# Basic settings
# =========================

sample_order = ["A1", "A2", "A3", "H1", "H2", "H3"]

condition_map = {
    "A1": "Aposymbiotic",
    "A2": "Aposymbiotic",
    "A3": "Aposymbiotic",
    "H1": "Symbiotic",
    "H2": "Symbiotic",
    "H3": "Symbiotic",
}

sample_conditions = [condition_map[s] for s in sample_order]

colors = [
    "#D8B88A" if cond == "Aposymbiotic" else "#F28E2B"
    for cond in sample_conditions
]

condition_palette = {
    "Aposymbiotic": "#D8B88A",
    "Symbiotic": "#F28E2B",
}

mod_order = ["m6A", "m1A", "A_I", "m5C", "hm5C", "m7G", "G_I", "psU"]

feature_map = {
    "m6A": "DRACH",
    "m1A": "NNANN",
    "A_I": "NNANN",
    "m5C": "NNCNN",
    "hm5C": "NNCNN",
    "m7G": "NNGNN",
    "G_I": "NNGNN",
    "psU": "NNTNN",
}


# =========================
# Helper functions
# =========================

def p_to_label(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return ""


def add_significance_labels(
    ax,
    data,
    x_col,
    y_col,
    group_col,
    x_order,
    group1="Aposymbiotic",
    group2="Symbiotic",
    y_offset_ratio=0.06,
):
    """
    Add significance labels for group1 vs group2 within each x category.
    Only significant labels are shown.
    P-values are BH-adjusted across modification types.
    """

    raw_pvals = []
    mods_for_test = []

    for mod in x_order:
        sub = data[data[x_col] == mod]

        vals1 = sub.loc[sub[group_col] == group1, y_col].dropna()
        vals2 = sub.loc[sub[group_col] == group2, y_col].dropna()

        if len(vals1) >= 2 and len(vals2) >= 2:
            try:
                stat, p = mannwhitneyu(vals1, vals2, alternative="two-sided")
            except Exception:
                p = np.nan
        else:
            p = np.nan

        raw_pvals.append(p)
        mods_for_test.append(mod)

    valid = ~pd.isna(raw_pvals)

    adjusted_pvals = np.full(len(raw_pvals), np.nan)

    if np.sum(valid) > 0:
        adjusted_pvals[valid] = multipletests(
            np.array(raw_pvals)[valid],
            method="fdr_bh"
        )[1]

    y_min = data[y_col].min()
    y_max = data[y_col].max()
    y_range = y_max - y_min if y_max > y_min else 1

    for i, mod in enumerate(mods_for_test):
        padj = adjusted_pvals[i]
        label = p_to_label(padj) if not np.isnan(padj) else ""

        if label == "":
            continue

        sub = data[data[x_col] == mod]
        local_ymax = sub[y_col].max()

        y = local_ymax + y_offset_ratio * y_range

        # bracket location around each x category
        x = x_order.index(mod)
        x1 = x - 0.22
        x2 = x + 0.22

        ax.plot(
            [x1, x1, x2, x2],
            [y, y + 0.01 * y_range, y + 0.01 * y_range, y],
            lw=1,
            color="black"
        )

        ax.text(
            x,
            y + 0.015 * y_range,
            label,
            ha="center",
            va="bottom",
            fontsize=12,
            color="black"
        )

    # add some space above the plot
    ax.set_ylim(y_min, y_max + 0.18 * y_range)


# =========================
# Load summary table
# =========================

df = pd.read_csv(summary_file, sep="\t")

df["condition"] = df["sample"].map(condition_map)

df["sample"] = pd.Categorical(
    df["sample"],
    categories=sample_order,
    ordered=True
)

df["condition"] = pd.Categorical(
    df["condition"],
    categories=["Aposymbiotic", "Symbiotic"],
    ordered=True
)

df["modification"] = pd.Categorical(
    df["modification"],
    categories=mod_order,
    ordered=True
)


# =========================
# Plot 1: number of read-level predictions
# =========================

plt.figure(figsize=(10, 5))

ax = sns.barplot(
    data=df,
    x="modification",
    y="n_rows",
    hue="condition",
    order=mod_order,
    palette=condition_palette,
    errorbar="se"
)

sns.stripplot(
    data=df,
    x="modification",
    y="n_rows",
    hue="condition",
    order=mod_order,
    dodge=True,
    color="black",
    size=3,
    jitter=True,
    ax=ax
)

plt.ylabel("Number of read-level predictions")
plt.xlabel("Modification type")
plt.xticks(rotation=45, ha="right")

# Remove duplicated legend from barplot + stripplot
handles, labels = ax.get_legend_handles_labels()
ax.legend(
    handles[:2],
    labels[:2],
    title="Condition",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    borderaxespad=0
)

plt.tight_layout()
plt.savefig(f"{outdir}/QC_n_read_level_predictions.png", dpi=300, bbox_inches="tight")
plt.savefig(f"{outdir}/QC_n_read_level_predictions.pdf", bbox_inches="tight")
plt.close()


# =========================
# Plot 2: unique transcriptomic sites
# =========================

plt.figure(figsize=(10, 5))

ax = sns.barplot(
    data=df,
    x="modification",
    y="n_unique_transcript_sites",
    hue="condition",
    order=mod_order,
    palette=condition_palette,
    errorbar="se"
)

sns.stripplot(
    data=df,
    x="modification",
    y="n_unique_transcript_sites",
    hue="condition",
    order=mod_order,
    dodge=True,
    color="black",
    size=3,
    jitter=True,
    ax=ax
)

plt.ylabel("Number of unique transcriptomic sites")
plt.xlabel("Modification type")
plt.xticks(rotation=45, ha="right")

handles, labels = ax.get_legend_handles_labels()
ax.legend(
    handles[:2],
    labels[:2],
    title="Condition",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    borderaxespad=0
)

plt.tight_layout()
plt.savefig(f"{outdir}/QC_n_unique_sites.png", dpi=300, bbox_inches="tight")
plt.savefig(f"{outdir}/QC_n_unique_sites.pdf", bbox_inches="tight")
plt.close()


# =========================
# Plot 3: high-confidence modified fraction with statistics
# =========================
# Change this if your column name is different.
# For probability >= 0.9, it should usually be:
metric = "modified_fraction_total_prob_ge_0p9"

plot_df = df.dropna(subset=[metric]).copy()

plt.figure(figsize=(10, 5))

ax = sns.boxplot(
    data=plot_df,
    x="modification",
    y=metric,
    hue="condition",
    order=mod_order,
    palette=condition_palette,
    showfliers=False,
    linewidth=1
)

sns.stripplot(
    data=plot_df,
    x="modification",
    y=metric,
    hue="condition",
    order=mod_order,
    dodge=True,
    color="black",
    size=4,
    jitter=True,
    ax=ax
)

add_significance_labels(
    ax=ax,
    data=plot_df,
    x_col="modification",
    y_col=metric,
    group_col="condition",
    x_order=mod_order
)

plt.ylabel("Fraction of high-confidence modified reads\n(probability ≥ 0.9)")
plt.xlabel("Modification type")
plt.xticks(rotation=45, ha="right")

handles, labels = ax.get_legend_handles_labels()
ax.legend(
    handles[:2],
    labels[:2],
    title="Condition",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    borderaxespad=0
)

plt.tight_layout()
plt.savefig(f"{outdir}/QC_modification_fraction_prob_ge_0p9.png", dpi=300, bbox_inches="tight")
plt.savefig(f"{outdir}/QC_modification_fraction_prob_ge_0p9.pdf", bbox_inches="tight")
plt.close()


# =========================
# Plot 4: heatmap of modified fraction
# =========================

heat = df.pivot_table(
    index="modification",
    columns="sample",
    values=metric
)

plt.figure(figsize=(8, 5))

ax = sns.heatmap(
    heat.loc[mod_order, sample_order],
    cmap="YlOrBr",
    annot=True,
    fmt=".3f",
    linewidths=0.4,
    linecolor="white"
)

plt.xlabel("Sample")
plt.ylabel("Modification type")
plt.tight_layout()
plt.savefig(f"{outdir}/QC_modification_fraction_heatmap_prob_ge_0p9.png", dpi=300, bbox_inches="tight")
plt.savefig(f"{outdir}/QC_modification_fraction_heatmap_prob_ge_0p9.pdf", bbox_inches="tight")
plt.close()

# =========================
# Plot 5: probability density by modification
# =========================

N_PER_FILE = 200000  # reduce to 50000 if files are very large

density_records = []

for sample in sample_order:
    condition = condition_map[sample]

    for mod in mod_order:
        feature_type = feature_map[mod]

        pred_file = os.path.join(
            pred_dir,
            f"{sample}_{mod}_{feature_type}_prediction.tsv"
        )

        if not os.path.exists(pred_file):
            print(f"Missing file, skip: {pred_file}")
            continue

        print(f"Reading probability values from: {pred_file}")

        try:
            # input has no header:
            # transcript_id, site, motif, read_id, label, probability
            tmp = pd.read_csv(
                pred_file,
                sep="\t",
                header=None,
                usecols=[5],
                names=["probability"]
            )

            tmp["probability"] = pd.to_numeric(
                tmp["probability"],
                errors="coerce"
            )

            tmp = tmp.dropna()

            if len(tmp) > N_PER_FILE:
                tmp = tmp.sample(
                    n=N_PER_FILE,
                    random_state=123
                )

            tmp["sample"] = sample
            tmp["condition"] = condition
            tmp["modification"] = mod

            density_records.append(tmp)

        except Exception as e:
            print(f"Failed to read {pred_file}")
            print(e)


density_df = pd.concat(density_records, ignore_index=True)

density_df["condition"] = pd.Categorical(
    density_df["condition"],
    categories=["Aposymbiotic", "Symbiotic"],
    ordered=True
)

density_df["modification"] = pd.Categorical(
    density_df["modification"],
    categories=mod_order,
    ordered=True
)

g = sns.FacetGrid(
    density_df,
    col="modification",
    col_order=mod_order,
    hue="condition",
    hue_order=["Aposymbiotic", "Symbiotic"],
    palette=condition_palette,
    col_wrap=4,
    height=2.4,
    aspect=1.25,
    sharex=True,
    sharey=False
)

g.map_dataframe(
    sns.kdeplot,
    x="probability",
    fill=True,
    alpha=0.35,
    linewidth=1.2,
    common_norm=False
)

g.set_axis_labels("Prediction probability", "Density")
g.set_titles("{col_name}")

for ax in g.axes.flatten():
    ax.set_xlim(0, 1)
    ax.axvline(0.9, color="black", linestyle="--", linewidth=0.8)
    ax.tick_params(axis="x", rotation=0)

g.add_legend(
    title="Condition",
    bbox_to_anchor=(1.02, 0.5),
    loc="center left",
    borderaxespad=0
)

plt.tight_layout()
plt.savefig(f"{outdir}/QC_prediction_probability_density_by_modification.png", dpi=300, bbox_inches="tight")
plt.savefig(f"{outdir}/QC_prediction_probability_density_by_modification.pdf", bbox_inches="tight")
plt.close()
