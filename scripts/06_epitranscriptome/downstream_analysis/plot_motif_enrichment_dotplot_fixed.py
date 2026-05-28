#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl


# ============================================================
# User settings
# ============================================================

INFILE = "../TandemMod_site_level_differential_analysis_p0p9/motif_analysis/Motif_enrichment_by_modification_and_direction.tsv"
OUTDIR = "./motif_plots"
os.makedirs(OUTDIR, exist_ok=True)

FDR_CUTOFF = 0.05
MIN_ODDS_RATIO = 1

# keep top N significant enriched motifs per direction
TOP_N_PER_DIRECTION = 20

DIRECTION_ORDER = [
    "Higher in aposymbiotic",
    "Higher in symbiotic",
]

MOD_LABEL_MAP = {
    "m6A": "m6A",
    "m1A": "m1A",
    "A_I": "A-to-I",
    "m5C": "m5C",
    "hm5C": "hm5C",
    "m7G": "m7G",
    "G_I": "G-to-I",
    "psU": "psU",
}


# ============================================================
# Load and prepare data
# ============================================================

df = pd.read_csv(INFILE, sep="\t")

plot_df = df[
    (df["fdr"] < FDR_CUTOFF) &
    (df["odds_ratio"] > MIN_ODDS_RATIO)
].copy()

plot_df["modification_label"] = plot_df["modification"].map(MOD_LABEL_MAP)
plot_df["modification_label"] = plot_df["modification_label"].fillna(plot_df["modification"])

plot_df["motif_label"] = (
    plot_df["modification_label"].astype(str) +
    " | " +
    plot_df["motif"].astype(str)
)

plot_df["minus_log10_fdr_plot"] = plot_df["minus_log10_fdr"].replace([np.inf, -np.inf], np.nan)
plot_df["minus_log10_fdr_plot"] = plot_df["minus_log10_fdr_plot"].clip(upper=50)

plot_df["log2_odds_ratio_plot"] = plot_df["log2_odds_ratio"].replace([np.inf, -np.inf], np.nan)

plot_df = plot_df.dropna(
    subset=[
        "direction_label",
        "motif_label",
        "log2_odds_ratio_plot",
        "minus_log10_fdr_plot",
        "target_fraction"
    ]
).copy()

if TOP_N_PER_DIRECTION is not None:
    plot_df = (
        plot_df
        .sort_values(
            ["direction_label", "fdr", "log2_odds_ratio_plot"],
            ascending=[True, True, False]
        )
        .groupby("direction_label", group_keys=False)
        .head(TOP_N_PER_DIRECTION)
        .copy()
    )

plot_df["direction_label"] = pd.Categorical(
    plot_df["direction_label"],
    categories=DIRECTION_ORDER,
    ordered=True
)

if plot_df.shape[0] == 0:
    raise ValueError("No significant enriched motifs found under current filters.")


# ============================================================
# Figure layout
# ============================================================

plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 16,
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 11,
})

fig = plt.figure(figsize=(15.5, 8.5))

# two plot panels + colorbar + size legend axis
gs = fig.add_gridspec(
    nrows=1,
    ncols=4,
    width_ratios=[1.35, 1.35, 0.10, 0.40],
    wspace=0.55
)

ax_left = fig.add_subplot(gs[0, 0])
ax_right = fig.add_subplot(gs[0, 1], sharex=ax_left)
cax = fig.add_subplot(gs[0, 2])
lax = fig.add_subplot(gs[0, 3])

axes = [ax_left, ax_right]


# ============================================================
# Color and size scaling
# ============================================================

vmin = plot_df["target_fraction"].min()
vmax = plot_df["target_fraction"].max()
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
cmap = plt.cm.YlOrBr

size_min = 40
size_max = 260
sig_min = plot_df["minus_log10_fdr_plot"].min()
sig_max = plot_df["minus_log10_fdr_plot"].max()

def scale_size(x):
    if pd.isna(x):
        return size_min
    if sig_max == sig_min:
        return (size_min + size_max) / 2
    return size_min + (x - sig_min) / (sig_max - sig_min) * (size_max - size_min)


# ============================================================
# Draw each panel
# ============================================================

xmax = plot_df["log2_odds_ratio_plot"].max()

for ax, direction in zip(axes, DIRECTION_ORDER):
    sub = plot_df[plot_df["direction_label"] == direction].copy()
    sub = sub.sort_values("log2_odds_ratio_plot", ascending=True).reset_index(drop=True)

    y_pos = np.arange(sub.shape[0])
    sizes = sub["minus_log10_fdr_plot"].apply(scale_size)

    ax.scatter(
        sub["log2_odds_ratio_plot"],
        y_pos,
        s=sizes,
        c=sub["target_fraction"],
        cmap=cmap,
        norm=norm,
        edgecolor="black",
        linewidth=0.4,
        alpha=0.95
    )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(sub["motif_label"])

    ax.axvline(0, color="black", linewidth=0.9)
    ax.grid(axis="x", alpha=0.25)
    ax.grid(axis="y", alpha=0.12)

    ax.set_title(direction)
    ax.set_xlabel("log2 odds ratio")
    ax.set_xlim(-0.2, xmax + 0.5)

    if ax is ax_left:
        ax.set_ylabel("Motif")
    else:
        ax.set_ylabel("")


# ============================================================
# Colorbar
# ============================================================

sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])

cbar = fig.colorbar(sm, cax=cax)
cbar.set_label("Fraction of target sites\nwith motif", rotation=90, labelpad=12)


# ============================================================
# Size legend
# ============================================================

lax.axis("off")

candidate_vals = [5, 10, 20, 30, 40, 50]
legend_vals = [v for v in candidate_vals if sig_min <= v <= sig_max]

if len(legend_vals) == 0:
    legend_vals = sorted(list(set([
        round(sig_min, 1),
        round((sig_min + sig_max) / 2, 1),
        round(sig_max, 1)
    ])))

handles = [
    plt.scatter(
        [], [], s=scale_size(v),
        facecolor="lightgray",
        edgecolor="black",
        linewidth=0.4
    )
    for v in legend_vals
]

labels = [f"{v:g}" for v in legend_vals]

lax.legend(
    handles,
    labels,
    title="-log10(FDR)",
    loc="center left",
    frameon=False,
    borderaxespad=0.0,
    handletextpad=1.0,
    labelspacing=1.2
)


# ============================================================
# Final touches
# ============================================================

fig.suptitle(
    "Enriched motif contexts among differentially modified sites",
    fontsize=18,
    y=0.98
)

fig.subplots_adjust(
    left=0.24,
    right=0.97,
    top=0.88,
    bottom=0.10
)

out_png = os.path.join(OUTDIR, "Top_enriched_motifs_dotplot_fixed.png")
out_pdf = os.path.join(OUTDIR, "Top_enriched_motifs_dotplot_fixed.pdf")

plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.savefig(out_pdf, bbox_inches="tight")
plt.close()

print("Saved:")
print(out_png)
print(out_pdf)
