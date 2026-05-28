#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl


# ============================================================
# User settings
# ============================================================

INFILE = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod/TandemMod_site_level_differential_analysis_p0p9/motif_analysis/Motif_enrichment_by_modification_and_direction.tsv"

OUTDIR = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod/TandemMod_site_level_differential_analysis_p0p9/motif_analysis/plots_by_modification"
os.makedirs(OUTDIR, exist_ok=True)

FDR_CUTOFF = 0.05
MIN_ODDS_RATIO = 1

# for all significant motifs, set to None.
TOP_N_PER_DIRECTION = None

MODS = ["m6A", "m1A", "A_I", "m5C", "hm5C", "m7G", "G_I", "psU"]

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

DIRECTION_ORDER = [
    "Higher in aposymbiotic",
    "Higher in symbiotic",
]


# ============================================================
# Helper functions
# ============================================================

def scale_size(x, sig_min, sig_max, size_min=45, size_max=260):
    if pd.isna(x):
        return size_min
    if sig_max == sig_min:
        return (size_min + size_max) / 2
    return size_min + (x - sig_min) / (sig_max - sig_min) * (size_max - size_min)


def plot_one_modification(mod, df, outdir):
    mod_label = MOD_LABEL_MAP.get(mod, mod)

    sub_all = df[df["modification"] == mod].copy()

    # Keep only significantly enriched motifs
    sub_all = sub_all[
        (sub_all["fdr"] < FDR_CUTOFF) &
        (sub_all["odds_ratio"] > MIN_ODDS_RATIO)
    ].copy()

    if sub_all.shape[0] == 0:
        print(f"No significant enriched motifs for {mod_label}; skipping plot.")
        return

    sub_all["motif_label"] = sub_all["motif"].astype(str)

    sub_all["minus_log10_fdr_plot"] = (
        sub_all["minus_log10_fdr"]
        .replace([np.inf, -np.inf], np.nan)
        .clip(upper=50)
    )

    sub_all["log2_odds_ratio_plot"] = (
        sub_all["log2_odds_ratio"]
        .replace([np.inf, -np.inf], np.nan)
    )

    sub_all = sub_all.dropna(
        subset=[
            "direction_label",
            "motif_label",
            "log2_odds_ratio_plot",
            "minus_log10_fdr_plot",
            "target_fraction",
        ]
    )

    if TOP_N_PER_DIRECTION is not None:
        sub_all = (
            sub_all
            .sort_values(
                ["direction_label", "fdr", "log2_odds_ratio_plot"],
                ascending=[True, True, False]
            )
            .groupby("direction_label", group_keys=False)
            .head(TOP_N_PER_DIRECTION)
            .copy()
        )

    if sub_all.shape[0] == 0:
        print(f"No plottable motifs for {mod_label}; skipping plot.")
        return

    # Dynamic height based on the larger number of motifs in either direction
    max_n_motifs = (
        sub_all
        .groupby("direction_label")
        .size()
        .reindex(DIRECTION_ORDER)
        .fillna(0)
        .max()
    )

    fig_height = max(4.5, 0.32 * max_n_motifs + 1.8)

    fig = plt.figure(figsize=(11.5, fig_height))

    gs = fig.add_gridspec(
        nrows=1,
        ncols=4,
        width_ratios=[1.35, 1.35, 0.10, 0.36],
        wspace=0.55
    )

    ax_left = fig.add_subplot(gs[0, 0])
    ax_right = fig.add_subplot(gs[0, 1], sharex=ax_left)
    cax = fig.add_subplot(gs[0, 2])
    lax = fig.add_subplot(gs[0, 3])

    axes = [ax_left, ax_right]

    # Color scale for target motif fraction
    vmin = sub_all["target_fraction"].min()
    vmax = sub_all["target_fraction"].max()

    if vmin == vmax:
        vmin = 0
        vmax = max(vmax, 0.01)

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.cm.YlOrBr

    # Size scale for -log10(FDR)
    sig_min = sub_all["minus_log10_fdr_plot"].min()
    sig_max = sub_all["minus_log10_fdr_plot"].max()

    xmax = sub_all["log2_odds_ratio_plot"].max()

    for ax, direction in zip(axes, DIRECTION_ORDER):
        sub = sub_all[sub_all["direction_label"] == direction].copy()

        if sub.shape[0] == 0:
            ax.text(
                0.5,
                0.5,
                "No significant\nenriched motifs",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=12
            )
            ax.set_title(direction)
            ax.set_xlabel("log2 odds ratio")
            ax.set_yticks([])
            ax.axvline(0, color="black", linewidth=0.9)
            ax.set_xlim(-0.2, max(1, xmax + 0.5))
            continue

        # strongest enrichment at top
        sub = sub.sort_values("log2_odds_ratio_plot", ascending=True).reset_index(drop=True)

        y_pos = np.arange(sub.shape[0])
        sizes = sub["minus_log10_fdr_plot"].apply(
            lambda x: scale_size(x, sig_min, sig_max)
        )

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

    # Colorbar
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label(
        "Fraction of target sites\nwith motif",
        rotation=90,
        labelpad=12
    )

    # Size legend
    lax.axis("off")

    candidate_vals = [5, 10, 20, 30, 40, 50]
    legend_vals = [v for v in candidate_vals if sig_min <= v <= sig_max]

    if len(legend_vals) == 0:
        legend_vals = sorted(list(set([
            round(sig_min, 1),
            round((sig_min + sig_max) / 2, 1),
            round(sig_max, 1),
        ])))

    handles = [
        plt.scatter(
            [],
            [],
            s=scale_size(v, sig_min, sig_max),
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

    fig.suptitle(
        f"{mod_label}: enriched motif contexts among differentially modified sites",
        fontsize=16,
        y=0.98
    )

    fig.subplots_adjust(
        left=0.18,
        right=0.96,
        top=0.88,
        bottom=0.12
    )

    out_png = os.path.join(
        outdir,
        f"{mod}_motif_enrichment_dotplot.png"
    )

    out_pdf = os.path.join(
        outdir,
        f"{mod}_motif_enrichment_dotplot.pdf"
    )

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()

    print(f"Saved {mod_label}:")
    print(out_png)
    print(out_pdf)


# ============================================================
# Main
# ============================================================

plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 13,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 10,
})

df = pd.read_csv(INFILE, sep="\t")

required_cols = [
    "modification",
    "direction_label",
    "motif",
    "target_fraction",
    "odds_ratio",
    "fdr",
    "minus_log10_fdr",
    "log2_odds_ratio",
]

missing = [c for c in required_cols if c not in df.columns]
if len(missing) > 0:
    raise ValueError(f"Missing columns in input file: {missing}")

for mod in MODS:
    plot_one_modification(mod, df, OUTDIR)

print("Finished generating per-modification motif enrichment plots.")
