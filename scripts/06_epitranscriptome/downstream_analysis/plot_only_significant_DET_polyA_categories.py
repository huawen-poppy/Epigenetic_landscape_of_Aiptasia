#!/usr/bin/env python

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# ============================================================
# User settings
# ============================================================

DET_FILE = "../TandemMod_site_level_differential_analysis_p0p9/integration_with_DET/Differential_modification_DET_merged_transcript_modification_level.tsv"
POLYA_FILE = "../TandemMod_site_level_differential_analysis_p0p9/integration_with_polyA/Differential_modification_polyA_merged_transcript_modification_level.tsv"

OUTDIR = "updated_overlap_category_plots"
os.makedirs(OUTDIR, exist_ok=True)

MODS = ["m6A", "m1A", "A_I", "m5C", "hm5C", "m7G", "G_I", "psU"]

COND_COLORS = {
    "DET up in aposymbiotic": "#D8B88A",
    "DET up in symbiotic": "#F28E2B",
    "Longer poly(A) in aposymbiotic": "#D8B88A",
    "Longer poly(A) in symbiotic": "#F28E2B",
}

sns.set_theme(style="whitegrid", font_scale=1.3)


# ============================================================
# Helper function
# ============================================================

def make_two_category_plot(
    infile,
    category_col,
    keep_categories,
    category_order,
    y_label_count,
    y_label_fraction,
    title_prefix,
    out_prefix,
):
    df = pd.read_csv(infile, sep="\t")

    required_cols = ["transcript_id", "modification", "n_dm_sites", category_col]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {infile}: {missing}")

    # Keep only differentially modified transcript-modification pairs
    df_dm = df[df["n_dm_sites"] > 0].copy()

    # Keep only the two biologically meaningful categories
    df_sig = df_dm[df_dm[category_col].isin(keep_categories)].copy()

    # Count unique transcripts per modification and category
    summary = (
        df_sig
        .groupby(["modification", category_col], observed=False)
        .agg(n_transcripts=("transcript_id", "nunique"))
        .reset_index()
    )

    # Add missing modification/category combinations as zeros
    full_index = pd.MultiIndex.from_product(
        [MODS, category_order],
        names=["modification", category_col]
    )

    summary = (
        summary
        .set_index(["modification", category_col])
        .reindex(full_index, fill_value=0)
        .reset_index()
    )

    # Fraction among the two significant categories only
    summary["total_significant_overlap"] = (
        summary
        .groupby("modification", observed=False)["n_transcripts"]
        .transform("sum")
    )

    summary["fraction_among_significant_overlap"] = (
        summary["n_transcripts"] / summary["total_significant_overlap"]
    )

    summary.loc[
        summary["total_significant_overlap"] == 0,
        "fraction_among_significant_overlap"
    ] = 0

    summary_out = os.path.join(OUTDIR, f"{out_prefix}_two_category_summary.tsv")
    summary.to_csv(summary_out, sep="\t", index=False)

    # ========================================================
    # Plot 1: counts
    # ========================================================

    plt.figure(figsize=(8, 5))

    ax = sns.barplot(
        data=summary,
        x="modification",
        y="n_transcripts",
        hue=category_col,
        order=MODS,
        hue_order=category_order,
        palette=COND_COLORS,
        edgecolor="white",
        linewidth=0.6,
    )

    ax.set_xlabel("Modification type")
    ax.set_ylabel(y_label_count)
    ax.set_title(f"{title_prefix}: transcript counts")
    ax.legend(
        title="Category",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
    )

    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    plt.savefig(
        os.path.join(OUTDIR, f"{out_prefix}_counts_only_two_categories.png"),
        dpi=300,
        bbox_inches="tight"
    )
    plt.savefig(
        os.path.join(OUTDIR, f"{out_prefix}_counts_only_two_categories.pdf"),
        bbox_inches="tight"
    )
    plt.close()

    # ========================================================
    # Plot 2: fraction among significant categories only
    # ========================================================

    plt.figure(figsize=(8, 5))

    ax = sns.barplot(
        data=summary,
        x="modification",
        y="fraction_among_significant_overlap",
        hue=category_col,
        order=MODS,
        hue_order=category_order,
        palette=COND_COLORS,
        edgecolor="white",
        linewidth=0.6,
    )

    ax.set_xlabel("Modification type")
    ax.set_ylabel(y_label_fraction)
    ax.set_title(f"{title_prefix}: fraction among significant overlap")
    ax.set_ylim(0, 1)

    ax.legend(
        title="Category",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
    )

    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    plt.savefig(
        os.path.join(OUTDIR, f"{out_prefix}_fraction_only_two_categories.png"),
        dpi=300,
        bbox_inches="tight"
    )
    plt.savefig(
        os.path.join(OUTDIR, f"{out_prefix}_fraction_only_two_categories.pdf"),
        bbox_inches="tight"
    )
    plt.close()

    print(f"Saved summary and plots for {out_prefix}")
    print(summary_out)


# ============================================================
# 1. DET category plot: only up in apo / up in symbiotic
# ============================================================

make_two_category_plot(
    infile=DET_FILE,
    category_col="DET_category",
    keep_categories=[
        "DET up in aposymbiotic",
        "DET up in symbiotic",
    ],
    category_order=[
        "DET up in aposymbiotic",
        "DET up in symbiotic",
    ],
    y_label_count="Number of differentially modified transcripts",
    y_label_fraction="Fraction among DM ∩ DET transcripts",
    title_prefix="Differentially modified transcripts overlapping DETs",
    out_prefix="DM_transcripts_DET_overlap",
)


# ============================================================
# 2. poly(A) category plot: only longer in apo / longer in symbiotic
# ============================================================

make_two_category_plot(
    infile=POLYA_FILE,
    category_col="polyA_category",
    keep_categories=[
        "Longer poly(A) in aposymbiotic",
        "Longer poly(A) in symbiotic",
    ],
    category_order=[
        "Longer poly(A) in aposymbiotic",
        "Longer poly(A) in symbiotic",
    ],
    y_label_count="Number of differentially modified transcripts",
    y_label_fraction="Fraction among DM ∩ differential-poly(A) transcripts",
    title_prefix="Differentially modified transcripts overlapping differential poly(A)",
    out_prefix="DM_transcripts_polyA_overlap",
)

print("Done.")
