#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import mannwhitneyu, ks_2samp, fisher_exact
from statsmodels.stats.multitest import multipletests


# ============================================================
# Settings
# ============================================================

DM_ANALYSIS_DIR = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod/TandemMod_site_level_differential_analysis_p0p9"

# From your previous position analysis
POSITION_FILE = f"{DM_ANALYSIS_DIR}/modification_position_analysis/Compact_modification_site_position_annotation.tsv.gz"

# From your motif enrichment analysis
MOTIF_ENRICHMENT_FILE = f"{DM_ANALYSIS_DIR}/motif_analysis/Motif_enrichment_by_modification_and_direction.tsv"

OUT_DIR = f"{DM_ANALYSIS_DIR}/enriched_motif_position_analysis"
os.makedirs(OUT_DIR, exist_ok=True)

PLOT_DIR = f"{OUT_DIR}/plots"
os.makedirs(PLOT_DIR, exist_ok=True)


# Enriched motif definition
ENRICH_FDR_CUTOFF = 0.05
MIN_ODDS_RATIO = 1.0

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
    "higher_in_aposymbiotic",
    "higher_in_symbiotic",
]

DIRECTION_LABELS = {
    "higher_in_aposymbiotic": "Higher in aposymbiotic",
    "higher_in_symbiotic": "Higher in symbiotic",
    "not_significant": "Not significant",
}

DIRECTION_COLORS = {
    "Higher in aposymbiotic": "#D8B88A",
    "Higher in symbiotic": "#F28E2B",
}

REGION_ORDER = [
    "five_prime_UTR",
    "CDS",
    "three_prime_UTR",
    "unannotated_region",
    "not_in_TransDecoder",
]

REGION_LABELS = {
    "five_prime_UTR": "5′ UTR",
    "CDS": "CDS",
    "three_prime_UTR": "3′ UTR",
    "unannotated_region": "Unannotated region",
    "not_in_TransDecoder": "Not in TransDecoder",
}

REGION_COLORS = {
    "5′ UTR": "#4E79A7",
    "CDS": "#59A14F",
    "3′ UTR": "#F28E2B",
    "Unannotated region": "#9C755F",
    "Not in TransDecoder": "#BDBDBD",
}

N_POSITION_BINS = 20


# ============================================================
# Helper functions
# ============================================================

def normalize_direction(x):
    if pd.isna(x):
        return np.nan

    x = str(x).strip()

    direction_map = {
        "higher_in_aposymbiotic": "higher_in_aposymbiotic",
        "Higher in aposymbiotic": "higher_in_aposymbiotic",
        "aposymbiotic": "higher_in_aposymbiotic",

        "higher_in_symbiotic": "higher_in_symbiotic",
        "Higher in symbiotic": "higher_in_symbiotic",
        "symbiotic": "higher_in_symbiotic",
    }

    return direction_map.get(x, x)


def add_position_bins(df, n_bins=20):
    df = df.copy()

    bins = np.linspace(0, 1, n_bins + 1)

    labels = [
        f"{bins[i]:.2f}-{bins[i + 1]:.2f}"
        for i in range(n_bins)
    ]

    df["relative_position_bin"] = pd.cut(
        df["relative_position"],
        bins=bins,
        labels=labels,
        include_lowest=True
    )

    df["relative_position_bin_mid"] = df["relative_position_bin"].apply(
        lambda x: np.nan if pd.isna(x)
        else (
            float(str(x).split("-")[0]) +
            float(str(x).split("-")[1])
        ) / 2
    )

    return df


def safe_mannwhitney(x, y):
    x = pd.Series(x).dropna()
    y = pd.Series(y).dropna()

    if len(x) < 5 or len(y) < 5:
        return np.nan, np.nan

    stat, p = mannwhitneyu(x, y, alternative="two-sided")
    return stat, p


def safe_ks(x, y):
    x = pd.Series(x).dropna()
    y = pd.Series(y).dropna()

    if len(x) < 5 or len(y) < 5:
        return np.nan, np.nan

    stat, p = ks_2samp(x, y, alternative="two-sided")
    return stat, p


def fisher_region_direction_test(df):
    """
    For each modification and transcript region, test whether that region is
    differently represented between enriched-motif sites higher in symbiotic
    and enriched-motif sites higher in aposymbiotic.

    Table:
                         in region    not in region
    higher symbiotic         a              b
    higher aposymbiotic      c              d
    """

    records = []

    for mod in MODS:
        mod_df = df[df["modification"] == mod].copy()

        if mod_df.shape[0] == 0:
            continue

        sym = mod_df[mod_df["direction_norm"] == "higher_in_symbiotic"]
        apo = mod_df[mod_df["direction_norm"] == "higher_in_aposymbiotic"]

        if sym.shape[0] == 0 or apo.shape[0] == 0:
            continue

        for region in REGION_ORDER:
            a = (sym["transcript_region"] == region).sum()
            b = (sym["transcript_region"] != region).sum()
            c = (apo["transcript_region"] == region).sum()
            d = (apo["transcript_region"] != region).sum()

            try:
                odds_ratio, p_value = fisher_exact(
                    [[a, b], [c, d]],
                    alternative="two-sided"
                )
            except Exception:
                odds_ratio, p_value = np.nan, np.nan

            records.append(
                {
                    "modification": mod,
                    "region": region,
                    "region_label": REGION_LABELS.get(region, region),

                    "symbiotic_in_region": int(a),
                    "symbiotic_not_region": int(b),
                    "aposymbiotic_in_region": int(c),
                    "aposymbiotic_not_region": int(d),

                    "symbiotic_region_fraction": a / (a + b)
                    if (a + b) > 0 else np.nan,

                    "aposymbiotic_region_fraction": c / (c + d)
                    if (c + d) > 0 else np.nan,

                    "odds_ratio_sym_vs_apo": odds_ratio,
                    "p_value": p_value,
                }
            )

    out = pd.DataFrame(records)

    out["fdr"] = np.nan

    valid = out["p_value"].notna()

    if valid.sum() > 0:
        out.loc[valid, "fdr"] = multipletests(
            out.loc[valid, "p_value"],
            method="fdr_bh"
        )[1]

    out["log2_odds_ratio_sym_vs_apo"] = np.log2(
        out["odds_ratio_sym_vs_apo"].replace(0, np.nan)
    )

    out["minus_log10_fdr"] = -np.log10(
        out["fdr"].replace(0, np.nan)
    )

    return out


# ============================================================
# Load data
# ============================================================

print("Loading position annotation:")
print(POSITION_FILE)

pos = pd.read_csv(POSITION_FILE, sep="\t")

print("Position table shape:", pos.shape)

print("Loading motif enrichment:")
print(MOTIF_ENRICHMENT_FILE)

enrich = pd.read_csv(MOTIF_ENRICHMENT_FILE, sep="\t")

print("Motif enrichment table shape:", enrich.shape)


# ============================================================
# Prepare enriched motif table
# ============================================================

if "direction" in enrich.columns:
    enrich["direction_norm"] = enrich["direction"].apply(normalize_direction)
elif "direction_label" in enrich.columns:
    enrich["direction_norm"] = enrich["direction_label"].apply(normalize_direction)
else:
    raise ValueError("Motif enrichment file must contain direction or direction_label column.")

required_enrich_cols = [
    "modification",
    "direction_norm",
    "motif",
    "fdr",
    "odds_ratio",
]

missing_enrich = [
    c for c in required_enrich_cols
    if c not in enrich.columns
]

if missing_enrich:
    raise ValueError(f"Missing required enrichment columns: {missing_enrich}")

enrich["motif"] = enrich["motif"].astype(str).str.upper().str.strip()
enrich["modification"] = enrich["modification"].astype(str).str.strip()

sig_enrich = enrich[
    (enrich["fdr"] < ENRICH_FDR_CUTOFF) &
    (enrich["odds_ratio"] > MIN_ODDS_RATIO) &
    (enrich["direction_norm"].isin(DIRECTION_ORDER))
].copy()

sig_enrich_keys = sig_enrich[
    ["modification", "direction_norm", "motif"]
].drop_duplicates()

sig_enrich_keys_out = os.path.join(
    OUT_DIR,
    "Significant_enriched_motif_keys.tsv"
)

sig_enrich_keys.to_csv(
    sig_enrich_keys_out,
    sep="\t",
    index=False
)

print("Number of significant enriched motif keys:", sig_enrich_keys.shape[0])
print("Saved enriched motif keys:", sig_enrich_keys_out)


# ============================================================
# Prepare position table
# ============================================================

required_pos_cols = [
    "site_id",
    "modification",
    "motif",
    "direction",
    "relative_position",
    "transcript_region",
    "transcript_region_label",
]

missing_pos = [
    c for c in required_pos_cols
    if c not in pos.columns
]

if missing_pos:
    raise ValueError(f"Missing required position columns: {missing_pos}")

pos["motif"] = pos["motif"].astype(str).str.upper().str.strip()
pos["modification"] = pos["modification"].astype(str).str.strip()
pos["direction_norm"] = pos["direction"].apply(normalize_direction)

# Keep differential sites only
pos_dm = pos[pos["direction_norm"].isin(DIRECTION_ORDER)].copy()

print("Differential sites in position table:", pos_dm.shape)


# ============================================================
# Subset to sites whose motif is significantly enriched
# for the same modification and direction
# ============================================================

enriched_pos = pos_dm.merge(
    sig_enrich_keys,
    on=["modification", "direction_norm", "motif"],
    how="inner"
)

enriched_pos = enriched_pos.drop_duplicates("site_id").copy()

enriched_pos["direction_label"] = enriched_pos["direction_norm"].map(
    DIRECTION_LABELS
)

enriched_pos["modification_label"] = enriched_pos["modification"].map(
    MOD_LABEL_MAP
)

enriched_pos["modification_label"] = enriched_pos["modification_label"].fillna(
    enriched_pos["modification"]
)

enriched_pos = add_position_bins(
    enriched_pos,
    n_bins=N_POSITION_BINS
)

enriched_pos_out = os.path.join(
    OUT_DIR,
    "Differential_modification_sites_with_significant_enriched_motifs_position.tsv.gz"
)

enriched_pos.to_csv(
    enriched_pos_out,
    sep="\t",
    index=False,
    compression="gzip"
)

print("Enriched-motif differential sites:", enriched_pos.shape)
print("Saved:", enriched_pos_out)


# ============================================================
# Summary 1. Counts by modification and direction
# ============================================================

counts = (
    enriched_pos
    .groupby(["modification", "direction_norm", "direction_label"], observed=False)
    .agg(
        n_sites=("site_id", "nunique"),
        n_transcripts=("transcriptome_id", "nunique"),
        n_motifs=("motif", "nunique"),
    )
    .reset_index()
)

counts_out = os.path.join(
    OUT_DIR,
    "Enriched_motif_site_counts_by_modification_direction.tsv"
)

counts.to_csv(counts_out, sep="\t", index=False)

print("Saved counts:", counts_out)


# ============================================================
# Summary 2. Region composition
# ============================================================

region_comp = (
    enriched_pos
    .groupby(
        [
            "modification",
            "direction_norm",
            "direction_label",
            "transcript_region",
            "transcript_region_label",
        ],
        observed=False
    )
    .size()
    .reset_index(name="n_sites")
)

region_comp["total_sites"] = (
    region_comp
    .groupby(["modification", "direction_norm"], observed=False)["n_sites"]
    .transform("sum")
)

region_comp["fraction"] = region_comp["n_sites"] / region_comp["total_sites"]

region_comp_out = os.path.join(
    OUT_DIR,
    "Enriched_motif_site_region_composition.tsv"
)

region_comp.to_csv(region_comp_out, sep="\t", index=False)

print("Saved region composition:", region_comp_out)


# ============================================================
# Summary 3. Relative-position comparison:
# higher symbiotic vs higher aposymbiotic
# ============================================================

relpos_records = []

for mod in MODS:
    sub = enriched_pos[
        (enriched_pos["modification"] == mod) &
        (enriched_pos["relative_position"].notna())
    ].copy()

    if sub.shape[0] == 0:
        continue

    apo = sub.loc[
        sub["direction_norm"] == "higher_in_aposymbiotic",
        "relative_position"
    ].dropna()

    sym = sub.loc[
        sub["direction_norm"] == "higher_in_symbiotic",
        "relative_position"
    ].dropna()

    mw_stat, mw_p = safe_mannwhitney(sym, apo)
    ks_stat, ks_p = safe_ks(sym, apo)

    relpos_records.append(
        {
            "modification": mod,
            "n_aposymbiotic": len(apo),
            "n_symbiotic": len(sym),

            "median_relative_position_aposymbiotic": apo.median()
            if len(apo) > 0 else np.nan,

            "median_relative_position_symbiotic": sym.median()
            if len(sym) > 0 else np.nan,

            "mean_relative_position_aposymbiotic": apo.mean()
            if len(apo) > 0 else np.nan,

            "mean_relative_position_symbiotic": sym.mean()
            if len(sym) > 0 else np.nan,

            "delta_median_sym_minus_apo": (
                sym.median() - apo.median()
                if len(sym) > 0 and len(apo) > 0 else np.nan
            ),

            "mannwhitney_u": mw_stat,
            "mannwhitney_p": mw_p,
            "ks_statistic": ks_stat,
            "ks_p": ks_p,
        }
    )

relpos_tests = pd.DataFrame(relpos_records)

for pcol in ["mannwhitney_p", "ks_p"]:
    fdr_col = pcol.replace("_p", "_fdr")
    relpos_tests[fdr_col] = np.nan

    valid = relpos_tests[pcol].notna()

    if valid.sum() > 0:
        relpos_tests.loc[valid, fdr_col] = multipletests(
            relpos_tests.loc[valid, pcol],
            method="fdr_bh"
        )[1]

relpos_tests_out = os.path.join(
    OUT_DIR,
    "Enriched_motif_relative_position_sym_vs_apo_tests.tsv"
)

relpos_tests.to_csv(
    relpos_tests_out,
    sep="\t",
    index=False
)

print("Saved relative-position tests:", relpos_tests_out)


# ============================================================
# Summary 4. Region difference between directions
# ============================================================

region_direction_tests = fisher_region_direction_test(enriched_pos)

region_direction_tests_out = os.path.join(
    OUT_DIR,
    "Enriched_motif_region_sym_vs_apo_Fisher_tests.tsv"
)

region_direction_tests.to_csv(
    region_direction_tests_out,
    sep="\t",
    index=False
)

print("Saved region direction tests:", region_direction_tests_out)


# ============================================================
# Summary 5. Binned relative-position distribution
# ============================================================

binned = (
    enriched_pos
    .groupby(
        [
            "modification",
            "direction_norm",
            "direction_label",
            "relative_position_bin",
            "relative_position_bin_mid",
        ],
        observed=False
    )
    .size()
    .reset_index(name="n_sites")
)

binned["total_sites"] = (
    binned
    .groupby(["modification", "direction_norm"], observed=False)["n_sites"]
    .transform("sum")
)

binned["fraction"] = binned["n_sites"] / binned["total_sites"]

binned_out = os.path.join(
    OUT_DIR,
    "Enriched_motif_binned_relative_position_distribution.tsv"
)

binned.to_csv(binned_out, sep="\t", index=False)

print("Saved binned distribution:", binned_out)


# ============================================================
# Plots
# ============================================================

sns.set_theme(style="whitegrid", font_scale=1.2)


# ------------------------------------------------------------
# Plot 1. Relative-position density, enriched motifs only
# ------------------------------------------------------------

density_df = enriched_pos[
    enriched_pos["relative_position"].notna()
].copy()

density_df["direction_label"] = pd.Categorical(
    density_df["direction_label"],
    categories=[
        "Higher in aposymbiotic",
        "Higher in symbiotic",
    ],
    ordered=True
)

g = sns.FacetGrid(
    density_df,
    col="modification",
    col_order=MODS,
    col_wrap=4,
    hue="direction_label",
    hue_order=[
        "Higher in aposymbiotic",
        "Higher in symbiotic",
    ],
    palette=DIRECTION_COLORS,
    height=3,
    aspect=1.1,
    sharex=True,
    sharey=False
)

g.map_dataframe(
    sns.kdeplot,
    x="relative_position",
    linewidth=1.5,
    fill=False,
    common_norm=False,
    clip=(0, 1),
    cut=0
)

g.set_axis_labels("Relative transcript position", "Density")
g.set_titles("{col_name}")

for ax in g.axes.flatten():
    ax.set_xlim(0, 1)
    ax.axvline(0.5, color="black", linestyle="--", linewidth=0.6, alpha=0.6)

g.add_legend(
    title="Enriched-motif site group",
    bbox_to_anchor=(1.02, 0.5),
    loc="center left",
    frameon=False
)

plt.tight_layout(rect=[0, 0, 0.88, 1])

density_png = os.path.join(
    PLOT_DIR,
    "Enriched_motif_relative_position_density_by_modification.png"
)

density_pdf = os.path.join(
    PLOT_DIR,
    "Enriched_motif_relative_position_density_by_modification.pdf"
)

plt.savefig(density_png, dpi=300, bbox_inches="tight")
plt.savefig(density_pdf, bbox_inches="tight")
plt.close()


# ------------------------------------------------------------
# Plot 2. Binned relative-position line plot
# ------------------------------------------------------------

binned_plot = binned.copy()

binned_plot["direction_label"] = pd.Categorical(
    binned_plot["direction_label"],
    categories=[
        "Higher in aposymbiotic",
        "Higher in symbiotic",
    ],
    ordered=True
)

g = sns.FacetGrid(
    binned_plot,
    col="modification",
    col_order=MODS,
    col_wrap=4,
    hue="direction_label",
    hue_order=[
        "Higher in aposymbiotic",
        "Higher in symbiotic",
    ],
    palette=DIRECTION_COLORS,
    height=3,
    aspect=1.1,
    sharex=True,
    sharey=False
)

g.map_dataframe(
    sns.lineplot,
    x="relative_position_bin_mid",
    y="fraction",
    marker="o",
    linewidth=1.5
)

g.set_axis_labels("Relative transcript position", "Fraction of enriched-motif sites")
g.set_titles("{col_name}")

for ax in g.axes.flatten():
    ax.set_xlim(0, 1)

g.add_legend(
    title="Enriched-motif site group",
    bbox_to_anchor=(1.02, 0.5),
    loc="center left",
    frameon=False
)

plt.tight_layout(rect=[0, 0, 0.88, 1])

binned_png = os.path.join(
    PLOT_DIR,
    "Enriched_motif_binned_relative_position_by_modification.png"
)

binned_pdf = os.path.join(
    PLOT_DIR,
    "Enriched_motif_binned_relative_position_by_modification.pdf"
)

plt.savefig(binned_png, dpi=300, bbox_inches="tight")
plt.savefig(binned_pdf, bbox_inches="tight")
plt.close()


# ------------------------------------------------------------
# Plot 3. Region composition stacked bar
# ------------------------------------------------------------

region_label_order = [REGION_LABELS[x] for x in REGION_ORDER]

plot_region = region_comp.copy()

plot_region["transcript_region_label"] = pd.Categorical(
    plot_region["transcript_region_label"],
    categories=region_label_order,
    ordered=True
)

plot_region["x_label"] = (
    plot_region["modification"].astype(str) +
    "\n" +
    plot_region["direction_label"].str.replace("Higher in ", "High in ")
)

pivot_region = plot_region.pivot_table(
    index=["modification", "direction_label", "x_label"],
    columns="transcript_region_label",
    values="fraction",
    aggfunc="sum",
    fill_value=0,
    observed=False
)

pivot_region = pivot_region.reindex(columns=region_label_order, fill_value=0)

fig, ax = plt.subplots(figsize=(12, 5.5))

x = np.arange(pivot_region.shape[0])
bottom = np.zeros(pivot_region.shape[0])

for region_label in region_label_order:
    vals = pivot_region[region_label].values

    ax.bar(
        x,
        vals,
        bottom=bottom,
        label=region_label,
        color=REGION_COLORS.get(region_label, "gray"),
        edgecolor="white",
        linewidth=0.4
    )

    bottom += vals

ax.set_xticks(x)
ax.set_xticklabels(
    [idx[2] for idx in pivot_region.index],
    rotation=45,
    ha="right"
)

ax.set_ylabel("Fraction of enriched-motif sites")
ax.set_xlabel("Modification type and direction")

ax.legend(
    title="Transcript region",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    frameon=False
)

plt.tight_layout()

region_png = os.path.join(
    PLOT_DIR,
    "Enriched_motif_region_composition_stacked_bar.png"
)

region_pdf = os.path.join(
    PLOT_DIR,
    "Enriched_motif_region_composition_stacked_bar.pdf"
)

plt.savefig(region_png, dpi=300, bbox_inches="tight")
plt.savefig(region_pdf, bbox_inches="tight")
plt.close()


# ------------------------------------------------------------
# Plot 4. Region difference between directions
# ------------------------------------------------------------

region_test_plot = region_direction_tests.copy()

region_test_plot = region_test_plot[
    (region_test_plot["fdr"] < 0.05) &
    (region_test_plot["odds_ratio_sym_vs_apo"] > 1)
].copy()

if region_test_plot.shape[0] > 0:
    region_test_plot["mod_region"] = (
        region_test_plot["modification"].astype(str) +
        " | " +
        region_test_plot["region_label"].astype(str)
    )

    region_test_plot["minus_log10_fdr_plot"] = (
        region_test_plot["minus_log10_fdr"]
        .replace([np.inf, -np.inf], np.nan)
        .clip(upper=50)
    )

    region_test_plot = region_test_plot.sort_values(
        "log2_odds_ratio_sym_vs_apo",
        ascending=True
    )

    plt.figure(figsize=(7, max(4, 0.35 * region_test_plot.shape[0] + 1.5)))

    ax = plt.scatter(
        region_test_plot["log2_odds_ratio_sym_vs_apo"],
        region_test_plot["mod_region"],
        s=region_test_plot["minus_log10_fdr_plot"] * 8 + 25,
        c=region_test_plot["symbiotic_region_fraction"],
        cmap="YlOrBr",
        edgecolor="black",
        linewidth=0.4
    )

    plt.axvline(0, color="black", linestyle="--", linewidth=0.8)

    cbar = plt.colorbar(ax)
    cbar.set_label("Fraction in region\namong symbiotic-higher sites")

    plt.xlabel("log2 odds ratio\nsymbiotic-higher vs aposymbiotic-higher")
    plt.ylabel("Modification | region")

    plt.tight_layout()

    region_test_png = os.path.join(
        PLOT_DIR,
        "Enriched_motif_region_sym_vs_apo_enrichment_dotplot.png"
    )

    region_test_pdf = os.path.join(
        PLOT_DIR,
        "Enriched_motif_region_sym_vs_apo_enrichment_dotplot.pdf"
    )

    plt.savefig(region_test_png, dpi=300, bbox_inches="tight")
    plt.savefig(region_test_pdf, bbox_inches="tight")
    plt.close()

else:
    print("No significant region differences between directions for plotting.")


# ============================================================
# Final summary
# ============================================================

summary = (
    enriched_pos
    .groupby("modification", observed=False)
    .agg(
        n_enriched_motif_sites=("site_id", "nunique"),
        n_transcripts=("transcriptome_id", "nunique"),
        n_enriched_motifs=("motif", "nunique"),
        median_relative_position=("relative_position", "median"),
        mean_relative_position=("relative_position", "mean"),
    )
    .reset_index()
)

summary_out = os.path.join(
    OUT_DIR,
    "Enriched_motif_position_summary_by_modification.tsv"
)

summary.to_csv(summary_out, sep="\t", index=False)

print("\nFinished enriched-motif position analysis.")
print("\nMain outputs:")
print(enriched_pos_out)
print(counts_out)
print(region_comp_out)
print(relpos_tests_out)
print(region_direction_tests_out)
print(binned_out)
print(summary_out)

print("\nPlots saved to:")
print(PLOT_DIR)
