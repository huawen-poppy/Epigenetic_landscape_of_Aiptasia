#!/usr/bin/env python

import os
import glob
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


# ============================================================
# Settings
# ============================================================

TANDEMMOD_DIR = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod"

DM_ANALYSIS_DIR = f"{TANDEMMOD_DIR}/TandemMod_site_level_differential_analysis_p0p9"

OUT_DIR = f"{DM_ANALYSIS_DIR}/integration_with_DET"
os.makedirs(OUT_DIR, exist_ok=True)

PLOT_DIR = f"{OUT_DIR}/plots"
os.makedirs(PLOT_DIR, exist_ok=True)

# Differential modification result file
DM_FILE = glob.glob(
    os.path.join(
        DM_ANALYSIS_DIR,
        "TandemMod_differential_modification_results_p0.9*.tsv"
    )
)[0]

# =========================
# EDIT THIS: your DET file
# =========================
DET_FILE = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/isoquant/outputs/final_all_ref/OUT/diffexpr-results-isoquant-genes-6samples_update_from_transcript_model_transcript_level_raw.csv"

# =========================
# EDIT THESE according to your DET table
# =========================
DET_TRANSCRIPT_COL = "transcript_id"
DET_LOG2FC_COL = "log2FoldChange"
DET_FDR_COL = "padj"

# If your DET file already contains only significant DETs, set DET_FILE_IS_ONLY_SIGNIFICANT = True
DET_FILE_IS_ONLY_SIGNIFICANT = False

DET_FDR_CUTOFF = 0.05

# Use 0.0 if you only want FDR-based DET definition.
# Use 1.0 if you want |log2FC| >= 1 as additional threshold.
DET_ABS_LOG2FC_CUTOFF = 1.0

MODS = ["m6A", "m1A", "A_I", "m5C", "hm5C", "m7G", "G_I", "psU"]

DIRECTION_ORDER = [
    "higher_in_aposymbiotic",
    "higher_in_symbiotic",
]

DIRECTION_LABELS = {
    "higher_in_aposymbiotic": "Higher modification in aposymbiotic",
    "higher_in_symbiotic": "Higher modification in symbiotic",
    "not_significant": "Not significant",
    "mixed": "Mixed direction",
}

DIRECTION_COLORS = {
    "Higher modification in aposymbiotic": "#D8B88A",
    "Higher modification in symbiotic": "#F28E2B",
    "Mixed direction": "#7F7F7F",
    "Not significant": "lightgray",
}

DET_COLORS = {
    "DET up in symbiotic": "#F28E2B",
    "DET up in aposymbiotic": "#D8B88A",
    "Not DET": "lightgray",
}


# ============================================================
# Helper functions
# ============================================================

def classify_det(row):
    if DET_FILE_IS_ONLY_SIGNIFICANT:
        is_det = True
    else:
        is_det = (
            row[DET_FDR_COL] < DET_FDR_CUTOFF and
            abs(row[DET_LOG2FC_COL]) >= DET_ABS_LOG2FC_CUTOFF
        )

    if not is_det:
        return "Not DET"

    if row[DET_LOG2FC_COL] > 0:
        return "DET up in symbiotic"
    elif row[DET_LOG2FC_COL] < 0:
        return "DET up in aposymbiotic"
    else:
        return "DET"


def fisher_overlap(table_df, target_col, target_value, det_col="is_DET"):
    """
    Test enrichment of DET among target transcripts.

    Table:
                    DET      non-DET
    target           a          b
    non-target       c          d
    """

    target = table_df[target_col] == target_value

    a = ((target) & (table_df[det_col])).sum()
    b = ((target) & (~table_df[det_col])).sum()
    c = ((~target) & (table_df[det_col])).sum()
    d = ((~target) & (~table_df[det_col])).sum()

    contingency = [[a, b], [c, d]]

    try:
        odds_ratio, p_value = fisher_exact(contingency, alternative="greater")
    except Exception:
        odds_ratio, p_value = np.nan, np.nan

    return {
        "a_target_and_DET": a,
        "b_target_not_DET": b,
        "c_non_target_DET": c,
        "d_non_target_not_DET": d,
        "odds_ratio": odds_ratio,
        "p_value": p_value,
        "target_total": a + b,
        "background_total": a + b + c + d,
        "target_DET_fraction": a / (a + b) if (a + b) > 0 else np.nan,
        "background_DET_fraction": (a + c) / (a + b + c + d) if (a + b + c + d) > 0 else np.nan,
    }


def build_transcript_modification_summary(dm):
    """
    Convert site-level DM results into transcript-modification-level summary.
    One row = transcriptome_id + modification.
    """

    dm = dm.copy()
    dm["is_dm_site"] = dm["direction"].isin(DIRECTION_ORDER)

    records = []

    for (tx, mod), sub in dm.groupby(["transcriptome_id", "modification"], observed=False):
        sig = sub[sub["is_dm_site"]].copy()

        n_sig = sig.shape[0]
        n_higher_sym = (sig["direction"] == "higher_in_symbiotic").sum()
        n_higher_apo = (sig["direction"] == "higher_in_aposymbiotic").sum()

        if n_sig == 0:
            dm_category = "not_significant"
            rep = sub.iloc[sub["delta_mean_mod_rate_sym_minus_apo"].abs().argmax()]
        elif n_higher_sym > 0 and n_higher_apo > 0:
            dm_category = "mixed"
            rep = sig.iloc[sig["delta_mean_mod_rate_sym_minus_apo"].abs().argmax()]
        elif n_higher_sym > 0:
            dm_category = "higher_in_symbiotic"
            rep = sig.iloc[sig["delta_mean_mod_rate_sym_minus_apo"].abs().argmax()]
        else:
            dm_category = "higher_in_aposymbiotic"
            rep = sig.iloc[sig["delta_mean_mod_rate_sym_minus_apo"].abs().argmax()]

        records.append({
            "transcript_id": tx,
            "modification": mod,
            "n_tested_sites": sub.shape[0],
            "n_dm_sites": n_sig,
            "n_dm_sites_higher_symbiotic": n_higher_sym,
            "n_dm_sites_higher_aposymbiotic": n_higher_apo,
            "dm_category": dm_category,
            "dm_category_label": DIRECTION_LABELS[dm_category],
            "representative_site_id": rep["site_id"],
            "representative_site": rep["site"],
            "representative_motif": rep["motif"],
            "representative_delta_mod_rate_sym_minus_apo": rep["delta_mean_mod_rate_sym_minus_apo"],
            "representative_fdr": rep["fdr_by_modification"],
            "min_fdr": sub["fdr_by_modification"].min(),
            "max_abs_delta_mod_rate": sub["delta_mean_mod_rate_sym_minus_apo"].abs().max(),
            "mean_delta_mod_rate_all_sites": sub["delta_mean_mod_rate_sym_minus_apo"].mean(),
        })

    return pd.DataFrame(records)


# ============================================================
# Load data
# ============================================================

print("Loading differential modification results:")
print(DM_FILE)

dm = pd.read_csv(DM_FILE, sep="\t")

print("Loading DET results:")
print(DET_FILE)

det = pd.read_csv(DET_FILE, sep=',',index_col=0)

# Check columns
for col in [DET_TRANSCRIPT_COL, DET_LOG2FC_COL, DET_FDR_COL]:
    if col not in det.columns:
        raise ValueError(f"DET file does not contain required column: {col}")

det = det.copy()
det = det.rename(columns={DET_TRANSCRIPT_COL: "transcript_id"})

det[DET_LOG2FC_COL] = pd.to_numeric(det[DET_LOG2FC_COL], errors="coerce")
det[DET_FDR_COL] = pd.to_numeric(det[DET_FDR_COL], errors="coerce")

det["DET_category"] = det.apply(classify_det, axis=1)
det["is_DET"] = det["DET_category"] != "Not DET"

# If duplicated transcripts exist, keep the row with smallest FDR
det = (
    det
    .sort_values(DET_FDR_COL)
    .drop_duplicates("transcript_id")
    .copy()
)


# ============================================================
# Build transcript-modification summary
# ============================================================

txmod = build_transcript_modification_summary(dm)

txmod_out = os.path.join(
    OUT_DIR,
    "Transcript_modification_summary_for_DET_integration.tsv"
)

txmod.to_csv(txmod_out, sep="\t", index=False)

print("Saved transcript-modification summary:")
print(txmod_out)


# ============================================================
# Merge with DET
# ============================================================

merged = txmod.merge(
    det,
    on="transcript_id",
    how="inner"
)

merged_out = os.path.join(
    OUT_DIR,
    "Differential_modification_DET_merged_transcript_modification_level.tsv"
)

merged.to_csv(merged_out, sep="\t", index=False)

print("Saved merged DM × DET table:")
print(merged_out)


# ============================================================
# Overall transcript-level summary across all modifications
# ============================================================

overall_tx = (
    txmod
    .groupby("transcript_id", as_index=False)
    .agg(
        n_tested_modification_classes=("modification", "nunique"),
        n_total_tested_sites=("n_tested_sites", "sum"),
        n_total_dm_sites=("n_dm_sites", "sum"),
        n_dm_sites_higher_symbiotic=("n_dm_sites_higher_symbiotic", "sum"),
        n_dm_sites_higher_aposymbiotic=("n_dm_sites_higher_aposymbiotic", "sum"),
        min_fdr=("min_fdr", "min"),
        max_abs_delta_mod_rate=("max_abs_delta_mod_rate", "max"),
    )
)

overall_tx["has_DM"] = overall_tx["n_total_dm_sites"] > 0

overall_merged = overall_tx.merge(
    det,
    on="transcript_id",
    how="inner"
)

overall_merged_out = os.path.join(
    OUT_DIR,
    "Differential_modification_DET_merged_transcript_overall.tsv"
)

overall_merged.to_csv(overall_merged_out, sep="\t", index=False)


# ============================================================
# Fisher enrichment tests
# ============================================================

enrichment_records = []

# Overall: any DM transcript vs DET
tmp = overall_merged.copy()
tmp["target_any_DM"] = tmp["has_DM"]

res = fisher_overlap(
    table_df=tmp,
    target_col="target_any_DM",
    target_value=True,
    det_col="is_DET"
)

res.update({
    "test_level": "any_modification",
    "modification": "any",
    "direction": "any_DM",
    "direction_label": "Any differential modification",
})

enrichment_records.append(res)

# Per modification and direction
for mod in MODS:
    sub = merged[merged["modification"] == mod].copy()

    if sub.shape[0] == 0:
        continue

    sub["is_dm_this_mod"] = sub["n_dm_sites"] > 0

    res = fisher_overlap(
        table_df=sub,
        target_col="is_dm_this_mod",
        target_value=True,
        det_col="is_DET"
    )

    res.update({
        "test_level": "modification",
        "modification": mod,
        "direction": "any_DM",
        "direction_label": "Any differential modification",
    })

    enrichment_records.append(res)

    for direction in ["higher_in_symbiotic", "higher_in_aposymbiotic", "mixed"]:
        sub["target_direction"] = sub["dm_category"] == direction

        if sub["target_direction"].sum() == 0:
            continue

        res = fisher_overlap(
            table_df=sub,
            target_col="target_direction",
            target_value=True,
            det_col="is_DET"
        )

        res.update({
            "test_level": "modification_direction",
            "modification": mod,
            "direction": direction,
            "direction_label": DIRECTION_LABELS[direction],
        })

        enrichment_records.append(res)

enrichment = pd.DataFrame(enrichment_records)

valid = enrichment["p_value"].notna()

enrichment["fdr"] = np.nan

if valid.sum() > 0:
    enrichment.loc[valid, "fdr"] = multipletests(
        enrichment.loc[valid, "p_value"],
        method="fdr_bh"
    )[1]

enrichment["log2_odds_ratio"] = np.log2(enrichment["odds_ratio"].replace(0, np.nan))

enrichment_out = os.path.join(
    OUT_DIR,
    "DM_transcript_enrichment_among_DET_Fisher_tests.tsv"
)

enrichment.to_csv(enrichment_out, sep="\t", index=False)

print("Saved enrichment table:")
print(enrichment_out)


# ============================================================
# Plot 1. DM transcripts overlapping DETs
# ============================================================

plot_df = merged.copy()

plot_df["is_dm_this_mod"] = plot_df["n_dm_sites"] > 0

count_df = (
    plot_df
    .groupby(["modification", "is_dm_this_mod", "DET_category"], observed=False)
    .size()
    .reset_index(name="n_transcripts")
)

count_df = count_df[count_df["is_dm_this_mod"]].copy()

count_df["DET_category"] = pd.Categorical(
    count_df["DET_category"],
    categories=[
        "DET up in aposymbiotic",
        "Not DET",
        "DET up in symbiotic",
    ],
    ordered=True
)

sns.set_theme(style="whitegrid", font_scale=1.3)

plt.figure(figsize=(9, 5))
ax = sns.barplot(
    data=count_df,
    x="modification",
    y="n_transcripts",
    hue="DET_category",
    order=MODS,
    palette=DET_COLORS
)

plt.xlabel("Modification type")
plt.ylabel("Number of differentially modified transcripts")
plt.xticks(rotation=45, ha="right")

ax.legend(
    title="DET status",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    frameon=False
)

plt.tight_layout()

plt.savefig(
    os.path.join(PLOT_DIR, "DM_transcripts_DET_overlap_counts_by_modification.png"),
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    os.path.join(PLOT_DIR, "DM_transcripts_DET_overlap_counts_by_modification.pdf"),
    bbox_inches="tight"
)

plt.close()


# ============================================================
# Plot 2. Fisher enrichment odds ratios
# ============================================================

forest_df = enrichment[
    enrichment["test_level"].isin(["modification"])
].copy()

forest_df = forest_df[forest_df["direction"] == "any_DM"].copy()
forest_df = forest_df[forest_df["modification"] != "any"].copy()

forest_df["minus_log10_fdr"] = -np.log10(forest_df["fdr"].replace(0, np.nan))

plt.figure(figsize=(7, 4.8))

ax = sns.scatterplot(
    data=forest_df,
    x="log2_odds_ratio",
    y="modification",
    size="minus_log10_fdr",
    hue="fdr",
    palette="viridis_r",
    sizes=(40, 250),
    edgecolor="black",
    linewidth=0.4
)

ax.axvline(0, color="black", linestyle="--", linewidth=0.8)

plt.xlabel("log2 odds ratio for DET enrichment")
plt.ylabel("Modification type")
plt.title("Enrichment of DETs among differentially modified transcripts")

ax.legend(
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    frameon=False
)

plt.tight_layout()

plt.savefig(
    os.path.join(PLOT_DIR, "DM_DET_enrichment_log2OR_by_modification.png"),
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    os.path.join(PLOT_DIR, "DM_DET_enrichment_log2OR_by_modification.pdf"),
    bbox_inches="tight"
)

plt.close()


# ============================================================
# Plot 3. DET log2FC vs representative modification-rate change
# ============================================================

scatter_df = merged.copy()

# downsample non-DM for readability
dm_rows = scatter_df[scatter_df["n_dm_sites"] > 0]
non_dm_rows = scatter_df[scatter_df["n_dm_sites"] == 0]

if non_dm_rows.shape[0] > 30000:
    non_dm_rows = non_dm_rows.sample(n=30000, random_state=123)

scatter_df = pd.concat([dm_rows, non_dm_rows], ignore_index=True)

scatter_df["plot_category"] = np.where(
    scatter_df["n_dm_sites"] > 0,
    scatter_df["dm_category_label"],
    "Not significant"
)

g = sns.FacetGrid(
    scatter_df,
    col="modification",
    col_order=MODS,
    col_wrap=4,
    height=3,
    aspect=1.1,
    sharex=True,
    sharey=True
)

def scatter_panel(data, **kwargs):
    ax = plt.gca()

    # plot non-significant first
    for cat in [
        "Not significant",
        "Higher modification in aposymbiotic",
        "Higher modification in symbiotic",
        "Mixed direction",
    ]:
        sub = data[data["plot_category"] == cat]

        if sub.shape[0] == 0:
            continue

        ax.scatter(
            sub[DET_LOG2FC_COL],
            sub["representative_delta_mod_rate_sym_minus_apo"],
            s=8 if cat == "Not significant" else 14,
            alpha=0.25 if cat == "Not significant" else 0.75,
            color=DIRECTION_COLORS.get(cat, "gray"),
            edgecolors="none",
        )

    ax.axhline(0, color="black", linewidth=0.7)
    ax.axvline(0, color="black", linewidth=0.7)

g.map_dataframe(scatter_panel)

g.set_axis_labels(
    "Transcript expression log2FC\nSymbiotic - Aposymbiotic",
    "Representative Δ modification rate\nSymbiotic - Aposymbiotic"
)

g.set_titles("{col_name}")

handles = [
    plt.Line2D([0], [0], marker="o", color="w", label=k,
               markerfacecolor=v, markersize=7)
    for k, v in DIRECTION_COLORS.items()
]

g.fig.legend(
    handles=handles,
    title="Modification status",
    bbox_to_anchor=(1.02, 0.5),
    loc="center left",
    frameon=False
)

plt.tight_layout(rect=[0, 0, 0.88, 1])

plt.savefig(
    os.path.join(PLOT_DIR, "DET_log2FC_vs_delta_modification_rate_by_modification.png"),
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    os.path.join(PLOT_DIR, "DET_log2FC_vs_delta_modification_rate_by_modification.pdf"),
    bbox_inches="tight"
)

plt.close()


# ============================================================
# Plot 4. Expression-change magnitude by DM category
# ============================================================

box_df = merged.copy()
box_df["abs_log2FC"] = box_df[DET_LOG2FC_COL].abs()

box_df["dm_category_label"] = pd.Categorical(
    box_df["dm_category_label"],
    categories=[
        "Not significant",
        "Higher modification in aposymbiotic",
        "Higher modification in symbiotic",
        "Mixed direction",
    ],
    ordered=True
)

plt.figure(figsize=(9, 5))

ax = sns.boxplot(
    data=box_df,
    x="modification",
    y="abs_log2FC",
    hue="dm_category_label",
    order=MODS,
    showfliers=False,
    palette=DIRECTION_COLORS
)

plt.xlabel("Modification type")
plt.ylabel("|Transcript expression log2FC|")
plt.xticks(rotation=45, ha="right")

ax.legend(
    title="DM category",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    frameon=False
)

plt.tight_layout()

plt.savefig(
    os.path.join(PLOT_DIR, "Expression_change_magnitude_by_DM_category.png"),
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    os.path.join(PLOT_DIR, "Expression_change_magnitude_by_DM_category.pdf"),
    bbox_inches="tight"
)

plt.close()

print("Finished DM × DET integration.")
print("Outputs saved to:")
print(OUT_DIR)
