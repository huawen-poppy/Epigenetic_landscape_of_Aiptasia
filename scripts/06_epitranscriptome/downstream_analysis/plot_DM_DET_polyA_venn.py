#!/usr/bin/env python

import os
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib_venn import venn2, venn3


# ============================================================
# User settings
# ============================================================

OUTDIR = "DM_DET_polyA_venn_plots"
os.makedirs(OUTDIR, exist_ok=True)

DET_OVERALL_FILE = "../TandemMod_site_level_differential_analysis_p0p9/integration_with_DET/Differential_modification_DET_merged_transcript_overall.tsv"
POLYA_OVERALL_FILE = "../TandemMod_site_level_differential_analysis_p0p9/integration_with_polyA/Differential_modification_polyA_merged_transcript_overall.tsv"

# Plot colors
COLOR_DM = "#7F7F7F"
COLOR_DET = "#F28E2B"
COLOR_DPA = "#D8B88A"


# ============================================================
# Helper functions
# ============================================================

def to_bool_series(s):
    """
    Robustly convert True/False-like columns to boolean.
    """
    if s.dtype == bool:
        return s

    return (
        s.astype(str)
        .str.strip()
        .str.lower()
        .isin(["true", "1", "yes", "y", "t"])
    )


def get_transcript_col(df):
    candidates = ["transcript_id", "transcriptome_id", "contig", "gene_id"]

    for c in candidates:
        if c in df.columns:
            return c

    raise ValueError(
        f"Could not find transcript ID column. Available columns: {df.columns.tolist()}"
    )


def get_bool_col(df, preferred_cols, fallback_category_cols=None, positive_values=None):
    """
    Return boolean series based on preferred boolean columns.
    If not found, use category columns.
    """
    for c in preferred_cols:
        if c in df.columns:
            return to_bool_series(df[c])

    if fallback_category_cols is not None:
        for c in fallback_category_cols:
            if c in df.columns:
                if positive_values is None:
                    raise ValueError("positive_values must be provided for category fallback.")

                return df[c].astype(str).isin(positive_values)

    raise ValueError(
        f"Could not find any of {preferred_cols} or fallback columns {fallback_category_cols}. "
        f"Available columns: {df.columns.tolist()}"
    )


def save_set_list(values, filename):
    pd.Series(sorted(values)).to_csv(
        os.path.join(OUTDIR, filename),
        index=False,
        header=False
    )


def format_venn_labels(v):
    """
    Increase font size for Venn labels.
    """
    if v is None:
        return

    for text in v.set_labels:
        if text is not None:
            text.set_fontsize(13)

    for text in v.subset_labels:
        if text is not None:
            text.set_fontsize(12)


# ============================================================
# Load files
# ============================================================

det_df = pd.read_csv(DET_OVERALL_FILE, sep="\t")
polya_df = pd.read_csv(POLYA_OVERALL_FILE, sep="\t")

det_tx_col = get_transcript_col(det_df)
polya_tx_col = get_transcript_col(polya_df)

det_df = det_df.rename(columns={det_tx_col: "transcript_id"})
polya_df = polya_df.rename(columns={polya_tx_col: "transcript_id"})

det_df["transcript_id"] = det_df["transcript_id"].astype(str)
polya_df["transcript_id"] = polya_df["transcript_id"].astype(str)


# ============================================================
# Define sets
# ============================================================

# DM and DET from the DET integration file
det_df["has_DM_bool"] = get_bool_col(
    det_df,
    preferred_cols=["has_DM", "is_DM", "has_differential_modification"]
)

det_df["is_DET_bool"] = get_bool_col(
    det_df,
    preferred_cols=["is_DET", "has_DET", "is_differentially_expressed"],
    fallback_category_cols=["DET_category", "det_category"],
    positive_values=[
        "DET up in symbiotic",
        "DET up in aposymbiotic",
        "up_in_symbiotic",
        "up_in_aposymbiotic",
        "up",
        "down",
    ]
)

DM_for_DET = set(det_df.loc[det_df["has_DM_bool"], "transcript_id"])
DET = set(det_df.loc[det_df["is_DET_bool"], "transcript_id"])


# DM and DPA from the polyA integration file
polya_df["has_DM_bool"] = get_bool_col(
    polya_df,
    preferred_cols=["has_DM", "is_DM", "has_differential_modification"]
)

polya_df["is_DPA_bool"] = get_bool_col(
    polya_df,
    preferred_cols=["is_DPA", "has_DPA", "is_differential_polyA", "is_differential_polya"],
    fallback_category_cols=["polyA_category", "polya_category"],
    positive_values=[
        "Longer poly(A) in aposymbiotic",
        "Longer poly(A) in symbiotic",
        "longer_in_aposymbiotic",
        "longer_in_symbiotic",
    ]
)

DM_for_DPA = set(polya_df.loc[polya_df["has_DM_bool"], "transcript_id"])
DPA = set(polya_df.loc[polya_df["is_DPA_bool"], "transcript_id"])


# For the 3-way comparison, restrict to transcripts present in both integration files
merged = det_df[
    ["transcript_id", "has_DM_bool", "is_DET_bool"]
].merge(
    polya_df[["transcript_id", "has_DM_bool", "is_DPA_bool"]],
    on="transcript_id",
    how="inner",
    suffixes=("_det", "_polya")
)

# Usually these two DM columns should agree. Use OR to be inclusive.
merged["has_DM_3way"] = merged["has_DM_bool_det"] | merged["has_DM_bool_polya"]

DM_3WAY = set(merged.loc[merged["has_DM_3way"], "transcript_id"])
DET_3WAY = set(merged.loc[merged["is_DET_bool"], "transcript_id"])
DPA_3WAY = set(merged.loc[merged["is_DPA_bool"], "transcript_id"])


# ============================================================
# Save transcript lists
# ============================================================

save_set_list(DM_for_DET, "transcripts_DM_for_DET_overlap.txt")
save_set_list(DET, "transcripts_DET.txt")

save_set_list(DM_for_DPA, "transcripts_DM_for_DPA_overlap.txt")
save_set_list(DPA, "transcripts_differential_polyA.txt")

save_set_list(DM_3WAY, "transcripts_DM_3way_background.txt")
save_set_list(DET_3WAY, "transcripts_DET_3way_background.txt")
save_set_list(DPA_3WAY, "transcripts_differential_polyA_3way_background.txt")

# Also save full merged 3-way table
merged.to_csv(
    os.path.join(OUTDIR, "DM_DET_DPA_3way_merged_background.tsv"),
    sep="\t",
    index=False
)


# ============================================================
# Print summary
# ============================================================

print("========== Set sizes ==========")
print(f"DM transcripts in DET background: {len(DM_for_DET)}")
print(f"DET transcripts: {len(DET)}")
print(f"DM ∩ DET: {len(DM_for_DET & DET)}")

print(f"DM transcripts in polyA background: {len(DM_for_DPA)}")
print(f"Differential-polyA transcripts: {len(DPA)}")
print(f"DM ∩ differential-polyA: {len(DM_for_DPA & DPA)}")

print(f"3-way background transcripts: {merged['transcript_id'].nunique()}")
print(f"DM transcripts in 3-way background: {len(DM_3WAY)}")
print(f"DET transcripts in 3-way background: {len(DET_3WAY)}")
print(f"Differential-polyA transcripts in 3-way background: {len(DPA_3WAY)}")
print(f"DM ∩ DET ∩ DPA: {len(DM_3WAY & DET_3WAY & DPA_3WAY)}")


# ============================================================
# Plot 1: DM vs DET
# ============================================================

plt.figure(figsize=(5.2, 4.8))

v = venn2(
    [DM_for_DET, DET],
    set_labels=("Differentially modified\ntranscripts", "DETs"),
    set_colors=(COLOR_DM, COLOR_DET),
    alpha=0.65
)

format_venn_labels(v)

plt.title("Differential modification and transcript expression", fontsize=14)
plt.tight_layout()

plt.savefig(
    os.path.join(OUTDIR, "Venn_DM_vs_DET.png"),
    dpi=300,
    bbox_inches="tight"
)
plt.savefig(
    os.path.join(OUTDIR, "Venn_DM_vs_DET.pdf"),
    bbox_inches="tight"
)
plt.close()


# ============================================================
# Plot 2: DM vs differential polyA
# ============================================================

plt.figure(figsize=(5.2, 4.8))

v = venn2(
    [DM_for_DPA, DPA],
    set_labels=("Differentially modified\ntranscripts", "Differential-poly(A)\ntranscripts"),
    set_colors=(COLOR_DM, COLOR_DPA),
    alpha=0.65
)

format_venn_labels(v)

plt.title("Differential modification and poly(A) tail length", fontsize=14)
plt.tight_layout()

plt.savefig(
    os.path.join(OUTDIR, "Venn_DM_vs_differential_polyA.png"),
    dpi=300,
    bbox_inches="tight"
)
plt.savefig(
    os.path.join(OUTDIR, "Venn_DM_vs_differential_polyA.pdf"),
    bbox_inches="tight"
)
plt.close()


# ============================================================
# Plot 3: DM vs DET vs differential polyA
# ============================================================

plt.figure(figsize=(6.2, 5.6))

v = venn3(
    [DM_3WAY, DET_3WAY, DPA_3WAY],
    set_labels=(
        "Differentially modified\ntranscripts",
        "DETs",
        "Differential-poly(A)\ntranscripts"
    ),
    set_colors=(COLOR_DM, COLOR_DET, COLOR_DPA),
    alpha=0.65
)

format_venn_labels(v)

plt.title("Overlap among RNA modification, expression and poly(A) changes", fontsize=14)
plt.tight_layout()

plt.savefig(
    os.path.join(OUTDIR, "Venn_DM_DET_differential_polyA_3way.png"),
    dpi=300,
    bbox_inches="tight"
)
plt.savefig(
    os.path.join(OUTDIR, "Venn_DM_DET_differential_polyA_3way.pdf"),
    bbox_inches="tight"
)
plt.close()


print("\nSaved Venn plots and transcript lists to:")
print(OUTDIR)
