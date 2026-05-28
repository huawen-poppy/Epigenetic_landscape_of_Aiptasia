#!/usr/bin/env python

import os
import re
import glob
import gzip
import warnings
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import fisher_exact, mannwhitneyu, ks_2samp
from statsmodels.stats.multitest import multipletests


# ============================================================
# User settings
# ============================================================

TANDEMMOD_DIR = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod"

DM_ANALYSIS_DIR = f"{TANDEMMOD_DIR}/TandemMod_site_level_differential_analysis_p0p9"

# =========================
# EDIT THIS
# =========================
TRANSD_GFF3 = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/rnabloom/busco_test/isoquant_input/isoquant_extend_transcript.fa.transdecoder.gff3"

OUT_DIR = f"{DM_ANALYSIS_DIR}/modification_position_analysis"
os.makedirs(OUT_DIR, exist_ok=True)

PLOT_DIR = f"{OUT_DIR}/plots"
os.makedirs(PLOT_DIR, exist_ok=True)

# Differential modification result file
DM_FILES = glob.glob(
    os.path.join(
        DM_ANALYSIS_DIR,
        "TandemMod_differential_modification_results_p0.9*.tsv"
    )
)

if len(DM_FILES) == 0:
    raise FileNotFoundError(
        f"No differential modification result file found in {DM_ANALYSIS_DIR}"
    )

DM_FILE = DM_FILES[0]

# Modification types
MODS = ["m6A", "m1A", "A_I", "m5C", "hm5C", "m7G", "G_I", "psU"]

# Differential modification directions
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
    "Not significant": "lightgray",
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

# Number of bins for relative-position metagene-like summaries
N_POSITION_BINS = 20

# Whether to save the full site-annotated table.
# This can be large, but useful.
SAVE_FULL_ANNOTATED_TABLE = True


# ============================================================
# Helper functions
# ============================================================

def open_maybe_gzip(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_gff3_attributes(attr_string):
    """
    Parse GFF3 attribute column into dictionary.
    Example:
        ID=AIPGENE10.p1;Parent=GENE.AIPGENE10~~AIPGENE10.p1
    """
    attrs = {}

    if pd.isna(attr_string):
        return attrs

    for item in str(attr_string).split(";"):
        item = item.strip()
        if item == "":
            continue

        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value.strip().strip('"')
        else:
            attrs[item] = ""

    return attrs


def read_transdecoder_gff3(gff3_path):
    """
    Read transcript-based TransDecoder GFF3.

    Expected format:
        seqid = transcript ID
        type = gene / mRNA / exon / CDS / five_prime_UTR / three_prime_UTR
        start/end = transcript-coordinate positions

    Returns:
        feature_df
        region_df
        transcript_length_df
        alias_map
    """

    records = []

    with open_maybe_gzip(gff3_path, "rt") as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attributes = parts

            try:
                start = int(start)
                end = int(end)
            except ValueError:
                continue

            attrs = parse_gff3_attributes(attributes)

            records.append(
                {
                    "transcript_id": seqid,
                    "source": source,
                    "feature_type": feature_type,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "phase": phase,
                    "attributes": attributes,
                    "ID": attrs.get("ID", np.nan),
                    "Parent": attrs.get("Parent", np.nan),
                    "Name": attrs.get("Name", np.nan),
                }
            )

    feature_df = pd.DataFrame(records)

    if feature_df.shape[0] == 0:
        raise ValueError(f"No valid GFF3 records were read from {gff3_path}")

    # Make sure start <= end
    feature_df["start2"] = feature_df[["start", "end"]].min(axis=1)
    feature_df["end2"] = feature_df[["start", "end"]].max(axis=1)
    feature_df["start"] = feature_df["start2"]
    feature_df["end"] = feature_df["end2"]
    feature_df = feature_df.drop(columns=["start2", "end2"])

    # Transcript length: because this GFF3 is transcript-based,
    # max end per seqid is a good transcript length estimate.
    # Prefer mRNA/gene/exon features if present, but max across all is robust.
    transcript_length_df = (
        feature_df
        .groupby("transcript_id", as_index=False)
        .agg(transcript_length=("end", "max"))
    )

    # Region features of interest
    region_df = feature_df[
        feature_df["feature_type"].isin(
            ["five_prime_UTR", "CDS", "three_prime_UTR"]
        )
    ].copy()

    # Alias map:
    # TandemMod IDs should usually match seqid, e.g. AIPGENE10.
    # But we also add mRNA IDs, e.g. AIPGENE10.p1 -> AIPGENE10,
    # just in case.
    alias_map = {}

    for tid in feature_df["transcript_id"].dropna().unique():
        alias_map[str(tid)] = str(tid)

    for _, row in feature_df.iterrows():
        seqid = str(row["transcript_id"])

        if pd.notna(row["ID"]):
            alias_map[str(row["ID"])] = seqid

        if pd.notna(row["Parent"]):
            # Parent may contain GENE.AIPGENE10~~AIPGENE10.p1
            alias_map[str(row["Parent"])] = seqid

        # Also map common p1/p2-style ORF IDs back to transcript ID if possible
        if pd.notna(row["ID"]):
            rid = str(row["ID"])
            rid_stripped = re.sub(r"\.p\d+$", "", rid)
            if rid_stripped != rid:
                alias_map[rid] = rid_stripped

    return feature_df, region_df, transcript_length_df, alias_map


def map_transcript_id(tid, alias_map):
    """
    Map TandemMod transcriptome_id to TransDecoder transcript_id.

    First exact match, then common fallbacks.
    """
    tid = str(tid)

    if tid in alias_map:
        return alias_map[tid]

    # If TandemMod ID has protein suffix, remove it
    stripped = re.sub(r"\.p\d+$", "", tid)
    if stripped in alias_map:
        return alias_map[stripped]

    return tid


def annotate_sites_by_region(site_df, region_df, transcript_length_df):
    """
    Efficiently annotate all tested sites with TransDecoder transcript regions.

    site_df must contain:
        transcript_id_mapped
        site_numeric

    Returns site_df with:
        transcript_length
        relative_position
        transcript_region
    """

    site_df = site_df.copy()

    # Add transcript length
    site_df = site_df.merge(
        transcript_length_df,
        left_on="transcript_id_mapped",
        right_on="transcript_id",
        how="left",
        suffixes=("", "_gff")
    )

    if "transcript_id_gff" in site_df.columns:
        site_df = site_df.drop(columns=["transcript_id_gff"])

    site_df["relative_position"] = (
        site_df["site_numeric"] / site_df["transcript_length"]
    )

    # Keep relative positions in valid range only
    site_df.loc[
        (site_df["relative_position"] < 0) |
        (site_df["relative_position"] > 1),
        "relative_position"
    ] = np.nan

    site_df["transcript_region"] = np.where(
        site_df["transcript_length"].notna(),
        "unannotated_region",
        "not_in_TransDecoder"
    )

    # Group intervals by transcript for faster annotation
    region_grouped = {
        tid: sub.sort_values(["start", "end"]).copy()
        for tid, sub in region_df.groupby("transcript_id")
    }

    # Priority: UTR/CDS intervals usually do not overlap, but if they do,
    # this priority prevents unannotated from overwriting real features.
    region_priority = [
        "five_prime_UTR",
        "CDS",
        "three_prime_UTR",
    ]

    annotated_parts = []

    n_groups = site_df["transcript_id_mapped"].nunique()
    print(f"Annotating sites by transcript region across {n_groups:,} transcripts...")

    for i, (tid, sub_sites) in enumerate(site_df.groupby("transcript_id_mapped"), start=1):
        if i % 5000 == 0:
            print(f"  annotated {i:,}/{n_groups:,} transcript groups")

        sub_sites = sub_sites.copy()

        if tid not in region_grouped:
            annotated_parts.append(sub_sites)
            continue

        sub_regions = region_grouped[tid]

        for region in region_priority:
            rr = sub_regions[sub_regions["feature_type"] == region]

            if rr.shape[0] == 0:
                continue

            mask_total = pd.Series(False, index=sub_sites.index)

            for _, r in rr.iterrows():
                mask = (
                    (sub_sites["site_numeric"] >= r["start"]) &
                    (sub_sites["site_numeric"] <= r["end"])
                )
                mask_total = mask_total | mask

            sub_sites.loc[mask_total, "transcript_region"] = region

        annotated_parts.append(sub_sites)

    annotated = pd.concat(annotated_parts, ignore_index=True)

    annotated["transcript_region_label"] = annotated["transcript_region"].map(
        REGION_LABELS
    )

    annotated["transcript_region_label"] = annotated["transcript_region_label"].fillna(
        annotated["transcript_region"]
    )

    return annotated


def fisher_region_enrichment(df):
    """
    Region enrichment of differential sites compared to all other tested sites
    within the same modification class.

    Tests:
        1. any differential site vs non-DM
        2. higher-in-symbiotic vs all other sites
        3. higher-in-aposymbiotic vs all other sites
    """

    records = []

    test_groups = [
        ("any_differential", "Any differential", lambda x: x["is_dm_site"]),
        (
            "higher_in_symbiotic",
            "Higher in symbiotic",
            lambda x: x["direction"] == "higher_in_symbiotic",
        ),
        (
            "higher_in_aposymbiotic",
            "Higher in aposymbiotic",
            lambda x: x["direction"] == "higher_in_aposymbiotic",
        ),
    ]

    for mod in MODS:
        mod_df = df[df["modification"] == mod].copy()

        if mod_df.shape[0] == 0:
            continue

        for test_name, test_label, target_func in test_groups:
            target_mask = target_func(mod_df)

            for region in REGION_ORDER:
                region_mask = mod_df["transcript_region"] == region

                a = (target_mask & region_mask).sum()
                b = (target_mask & ~region_mask).sum()
                c = (~target_mask & region_mask).sum()
                d = (~target_mask & ~region_mask).sum()

                if (a + b) == 0 or (c + d) == 0:
                    odds_ratio, p_value = np.nan, np.nan
                else:
                    try:
                        odds_ratio, p_value = fisher_exact(
                            [[a, b], [c, d]],
                            alternative="greater"
                        )
                    except Exception:
                        odds_ratio, p_value = np.nan, np.nan

                records.append(
                    {
                        "modification": mod,
                        "test_name": test_name,
                        "test_label": test_label,
                        "region": region,
                        "region_label": REGION_LABELS[region],
                        "target_in_region": int(a),
                        "target_not_region": int(b),
                        "background_in_region": int(c),
                        "background_not_region": int(d),
                        "target_total": int(a + b),
                        "background_total": int(c + d),
                        "target_region_fraction": a / (a + b)
                        if (a + b) > 0 else np.nan,
                        "background_region_fraction": c / (c + d)
                        if (c + d) > 0 else np.nan,
                        "odds_ratio": odds_ratio,
                        "p_value": p_value,
                    }
                )

    out = pd.DataFrame(records)

    out["fdr"] = np.nan

    # BH correction within each test type across all modification-region combinations
    for test_name, idx in out.groupby("test_name").groups.items():
        idx = list(idx)
        pvals = out.loc[idx, "p_value"]
        valid = pvals.notna()

        if valid.sum() > 0:
            valid_idx = pvals[valid].index
            out.loc[valid_idx, "fdr"] = multipletests(
                pvals[valid],
                method="fdr_bh"
            )[1]

    out["log2_odds_ratio"] = np.log2(out["odds_ratio"].replace(0, np.nan))
    out["minus_log10_fdr"] = -np.log10(out["fdr"].replace(0, np.nan))

    return out


def relative_position_tests(df):
    """
    Compare relative position distributions:
        any DM vs non-DM
        higher in symbiotic vs all others
        higher in aposymbiotic vs all others

    within each modification type.

    Uses Mann-Whitney U and Kolmogorov-Smirnov tests.
    """

    records = []

    test_groups = [
        ("any_differential", "Any differential", lambda x: x["is_dm_site"]),
        (
            "higher_in_symbiotic",
            "Higher in symbiotic",
            lambda x: x["direction"] == "higher_in_symbiotic",
        ),
        (
            "higher_in_aposymbiotic",
            "Higher in aposymbiotic",
            lambda x: x["direction"] == "higher_in_aposymbiotic",
        ),
    ]

    for mod in MODS:
        mod_df = df[
            (df["modification"] == mod) &
            (df["relative_position"].notna())
        ].copy()

        if mod_df.shape[0] == 0:
            continue

        for test_name, test_label, target_func in test_groups:
            target_mask = target_func(mod_df)

            target = mod_df.loc[target_mask, "relative_position"].dropna()
            background = mod_df.loc[~target_mask, "relative_position"].dropna()

            if len(target) < 5 or len(background) < 5:
                mw_stat, mw_p = np.nan, np.nan
                ks_stat, ks_p = np.nan, np.nan
            else:
                mw_stat, mw_p = mannwhitneyu(
                    target,
                    background,
                    alternative="two-sided"
                )

                ks_stat, ks_p = ks_2samp(
                    target,
                    background,
                    alternative="two-sided"
                )

            records.append(
                {
                    "modification": mod,
                    "test_name": test_name,
                    "test_label": test_label,
                    "n_target": len(target),
                    "n_background": len(background),
                    "median_relative_position_target": target.median()
                    if len(target) > 0 else np.nan,
                    "median_relative_position_background": background.median()
                    if len(background) > 0 else np.nan,
                    "mean_relative_position_target": target.mean()
                    if len(target) > 0 else np.nan,
                    "mean_relative_position_background": background.mean()
                    if len(background) > 0 else np.nan,
                    "mannwhitney_u": mw_stat,
                    "mannwhitney_p": mw_p,
                    "ks_statistic": ks_stat,
                    "ks_p": ks_p,
                }
            )

    out = pd.DataFrame(records)

    for pcol in ["mannwhitney_p", "ks_p"]:
        fdr_col = pcol.replace("_p", "_fdr")
        out[fdr_col] = np.nan

        valid = out[pcol].notna()

        if valid.sum() > 0:
            out.loc[valid, fdr_col] = multipletests(
                out.loc[valid, pcol],
                method="fdr_bh"
            )[1]

    return out


def add_position_bins(df, n_bins=20):
    df = df.copy()

    bins = np.linspace(0, 1, n_bins + 1)
    labels = [
        f"{bins[i]:.2f}-{bins[i+1]:.2f}"
        for i in range(n_bins)
    ]

    df["relative_position_bin"] = pd.cut(
        df["relative_position"],
        bins=bins,
        labels=labels,
        include_lowest=True
    )

    df["relative_position_bin_mid"] = df["relative_position_bin"].apply(
        lambda x: np.nan if pd.isna(x) else (
            float(str(x).split("-")[0]) + float(str(x).split("-")[1])
        ) / 2
    )

    return df


# ============================================================
# Main analysis
# ============================================================

print("Using differential modification result file:")
print(DM_FILE)

print("Using TransDecoder GFF3:")
print(TRANSD_GFF3)

# -----------------------------
# Read TransDecoder GFF3
# -----------------------------

feature_df, region_df, transcript_length_df, alias_map = read_transdecoder_gff3(
    TRANSD_GFF3
)

feature_out = os.path.join(OUT_DIR, "Parsed_TransDecoder_GFF3_features.tsv")
region_out = os.path.join(OUT_DIR, "Parsed_TransDecoder_regions.tsv")
length_out = os.path.join(OUT_DIR, "Parsed_TransDecoder_transcript_lengths.tsv")

feature_df.to_csv(feature_out, sep="\t", index=False)
region_df.to_csv(region_out, sep="\t", index=False)
transcript_length_df.to_csv(length_out, sep="\t", index=False)

print(f"Parsed GFF3 features: {feature_df.shape}")
print(f"Region features: {region_df.shape}")
print(f"Transcript lengths: {transcript_length_df.shape}")


# -----------------------------
# Read differential modification results
# -----------------------------

dm = pd.read_csv(DM_FILE, sep="\t")

required_cols = [
    "site_id",
    "modification",
    "transcriptome_id",
    "site",
    "motif",
    "direction",
    "fdr_by_modification",
    "delta_mean_mod_rate_sym_minus_apo",
]

missing = [c for c in required_cols if c not in dm.columns]

if missing:
    raise ValueError(f"Differential modification result missing columns: {missing}")

# One row per tested site
dm = dm.drop_duplicates("site_id").copy()

dm["site_numeric"] = pd.to_numeric(dm["site"], errors="coerce")

dm = dm.dropna(subset=["site_numeric"]).copy()

dm["site_numeric"] = dm["site_numeric"].astype(int)

dm["transcript_id_mapped"] = dm["transcriptome_id"].apply(
    lambda x: map_transcript_id(x, alias_map)
)

dm["is_dm_site"] = dm["direction"].isin(DIRECTION_ORDER)

dm["direction_label"] = dm["direction"].map(DIRECTION_LABELS)
dm["direction_label"] = dm["direction_label"].fillna("Not significant")

dm["modification"] = pd.Categorical(
    dm["modification"],
    categories=MODS,
    ordered=True
)

print("DM table after preprocessing:")
print(dm.shape)

# -----------------------------
# Annotate sites by transcript region
# -----------------------------

annotated = annotate_sites_by_region(
    site_df=dm,
    region_df=region_df,
    transcript_length_df=transcript_length_df
)

annotated = add_position_bins(
    annotated,
    n_bins=N_POSITION_BINS
)

# Save full annotated table
if SAVE_FULL_ANNOTATED_TABLE:
    annotated_out = os.path.join(
        OUT_DIR,
        "All_tested_modification_sites_with_transcript_position_annotation.tsv.gz"
    )

    annotated.to_csv(
        annotated_out,
        sep="\t",
        index=False,
        compression="gzip"
    )

    print("Saved full annotated site table:")
    print(annotated_out)

# Smaller compact table
compact_cols = [
    "site_id",
    "modification",
    "transcriptome_id",
    "transcript_id_mapped",
    "site",
    "motif",
    "direction",
    "direction_label",
    "is_dm_site",
    "fdr_by_modification",
    "delta_mean_mod_rate_sym_minus_apo",
    "transcript_length",
    "relative_position",
    "relative_position_bin",
    "relative_position_bin_mid",
    "transcript_region",
    "transcript_region_label",
]

compact_out = os.path.join(
    OUT_DIR,
    "Compact_modification_site_position_annotation.tsv.gz"
)

annotated[compact_cols].to_csv(
    compact_out,
    sep="\t",
    index=False,
    compression="gzip"
)

print("Saved compact annotated table:")
print(compact_out)


# ============================================================
# Summary 1. Region composition
# ============================================================

region_counts = (
    annotated
    .groupby(
        [
            "modification",
            "direction",
            "direction_label",
            "transcript_region",
            "transcript_region_label",
        ],
        observed=False
    )
    .size()
    .reset_index(name="n_sites")
)

region_counts["total_sites_in_mod_direction"] = (
    region_counts
    .groupby(["modification", "direction"], observed=False)["n_sites"]
    .transform("sum")
)

region_counts["fraction"] = (
    region_counts["n_sites"] /
    region_counts["total_sites_in_mod_direction"]
)

region_counts_out = os.path.join(
    OUT_DIR,
    "Region_composition_by_modification_and_direction.tsv"
)

region_counts.to_csv(region_counts_out, sep="\t", index=False)

print("Saved region composition table:")
print(region_counts_out)


# ============================================================
# Summary 2. Region enrichment Fisher tests
# ============================================================

region_enrichment = fisher_region_enrichment(annotated)

region_enrichment_out = os.path.join(
    OUT_DIR,
    "Region_enrichment_Fisher_tests.tsv"
)

region_enrichment.to_csv(
    region_enrichment_out,
    sep="\t",
    index=False
)

print("Saved region enrichment table:")
print(region_enrichment_out)


# ============================================================
# Summary 3. Relative position tests
# ============================================================

relpos_tests = relative_position_tests(annotated)

relpos_tests_out = os.path.join(
    OUT_DIR,
    "Relative_position_distribution_tests.tsv"
)

relpos_tests.to_csv(
    relpos_tests_out,
    sep="\t",
    index=False
)

print("Saved relative-position tests:")
print(relpos_tests_out)


# ============================================================
# Summary 4. Binned relative position distribution
# ============================================================

# For plotting, collapse directions into:
# all tested, higher in apo, higher in sym
plot_parts = []

all_tested = annotated[annotated["relative_position"].notna()].copy()
all_tested["position_group"] = "All tested sites"
plot_parts.append(all_tested)

higher_apo = annotated[
    (annotated["direction"] == "higher_in_aposymbiotic") &
    (annotated["relative_position"].notna())
].copy()
higher_apo["position_group"] = "Higher in aposymbiotic"
plot_parts.append(higher_apo)

higher_sym = annotated[
    (annotated["direction"] == "higher_in_symbiotic") &
    (annotated["relative_position"].notna())
].copy()
higher_sym["position_group"] = "Higher in symbiotic"
plot_parts.append(higher_sym)

position_plot_df = pd.concat(plot_parts, ignore_index=True)

binned = (
    position_plot_df
    .groupby(
        [
            "modification",
            "position_group",
            "relative_position_bin",
            "relative_position_bin_mid",
        ],
        observed=False
    )
    .size()
    .reset_index(name="n_sites")
)

binned["total_sites_in_group"] = (
    binned
    .groupby(["modification", "position_group"], observed=False)["n_sites"]
    .transform("sum")
)

binned["fraction"] = binned["n_sites"] / binned["total_sites_in_group"]

binned_out = os.path.join(
    OUT_DIR,
    "Binned_relative_position_distribution.tsv"
)

binned.to_csv(
    binned_out,
    sep="\t",
    index=False
)

print("Saved binned relative-position distribution:")
print(binned_out)


# ============================================================
# Plots
# ============================================================

sns.set_theme(style="whitegrid", font_scale=1.2)

# ------------------------------------------------------------
# Plot 1. Region composition stacked bars
# ------------------------------------------------------------

comp = region_counts.copy()

# Keep key directions only for clarity
comp = comp[
    comp["direction"].isin(
        ["not_significant", "higher_in_aposymbiotic", "higher_in_symbiotic"]
    )
].copy()

comp["group_label"] = comp["direction"].map(DIRECTION_LABELS)
comp["group_label"] = comp["group_label"].fillna("Not significant")

# Plot only DM directions in main plot
comp_dm = comp[
    comp["direction"].isin(["higher_in_aposymbiotic", "higher_in_symbiotic"])
].copy()

# Ensure region label order
region_label_order = [REGION_LABELS[x] for x in REGION_ORDER]

comp_dm["transcript_region_label"] = pd.Categorical(
    comp_dm["transcript_region_label"],
    categories=region_label_order,
    ordered=True
)

# Build stacked bar manually
for mod in MODS:
    sub = comp_dm[comp_dm["modification"] == mod].copy()

    if sub.shape[0] == 0:
        continue

# Combined stacked bar for all modifications and directions
pivot_comp = comp_dm.pivot_table(
    index=["modification", "group_label"],
    columns="transcript_region_label",
    values="fraction",
    aggfunc="sum",
    fill_value=0,
    observed=False
)

pivot_comp = pivot_comp.reindex(columns=region_label_order, fill_value=0)

fig, ax = plt.subplots(figsize=(12, 5.5))

bottom = np.zeros(pivot_comp.shape[0])
x = np.arange(pivot_comp.shape[0])

for region_label in region_label_order:
    vals = pivot_comp[region_label].values

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

xlabels = [
    f"{idx[0]}\n{idx[1].replace('Higher in ', 'High in ')}"
    for idx in pivot_comp.index
]

ax.set_xticks(x)
ax.set_xticklabels(xlabels, rotation=45, ha="right")
ax.set_ylabel("Fraction of differential sites")
ax.set_xlabel("Modification type and direction")
ax.legend(
    title="Transcript region",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    frameon=False
)

plt.tight_layout()

plt.savefig(
    os.path.join(PLOT_DIR, "Region_composition_differential_sites_stacked_bar.png"),
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    os.path.join(PLOT_DIR, "Region_composition_differential_sites_stacked_bar.pdf"),
    bbox_inches="tight"
)

plt.close()


# ------------------------------------------------------------
# Plot 2. Region enrichment dotplot
# ------------------------------------------------------------

enrich_plot = region_enrichment.copy()

enrich_plot = enrich_plot[
    enrich_plot["test_name"].isin(
        ["higher_in_aposymbiotic", "higher_in_symbiotic"]
    )
].copy()

# Avoid infinite plotting
enrich_plot["minus_log10_fdr_plot"] = enrich_plot["minus_log10_fdr"].clip(upper=50)

enrich_plot["region_label"] = pd.Categorical(
    enrich_plot["region_label"],
    categories=region_label_order,
    ordered=True
)

enrich_plot["mod_region"] = (
    enrich_plot["modification"].astype(str) +
    " | " +
    enrich_plot["region_label"].astype(str)
)

# Keep only results with OR > 1 for cleaner enrichment plot
enrich_plot2 = enrich_plot[
    (enrich_plot["odds_ratio"] > 1) &
    (enrich_plot["fdr"] < 0.05)
].copy()

if enrich_plot2.shape[0] > 0:
    g = sns.FacetGrid(
        enrich_plot2,
        col="test_label",
        col_order=["Higher in aposymbiotic", "Higher in symbiotic"],
        height=6,
        aspect=0.9,
        sharey=False,
        sharex=True
    )

    def dot_panel(data, **kwargs):
        ax = plt.gca()

        data = data.sort_values("log2_odds_ratio", ascending=True)

        sc = ax.scatter(
            data["log2_odds_ratio"],
            data["mod_region"],
            s=data["minus_log10_fdr_plot"] * 8 + 20,
            c=data["target_region_fraction"],
            cmap="YlOrBr",
            edgecolor="black",
            linewidth=0.3
        )

        ax.axvline(0, color="black", linestyle="--", linewidth=0.8)
        ax.set_xlabel("log2 odds ratio")
        ax.set_ylabel("Modification | region")

    g.map_dataframe(dot_panel)
    g.set_titles("{col_name}")

    plt.tight_layout()

    plt.savefig(
        os.path.join(PLOT_DIR, "Region_enrichment_dotplot_significant.png"),
        dpi=300,
        bbox_inches="tight"
    )

    plt.savefig(
        os.path.join(PLOT_DIR, "Region_enrichment_dotplot_significant.pdf"),
        bbox_inches="tight"
    )

    plt.close()

else:
    print("No significant region enrichment with OR > 1 for dotplot.")


# ------------------------------------------------------------
# Plot 3. Relative position density
# ------------------------------------------------------------

density_df = annotated[
    annotated["relative_position"].notna()
].copy()

# Downsample all tested/not significant for plotting readability
plot_density_parts = []

for mod in MODS:
    mod_df = density_df[density_df["modification"] == mod].copy()

    # All tested background
    bg = mod_df.copy()
    bg["position_group"] = "All tested sites"

    if bg.shape[0] > 100000:
        bg = bg.sample(n=100000, random_state=123)

    plot_density_parts.append(bg)

    # DM directions
    for direction, label in [
        ("higher_in_aposymbiotic", "Higher in aposymbiotic"),
        ("higher_in_symbiotic", "Higher in symbiotic"),
    ]:
        target = mod_df[mod_df["direction"] == direction].copy()
        target["position_group"] = label

        if target.shape[0] > 100000:
            target = target.sample(n=100000, random_state=123)

        plot_density_parts.append(target)

density_plot_df = pd.concat(plot_density_parts, ignore_index=True)

position_group_palette = {
    "All tested sites": "lightgray",
    "Higher in aposymbiotic": "#D8B88A",
    "Higher in symbiotic": "#F28E2B",
}

g = sns.FacetGrid(
    density_plot_df,
    col="modification",
    col_order=MODS,
    col_wrap=4,
    hue="position_group",
    hue_order=[
        "All tested sites",
        "Higher in aposymbiotic",
        "Higher in symbiotic",
    ],
    palette=position_group_palette,
    height=3,
    aspect=1.1,
    sharex=True,
    sharey=False
)

g.map_dataframe(
    sns.kdeplot,
    x="relative_position",
    linewidth=1.4,
    fill=False,
    common_norm=False,
    clip=(0, 1),
    cut=0
)

g.set_axis_labels("Relative transcript position", "Density")
g.set_titles("{col_name}")

for ax in g.axes.flatten():
    ax.set_xlim(0, 1)
    ax.axvline(0.5, color="black", linewidth=0.5, linestyle="--", alpha=0.5)

g.add_legend(
    title="Site group",
    bbox_to_anchor=(1.02, 0.5),
    loc="center left",
    frameon=False
)

plt.tight_layout(rect=[0, 0, 0.88, 1])

plt.savefig(
    os.path.join(PLOT_DIR, "Relative_position_density_by_modification.png"),
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    os.path.join(PLOT_DIR, "Relative_position_density_by_modification.pdf"),
    bbox_inches="tight"
)

plt.close()


# ------------------------------------------------------------
# Plot 4. Binned relative position metagene-like line plot
# ------------------------------------------------------------

binned_plot = binned.copy()

binned_plot["position_group"] = pd.Categorical(
    binned_plot["position_group"],
    categories=[
        "All tested sites",
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
    hue="position_group",
    hue_order=[
        "All tested sites",
        "Higher in aposymbiotic",
        "Higher in symbiotic",
    ],
    palette=position_group_palette,
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

g.set_axis_labels("Relative transcript position", "Fraction of sites")
g.set_titles("{col_name}")

for ax in g.axes.flatten():
    ax.set_xlim(0, 1)

g.add_legend(
    title="Site group",
    bbox_to_anchor=(1.02, 0.5),
    loc="center left",
    frameon=False
)

plt.tight_layout(rect=[0, 0, 0.88, 1])

plt.savefig(
    os.path.join(PLOT_DIR, "Binned_relative_position_distribution_by_modification.png"),
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    os.path.join(PLOT_DIR, "Binned_relative_position_distribution_by_modification.pdf"),
    bbox_inches="tight"
)

plt.close()


# ============================================================
# Final summary
# ============================================================

summary_records = []

for mod in MODS:
    sub = annotated[annotated["modification"] == mod].copy()

    if sub.shape[0] == 0:
        continue

    dm_sub = sub[sub["is_dm_site"]]
    sym_sub = sub[sub["direction"] == "higher_in_symbiotic"]
    apo_sub = sub[sub["direction"] == "higher_in_aposymbiotic"]

    summary_records.append(
        {
            "modification": mod,
            "n_tested_sites": sub.shape[0],
            "n_sites_with_transdecoder_annotation": (
                sub["transcript_region"] != "not_in_TransDecoder"
            ).sum(),
            "percent_sites_with_transdecoder_annotation": (
                (sub["transcript_region"] != "not_in_TransDecoder").mean() * 100
            ),
            "n_differential_sites": dm_sub.shape[0],
            "n_higher_in_symbiotic": sym_sub.shape[0],
            "n_higher_in_aposymbiotic": apo_sub.shape[0],
            "median_relative_position_all_tested": sub["relative_position"].median(),
            "median_relative_position_differential": dm_sub["relative_position"].median(),
            "median_relative_position_higher_symbiotic": sym_sub["relative_position"].median(),
            "median_relative_position_higher_aposymbiotic": apo_sub["relative_position"].median(),
        }
    )

summary_df = pd.DataFrame(summary_records)

summary_out = os.path.join(
    OUT_DIR,
    "Modification_position_analysis_summary_by_modification.tsv"
)

summary_df.to_csv(summary_out, sep="\t", index=False)

print("\nSaved final summary:")
print(summary_out)

print("\nMain output tables:")
print(compact_out)
print(region_counts_out)
print(region_enrichment_out)
print(relpos_tests_out)
print(binned_out)
print(summary_out)

print("\nPlots saved to:")
print(PLOT_DIR)

print("\nAnalysis finished successfully.")
