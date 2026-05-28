#!/usr/bin/env python

import os
import glob
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import fisher_exact, spearmanr, mannwhitneyu
from statsmodels.stats.multitest import multipletests


# ============================================================
# User settings
# ============================================================

TANDEMMOD_DIR = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod"

DM_ANALYSIS_DIR = f"{TANDEMMOD_DIR}/TandemMod_site_level_differential_analysis_p0p9"

OUT_DIR = f"{DM_ANALYSIS_DIR}/integration_with_polyA"
os.makedirs(OUT_DIR, exist_ok=True)

PLOT_DIR = f"{OUT_DIR}/plots"
os.makedirs(PLOT_DIR, exist_ok=True)


POLYA_FILE = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/pipeline-polya-diff/report/polya_diff_per_transcript.tsv"


# ============================================================
# Important direction setting
# ============================================================
# The polyA table has:
# median diff = treatment median - control median
#
# Now:
#   control = aposymbiotic
#   treatment = symbiotic
# then:
#   delta_polyA_sym_minus_apo = median diff
# so use:
#   POLYA_DELTA_MULTIPLIER = 1
#

POLYA_DELTA_MULTIPLIER = 1


# ============================================================
# Differential polyA definition
# ============================================================

POLYA_FDR_CUTOFF = 0.05

# Minimum absolute median tail-length difference in nt.
POLYA_ABS_DELTA_CUTOFF = 5


# ============================================================
# Modification settings
# ============================================================

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

POLYA_COLORS = {
    "Longer poly(A) in aposymbiotic": "#D8B88A",
    "Longer poly(A) in symbiotic": "#F28E2B",
    "Not differential poly(A)": "lightgray",
}


# ============================================================
# Locate differential modification result file
# ============================================================

dm_files = glob.glob(
    os.path.join(
        DM_ANALYSIS_DIR,
        "TandemMod_differential_modification_results_p0.9*.tsv"
    )
)

if len(dm_files) == 0:
    raise FileNotFoundError(
        f"Cannot find differential modification result file in {DM_ANALYSIS_DIR}"
    )

DM_FILE = dm_files[0]

print("Using differential modification file:")
print(DM_FILE)


# ============================================================
# Helper functions
# ============================================================

def read_polya_result_file(path, delta_multiplier=1):
    """
    Read your differential polyA result table.

    Expected structure, based on your example:

                    count     count      median      median      median      u       p-value     FDR
    group           control   treatment  control     diff        treatment
    contig
    AIPGENE20026    89813     47641      59.0        -5.02       53.98       ...

    This function returns a flattened table with:

        transcript_id
        count_control
        count_treatment
        median_polyA_control
        median_polyA_treatment
        median_diff_treatment_minus_control
        delta_polyA_sym_minus_apo
        p_value
        FDR
    """

    print("Reading polyA file:")
    print(path)

    # First attempt: pandas MultiIndex-style header
    try:
        df = pd.read_csv(
            path,
            sep="\t",
            header=[0, 1],
            index_col=0
        )

        df = df.reset_index()

        flat_cols = []
        for c in df.columns:
            if isinstance(c, tuple):
                c1 = str(c[0]).strip()
                c2 = str(c[1]).strip()

                if c1.startswith("Unnamed") and c2.startswith("Unnamed"):
                    flat_cols.append("transcript_id")
                elif c2.startswith("Unnamed") or c2 == "":
                    flat_cols.append(c1)
                elif c1.startswith("Unnamed") or c1 == "":
                    flat_cols.append(c2)
                else:
                    flat_cols.append(f"{c1}_{c2}")
            else:
                flat_cols.append(str(c).strip())

        df.columns = flat_cols

    except Exception as e:
        print("MultiIndex read failed; trying regular read.")
        print(e)

        df = pd.read_csv(path, sep="\t")
        df.columns = [str(c).strip() for c in df.columns]

    # Clean column names
    df.columns = [
        str(c)
        .strip()
        .replace(" ", "_")
        .replace("/", "_")
        .replace("-", "_")
        for c in df.columns
    ]

    # Drop possible header rows that got read as data
    first_col = df.columns[0]
    df = df[~df[first_col].astype(str).isin(["contig", "group", "nan"])].copy()

    # Helper: robust column finder
    def find_col(possible_names, exact=True):
        lower_map = {c.lower(): c for c in df.columns}

        for name in possible_names:
            if name.lower() in lower_map:
                return lower_map[name.lower()]

        if not exact:
            for col in df.columns:
                col_low = col.lower()
                for name in possible_names:
                    if name.lower() in col_low:
                        return col

        return None

    transcript_col = find_col(
        ["transcript_id", "contig", "index"],
        exact=True
    )

    if transcript_col is None:
        transcript_col = df.columns[0]

    count_control_col = find_col(
        ["count_control", "count.control"],
        exact=False
    )

    count_treatment_col = find_col(
        ["count_treatment", "count.treatment"],
        exact=False
    )

    median_control_col = find_col(
        ["median_control", "median.control"],
        exact=False
    )

    median_treatment_col = find_col(
        ["median_treatment", "median.treatment"],
        exact=False
    )

    median_diff_col = find_col(
        ["median_diff", "median.diff", "diff"],
        exact=False
    )

    p_col = find_col(
        ["p_value", "p.value", "pvalue", "p"],
        exact=False
    )

    fdr_col = find_col(
        ["FDR", "fdr", "padj", "qvalue", "q_value"],
        exact=False
    )

    if median_diff_col is None:
        raise ValueError(
            "Could not find the median diff column. "
            f"Detected columns are: {df.columns.tolist()}"
        )

    if fdr_col is None:
        raise ValueError(
            "Could not find the FDR column. "
            f"Detected columns are: {df.columns.tolist()}"
        )

    out = pd.DataFrame()
    out["transcript_id"] = df[transcript_col].astype(str)

    if count_control_col is not None:
        out["count_control"] = pd.to_numeric(
            df[count_control_col],
            errors="coerce"
        )

    if count_treatment_col is not None:
        out["count_treatment"] = pd.to_numeric(
            df[count_treatment_col],
            errors="coerce"
        )

    if median_control_col is not None:
        out["median_polyA_control"] = pd.to_numeric(
            df[median_control_col],
            errors="coerce"
        )

    if median_treatment_col is not None:
        out["median_polyA_treatment"] = pd.to_numeric(
            df[median_treatment_col],
            errors="coerce"
        )

    out["median_diff_treatment_minus_control"] = pd.to_numeric(
        df[median_diff_col],
        errors="coerce"
    )

    out["delta_polyA_sym_minus_apo"] = (
        out["median_diff_treatment_minus_control"] * delta_multiplier
    )

    if p_col is not None:
        out["p_value"] = pd.to_numeric(
            df[p_col],
            errors="coerce"
        )

    out["FDR"] = pd.to_numeric(
        df[fdr_col],
        errors="coerce"
    )

    out = out.dropna(
        subset=[
            "transcript_id",
            "delta_polyA_sym_minus_apo",
            "FDR"
        ]
    )

    out = out.drop_duplicates("transcript_id").copy()

    return out


def classify_polyA(row):
    """
    Classify differential polyA direction.
    Positive delta means longer poly(A) in symbiotic.
    Negative delta means longer poly(A) in aposymbiotic.
    """

    is_sig = (
        row["FDR"] < POLYA_FDR_CUTOFF and
        abs(row["delta_polyA_sym_minus_apo"]) >= POLYA_ABS_DELTA_CUTOFF
    )

    if not is_sig:
        return "Not differential poly(A)"

    if row["delta_polyA_sym_minus_apo"] > 0:
        return "Longer poly(A) in symbiotic"

    if row["delta_polyA_sym_minus_apo"] < 0:
        return "Longer poly(A) in aposymbiotic"

    return "Not differential poly(A)"


def fisher_overlap(table_df, target_col, target_value, dpa_col="is_DPA"):
    """
    Fisher's exact test for enrichment of differential polyA among target transcripts.

                    DPA      non-DPA
    target           a          b
    non-target       c          d
    """

    target = table_df[target_col] == target_value

    a = ((target) & (table_df[dpa_col])).sum()
    b = ((target) & (~table_df[dpa_col])).sum()
    c = ((~target) & (table_df[dpa_col])).sum()
    d = ((~target) & (~table_df[dpa_col])).sum()

    try:
        odds_ratio, p_value = fisher_exact(
            [[a, b], [c, d]],
            alternative="greater"
        )
    except Exception:
        odds_ratio, p_value = np.nan, np.nan

    return {
        "a_target_and_DPA": int(a),
        "b_target_not_DPA": int(b),
        "c_non_target_DPA": int(c),
        "d_non_target_not_DPA": int(d),
        "odds_ratio": odds_ratio,
        "p_value": p_value,
        "target_total": int(a + b),
        "background_total": int(a + b + c + d),
        "target_DPA_fraction": a / (a + b) if (a + b) > 0 else np.nan,
        "background_DPA_fraction": (a + c) / (a + b + c + d)
        if (a + b + c + d) > 0 else np.nan,
    }


def safe_spearman(x, y):
    valid = pd.notna(x) & pd.notna(y)

    if valid.sum() < 5:
        return np.nan, np.nan, int(valid.sum())

    rho, pval = spearmanr(x[valid], y[valid])
    return rho, pval, int(valid.sum())


def build_transcript_modification_summary(dm):
    """
    Convert site-level differential modification results into transcript-modification-level summary.

    One row = transcript_id + modification.

    Important columns:
        n_dm_sites
        dm_category
        representative_delta_mod_rate_sym_minus_apo
        min_fdr
    """

    dm = dm.copy()
    dm["is_dm_site"] = dm["direction"].isin(DIRECTION_ORDER)

    records = []

    group_cols = ["transcriptome_id", "modification"]

    for (tx, mod), sub in dm.groupby(group_cols, observed=False):
        sub = sub.copy()
        sig = sub[sub["is_dm_site"]].copy()

        n_sig = sig.shape[0]
        n_higher_sym = (sig["direction"] == "higher_in_symbiotic").sum()
        n_higher_apo = (sig["direction"] == "higher_in_aposymbiotic").sum()

        if n_sig == 0:
            dm_category = "not_significant"
            choose_df = sub.copy()
        elif n_higher_sym > 0 and n_higher_apo > 0:
            dm_category = "mixed"
            choose_df = sig.copy()
        elif n_higher_sym > 0:
            dm_category = "higher_in_symbiotic"
            choose_df = sig.copy()
        else:
            dm_category = "higher_in_aposymbiotic"
            choose_df = sig.copy()

        # Representative site: largest absolute modification-rate change
        choose_df = choose_df.dropna(
            subset=["delta_mean_mod_rate_sym_minus_apo"]
        )

        if choose_df.shape[0] == 0:
            rep = sub.iloc[0]
        else:
            rep_idx = choose_df["delta_mean_mod_rate_sym_minus_apo"].abs().idxmax()
            rep = choose_df.loc[rep_idx]

        records.append(
            {
                "transcript_id": tx,
                "modification": mod,
                "n_tested_sites": sub.shape[0],
                "n_dm_sites": n_sig,
                "n_dm_sites_higher_symbiotic": int(n_higher_sym),
                "n_dm_sites_higher_aposymbiotic": int(n_higher_apo),
                "dm_category": dm_category,
                "dm_category_label": DIRECTION_LABELS[dm_category],
                "representative_site_id": rep.get("site_id", np.nan),
                "representative_site": rep.get("site", np.nan),
                "representative_motif": rep.get("motif", np.nan),
                "representative_delta_mod_rate_sym_minus_apo": rep.get(
                    "delta_mean_mod_rate_sym_minus_apo",
                    np.nan
                ),
                "representative_fdr": rep.get("fdr_by_modification", np.nan),
                "min_fdr": sub["fdr_by_modification"].min(),
                "max_abs_delta_mod_rate": sub[
                    "delta_mean_mod_rate_sym_minus_apo"
                ].abs().max(),
                "mean_delta_mod_rate_all_sites": sub[
                    "delta_mean_mod_rate_sym_minus_apo"
                ].mean(),
            }
        )

    return pd.DataFrame(records)


# ============================================================
# Load data
# ============================================================

print("Loading differential modification results:")
dm = pd.read_csv(DM_FILE, sep="\t")

required_dm_cols = [
    "site_id",
    "transcriptome_id",
    "modification",
    "site",
    "motif",
    "direction",
    "fdr_by_modification",
    "delta_mean_mod_rate_sym_minus_apo",
]

missing_dm_cols = [c for c in required_dm_cols if c not in dm.columns]

if missing_dm_cols:
    raise ValueError(
        f"Differential modification file is missing columns: {missing_dm_cols}"
    )

print("DM shape:", dm.shape)

print("Loading and formatting polyA results:")
polyA = read_polya_result_file(
    POLYA_FILE,
    delta_multiplier=POLYA_DELTA_MULTIPLIER
)

polyA["polyA_category"] = polyA.apply(classify_polyA, axis=1)
polyA["is_DPA"] = polyA["polyA_category"] != "Not differential poly(A)"

polyA_out = os.path.join(
    OUT_DIR,
    "Formatted_differential_polyA_table_for_integration.tsv"
)

polyA.to_csv(polyA_out, sep="\t", index=False)

print("Formatted polyA table saved:")
print(polyA_out)

print("\nPolyA category counts:")
print(polyA["polyA_category"].value_counts())


# ============================================================
# Build transcript-modification summary
# ============================================================

print("\nBuilding transcript-modification summary...")

txmod = build_transcript_modification_summary(dm)

txmod_out = os.path.join(
    OUT_DIR,
    "Transcript_modification_summary_for_polyA_integration.tsv"
)

txmod.to_csv(txmod_out, sep="\t", index=False)

print("Saved transcript-modification summary:")
print(txmod_out)
print("txmod shape:", txmod.shape)


# ============================================================
# Merge with polyA
# ============================================================

merged = txmod.merge(
    polyA,
    on="transcript_id",
    how="inner"
)

merged_out = os.path.join(
    OUT_DIR,
    "Differential_modification_polyA_merged_transcript_modification_level.tsv"
)

merged.to_csv(merged_out, sep="\t", index=False)

print("\nSaved merged DM × polyA table:")
print(merged_out)
print("Merged shape:", merged.shape)


# ============================================================
# Overall transcript-level table across all modifications
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
    polyA,
    on="transcript_id",
    how="inner"
)

overall_merged_out = os.path.join(
    OUT_DIR,
    "Differential_modification_polyA_merged_transcript_overall.tsv"
)

overall_merged.to_csv(overall_merged_out, sep="\t", index=False)

print("Saved overall transcript-level DM × polyA table:")
print(overall_merged_out)


# ============================================================
# Spearman correlation:
# delta modification rate vs delta polyA
# ============================================================

corr_records = []

for mod in MODS:
    sub = merged[merged["modification"] == mod].copy()

    if sub.shape[0] == 0:
        continue

    for subset_name, subset_df in [
        ("all_tested_transcripts", sub),
        ("differentially_modified_transcripts", sub[sub["n_dm_sites"] > 0]),
    ]:
        rho, p_value, n = safe_spearman(
            subset_df["representative_delta_mod_rate_sym_minus_apo"],
            subset_df["delta_polyA_sym_minus_apo"]
        )

        corr_records.append(
            {
                "modification": mod,
                "subset": subset_name,
                "n_transcripts": n,
                "spearman_rho": rho,
                "p_value": p_value,
            }
        )

corr = pd.DataFrame(corr_records)

corr["fdr"] = np.nan

for subset_name, idx in corr.groupby("subset").groups.items():
    idx = list(idx)
    pvals = corr.loc[idx, "p_value"]
    valid = pvals.notna()

    if valid.sum() > 0:
        valid_idx = pvals[valid].index
        corr.loc[valid_idx, "fdr"] = multipletests(
            pvals[valid],
            method="fdr_bh"
        )[1]

corr_out = os.path.join(
    OUT_DIR,
    "Spearman_correlation_delta_modification_vs_delta_polyA.tsv"
)

corr.to_csv(corr_out, sep="\t", index=False)

print("\nSaved Spearman correlation table:")
print(corr_out)


# ============================================================
# Fisher enrichment:
# DM transcripts enriched among differential-polyA transcripts
# ============================================================

enrichment_records = []

# Overall: any modification
tmp = overall_merged.copy()
tmp["target_any_DM"] = tmp["has_DM"]

res = fisher_overlap(
    table_df=tmp,
    target_col="target_any_DM",
    target_value=True,
    dpa_col="is_DPA"
)

res.update(
    {
        "test_level": "any_modification",
        "modification": "any",
        "direction": "any_DM",
        "direction_label": "Any differential modification",
    }
)

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
        dpa_col="is_DPA"
    )

    res.update(
        {
            "test_level": "modification",
            "modification": mod,
            "direction": "any_DM",
            "direction_label": "Any differential modification",
        }
    )

    enrichment_records.append(res)

    for direction in ["higher_in_symbiotic", "higher_in_aposymbiotic", "mixed"]:
        sub["target_direction"] = sub["dm_category"] == direction

        if sub["target_direction"].sum() == 0:
            continue

        res = fisher_overlap(
            table_df=sub,
            target_col="target_direction",
            target_value=True,
            dpa_col="is_DPA"
        )

        res.update(
            {
                "test_level": "modification_direction",
                "modification": mod,
                "direction": direction,
                "direction_label": DIRECTION_LABELS[direction],
            }
        )

        enrichment_records.append(res)

enrichment = pd.DataFrame(enrichment_records)

enrichment["fdr"] = np.nan

valid = enrichment["p_value"].notna()

if valid.sum() > 0:
    enrichment.loc[valid, "fdr"] = multipletests(
        enrichment.loc[valid, "p_value"],
        method="fdr_bh"
    )[1]

enrichment["log2_odds_ratio"] = np.log2(
    enrichment["odds_ratio"].replace(0, np.nan)
)

enrichment["minus_log10_fdr"] = -np.log10(
    enrichment["fdr"].replace(0, np.nan)
)

enrichment_out = os.path.join(
    OUT_DIR,
    "DM_transcript_enrichment_among_differential_polyA_Fisher_tests.tsv"
)

enrichment.to_csv(enrichment_out, sep="\t", index=False)

print("\nSaved enrichment table:")
print(enrichment_out)


# ============================================================
# Mann-Whitney tests:
# delta polyA distribution by DM category
# ============================================================

mw_records = []

for mod in MODS:
    sub = merged[merged["modification"] == mod].copy()

    if sub.shape[0] == 0:
        continue

    baseline = sub[sub["dm_category"] == "not_significant"][
        "delta_polyA_sym_minus_apo"
    ].dropna()

    for category in [
        "higher_in_symbiotic",
        "higher_in_aposymbiotic",
        "mixed",
    ]:
        target = sub[sub["dm_category"] == category][
            "delta_polyA_sym_minus_apo"
        ].dropna()

        if len(target) >= 3 and len(baseline) >= 3:
            stat, p_value = mannwhitneyu(
                target,
                baseline,
                alternative="two-sided"
            )
        else:
            stat, p_value = np.nan, np.nan

        mw_records.append(
            {
                "modification": mod,
                "dm_category": category,
                "dm_category_label": DIRECTION_LABELS[category],
                "n_target": len(target),
                "n_baseline": len(baseline),
                "median_delta_polyA_target": target.median()
                if len(target) > 0 else np.nan,
                "median_delta_polyA_baseline": baseline.median()
                if len(baseline) > 0 else np.nan,
                "mannwhitney_u": stat,
                "p_value": p_value,
            }
        )

mw = pd.DataFrame(mw_records)

mw["fdr"] = np.nan

valid = mw["p_value"].notna()

if valid.sum() > 0:
    mw.loc[valid, "fdr"] = multipletests(
        mw.loc[valid, "p_value"],
        method="fdr_bh"
    )[1]

mw_out = os.path.join(
    OUT_DIR,
    "MannWhitney_delta_polyA_by_DM_category.tsv"
)

mw.to_csv(mw_out, sep="\t", index=False)

print("Saved Mann-Whitney test table:")
print(mw_out)


# ============================================================
# Plot 1. Delta modification rate vs delta polyA
# ============================================================

sns.set_theme(style="whitegrid", font_scale=1.1)

scatter_df = merged.copy()

dm_rows = scatter_df[scatter_df["n_dm_sites"] > 0]
non_dm_rows = scatter_df[scatter_df["n_dm_sites"] == 0]

# Downsample non-DM rows for plotting only
if non_dm_rows.shape[0] > 30000:
    non_dm_rows = non_dm_rows.sample(n=30000, random_state=123)

scatter_df = pd.concat([dm_rows, non_dm_rows], ignore_index=True)

scatter_df["plot_category"] = np.where(
    scatter_df["n_dm_sites"] > 0,
    scatter_df["dm_category_label"],
    "Not significant"
)

plot_order = [
    "Not significant",
    "Higher modification in aposymbiotic",
    "Higher modification in symbiotic",
    "Mixed direction",
]

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

    for cat in plot_order:
        sub = data[data["plot_category"] == cat]

        if sub.shape[0] == 0:
            continue

        ax.scatter(
            sub["representative_delta_mod_rate_sym_minus_apo"],
            sub["delta_polyA_sym_minus_apo"],
            s=8 if cat == "Not significant" else 14,
            alpha=0.25 if cat == "Not significant" else 0.75,
            color=DIRECTION_COLORS.get(cat, "gray"),
            edgecolors="none",
        )

    ax.axhline(0, color="black", linewidth=0.7)
    ax.axvline(0, color="black", linewidth=0.7)

g.map_dataframe(scatter_panel)

g.set_axis_labels(
    "Representative Δ modification rate\nSymbiotic - Aposymbiotic",
    "Δ poly(A) tail length\nSymbiotic - Aposymbiotic"
)

g.set_titles("{col_name}")

handles = [
    plt.Line2D(
        [0],
        [0],
        marker="o",
        color="w",
        label=k,
        markerfacecolor=v,
        markersize=7
    )
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

scatter_png = os.path.join(
    PLOT_DIR,
    "Delta_modification_rate_vs_delta_polyA_by_modification.png"
)

scatter_pdf = os.path.join(
    PLOT_DIR,
    "Delta_modification_rate_vs_delta_polyA_by_modification.pdf"
)

plt.savefig(scatter_png, dpi=300, bbox_inches="tight")
plt.savefig(scatter_pdf, bbox_inches="tight")
plt.close()


# ============================================================
# Plot 2. Delta polyA by DM category
# ============================================================

box_df = merged.copy()

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

plt.figure(figsize=(10, 5))

ax = sns.boxplot(
    data=box_df,
    x="modification",
    y="delta_polyA_sym_minus_apo",
    hue="dm_category_label",
    order=MODS,
    showfliers=False,
    palette=DIRECTION_COLORS
)

ax.axhline(0, color="black", linewidth=0.8)

plt.xlabel("Modification type")
plt.ylabel("Δ poly(A) tail length\nSymbiotic - Aposymbiotic")
plt.xticks(rotation=45, ha="right")

ax.legend(
    title="DM category",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    frameon=False
)

plt.tight_layout()

box_png = os.path.join(
    PLOT_DIR,
    "Delta_polyA_by_DM_category.png"
)

box_pdf = os.path.join(
    PLOT_DIR,
    "Delta_polyA_by_DM_category.pdf"
)

plt.savefig(box_png, dpi=300, bbox_inches="tight")
plt.savefig(box_pdf, bbox_inches="tight")
plt.close()


# ============================================================
# Plot 3. Fisher enrichment odds ratios
# ============================================================

forest_df = enrichment[
    (enrichment["test_level"] == "modification") &
    (enrichment["direction"] == "any_DM")
].copy()

forest_df["modification"] = pd.Categorical(
    forest_df["modification"],
    categories=MODS,
    ordered=True
)

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

plt.xlabel("log2 odds ratio for differential-poly(A) enrichment")
plt.ylabel("Modification type")
plt.title(
    "Enrichment of differential-poly(A) transcripts\n"
    "among differentially modified transcripts"
)

ax.legend(
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    frameon=False
)

plt.tight_layout()

forest_png = os.path.join(
    PLOT_DIR,
    "DM_polyA_enrichment_log2OR_by_modification.png"
)

forest_pdf = os.path.join(
    PLOT_DIR,
    "DM_polyA_enrichment_log2OR_by_modification.pdf"
)

plt.savefig(forest_png, dpi=300, bbox_inches="tight")
plt.savefig(forest_pdf, bbox_inches="tight")
plt.close()


# ============================================================
# Plot 4. PolyA category composition among DM transcripts
# ============================================================

composition_df = merged[merged["n_dm_sites"] > 0].copy()

composition = (
    composition_df
    .groupby(["modification", "polyA_category"], observed=False)
    .size()
    .reset_index(name="n_transcripts")
)

composition["total"] = (
    composition
    .groupby("modification", observed=False)["n_transcripts"]
    .transform("sum")
)

composition["fraction"] = composition["n_transcripts"] / composition["total"]

composition["polyA_category"] = pd.Categorical(
    composition["polyA_category"],
    categories=[
        "Longer poly(A) in aposymbiotic",
        "Not differential poly(A)",
        "Longer poly(A) in symbiotic",
    ],
    ordered=True
)

composition_out = os.path.join(
    OUT_DIR,
    "PolyA_category_composition_among_DM_transcripts.tsv"
)

composition.to_csv(composition_out, sep="\t", index=False)

plt.figure(figsize=(9, 5))

ax = sns.barplot(
    data=composition,
    x="modification",
    y="fraction",
    hue="polyA_category",
    order=MODS,
    palette=POLYA_COLORS
)

plt.xlabel("Modification type")
plt.ylabel("Fraction of differentially modified transcripts")
plt.xticks(rotation=45, ha="right")

ax.legend(
    title="Poly(A) category",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    frameon=False
)

plt.tight_layout()

composition_png = os.path.join(
    PLOT_DIR,
    "PolyA_category_composition_among_DM_transcripts.png"
)

composition_pdf = os.path.join(
    PLOT_DIR,
    "PolyA_category_composition_among_DM_transcripts.pdf"
)

plt.savefig(composition_png, dpi=300, bbox_inches="tight")
plt.savefig(composition_pdf, bbox_inches="tight")
plt.close()


# ============================================================
# Save compact manuscript summary numbers
# ============================================================

summary_records = []

for mod in MODS:
    sub = merged[merged["modification"] == mod].copy()

    dm_sub = sub[sub["n_dm_sites"] > 0]
    non_dm_sub = sub[sub["n_dm_sites"] == 0]

    summary_records.append(
        {
            "modification": mod,
            "n_transcripts_in_overlap_background": sub["transcript_id"].nunique(),
            "n_differentially_modified_transcripts": dm_sub["transcript_id"].nunique(),
            "n_non_differentially_modified_transcripts": non_dm_sub["transcript_id"].nunique(),
            "n_DPA_transcripts_among_DM": dm_sub["is_DPA"].sum(),
            "fraction_DPA_transcripts_among_DM": dm_sub["is_DPA"].mean()
            if dm_sub.shape[0] > 0 else np.nan,
            "median_delta_polyA_DM": dm_sub["delta_polyA_sym_minus_apo"].median(),
            "median_delta_polyA_nonDM": non_dm_sub["delta_polyA_sym_minus_apo"].median(),
        }
    )

summary_df = pd.DataFrame(summary_records)

summary_out = os.path.join(
    OUT_DIR,
    "DM_polyA_integration_summary_by_modification.tsv"
)

summary_df.to_csv(summary_out, sep="\t", index=False)


# ============================================================
# Done
# ============================================================

print("\nFinished DM × polyA integration.")
print("\nMain outputs:")
print(merged_out)
print(overall_merged_out)
print(corr_out)
print(enrichment_out)
print(mw_out)
print(summary_out)

print("\nPlots saved to:")
print(PLOT_DIR)

print("\nPlease check direction:")
print("delta_polyA_sym_minus_apo > 0 means longer poly(A) in symbiotic.")
print("delta_polyA_sym_minus_apo < 0 means longer poly(A) in aposymbiotic.")
print(f"Current POLYA_DELTA_MULTIPLIER = {POLYA_DELTA_MULTIPLIER}")
