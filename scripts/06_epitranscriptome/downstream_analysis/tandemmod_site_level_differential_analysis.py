#!/usr/bin/env python

import os
import warnings
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm


# ============================================================
# User settings
# ============================================================

TANDEMMOD_DIR = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod"

SITE_DIR = f"{TANDEMMOD_DIR}/TandemMod_site_level_prediction"

OUT_DIR = f"{TANDEMMOD_DIR}/TandemMod_site_level_differential_analysis_p0p9"

os.makedirs(OUT_DIR, exist_ok=True)

PLOT_DIR = f"{OUT_DIR}/plots"
os.makedirs(PLOT_DIR, exist_ok=True)


# -----------------------------
# Samples and conditions
# -----------------------------

SAMPLES = ["A1", "A2", "A3", "H1", "H2", "H3"]

CONDITION_MAP = {
    "A1": "Aposymbiotic",
    "A2": "Aposymbiotic",
    "A3": "Aposymbiotic",
    "H1": "Symbiotic",
    "H2": "Symbiotic",
    "H3": "Symbiotic",
}

CONDITION_ORDER = ["Aposymbiotic", "Symbiotic"]


# -----------------------------
# Modifications and feature types
# -----------------------------

MODS = ["m6A", "m1A", "A_I", "m5C", "hm5C", "m7G", "G_I", "psU"]

FEATURE_TYPES = {
    "m6A": "DRACH",
    "m1A": "NNANN",
    "A_I": "NNANN",
    "m5C": "NNCNN",
    "hm5C": "NNCNN",
    "m7G": "NNGNN",
    "G_I": "NNGNN",
    "psU": "NNTNN",
}


# -----------------------------
# Analysis thresholds
# -----------------------------

# Main probability cutoff to use
P_CUTOFF = "0.9"

MOD_READS_COL = f"mod_reads_p{P_CUTOFF}"
MOD_RATE_COL = f"mod_rate_p{P_CUTOFF}"

# Coverage filter
MIN_READS_PER_SAMPLE = 10

# Site must have enough coverage in at least this many samples per condition
MIN_SAMPLES_PER_CONDITION = 2

# Differential modification cutoff
FDR_CUTOFF = 0.05
DELTA_RATE_CUTOFF = 0.20


# -----------------------------
# Colors
# -----------------------------

CONDITION_PALETTE = {
    "Aposymbiotic": "#D8B88A",
    "Symbiotic": "#F28E2B",
}

MOD_PALETTE = {
    "m6A":  "#4E79A7",
    "m1A":  "#F28E2B",
    "A_I":  "#E15759",
    "m5C":  "#76B7B2",
    "hm5C": "#59A14F",
    "m7G":  "#EDC948",
    "G_I":  "#B07AA1",
    "psU":  "#9C755F",
}


# ============================================================
# Helper functions
# ============================================================

def p_to_stars(p):
    if pd.isna(p):
        return ""
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return ""


def safe_neglog10(x):
    if pd.isna(x):
        return np.nan
    if x <= 0:
        return 300
    return -np.log10(x)


def make_site_id(row):
    return (
        str(row["modification"]) + "|" +
        str(row["transcriptome_id"]) + "|" +
        str(row["site"]) + "|" +
        str(row["motif"])
    )


def extract_gene_like_id(transcriptome_id):
    """
    Simple placeholder gene/transcript ID extraction.

    For AIPGENE-style IDs:
        AIPGENE1900 -> AIPGENE1900

    For IsoQuant novel transcript IDs:
        transcript1412.scaffold93size630182.nic -> transcript1412.scaffold93size630182.nic

    Later, replace this with a real transcript-to-gene mapping if available.
    """
    return str(transcriptome_id)


def read_one_site_file(sample, mod):
    feature = FEATURE_TYPES[mod]

    path = os.path.join(
        SITE_DIR,
        f"{sample}_{mod}_{feature}_site_level_prediction.tsv"
    )

    if not os.path.exists(path):
        print(f"WARNING: missing file: {path}")
        return None

    print(f"Reading {sample} {mod}: {path}")

    df = pd.read_csv(path, sep="\t")

    required_cols = [
        "transcriptome_id",
        "site",
        "motif",
        MOD_READS_COL,
        "total_reads",
        MOD_RATE_COL,
    ]

    missing = [c for c in required_cols if c not in df.columns]
    if len(missing) > 0:
        raise ValueError(
            f"File {path} is missing required columns: {missing}"
        )

    df = df[required_cols].copy()

    df["sample"] = sample
    df["condition"] = CONDITION_MAP[sample]
    df["modification"] = mod
    df["feature_type"] = feature

    df["site"] = df["site"].astype(str)
    df["motif"] = df["motif"].astype(str)

    df[MOD_READS_COL] = pd.to_numeric(df[MOD_READS_COL], errors="coerce")
    df["total_reads"] = pd.to_numeric(df["total_reads"], errors="coerce")
    df[MOD_RATE_COL] = pd.to_numeric(df[MOD_RATE_COL], errors="coerce")

    df = df.dropna(subset=[MOD_READS_COL, "total_reads", MOD_RATE_COL])

    df[MOD_READS_COL] = df[MOD_READS_COL].astype(int)
    df["total_reads"] = df["total_reads"].astype(int)

    # In case of rare duplicated rows, aggregate them.
    df = (
        df
        .groupby(
            [
                "sample",
                "condition",
                "modification",
                "feature_type",
                "transcriptome_id",
                "site",
                "motif",
            ],
            as_index=False
        )
        .agg(
            {
                MOD_READS_COL: "sum",
                "total_reads": "sum",
            }
        )
    )

    df[MOD_RATE_COL] = df[MOD_READS_COL] / df["total_reads"]

    df["unmod_reads"] = df["total_reads"] - df[MOD_READS_COL]

    df["site_id"] = df.apply(make_site_id, axis=1)
    df["gene_like_id"] = df["transcriptome_id"].apply(extract_gene_like_id)

    return df


def run_quasibinomial_test(site_df):
    """
    Differential modification test for one site.

    Model:
        modified_reads / total_reads ~ condition

    Implemented as binomial GLM with total_reads as frequency weights.
    Scale is estimated using Pearson chi-square to approximate quasibinomial behavior.

    Returns:
        coefficient, p-value, model status
    """

    d = site_df.copy()

    d = d[d["total_reads"] > 0].copy()

    if d["condition"].nunique() < 2:
        return np.nan, np.nan, "only_one_condition"

    # Need at least two samples per condition after filtering
    n_by_cond = d.groupby("condition")["sample"].nunique().to_dict()

    if n_by_cond.get("Aposymbiotic", 0) < 2 or n_by_cond.get("Symbiotic", 0) < 2:
        return np.nan, np.nan, "too_few_samples"

    d["condition_binary"] = (d["condition"] == "Symbiotic").astype(int)

    y = d[MOD_READS_COL] / d["total_reads"]

    X = sm.add_constant(d["condition_binary"])

    try:
        model = sm.GLM(
            y,
            X,
            family=sm.families.Binomial(),
            freq_weights=d["total_reads"]
        )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fit = model.fit(scale="X2", maxiter=100)

        coef = fit.params["condition_binary"]
        pval = fit.pvalues["condition_binary"]

        return coef, pval, "quasibinomial_glm"

    except Exception:
        return np.nan, np.nan, "glm_failed"


def run_fisher_pooled(site_df):
    """
    Fallback / supporting test using pooled read counts.

    This does not model replicate structure, so use mainly as backup.
    """

    d = site_df.copy()

    apo = d[d["condition"] == "Aposymbiotic"]
    sym = d[d["condition"] == "Symbiotic"]

    apo_mod = apo[MOD_READS_COL].sum()
    apo_unmod = apo["unmod_reads"].sum()

    sym_mod = sym[MOD_READS_COL].sum()
    sym_unmod = sym["unmod_reads"].sum()

    table = np.array([
        [sym_mod, sym_unmod],
        [apo_mod, apo_unmod],
    ])

    try:
        odds, pval = fisher_exact(table, alternative="two-sided")
    except Exception:
        odds, pval = np.nan, np.nan

    return odds, pval


# ============================================================
# Step 1. Combine all site-level files
# ============================================================

all_site_tables = []

for sample in SAMPLES:
    for mod in MODS:
        one = read_one_site_file(sample, mod)
        if one is not None:
            all_site_tables.append(one)

combined = pd.concat(all_site_tables, ignore_index=True)

combined["sample"] = pd.Categorical(
    combined["sample"],
    categories=SAMPLES,
    ordered=True
)

combined["condition"] = pd.Categorical(
    combined["condition"],
    categories=CONDITION_ORDER,
    ordered=True
)

combined["modification"] = pd.Categorical(
    combined["modification"],
    categories=MODS,
    ordered=True
)

combined_out = os.path.join(
    OUT_DIR,
    f"TandemMod_combined_site_level_predictions_p{P_CUTOFF}.tsv"
)

combined.to_csv(combined_out, sep="\t", index=False)

print("\nSaved combined site-level table:")
print(combined_out)
print(combined.shape)


# ============================================================
# Step 2. Coverage filter
# ============================================================

combined["covered"] = combined["total_reads"] >= MIN_READS_PER_SAMPLE

coverage_summary = (
    combined
    .groupby(["site_id", "condition"], observed=False)
    .agg(
        n_samples_covered=("covered", "sum"),
        n_samples_total=("sample", "nunique"),
        mean_total_reads=("total_reads", "mean"),
        median_total_reads=("total_reads", "median")
    )
    .reset_index()
)

coverage_pivot = coverage_summary.pivot(
    index="site_id",
    columns="condition",
    values="n_samples_covered"
).fillna(0)

coverage_pivot = coverage_pivot.reset_index()

coverage_pivot["pass_coverage"] = (
    (coverage_pivot.get("Aposymbiotic", 0) >= MIN_SAMPLES_PER_CONDITION) &
    (coverage_pivot.get("Symbiotic", 0) >= MIN_SAMPLES_PER_CONDITION)
)

passing_site_ids = set(
    coverage_pivot.loc[
        coverage_pivot["pass_coverage"],
        "site_id"
    ]
)

filtered = combined[combined["site_id"].isin(passing_site_ids)].copy()

filtered_out = os.path.join(
    OUT_DIR,
    f"TandemMod_filtered_site_level_predictions_p{P_CUTOFF}_minReads{MIN_READS_PER_SAMPLE}_minSamples{MIN_SAMPLES_PER_CONDITION}.tsv"
)

filtered.to_csv(filtered_out, sep="\t", index=False)

coverage_out = os.path.join(
    OUT_DIR,
    f"TandemMod_site_coverage_summary_p{P_CUTOFF}.tsv"
)

coverage_pivot.to_csv(coverage_out, sep="\t", index=False)

print("\nSaved filtered site-level table:")
print(filtered_out)
print(filtered.shape)

print("\nNumber of sites before filtering:", combined["site_id"].nunique())
print("Number of sites after filtering:", filtered["site_id"].nunique())


# ============================================================
# Step 3. Differential modification test
# ============================================================

results = []

site_meta = (
    filtered
    .drop_duplicates("site_id")
    [
        [
            "site_id",
            "modification",
            "feature_type",
            "transcriptome_id",
            "gene_like_id",
            "site",
            "motif",
        ]
    ]
)

n_sites = filtered["site_id"].nunique()

print(f"\nRunning differential modification test for {n_sites} sites...")

for idx, (site_id, d) in enumerate(filtered.groupby("site_id"), start=1):

    if idx % 1000 == 0:
        print(f"Processed {idx}/{n_sites} sites")

    apo = d[d["condition"] == "Aposymbiotic"]
    sym = d[d["condition"] == "Symbiotic"]

    apo_mean_rate = apo[MOD_RATE_COL].mean()
    sym_mean_rate = sym[MOD_RATE_COL].mean()

    apo_median_rate = apo[MOD_RATE_COL].median()
    sym_median_rate = sym[MOD_RATE_COL].median()

    apo_pooled_mod = apo[MOD_READS_COL].sum()
    apo_pooled_total = apo["total_reads"].sum()
    sym_pooled_mod = sym[MOD_READS_COL].sum()
    sym_pooled_total = sym["total_reads"].sum()

    apo_pooled_rate = (
        apo_pooled_mod / apo_pooled_total
        if apo_pooled_total > 0 else np.nan
    )

    sym_pooled_rate = (
        sym_pooled_mod / sym_pooled_total
        if sym_pooled_total > 0 else np.nan
    )

    delta_mean_rate = sym_mean_rate - apo_mean_rate
    delta_pooled_rate = sym_pooled_rate - apo_pooled_rate

    coef, p_glm, test_status = run_quasibinomial_test(d)

    odds_fisher, p_fisher = run_fisher_pooled(d)

    # Main p-value: GLM if successful, otherwise Fisher fallback
    if not pd.isna(p_glm):
        p_main = p_glm
        main_test = test_status
    else:
        p_main = p_fisher
        main_test = "fisher_pooled_fallback"

    results.append(
        {
            "site_id": site_id,

            "n_samples_apo": apo["sample"].nunique(),
            "n_samples_symbiotic": sym["sample"].nunique(),

            "n_samples_apo_cov": (apo["total_reads"] >= MIN_READS_PER_SAMPLE).sum(),
            "n_samples_symbiotic_cov": (sym["total_reads"] >= MIN_READS_PER_SAMPLE).sum(),

            "apo_mean_mod_rate": apo_mean_rate,
            "symbiotic_mean_mod_rate": sym_mean_rate,
            "delta_mean_mod_rate_sym_minus_apo": delta_mean_rate,

            "apo_median_mod_rate": apo_median_rate,
            "symbiotic_median_mod_rate": sym_median_rate,

            "apo_pooled_mod_reads": apo_pooled_mod,
            "apo_pooled_total_reads": apo_pooled_total,
            "apo_pooled_mod_rate": apo_pooled_rate,

            "symbiotic_pooled_mod_reads": sym_pooled_mod,
            "symbiotic_pooled_total_reads": sym_pooled_total,
            "symbiotic_pooled_mod_rate": sym_pooled_rate,

            "delta_pooled_mod_rate_sym_minus_apo": delta_pooled_rate,

            "glm_coef_condition_symbiotic": coef,
            "p_glm_quasibinomial": p_glm,

            "fisher_odds_ratio_pooled": odds_fisher,
            "p_fisher_pooled": p_fisher,

            "p_value": p_main,
            "main_test": main_test,
        }
    )

results_df = pd.DataFrame(results)

results_df = results_df.merge(site_meta, on="site_id", how="left")

# FDR correction within each modification type
results_df["fdr_by_modification"] = np.nan

for mod in MODS:
    mask = results_df["modification"] == mod
    pvals = results_df.loc[mask, "p_value"]

    valid = pvals.notna()

    if valid.sum() > 0:
        results_df.loc[mask & valid, "fdr_by_modification"] = multipletests(
            pvals[valid],
            method="fdr_bh"
        )[1]

# Also global FDR across all tested sites
valid_global = results_df["p_value"].notna()

results_df["fdr_global"] = np.nan

if valid_global.sum() > 0:
    results_df.loc[valid_global, "fdr_global"] = multipletests(
        results_df.loc[valid_global, "p_value"],
        method="fdr_bh"
    )[1]

# Significance calls
results_df["direction"] = "not_significant"

results_df.loc[
    (
        (results_df["fdr_by_modification"] < FDR_CUTOFF) &
        (results_df["delta_mean_mod_rate_sym_minus_apo"] >= DELTA_RATE_CUTOFF)
    ),
    "direction"
] = "higher_in_symbiotic"

results_df.loc[
    (
        (results_df["fdr_by_modification"] < FDR_CUTOFF) &
        (results_df["delta_mean_mod_rate_sym_minus_apo"] <= -DELTA_RATE_CUTOFF)
    ),
    "direction"
] = "higher_in_aposymbiotic"

results_df["significant"] = results_df["direction"] != "not_significant"

results_df["neglog10_fdr"] = results_df["fdr_by_modification"].apply(safe_neglog10)

results_out = os.path.join(
    OUT_DIR,
    f"TandemMod_differential_modification_results_p{P_CUTOFF}_minReads{MIN_READS_PER_SAMPLE}_minSamples{MIN_SAMPLES_PER_CONDITION}.tsv"
)

results_df.to_csv(results_out, sep="\t", index=False)

sig_out = os.path.join(
    OUT_DIR,
    f"TandemMod_significant_differential_modification_sites_p{P_CUTOFF}_FDR{FDR_CUTOFF}_delta{DELTA_RATE_CUTOFF}.tsv"
)

results_df[results_df["significant"]].to_csv(sig_out, sep="\t", index=False)

print("\nSaved differential modification results:")
print(results_out)

print("\nSaved significant differential modification sites:")
print(sig_out)

print("\nSignificant site counts:")
print(results_df["direction"].value_counts())


# ============================================================
# Step 4. Summary tables
# ============================================================

summary_by_mod = (
    results_df
    .groupby(["modification", "direction"], observed=False)
    .size()
    .reset_index(name="n_sites")
)

summary_by_mod_out = os.path.join(
    OUT_DIR,
    f"TandemMod_differential_site_count_by_modification_p{P_CUTOFF}.tsv"
)

summary_by_mod.to_csv(summary_by_mod_out, sep="\t", index=False)

gene_summary = (
    results_df[results_df["significant"]]
    .groupby(["modification", "direction"], observed=False)
    .agg(
        n_sites=("site_id", "nunique"),
        n_transcripts=("transcriptome_id", "nunique"),
        n_gene_like_ids=("gene_like_id", "nunique")
    )
    .reset_index()
)

gene_summary_out = os.path.join(
    OUT_DIR,
    f"TandemMod_significant_site_transcript_gene_summary_p{P_CUTOFF}.tsv"
)

gene_summary.to_csv(gene_summary_out, sep="\t", index=False)

print("\nSaved summary tables.")


# ============================================================
# Step 5. Volcano plots
# ============================================================
results_out = os.path.join(
    OUT_DIR,
    f"TandemMod_differential_modification_results_p{P_CUTOFF}_minReads{MIN_READS_PER_SAMPLE}_minSamples{MIN_SAMPLES_PER_CONDITION}.tsv"
)

results_df=pd.read_csv(results_out, sep="\t")
sns.set_theme(style="whitegrid", font_scale=1.1)

volcano_df = results_df.copy()

volcano_df["direction_plot"] = volcano_df["direction"].map(
    {
        "higher_in_symbiotic": "Higher in symbiotic",
        "higher_in_aposymbiotic": "Higher in aposymbiotic",
        "not_significant": "Not significant",
    }
)

volcano_palette = {
    "Higher in symbiotic": "#F28E2B",
    "Higher in aposymbiotic": "#D8B88A",
    "Not significant": "lightgray",
}

g = sns.FacetGrid(
    volcano_df,
    col="modification",
    col_order=MODS,
    col_wrap=4,
    height=3.0,
    aspect=1.2,
    sharex=True,
    sharey=False
)

def volcano_panel(data, **kwargs):
    ax = plt.gca()

    sns.scatterplot(
        data=data,
        x="delta_mean_mod_rate_sym_minus_apo",
        y="neglog10_fdr",
        hue="direction_plot",
        palette=volcano_palette,
        hue_order=[
            "Higher in aposymbiotic",
            "Not significant",
            "Higher in symbiotic",
        ],
        s=16,
        linewidth=0,
        alpha=0.75,
        ax=ax,
        legend=False
    )

    ax.axvline(0, color="black", linewidth=0.8)
    ax.axvline(DELTA_RATE_CUTOFF, color="black", linestyle="--", linewidth=0.7)
    ax.axvline(-DELTA_RATE_CUTOFF, color="black", linestyle="--", linewidth=0.7)

    if FDR_CUTOFF > 0:
        ax.axhline(
            -np.log10(FDR_CUTOFF),
            color="black",
            linestyle="--",
            linewidth=0.7
        )

    ax.set_xlim(-1, 1)

g.map_dataframe(volcano_panel)

g.set_axis_labels(
    "Δ modification rate\nSymbiotic - Aposymbiotic",
    "-log10(FDR)"
)

g.set_titles("{col_name}")

# Custom legend
handles = [
    plt.Line2D(
        [0], [0],
        marker="o",
        color="w",
        label="Higher in aposymbiotic",
        markerfacecolor="#D8B88A",
        markersize=7
    ),
    plt.Line2D(
        [0], [0],
        marker="o",
        color="w",
        label="Not significant",
        markerfacecolor="lightgray",
        markersize=7
    ),
    plt.Line2D(
        [0], [0],
        marker="o",
        color="w",
        label="Higher in symbiotic",
        markerfacecolor="#F28E2B",
        markersize=7
    ),
]

g.fig.legend(
    handles=handles,
    title="Differential modification",
    bbox_to_anchor=(1.02, 0.5),
    loc="center left",
    frameon=False
)

plt.tight_layout(rect=[0, 0, 0.90, 1])

volcano_png = os.path.join(
    PLOT_DIR,
    f"Volcano_differential_modification_all_mods_p{P_CUTOFF}.png"
)

volcano_pdf = os.path.join(
    PLOT_DIR,
    f"Volcano_differential_modification_all_mods_p{P_CUTOFF}.pdf"
)

plt.savefig(volcano_png, dpi=300, bbox_inches="tight")
plt.savefig(volcano_pdf, bbox_inches="tight")
plt.close()

print("\nSaved volcano plots.")


# ============================================================
# Step 6. Bar plot: number of significant sites
# ============================================================

bar_df = (
    results_df[results_df["significant"]]
    .groupby(["modification", "direction"], observed=False)
    .size()
    .reset_index(name="n_sites")
)

# Make sure all combinations exist
all_combo = pd.MultiIndex.from_product(
    [
        MODS,
        ["higher_in_aposymbiotic", "higher_in_symbiotic"]
    ],
    names=["modification", "direction"]
).to_frame(index=False)

bar_df = all_combo.merge(
    bar_df,
    on=["modification", "direction"],
    how="left"
)

bar_df["n_sites"] = bar_df["n_sites"].fillna(0)

bar_df["direction_label"] = bar_df["direction"].map(
    {
        "higher_in_aposymbiotic": "Higher in aposymbiotic",
        "higher_in_symbiotic": "Higher in symbiotic",
    }
)

plt.figure(figsize=(8, 4.5))

ax = sns.barplot(
    data=bar_df,
    x="modification",
    y="n_sites",
    hue="direction_label",
    order=MODS,
    palette={
        "Higher in aposymbiotic": "#D8B88A",
        "Higher in symbiotic": "#F28E2B",
    }
)

plt.xlabel("Modification type")
plt.ylabel("Number of differential modification sites")
plt.xticks(rotation=45, ha="right")

ax.legend(
    title="Direction",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    frameon=False
)

plt.tight_layout()

bar_png = os.path.join(
    PLOT_DIR,
    f"Barplot_number_differential_modification_sites_p{P_CUTOFF}.png"
)

bar_pdf = os.path.join(
    PLOT_DIR,
    f"Barplot_number_differential_modification_sites_p{P_CUTOFF}.pdf"
)

plt.savefig(bar_png, dpi=300, bbox_inches="tight")
plt.savefig(bar_pdf, bbox_inches="tight")
plt.close()

print("Saved differential site count bar plot.")


# ============================================================
# Step 7. Heatmap of top differential sites
# ============================================================

top_n = 50

sig_sites = results_df[results_df["significant"]].copy()
filtered_out = os.path.join(
    OUT_DIR,
    f"TandemMod_filtered_site_level_predictions_p{P_CUTOFF}_minReads{MIN_READS_PER_SAMPLE}_minSamples{MIN_SAMPLES_PER_CONDITION}.tsv"
)

filtered=pd.read_csv(filtered_out, sep="\t")


if sig_sites.shape[0] > 0:

    sig_sites = sig_sites.sort_values(
        ["fdr_by_modification", "delta_mean_mod_rate_sym_minus_apo"],
        ascending=[True, False]
    )

    top_sites = sig_sites.head(top_n)["site_id"].tolist()

    heat_df = filtered[filtered["site_id"].isin(top_sites)].copy()

    heat_df = heat_df.merge(
        results_df[
            [
                "site_id",
                "modification",
                "transcriptome_id",
                "site",
                "motif",
                "direction",
                "fdr_by_modification",
                "delta_mean_mod_rate_sym_minus_apo",
            ]
        ].drop_duplicates("site_id"),
        on=["site_id", "modification", "transcriptome_id", "site", "motif"],
        how="left"
    )

    heat_df["site_label"] = (
        heat_df["modification"].astype(str) + " | " +
        heat_df["transcriptome_id"].astype(str) + ":" +
        heat_df["site"].astype(str) + " | " +
        heat_df["motif"].astype(str)
    )

    # Keep row order according to significance
    site_label_order = (
        heat_df[["site_id", "site_label"]]
        .drop_duplicates()
        .set_index("site_id")
        .loc[top_sites, "site_label"]
        .tolist()
    )

    heat_mat = heat_df.pivot_table(
        index="site_label",
        columns="sample",
        values=MOD_RATE_COL,
        aggfunc="mean"
    )

    heat_mat = heat_mat.reindex(site_label_order)
    heat_mat = heat_mat[SAMPLES]

    plt.figure(figsize=(8, max(6, top_n * 0.18)))

    ax = sns.heatmap(
        heat_mat,
        cmap="YlOrBr",
        vmin=0,
        vmax=1,
        linewidths=0.2,
        linecolor="white",
        mask=heat_mat.isna(),
        cbar_kws={"label": f"Modification rate, p ≥ {P_CUTOFF.replace('p', '.')}"}
    )

    plt.xlabel("Sample")
    plt.ylabel("Top differential modification sites")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    heat_png = os.path.join(
        PLOT_DIR,
        f"Heatmap_top{top_n}_differential_modification_sites_p{P_CUTOFF}.png"
    )

    heat_pdf = os.path.join(
        PLOT_DIR,
        f"Heatmap_top{top_n}_differential_modification_sites_p{P_CUTOFF}.pdf"
    )

    plt.savefig(heat_png, dpi=300, bbox_inches="tight")
    plt.savefig(heat_pdf, bbox_inches="tight")
    plt.close()

    print("Saved heatmap of top differential sites.")

else:
    print("No significant sites found. Skipping heatmap.")


# ============================================================
# Step 8. Boxplots for top example sites
# ============================================================

example_n = 12

if sig_sites.shape[0] > 0:

    example_sites = sig_sites.head(example_n)["site_id"].tolist()

    example_df = filtered[filtered["site_id"].isin(example_sites)].copy()

    meta_for_label = results_df[
        [
            "site_id",
            "modification",
            "transcriptome_id",
            "site",
            "motif",
            "direction",
            "fdr_by_modification",
            "delta_mean_mod_rate_sym_minus_apo",
        ]
    ].drop_duplicates("site_id")

    example_df = example_df.merge(
        meta_for_label,
        on=["site_id", "modification", "transcriptome_id", "site", "motif"],
        how="left"
    )

    example_df["site_label"] = (
        example_df["modification"].astype(str) + "\n" +
        example_df["transcriptome_id"].astype(str) + ":" +
        example_df["site"].astype(str)
    )

    label_order = (
        example_df[["site_id", "site_label"]]
        .drop_duplicates()
        .set_index("site_id")
        .loc[example_sites, "site_label"]
        .tolist()
    )

    example_df["site_label"] = pd.Categorical(
        example_df["site_label"],
        categories=label_order,
        ordered=True
    )

    g = sns.catplot(
        data=example_df,
        x="condition",
        y=MOD_RATE_COL,
        col="site_label",
        col_wrap=4,
        kind="box",
        hue="condition",
        order=CONDITION_ORDER,
        hue_order=CONDITION_ORDER,
        palette=CONDITION_PALETTE,
        showfliers=False,
        height=3,
        aspect=0.9,
        legend=False
    )

    # Add sample points
    for ax, site_label in zip(g.axes.flatten(), label_order):
        sub = example_df[example_df["site_label"] == site_label]

        sns.stripplot(
            data=sub,
            x="condition",
            y=MOD_RATE_COL,
            order=CONDITION_ORDER,
            hue="condition",
            hue_order=CONDITION_ORDER,
            palette=CONDITION_PALETTE,
            dodge=False,
            jitter=True,
            size=4,
            edgecolor="black",
            linewidth=0.4,
            ax=ax,
            legend=False
        )

        ax.set_xlabel("")
        ax.set_ylabel("Modification rate")
        ax.tick_params(axis="x", rotation=30)

    g.set_titles("{col_name}")

    plt.tight_layout()

    box_png = os.path.join(
        PLOT_DIR,
        f"Boxplot_top{example_n}_differential_modification_sites_p{P_CUTOFF}.png"
    )

    box_pdf = os.path.join(
        PLOT_DIR,
        f"Boxplot_top{example_n}_differential_modification_sites_p{P_CUTOFF}.pdf"
    )

    plt.savefig(box_png, dpi=300, bbox_inches="tight")
    plt.savefig(box_pdf, bbox_inches="tight")
    plt.close()

    print("Saved example site boxplots.")

else:
    print("No significant sites found. Skipping boxplots.")


# ============================================================
# Step 9. Save compact final summary
# ============================================================

compact_cols = [
    "modification",
    "transcriptome_id",
    "site",
    "motif",
    "gene_like_id",
    "apo_mean_mod_rate",
    "symbiotic_mean_mod_rate",
    "delta_mean_mod_rate_sym_minus_apo",
    "apo_pooled_mod_rate",
    "symbiotic_pooled_mod_rate",
    "delta_pooled_mod_rate_sym_minus_apo",
    "p_value",
    "fdr_by_modification",
    "fdr_global",
    "direction",
    "main_test",
    "n_samples_apo_cov",
    "n_samples_symbiotic_cov",
    "apo_pooled_total_reads",
    "symbiotic_pooled_total_reads",
]

compact_out = os.path.join(
    OUT_DIR,
    f"TandemMod_differential_modification_compact_results_p{P_CUTOFF}.tsv"
)

results_df[compact_cols].sort_values(
    ["modification", "fdr_by_modification"]
).to_csv(compact_out, sep="\t", index=False)

print("\nSaved compact results:")
print(compact_out)

print("\nAnalysis finished successfully.")
