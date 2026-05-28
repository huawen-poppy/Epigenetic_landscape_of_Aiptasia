#!/usr/bin/env python

import os
import glob
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# ============================================================
# User settings
# ============================================================

TANDEMMOD_DIR = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod"

ANALYSIS_DIR = f"{TANDEMMOD_DIR}/TandemMod_site_level_differential_analysis_p0p9"

PLOT_DIR = f"{ANALYSIS_DIR}/plots_publication"
os.makedirs(PLOT_DIR, exist_ok=True)

SUMMARY_DIR = f"{ANALYSIS_DIR}/summary_for_manuscript"
os.makedirs(SUMMARY_DIR, exist_ok=True)

P_CUTOFF_TEXT = "0.9"
FDR_CUTOFF = 0.05
DELTA_RATE_CUTOFF = 0.20

MODS = ["m6A", "m1A", "A_I", "m5C", "hm5C", "m7G", "G_I", "psU"]

SAMPLES = ["A1", "A2", "A3", "H1", "H2", "H3"]

CONDITION_MAP = {
    "A1": "Aposymbiotic",
    "A2": "Aposymbiotic",
    "A3": "Aposymbiotic",
    "H1": "Symbiotic",
    "H2": "Symbiotic",
    "H3": "Symbiotic",
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
    "Not significant": "lightgray",
}

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
# Locate input files
# ============================================================

result_candidates = glob.glob(
    os.path.join(
        ANALYSIS_DIR,
        "TandemMod_differential_modification_results_p0.9*.tsv"
    )
)

if len(result_candidates) == 0:
    raise FileNotFoundError("Could not find differential modification result file.")

RESULT_FILE = result_candidates[0]

filtered_candidates = glob.glob(
    os.path.join(
        ANALYSIS_DIR,
        "TandemMod_filtered_site_level_predictions_p0.9*.tsv"
    )
)

if len(filtered_candidates) == 0:
    raise FileNotFoundError("Could not find filtered site-level prediction file.")

FILTERED_FILE = filtered_candidates[0]

print("Using result file:")
print(RESULT_FILE)

print("Using filtered site-level file:")
print(FILTERED_FILE)


# ============================================================
# Load data
# ============================================================

results = pd.read_csv(RESULT_FILE, sep="\t")
filtered = pd.read_csv(FILTERED_FILE, sep="\t")

# Make sure key columns exist
required_result_cols = [
    "site_id",
    "modification",
    "transcriptome_id",
    "site",
    "motif",
    "direction",
    "fdr_by_modification",
    "delta_mean_mod_rate_sym_minus_apo",
]

missing = [c for c in required_result_cols if c not in results.columns]
if missing:
    raise ValueError(f"Missing columns from result file: {missing}")

# Identify modification rate column in filtered table
rate_candidates = [
    c for c in filtered.columns
    if c.startswith("mod_rate_p0p9") or c == "mod_rate_p0.9"
]

if len(rate_candidates) > 0:
    MOD_RATE_COL = rate_candidates[0]
else:
    # fallback from previous script
    if "mod_rate_p0p9" in filtered.columns:
        MOD_RATE_COL = "mod_rate_p0p9"
    elif "mod_rate_p0.9" in filtered.columns:
        MOD_RATE_COL = "mod_rate_p0.9"
    else:
        MOD_RATE_COL = "mod_rate_p0.9"

print(f"Using modification rate column: {MOD_RATE_COL}")

results["modification"] = pd.Categorical(
    results["modification"],
    categories=MODS,
    ordered=True
)

filtered["sample"] = pd.Categorical(
    filtered["sample"],
    categories=SAMPLES,
    ordered=True
)

filtered["condition"] = filtered["sample"].map(CONDITION_MAP)

filtered["condition"] = pd.Categorical(
    filtered["condition"],
    categories=["Aposymbiotic", "Symbiotic"],
    ordered=True
)


# ============================================================
# Summary counts and percentages
# ============================================================

tested_by_mod = (
    results
    .groupby("modification", observed=False)
    .agg(
        n_tested_sites=("site_id", "nunique"),
        n_tested_transcripts=("transcriptome_id", "nunique"),
    )
    .reset_index()
)

sig = results[results["direction"].isin(DIRECTION_ORDER)].copy()

sig_by_mod_dir = (
    sig
    .groupby(["modification", "direction"], observed=False)
    .agg(
        n_differential_sites=("site_id", "nunique"),
        n_differential_transcripts=("transcriptome_id", "nunique"),
    )
    .reset_index()
)

all_combo = pd.MultiIndex.from_product(
    [MODS, DIRECTION_ORDER],
    names=["modification", "direction"]
).to_frame(index=False)

sig_by_mod_dir = all_combo.merge(
    sig_by_mod_dir,
    on=["modification", "direction"],
    how="left"
)

sig_by_mod_dir[["n_differential_sites", "n_differential_transcripts"]] = (
    sig_by_mod_dir[["n_differential_sites", "n_differential_transcripts"]]
    .fillna(0)
    .astype(int)
)

summary = sig_by_mod_dir.merge(
    tested_by_mod,
    on="modification",
    how="left"
)

summary["percent_of_tested_sites"] = (
    summary["n_differential_sites"] /
    summary["n_tested_sites"] * 100
)

summary["direction_label"] = summary["direction"].map(DIRECTION_LABELS)

summary_out = os.path.join(
    SUMMARY_DIR,
    "TandemMod_differential_modification_counts_and_percentages.tsv"
)

summary.to_csv(summary_out, sep="\t", index=False)

# Total summary per modification
summary_total_mod = (
    summary
    .groupby("modification", observed=False)
    .agg(
        n_tested_sites=("n_tested_sites", "first"),
        n_tested_transcripts=("n_tested_transcripts", "first"),
        n_differential_sites=("n_differential_sites", "sum"),
        n_differential_transcripts=("n_differential_transcripts", "sum"),
    )
    .reset_index()
)

summary_total_mod["percent_differential_sites"] = (
    summary_total_mod["n_differential_sites"] /
    summary_total_mod["n_tested_sites"] * 100
)

summary_total_mod_out = os.path.join(
    SUMMARY_DIR,
    "TandemMod_differential_modification_total_by_modification.tsv"
)

summary_total_mod.to_csv(summary_total_mod_out, sep="\t", index=False)

# Overall numbers
n_tested_total = results["site_id"].nunique()
n_sig_total = sig["site_id"].nunique()
percent_sig_total = n_sig_total / n_tested_total * 100 if n_tested_total > 0 else np.nan

n_higher_sym = (sig["direction"] == "higher_in_symbiotic").sum()
n_higher_apo = (sig["direction"] == "higher_in_aposymbiotic").sum()

n_sig_transcripts = sig["transcriptome_id"].nunique()

if "gene_like_id" in sig.columns:
    n_sig_genes = sig["gene_like_id"].nunique()
else:
    n_sig_genes = np.nan

top_count_row = summary_total_mod.sort_values(
    "n_differential_sites",
    ascending=False
).iloc[0]

top_percent_row = summary_total_mod.sort_values(
    "percent_differential_sites",
    ascending=False
).iloc[0]

psu_row = summary_total_mod[
    summary_total_mod["modification"].astype(str) == "psU"
]

if psu_row.shape[0] > 0:
    psu_sig = int(psu_row["n_differential_sites"].iloc[0])
    psu_percent = psu_row["percent_differential_sites"].iloc[0]
else:
    psu_sig = np.nan
    psu_percent = np.nan

numbers = {
    "N_tested_total": n_tested_total,
    "N_sig_total": n_sig_total,
    "percent_sig_total": percent_sig_total,
    "N_higher_symbiotic": n_higher_sym,
    "N_higher_aposymbiotic": n_higher_apo,
    "N_sig_transcripts": n_sig_transcripts,
    "N_sig_genes": n_sig_genes,
    "top_count_mod": top_count_row["modification"],
    "top_count_n": int(top_count_row["n_differential_sites"]),
    "top_count_percent": top_count_row["percent_differential_sites"],
    "top_percent_mod": top_percent_row["modification"],
    "top_percent_value": top_percent_row["percent_differential_sites"],
    "top_percent_n": int(top_percent_row["n_differential_sites"]),
    "psU_sig": psu_sig,
    "psU_percent_tested": psu_percent,
}

numbers_df = pd.DataFrame([numbers])

numbers_out = os.path.join(
    SUMMARY_DIR,
    "TandemMod_key_numbers_for_results_section.tsv"
)

numbers_df.to_csv(numbers_out, sep="\t", index=False)

print("\nSaved summary tables:")
print(summary_out)
print(summary_total_mod_out)
print(numbers_out)


# ============================================================
# Generate ready-to-copy result text with exact numbers
# ============================================================

def fmt_int(x):
    if pd.isna(x):
        return "NA"
    return f"{int(x):,}"

def fmt_pct(x):
    if pd.isna(x):
        return "NA"
    return f"{x:.2f}"

auto_text = f"""
Site-level differential RNA modification analysis was performed using high-confidence TandemMod calls defined by a prediction probability threshold of ≥{P_CUTOFF_TEXT}. After coverage filtering, {fmt_int(n_tested_total)} candidate modification sites were tested across eight RNA modification classes. Using an FDR threshold of {FDR_CUTOFF} and an absolute mean modification-rate difference threshold of {DELTA_RATE_CUTOFF}, {fmt_int(n_sig_total)} sites were identified as differentially modified, corresponding to {fmt_pct(percent_sig_total)}% of all tested sites. Among these, {fmt_int(n_higher_sym)} sites showed higher modification rates in symbiotic samples, whereas {fmt_int(n_higher_apo)} sites showed higher modification rates in aposymbiotic samples.

The number of differentially modified sites varied substantially among modification types. {top_count_row['modification']} showed the largest number of differential sites ({fmt_int(top_count_row['n_differential_sites'])} sites; {fmt_pct(top_count_row['percent_differential_sites'])}% of tested {top_count_row['modification']} sites). When normalized by the number of tested sites per modification class, {top_percent_row['modification']} showed the highest relative proportion of differential modification ({fmt_pct(top_percent_row['percent_differential_sites'])}% of tested sites; {fmt_int(top_percent_row['n_differential_sites'])} sites). Pseudouridine accounted for {fmt_int(psu_sig)} differential sites, representing {fmt_pct(psu_percent)}% of tested psU sites.
"""

text_out = os.path.join(
    SUMMARY_DIR,
    "TandemMod_results_section_numbers_paragraph.txt"
)

with open(text_out, "w") as f:
    f.write(auto_text)

print("\nSaved auto-filled text paragraph:")
print(text_out)
print(auto_text)


# ============================================================
# Plot 1. Improved bar plot: counts and percentages
# ============================================================

sns.set_theme(style="whitegrid", font_scale=1.4)

fig, axes = plt.subplots(
    nrows=2,
    ncols=1,
    figsize=(9, 8),
    sharex=True,
    gridspec_kw={"height_ratios": [2.2, 1.5]}
)

# Panel 1: raw counts
ax = axes[0]

sns.barplot(
    data=summary,
    x="modification",
    y="n_differential_sites",
    hue="direction_label",
    order=MODS,
    hue_order=[
        "Higher in aposymbiotic",
        "Higher in symbiotic",
    ],
    palette={
        "Higher in aposymbiotic": "#D8B88A",
        "Higher in symbiotic": "#F28E2B",
    },
    ax=ax
)

ax.set_ylabel("Number of differential sites")
ax.set_xlabel("")
ax.set_title("Differential modification sites by modification type")

# Add labels on bars
for container in ax.containers:
    ax.bar_label(container, fmt="%.0f", fontsize=8, padding=2)

ax.legend(
    title="Direction",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    frameon=False
)

# Panel 2: percentages
ax2 = axes[1]

sns.barplot(
    data=summary,
    x="modification",
    y="percent_of_tested_sites",
    hue="direction_label",
    order=MODS,
    hue_order=[
        "Higher in aposymbiotic",
        "Higher in symbiotic",
    ],
    palette={
        "Higher in aposymbiotic": "#D8B88A",
        "Higher in symbiotic": "#F28E2B",
    },
    ax=ax2
)

ax2.set_ylabel("% of tested sites")
ax2.set_xlabel("Modification type")
ax2.tick_params(axis="x", rotation=45)

for label in ax2.get_xticklabels():
    label.set_horizontalalignment("right")

# Remove duplicated legend in second panel
leg = ax2.get_legend()
if leg is not None:
    leg.remove()

for container in ax2.containers:
    ax2.bar_label(container, fmt="%.1f", fontsize=8, padding=2)

plt.tight_layout()

bar_out_png = os.path.join(
    PLOT_DIR,
    "Improved_barplot_counts_and_percentages_differential_modification_sites.png"
)

bar_out_pdf = os.path.join(
    PLOT_DIR,
    "Improved_barplot_counts_and_percentages_differential_modification_sites.pdf"
)

plt.savefig(bar_out_png, dpi=300, bbox_inches="tight")
plt.savefig(bar_out_pdf, bbox_inches="tight")
plt.close()


# ============================================================
# Plot 2. Improved volcano plot with capped y-axis and counts
# ============================================================

volcano_df = results.copy()

volcano_df["direction_label"] = volcano_df["direction"].map(DIRECTION_LABELS)
volcano_df["direction_label"] = volcano_df["direction_label"].fillna("Not significant")

# cap -log10 FDR for visualization
if "neglog10_fdr" in volcano_df.columns:
    volcano_df["neglog10_fdr_raw"] = volcano_df["neglog10_fdr"]
else:
    volcano_df["neglog10_fdr_raw"] = -np.log10(volcano_df["fdr_by_modification"].replace(0, np.nan))

Y_CAP = 100

volcano_df["neglog10_fdr_capped"] = volcano_df["neglog10_fdr_raw"].clip(upper=Y_CAP)

# Add counts to titles
title_counts = (
    summary_total_mod
    .set_index("modification")
    ["n_differential_sites"]
    .to_dict()
)

sns.set_theme(style="whitegrid", font_scale=1.1)

g = sns.FacetGrid(
    volcano_df,
    col="modification",
    col_order=MODS,
    col_wrap=4,
    height=3.0,
    aspect=1.15,
    sharex=True,
    sharey=True
)

def volcano_panel(data, **kwargs):
    ax = plt.gca()

    # Plot not significant first
    for label in [
        "Not significant",
        "Higher in aposymbiotic",
        "Higher in symbiotic",
    ]:
        sub = data[data["direction_label"] == label]

        if sub.shape[0] == 0:
            continue

        ax.scatter(
            sub["delta_mean_mod_rate_sym_minus_apo"],
            sub["neglog10_fdr_capped"],
            s=10 if label == "Not significant" else 14,
            c=DIRECTION_COLORS[label],
            alpha=0.35 if label == "Not significant" else 0.8,
            edgecolors="none",
            label=label,
        )

    ax.axvline(0, color="black", linewidth=0.8)
    ax.axvline(DELTA_RATE_CUTOFF, color="black", linestyle="--", linewidth=0.7)
    ax.axvline(-DELTA_RATE_CUTOFF, color="black", linestyle="--", linewidth=0.7)
    ax.axhline(-np.log10(FDR_CUTOFF), color="black", linestyle="--", linewidth=0.7)

    ax.set_xlim(-1, 1)
    ax.set_ylim(0, Y_CAP)

g.map_dataframe(volcano_panel)

g.set_axis_labels(
    "Δ modification rate\nSymbiotic - Aposymbiotic",
    f"-log10(FDR), capped at {Y_CAP}"
)

# Custom titles with counts
for ax, mod in zip(g.axes.flatten(), MODS):
    n_sig = title_counts.get(mod, 0)
    ax.set_title(f"{mod}\n{int(n_sig):,} differential sites", fontsize=11)

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

plt.tight_layout(rect=[0, 0, 0.88, 1])

volcano_out_png = os.path.join(
    PLOT_DIR,
    "Improved_volcano_differential_modification_all_mods_capped.png"
)

volcano_out_pdf = os.path.join(
    PLOT_DIR,
    "Improved_volcano_differential_modification_all_mods_capped.pdf"
)

plt.savefig(volcano_out_png, dpi=300, bbox_inches="tight")
plt.savefig(volcano_out_pdf, bbox_inches="tight")
plt.close()


# ============================================================
# Plot 3. Improved heatmap with short labels
# ============================================================

TOP_N_HEATMAP = 25

sig_sorted = sig.copy()

sig_sorted = sig_sorted.sort_values(
    ["fdr_by_modification", "delta_mean_mod_rate_sym_minus_apo"],
    ascending=[True, False]
)

top_sites = sig_sorted.head(TOP_N_HEATMAP)["site_id"].tolist()

heat_df = filtered[filtered["site_id"].isin(top_sites)].copy()

heat_meta = results[
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

heat_df = heat_df.merge(
    heat_meta,
    on=["site_id", "modification", "transcriptome_id", "site", "motif"],
    how="left"
)

# Short labels for plot
rank_map = {site: i + 1 for i, site in enumerate(top_sites)}

heat_df["site_rank"] = heat_df["site_id"].map(rank_map)

heat_df["short_label"] = (
    "S" + heat_df["site_rank"].astype(str) +
    " | " +
    heat_df["modification"].astype(str) +
    " | " +
    heat_df["transcriptome_id"].astype(str).str.replace(r"\.scaffold.*", "", regex=True) +
    ":" +
    heat_df["site"].astype(str)
)

# Save full mapping table
label_map = (
    heat_df[
        [
            "site_rank",
            "short_label",
            "site_id",
            "modification",
            "transcriptome_id",
            "site",
            "motif",
            "direction",
            "fdr_by_modification",
            "delta_mean_mod_rate_sym_minus_apo",
        ]
    ]
    .drop_duplicates()
    .sort_values("site_rank")
)

label_map_out = os.path.join(
    SUMMARY_DIR,
    "Top_heatmap_site_label_mapping.tsv"
)

label_map.to_csv(label_map_out, sep="\t", index=False)

short_label_order = (
    label_map
    .set_index("site_rank")
    .loc[list(range(1, len(top_sites) + 1)), "short_label"]
    .tolist()
)

heat_mat = heat_df.pivot_table(
    index="short_label",
    columns="sample",
    values=MOD_RATE_COL,
    aggfunc="mean"
)

heat_mat = heat_mat.reindex(short_label_order)
heat_mat = heat_mat[SAMPLES]

plt.figure(figsize=(8, max(6, TOP_N_HEATMAP * 0.28)))

ax = sns.heatmap(
    heat_mat,
    cmap="YlOrBr",
    vmin=0,
    vmax=1,
    linewidths=0.3,
    linecolor="white",
    mask=heat_mat.isna(),
    cbar_kws={"label": f"Modification rate, p ≥ {P_CUTOFF_TEXT}"}
)

plt.xlabel("Sample")
plt.ylabel(f"Top {TOP_N_HEATMAP} differential modification sites")
plt.xticks(rotation=45, ha="right")
plt.yticks(rotation=0)

plt.tight_layout()

heat_out_png = os.path.join(
    PLOT_DIR,
    f"Improved_heatmap_top{TOP_N_HEATMAP}_differential_modification_sites_short_labels.png"
)

heat_out_pdf = os.path.join(
    PLOT_DIR,
    f"Improved_heatmap_top{TOP_N_HEATMAP}_differential_modification_sites_short_labels.pdf"
)

plt.savefig(heat_out_png, dpi=300, bbox_inches="tight")
plt.savefig(heat_out_pdf, bbox_inches="tight")
plt.close()


# ============================================================
# Plot 4. Improved boxplots with short labels
# ============================================================

TOP_N_BOXPLOT = 8

example_sites = sig_sorted.head(TOP_N_BOXPLOT)["site_id"].tolist()

example_df = filtered[filtered["site_id"].isin(example_sites)].copy()

example_meta = results[
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
    example_meta,
    on=["site_id", "modification", "transcriptome_id", "site", "motif"],
    how="left"
)

example_rank_map = {
    site: i + 1
    for i, site in enumerate(example_sites)
}

example_df["site_rank"] = example_df["site_id"].map(example_rank_map)

example_df["short_label"] = (
    "Site " + example_df["site_rank"].astype(str) +
    "\n" +
    example_df["modification"].astype(str) +
    " | " +
    example_df["transcriptome_id"].astype(str).str.replace(r"\.scaffold.*", "", regex=True) +
    ":" +
    example_df["site"].astype(str)
)

example_label_map = (
    example_df[
        [
            "site_rank",
            "short_label",
            "site_id",
            "modification",
            "transcriptome_id",
            "site",
            "motif",
            "direction",
            "fdr_by_modification",
            "delta_mean_mod_rate_sym_minus_apo",
        ]
    ]
    .drop_duplicates()
    .sort_values("site_rank")
)

example_label_map_out = os.path.join(
    SUMMARY_DIR,
    "Top_boxplot_site_label_mapping.tsv"
)

example_label_map.to_csv(example_label_map_out, sep="\t", index=False)

label_order = (
    example_label_map
    .sort_values("site_rank")
    ["short_label"]
    .tolist()
)

example_df["short_label"] = pd.Categorical(
    example_df["short_label"],
    categories=label_order,
    ordered=True
)

sns.set_theme(style="whitegrid", font_scale=1.2)

g = sns.catplot(
    data=example_df,
    x="condition",
    y=MOD_RATE_COL,
    col="short_label",
    col_wrap=4,
    kind="box",
    hue="condition",
    order=["Aposymbiotic", "Symbiotic"],
    hue_order=["Aposymbiotic", "Symbiotic"],
    palette=CONDITION_PALETTE,
    showfliers=False,
    height=3.1,
    aspect=0.9,
    legend=False
)

for ax, label in zip(g.axes.flatten(), label_order):
    sub = example_df[example_df["short_label"] == label]

    sns.stripplot(
        data=sub,
        x="condition",
        y=MOD_RATE_COL,
        order=["Aposymbiotic", "Symbiotic"],
        hue="condition",
        hue_order=["Aposymbiotic", "Symbiotic"],
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
    ax.set_ylim(-0.05, 1.05)

g.set_titles("{col_name}")

plt.tight_layout()

box_out_png = os.path.join(
    PLOT_DIR,
    f"Improved_boxplot_top{TOP_N_BOXPLOT}_differential_modification_sites_short_labels.png"
)

box_out_pdf = os.path.join(
    PLOT_DIR,
    f"Improved_boxplot_top{TOP_N_BOXPLOT}_differential_modification_sites_short_labels.pdf"
)

plt.savefig(box_out_png, dpi=300, bbox_inches="tight")
plt.savefig(box_out_pdf, bbox_inches="tight")
plt.close()


print("\nImproved plots saved to:")
print(PLOT_DIR)

print("\nLabel mapping tables saved to:")
print(label_map_out)
print(example_label_map_out)
