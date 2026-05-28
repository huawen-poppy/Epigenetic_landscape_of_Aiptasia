#!/usr/bin/env python

import os
import glob
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


# =========================
# User settings
# =========================

TANDEMMOD_DIR = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod"

ANALYSIS_DIR = f"{TANDEMMOD_DIR}/TandemMod_site_level_differential_analysis_p0p9"

OUT_DIR = f"{ANALYSIS_DIR}/motif_analysis"
os.makedirs(OUT_DIR, exist_ok=True)

PLOT_DIR = f"{OUT_DIR}/plots"
os.makedirs(PLOT_DIR, exist_ok=True)

MODS = ["m6A", "m1A", "A_I", "m5C", "hm5C", "m7G", "G_I", "psU"]

DIRECTIONS = [
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


# =========================
# Load differential result table
# =========================

result_files = glob.glob(
    os.path.join(
        ANALYSIS_DIR,
        "TandemMod_differential_modification_results_p0.9*.tsv"
    )
)

if len(result_files) == 0:
    raise FileNotFoundError("Cannot find differential modification result file.")

result_file = result_files[0]

print("Using result file:")
print(result_file)

df = pd.read_csv(result_file, sep="\t")

required_cols = [
    "site_id",
    "modification",
    "motif",
    "direction",
    "fdr_by_modification",
    "delta_mean_mod_rate_sym_minus_apo",
]

missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns: {missing}")

# one row per tested site
df = df.drop_duplicates("site_id").copy()

df["modification"] = pd.Categorical(
    df["modification"],
    categories=MODS,
    ordered=True
)

df["direction_label"] = df["direction"].map(DIRECTION_LABELS)
df["direction_label"] = df["direction_label"].fillna("Not significant")


# =========================
# 1. Motif composition table
# =========================

composition = (
    df
    .groupby(["modification", "direction", "motif"], observed=False)
    .size()
    .reset_index(name="n_sites")
)

composition["total_in_group"] = (
    composition
    .groupby(["modification", "direction"], observed=False)["n_sites"]
    .transform("sum")
)

composition["fraction"] = composition["n_sites"] / composition["total_in_group"]

composition_out = os.path.join(
    OUT_DIR,
    "Motif_composition_by_modification_and_direction.tsv"
)

composition.to_csv(composition_out, sep="\t", index=False)


# =========================
# 2. Motif enrichment
# =========================

enrichment_records = []

for mod in MODS:
    mod_df = df[df["modification"] == mod].copy()

    if mod_df.shape[0] == 0:
        continue

    motifs = sorted(mod_df["motif"].dropna().unique())

    for direction in DIRECTIONS:
        target = mod_df[mod_df["direction"] == direction]
        background = mod_df[mod_df["direction"] != direction]

        if target.shape[0] == 0:
            continue

        for motif in motifs:
            a = ((target["motif"] == motif)).sum()
            b = ((target["motif"] != motif)).sum()
            c = ((background["motif"] == motif)).sum()
            d = ((background["motif"] != motif)).sum()

            table = [[a, b], [c, d]]

            try:
                odds_ratio, p_value = fisher_exact(table, alternative="greater")
            except Exception:
                odds_ratio, p_value = np.nan, np.nan

            target_fraction = a / (a + b) if (a + b) > 0 else np.nan
            background_fraction = c / (c + d) if (c + d) > 0 else np.nan

            enrichment_records.append(
                {
                    "modification": mod,
                    "direction": direction,
                    "direction_label": DIRECTION_LABELS[direction],
                    "motif": motif,
                    "n_target_with_motif": a,
                    "n_target_total": a + b,
                    "n_background_with_motif": c,
                    "n_background_total": c + d,
                    "target_fraction": target_fraction,
                    "background_fraction": background_fraction,
                    "odds_ratio": odds_ratio,
                    "p_value": p_value,
                }
            )

enrich = pd.DataFrame(enrichment_records)

# BH correction within each modification + direction
enrich["fdr"] = np.nan

for (mod, direction), sub_idx in enrich.groupby(
    ["modification", "direction"]
).groups.items():
    idx = list(sub_idx)
    pvals = enrich.loc[idx, "p_value"]
    valid = pvals.notna()

    if valid.sum() > 0:
        enrich.loc[np.array(idx)[valid], "fdr"] = multipletests(
            pvals[valid],
            method="fdr_bh"
        )[1]

enrich["minus_log10_fdr"] = -np.log10(enrich["fdr"].replace(0, np.nan))
enrich["log2_odds_ratio"] = np.log2(enrich["odds_ratio"].replace(0, np.nan))

enrichment_out = os.path.join(
    OUT_DIR,
    "Motif_enrichment_by_modification_and_direction.tsv"
)

enrich.to_csv(enrichment_out, sep="\t", index=False)


# =========================
# 3. Plot top enriched motifs
# =========================

plot_df = enrich.copy()

plot_df = plot_df[
    (plot_df["fdr"] < 0.05) &
    (plot_df["odds_ratio"] > 1)
].copy()

# keep top 8 motifs per modification + direction
plot_df = (
    plot_df
    .sort_values(["modification", "direction", "fdr", "odds_ratio"],
                 ascending=[True, True, True, False])
    .groupby(["modification", "direction"], observed=False)
    .head(8)
    .copy()
)

if plot_df.shape[0] > 0:

    plot_df["motif_label"] = (
        plot_df["modification"].astype(str) +
        " | " +
        plot_df["motif"].astype(str)
    )

    sns.set_theme(style="whitegrid", font_scale=1.1)

    g = sns.FacetGrid(
        plot_df,
        col="direction_label",
        col_order=[
            "Higher in aposymbiotic",
            "Higher in symbiotic",
        ],
        height=6,
        aspect=0.9,
        sharey=False,
        sharex=True
    )

    def dotplot_panel(data, **kwargs):
        ax = plt.gca()

        data = data.sort_values("log2_odds_ratio", ascending=True)

        ax.scatter(
            data["log2_odds_ratio"],
            data["motif_label"],
            s=data["minus_log10_fdr"] * 15,
            c=data["target_fraction"],
            cmap="YlOrBr",
            edgecolor="black",
            linewidth=0.3
        )

        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_xlabel("log2 odds ratio")
        ax.set_ylabel("Motif")
        ax.grid(True, axis="x", alpha=0.3)

    g.map_dataframe(dotplot_panel)

    g.set_titles("{col_name}")

    plt.tight_layout()

    motif_plot_png = os.path.join(
        PLOT_DIR,
        "Top_enriched_motifs_dotplot.png"
    )

    motif_plot_pdf = os.path.join(
        PLOT_DIR,
        "Top_enriched_motifs_dotplot.pdf"
    )

    plt.savefig(motif_plot_png, dpi=300, bbox_inches="tight")
    plt.savefig(motif_plot_pdf, bbox_inches="tight")
    plt.close()

else:
    print("No significantly enriched motifs found for plotting.")


# =========================
# 4. Motif composition stacked bar plot
# =========================

# For readability, keep top motifs per modification by total frequency
top_motifs = (
    df
    .groupby(["modification", "motif"], observed=False)
    .size()
    .reset_index(name="n")
    .sort_values(["modification", "n"], ascending=[True, False])
    .groupby("modification", observed=False)
    .head(8)
)

top_motif_set = set(
    zip(top_motifs["modification"].astype(str), top_motifs["motif"].astype(str))
)

comp_plot = composition[
    composition["direction"].isin(DIRECTIONS)
].copy()

comp_plot["motif_plot"] = comp_plot.apply(
    lambda r: r["motif"]
    if (str(r["modification"]), str(r["motif"])) in top_motif_set
    else "Other",
    axis=1
)

comp_plot2 = (
    comp_plot
    .groupby(["modification", "direction", "motif_plot"], observed=False)
    ["n_sites"]
    .sum()
    .reset_index()
)

comp_plot2["total"] = (
    comp_plot2
    .groupby(["modification", "direction"], observed=False)["n_sites"]
    .transform("sum")
)

comp_plot2["fraction"] = comp_plot2["n_sites"] / comp_plot2["total"]

comp_plot2["group"] = (
    comp_plot2["modification"].astype(str) +
    "\n" +
    comp_plot2["direction"].map(DIRECTION_LABELS)
)

composition_plot_out = os.path.join(
    OUT_DIR,
    "Motif_composition_for_plot.tsv"
)

comp_plot2.to_csv(composition_plot_out, sep="\t", index=False)

print("Saved motif composition and enrichment:")
print(composition_out)
print(enrichment_out)
print(composition_plot_out)
