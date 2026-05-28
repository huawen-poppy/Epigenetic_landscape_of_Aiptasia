#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker


# ============================================================
# User settings
# ============================================================

ENRICHMENT_FILE = "../TandemMod_site_level_differential_analysis_p0p9/motif_analysis/Motif_enrichment_by_modification_and_direction.tsv"
COMPOSITION_FILE = "../TandemMod_site_level_differential_analysis_p0p9/motif_analysis/Motif_composition_by_modification_and_direction.tsv"

OUTDIR = "./enriched_motif_logo_plots"
os.makedirs(OUTDIR, exist_ok=True)

FDR_CUTOFF = 0.05
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

DIRECTION_MAP = {
    "higher_in_aposymbiotic": "higher_in_aposymbiotic",
    "Higher in aposymbiotic": "higher_in_aposymbiotic",
    "aposymbiotic": "higher_in_aposymbiotic",

    "higher_in_symbiotic": "higher_in_symbiotic",
    "Higher in symbiotic": "higher_in_symbiotic",
    "symbiotic": "higher_in_symbiotic",
}

DIRECTION_LABEL_MAP = {
    "higher_in_aposymbiotic": "Higher in aposymbiotic",
    "higher_in_symbiotic": "Higher in symbiotic",
}


# ============================================================
# Helper functions
# ============================================================

def normalize_direction(x):
    if pd.isna(x):
        return np.nan
    x = str(x).strip()
    return DIRECTION_MAP.get(x, x)


def find_count_column(df):
    """
    Try to detect the count column in the motif composition table.
    Expected names may vary.
    """
    candidates = [
        "n_sites",
        "count",
        "n",
        "target_count",
        "n_target_sites",
        "motif_count",
        "site_count",
    ]
    for c in candidates:
        if c in df.columns:
            return c

    # fallback: choose first numeric column that's not obviously a fraction
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    numeric_cols = [c for c in numeric_cols if "fraction" not in c.lower()]
    if len(numeric_cols) > 0:
        return numeric_cols[0]

    raise ValueError(
        "Could not identify count column in composition table. "
        f"Available columns: {df.columns.tolist()}"
    )


def motifs_to_weighted_frequency_matrix(motif_weights, alphabet=("A", "C", "G", "T")):
    """
    motif_weights: list of tuples [(motif, weight), ...]
    Creates a weighted base-frequency matrix across 5-mer positions.
    """
    motif_weights = [
        (str(m).strip().upper(), float(w))
        for m, w in motif_weights
        if pd.notna(m) and pd.notna(w)
    ]

    motif_weights = [
        (m, w) for m, w in motif_weights
        if len(m) == 5 and w > 0
    ]

    if len(motif_weights) == 0:
        return None

    L = 5
    mat = pd.DataFrame(0.0, index=range(L), columns=list(alphabet))

    for motif, weight in motif_weights:
        for i, base in enumerate(motif):
            if base in mat.columns:
                mat.loc[i, base] += weight

    row_sums = mat.sum(axis=1)
    mat = mat.div(row_sums, axis=0)

    return mat


def draw_logo(ax, freq_df, title):
    logomaker.Logo(freq_df, ax=ax)
    ax.set_title(title, fontsize=14)
    ax.set_xlabel("Motif position")
    ax.set_ylabel("Base frequency")
    ax.set_ylim(0, 1.05)
    ax.set_xticks(range(len(freq_df.index)))
    ax.set_xticklabels([str(i + 1) for i in range(len(freq_df.index))])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def save_logo_grid(df, mods, out_png, out_pdf, title_prefix):
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()

    for ax, mod in zip(axes, mods):
        sub = df[df["modification"] == mod].copy()

        if sub.shape[0] == 0:
            ax.axis("off")
            ax.set_title(f"{MOD_LABEL_MAP.get(mod, mod)}\n(no enriched motifs)")
            continue

        motif_weights = list(zip(sub["motif"], sub["weight"]))
        freq_df = motifs_to_weighted_frequency_matrix(motif_weights)

        if freq_df is None:
            ax.axis("off")
            ax.set_title(f"{MOD_LABEL_MAP.get(mod, mod)}\n(no enriched motifs)")
            continue

        n_sites = int(sub["weight"].sum())
        n_unique_motifs = sub["motif"].nunique()

        freq_out = os.path.join(
            OUTDIR,
            f"{mod}_{title_prefix}_enriched_frequency_matrix.tsv"
        )
        freq_df.to_csv(freq_out, sep="\t", index=True)

        title = (
            f"{MOD_LABEL_MAP.get(mod, mod)}\n"
            f"{title_prefix.replace('_', ' ')} (n={n_sites}, motifs={n_unique_motifs})"
        )
        draw_logo(ax, freq_df, title)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()


# ============================================================
# Load data
# ============================================================

enrich = pd.read_csv(ENRICHMENT_FILE, sep="\t")
comp = pd.read_csv(COMPOSITION_FILE, sep="\t")

print("Enrichment table shape:", enrich.shape)
print("Composition table shape:", comp.shape)

# Normalize key columns
for df in [enrich, comp]:
    if "direction" in df.columns:
        df["direction_norm"] = df["direction"].apply(normalize_direction)
    elif "direction_label" in df.columns:
        df["direction_norm"] = df["direction_label"].apply(normalize_direction)
    else:
        raise ValueError(
            f"Could not find direction column in table with columns: {df.columns.tolist()}"
        )

    if "motif" not in df.columns:
        raise ValueError(f"Missing 'motif' column in table: {df.columns.tolist()}")

    if "modification" not in df.columns:
        raise ValueError(f"Missing 'modification' column in table: {df.columns.tolist()}")

    df["motif"] = df["motif"].astype(str).str.upper().str.strip()
    df["modification"] = df["modification"].astype(str).str.strip()

# Identify count column in composition table
count_col = find_count_column(comp)
print("Using count column from composition table:", count_col)

# Filter enriched motifs
required_enrich_cols = ["modification", "motif", "direction_norm", "fdr", "odds_ratio"]
missing_enrich = [c for c in required_enrich_cols if c not in enrich.columns]
if missing_enrich:
    raise ValueError(f"Missing required enrichment columns: {missing_enrich}")

enriched = enrich[
    (enrich["fdr"] < FDR_CUTOFF) &
    (enrich["odds_ratio"] > MIN_ODDS_RATIO)
].copy()

print("Number of significant enriched motif rows:", enriched.shape[0])

if enriched.shape[0] == 0:
    raise ValueError("No enriched motifs found under the current thresholds.")

# Keep only the keys needed for merging
enriched_keys = enriched[["modification", "direction_norm", "motif"]].drop_duplicates()

# Merge with composition table to get weights
merge_cols = ["modification", "direction_norm", "motif"]

comp_sub = comp.copy()
comp_sub["weight"] = pd.to_numeric(comp_sub[count_col], errors="coerce")
comp_sub = comp_sub.dropna(subset=["weight"])

enriched_weighted = comp_sub.merge(
    enriched_keys,
    on=merge_cols,
    how="inner"
)

print("Enriched motifs with weights shape:", enriched_weighted.shape)

if enriched_weighted.shape[0] == 0:
    raise ValueError(
        "After merging enrichment and composition tables, no rows remain. "
        "Please check that modification / direction / motif columns match between the two files."
    )

# Save the merged weighted table
merged_out = os.path.join(OUTDIR, "Enriched_motifs_weighted_table.tsv")
enriched_weighted.to_csv(merged_out, sep="\t", index=False)
print("Saved weighted enriched motif table:", merged_out)


# ============================================================
# 1. Pooled enriched-motif logos per modification
# ============================================================

pooled = (
    enriched_weighted
    .groupby(["modification", "motif"], as_index=False)["weight"]
    .sum()
)

pooled_png = os.path.join(OUTDIR, "Enriched_motif_logo_pooled.png")
pooled_pdf = os.path.join(OUTDIR, "Enriched_motif_logo_pooled.pdf")

save_logo_grid(
    pooled,
    MODS,
    pooled_png,
    pooled_pdf,
    title_prefix="all_enriched_motifs"
)

print("Saved pooled enriched-motif logo figure:")
print(pooled_png)
print(pooled_pdf)


# ============================================================
# 2. Direction-specific enriched-motif logos
# ============================================================

for direction in ["higher_in_aposymbiotic", "higher_in_symbiotic"]:
    sub = enriched_weighted[
        enriched_weighted["direction_norm"] == direction
    ].copy()

    if sub.shape[0] == 0:
        print(f"No enriched motifs for direction {direction}; skipping.")
        continue

    grouped = (
        sub
        .groupby(["modification", "motif"], as_index=False)["weight"]
        .sum()
    )

    out_png = os.path.join(OUTDIR, f"Enriched_motif_logo_{direction}.png")
    out_pdf = os.path.join(OUTDIR, f"Enriched_motif_logo_{direction}.pdf")

    save_logo_grid(
        grouped,
        MODS,
        out_png,
        out_pdf,
        title_prefix=direction
    )

    print("Saved direction-specific enriched-motif logo figure:")
    print(out_png)
    print(out_pdf)


# ============================================================
# 3. Also save per-modification single logos (optional but useful)
# ============================================================

single_dir = os.path.join(OUTDIR, "single_modification_logos")
os.makedirs(single_dir, exist_ok=True)

# pooled single-modification logos
for mod in MODS:
    sub = pooled[pooled["modification"] == mod].copy()

    if sub.shape[0] == 0:
        continue

    motif_weights = list(zip(sub["motif"], sub["weight"]))
    freq_df = motifs_to_weighted_frequency_matrix(motif_weights)

    if freq_df is None:
        continue

    fig, ax = plt.subplots(figsize=(5, 3.8))
    n_sites = int(sub["weight"].sum())
    n_unique_motifs = sub["motif"].nunique()

    draw_logo(
        ax,
        freq_df,
        f"{MOD_LABEL_MAP.get(mod, mod)}\nall enriched motifs (n={n_sites}, motifs={n_unique_motifs})"
    )

    out_png = os.path.join(single_dir, f"{mod}_all_enriched_motifs_logo.png")
    out_pdf = os.path.join(single_dir, f"{mod}_all_enriched_motifs_logo.pdf")

    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()

# direction-specific single-modification logos
for direction in ["higher_in_aposymbiotic", "higher_in_symbiotic"]:
    subdir = os.path.join(single_dir, direction)
    os.makedirs(subdir, exist_ok=True)

    sub_all = enriched_weighted[
        enriched_weighted["direction_norm"] == direction
    ].copy()

    for mod in MODS:
        sub = (
            sub_all[sub_all["modification"] == mod]
            .groupby(["modification", "motif"], as_index=False)["weight"]
            .sum()
        )

        if sub.shape[0] == 0:
            continue

        motif_weights = list(zip(sub["motif"], sub["weight"]))
        freq_df = motifs_to_weighted_frequency_matrix(motif_weights)

        if freq_df is None:
            continue

        fig, ax = plt.subplots(figsize=(5, 3.8))
        n_sites = int(sub["weight"].sum())
        n_unique_motifs = sub["motif"].nunique()

        draw_logo(
            ax,
            freq_df,
            f"{MOD_LABEL_MAP.get(mod, mod)}\n{DIRECTION_LABEL_MAP[direction]} (n={n_sites}, motifs={n_unique_motifs})"
        )

        out_png = os.path.join(subdir, f"{mod}_{direction}_enriched_logo.png")
        out_pdf = os.path.join(subdir, f"{mod}_{direction}_enriched_logo.pdf")

        plt.tight_layout()
        plt.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.savefig(out_pdf, bbox_inches="tight")
        plt.close()

print("Done.")
