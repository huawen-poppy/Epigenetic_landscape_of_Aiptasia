#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import logomaker


INFILE = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod/TandemMod_site_level_differential_analysis_p0p9/TandemMod_differential_modification_results_p0.9_minReads10_minSamples2.tsv"

OUTDIR = "./motif_logo_plots"
os.makedirs(OUTDIR, exist_ok=True)

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

DIFF_DIRECTIONS = ["higher_in_aposymbiotic", "higher_in_symbiotic"]

# If True: create separate logos by direction as well
MAKE_DIRECTION_SPECIFIC = True


# ============================================================
# Helper functions
# ============================================================

def motifs_to_frequency_matrix(motif_list, alphabet=("A", "C", "G", "T")):
    motif_list = [m for m in motif_list if isinstance(m, str)]
    motif_list = [m.strip().upper() for m in motif_list if len(m.strip()) == 5]

    if len(motif_list) == 0:
        return None

    L = 5
    mat = pd.DataFrame(0.0, index=range(L), columns=list(alphabet))

    for motif in motif_list:
        for i, base in enumerate(motif):
            if base in mat.columns:
                mat.loc[i, base] += 1

    mat = mat.div(mat.sum(axis=1), axis=0)
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


# ============================================================
# Load differential modification table
# ============================================================

df = pd.read_csv(INFILE, sep="\t")

required_cols = ["modification", "direction", "motif"]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

# Keep only differentially modified sites
df = df[df["direction"].isin(DIFF_DIRECTIONS)].copy()

# Clean motif
df["motif"] = df["motif"].astype(str).str.upper().str.strip()
df = df[df["motif"].str.len() == 5].copy()

print("Differentially modified rows retained:", df.shape[0])


# ============================================================
# 1. One pooled logo per modification
# ============================================================

fig, axes = plt.subplots(2, 4, figsize=(16, 8))
axes = axes.flatten()

for ax, mod in zip(axes, MODS):
    sub = df[df["modification"] == mod].copy()
    motifs = sub["motif"].tolist()

    freq_df = motifs_to_frequency_matrix(motifs)

    if freq_df is None:
        ax.axis("off")
        ax.set_title(f"{MOD_LABEL_MAP.get(mod, mod)}\n(no sites)")
        continue

    freq_out = os.path.join(OUTDIR, f"{mod}_pooled_frequency_matrix.tsv")
    freq_df.to_csv(freq_out, sep="\t", index=True)

    title = f"{MOD_LABEL_MAP.get(mod, mod)}\nall differential sites (n={len(motifs)})"
    draw_logo(ax, freq_df, title)

plt.tight_layout()
pooled_png = os.path.join(OUTDIR, "Modification_5mer_logo_pooled.png")
pooled_pdf = os.path.join(OUTDIR, "Modification_5mer_logo_pooled.pdf")
plt.savefig(pooled_png, dpi=300, bbox_inches="tight")
plt.savefig(pooled_pdf, bbox_inches="tight")
plt.close()

print("Saved pooled logo figure:")
print(pooled_png)


# ============================================================
# 2. Direction-specific logos
# ============================================================

if MAKE_DIRECTION_SPECIFIC:
    direction_label_map = {
        "higher_in_aposymbiotic": "Higher in aposymbiotic",
        "higher_in_symbiotic": "Higher in symbiotic",
    }

    for direction in DIFF_DIRECTIONS:
        fig, axes = plt.subplots(2, 4, figsize=(16, 8))
        axes = axes.flatten()

        for ax, mod in zip(axes, MODS):
            sub = df[
                (df["modification"] == mod) &
                (df["direction"] == direction)
            ].copy()

            motifs = sub["motif"].tolist()
            freq_df = motifs_to_frequency_matrix(motifs)

            if freq_df is None:
                ax.axis("off")
                ax.set_title(f"{MOD_LABEL_MAP.get(mod, mod)}\n(no sites)")
                continue

            freq_out = os.path.join(
                OUTDIR,
                f"{mod}_{direction}_frequency_matrix.tsv"
            )
            freq_df.to_csv(freq_out, sep="\t", index=True)

            title = (
                f"{MOD_LABEL_MAP.get(mod, mod)}\n"
                f"{direction_label_map[direction]} (n={len(motifs)})"
            )
            draw_logo(ax, freq_df, title)

        plt.tight_layout()
        out_png = os.path.join(OUTDIR, f"Modification_5mer_logo_{direction}.png")
        out_pdf = os.path.join(OUTDIR, f"Modification_5mer_logo_{direction}.pdf")
        plt.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.savefig(out_pdf, bbox_inches="tight")
        plt.close()

        print("Saved:", out_png)

print("Done.")
