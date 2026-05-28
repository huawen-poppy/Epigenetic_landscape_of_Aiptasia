import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# =========================
# Input / output paths
# =========================

pred_dir = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod/TandemMod_prediction"

outdir = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod/TandemMod_prediction_summary/plots_updated"
os.makedirs(outdir, exist_ok=True)


# =========================
# Sample / modification setup
# =========================

sample_order = ["A1", "A2", "A3", "H1", "H2", "H3"]

condition_map = {
    "A1": "Aposymbiotic",
    "A2": "Aposymbiotic",
    "A3": "Aposymbiotic",
    "H1": "Symbiotic",
    "H2": "Symbiotic",
    "H3": "Symbiotic",
}

mod_order = ["m6A", "m1A", "A_I", "m5C", "hm5C", "m7G", "G_I", "psU"]

feature_map = {
    "m6A": "DRACH",
    "m1A": "NNANN",
    "A_I": "NNANN",
    "m5C": "NNCNN",
    "hm5C": "NNCNN",
    "m7G": "NNGNN",
    "G_I": "NNGNN",
    "psU": "NNTNN",
}

mod_palette = {
    "m6A":  "#4E79A7",
    "m1A":  "#F28E2B",
    "A_I":  "#E15759",
    "m5C":  "#76B7B2",
    "hm5C": "#59A14F",
    "m7G":  "#EDC948",
    "G_I":  "#B07AA1",
    "psU":  "#9C755F",
}


# =========================
# Memory-safe probability reader
# =========================

def file_has_header(path):
    """
    Your TandemMod prediction files usually have no header:
    transcript_id site motif read_id label probability

    This function checks whether the first line has a numeric 6th column.
    """
    with open(path, "r") as f:
        first = f.readline().strip().split("\t")

    if len(first) < 6:
        return False

    try:
        float(first[5])
        return False
    except ValueError:
        return True


def sample_probabilities(
    path,
    prob_col=5,
    n_max=200000,
    chunksize=1000000,
    seed=123,
):
    """
    Read probability column from a large prediction file.
    Uses random-key sampling to avoid loading the whole file into memory.

    Input columns:
    0 transcript_id
    1 site
    2 motif
    3 read_id
    4 label
    5 probability
    """

    rng = np.random.default_rng(seed)

    sampled_probs = np.array([], dtype=float)
    sampled_keys = np.array([], dtype=float)

    header = 0 if file_has_header(path) else None

    reader = pd.read_csv(
        path,
        sep="\t",
        header=header,
        usecols=[prob_col],
        chunksize=chunksize,
    )

    for chunk in reader:
        probs = pd.to_numeric(chunk.iloc[:, 0], errors="coerce")
        probs = probs.dropna().to_numpy(dtype=float)

        # keep valid probability values only
        probs = probs[(probs >= 0) & (probs <= 1)]

        if len(probs) == 0:
            continue

        keys = rng.random(len(probs))

        sampled_probs = np.concatenate([sampled_probs, probs])
        sampled_keys = np.concatenate([sampled_keys, keys])

        # keep only the largest random keys, equivalent to random sampling
        if len(sampled_probs) > n_max * 2:
            keep_idx = np.argpartition(sampled_keys, -n_max)[-n_max:]
            sampled_probs = sampled_probs[keep_idx]
            sampled_keys = sampled_keys[keep_idx]

    if len(sampled_probs) > n_max:
        keep_idx = np.argpartition(sampled_keys, -n_max)[-n_max:]
        sampled_probs = sampled_probs[keep_idx]

    return sampled_probs


# =========================
# Plot: 8 sample panels
# =========================

sns.set_theme(style="whitegrid")

fig, axes = plt.subplots(
    nrows=2,
    ncols=4,
    figsize=(16, 7),
    sharex=True,
    sharey=False
)

axes = axes.flatten()

for ax, sample in zip(axes, sample_order):

    condition = condition_map[sample]

    for mod in mod_order:
        feature_type = feature_map[mod]

        pred_file = os.path.join(
            pred_dir,
            f"{sample}_{mod}_{feature_type}_prediction.tsv"
        )

        if not os.path.exists(pred_file):
            print(f"Missing file, skip: {pred_file}")
            continue

        print(f"Reading: {sample} {mod}")

        probs = sample_probabilities(
            pred_file,
            prob_col=5,
            n_max=200000,
            chunksize=1000000,
            seed=123,
        )

        if len(probs) < 10:
            print(f"Too few probability values, skip: {pred_file}")
            continue

        # KDE density curve
        sns.kdeplot(
            x=probs,
            ax=ax,
            color=mod_palette[mod],
            label=mod,
            linewidth=1.5,
            fill=False,
            clip=(0, 1),
            cut=0,
            bw_adjust=0.8,
        )

    ax.axvline(
        0.9,
        color="black",
        linestyle="--",
        linewidth=0.8
    )

    ax.set_title(f"{sample} ({condition})", fontsize=12)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Prediction probability")
    ax.set_ylabel("Density")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# =========================
# Shared legend outside
# =========================

handles, labels = axes[0].get_legend_handles_labels()

# remove legends from individual panels
for ax in axes:
    leg = ax.get_legend()
    if leg is not None:
        leg.remove()

fig.legend(
    handles,
    labels,
    title="Modification",
    bbox_to_anchor=(1.01, 0.5),
    loc="center left",
    frameon=False
)

plt.tight_layout(rect=[0, 0, 0.90, 1])

plt.savefig(
    os.path.join(outdir, "QC_prediction_probability_density_6samples_by_modification.png"),
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    os.path.join(outdir, "QC_prediction_probability_density_6samples_by_modification.pdf"),
    bbox_inches="tight"
)

plt.close()
