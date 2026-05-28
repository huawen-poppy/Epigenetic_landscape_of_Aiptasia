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
# Settings
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

N_PER_FILE = 200000


# =========================
# Read probability values
# =========================

def read_probability_and_label(path, n_max=200000, random_state=123):
    """
    TandemMod prediction format:
    0 transcript_id
    1 site
    2 motif
    3 read_id
    4 label: mod / unmod
    5 probability
    """

    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        usecols=[4, 5],
        names=["label", "probability"]
    )

    df["probability"] = pd.to_numeric(df["probability"], errors="coerce")
    df = df.dropna(subset=["probability"])

    df = df[
        (df["probability"] >= 0) &
        (df["probability"] <= 1)
    ]

    if len(df) > n_max:
        df = df.sample(n=n_max, random_state=random_state)

    return df


# =========================
# Build dataframe
# =========================

records = []

for sample in sample_order:
    for mod in mod_order:
        feature_type = feature_map[mod]

        pred_file = os.path.join(
            pred_dir,
            f"{sample}_{mod}_{feature_type}_prediction.tsv"
        )

        if not os.path.exists(pred_file):
            print(f"Missing file, skip: {pred_file}")
            continue

        print(f"Reading {sample} {mod}")

        tmp = read_probability_and_label(
            pred_file,
            n_max=N_PER_FILE,
            random_state=123
        )

        tmp["sample"] = sample
        tmp["condition"] = condition_map[sample]
        tmp["modification"] = mod

        records.append(tmp)

density_df = pd.concat(records, ignore_index=True)


# =========================
# Plot: full range with log y-axis
# =========================

sns.set_theme(style="whitegrid", font_scale=1.2)

fig, axes = plt.subplots(
    nrows=2,
    ncols=4,
    figsize=(18, 7),
    sharex=True,
    sharey=True
)

axes = axes.flatten()

bins = np.linspace(0, 1, 101)

for ax, sample in zip(axes, sample_order):
    sub_sample = density_df[density_df["sample"] == sample]

    for mod in mod_order:
        sub = sub_sample[sub_sample["modification"] == mod]

        if len(sub) == 0:
            continue

        sns.histplot(
            data=sub,
            x="probability",
            bins=bins,
            stat="density",
            element="step",
            fill=False,
            common_norm=False,
            color=mod_palette[mod],
            linewidth=1.4,
            ax=ax,
            label=mod
        )

    ax.axvline(
        0.9,
        color="black",
        linestyle="--",
        linewidth=0.9
    )

    ax.set_yscale("log")

    ax.set_title(
        f"{sample} ({condition_map[sample]})",
        fontsize=13
    )

    ax.set_xlim(0, 1)
    ax.set_xlabel("Prediction probability")
    ax.set_ylabel("Density, log scale")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# shared legend outside
handles, labels = axes[0].get_legend_handles_labels()

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
    os.path.join(outdir, "QC_prediction_probability_density_6samples_logY.png"),
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    os.path.join(outdir, "QC_prediction_probability_density_6samples_logY.pdf"),
    bbox_inches="tight"
)

plt.close()
