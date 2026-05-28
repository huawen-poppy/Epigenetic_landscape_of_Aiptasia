import os
import glob
import numpy as np
import pandas as pd

# =========================
# User settings
# =========================

PRED_DIR = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod/TandemMod_prediction"
OUT_DIR = "/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod/TandemMod_prediction_summary"

os.makedirs(OUT_DIR, exist_ok=True)

SAMPLES = ["A1", "A2", "A3", "H1", "H2", "H3"]

CONDITION_MAP = {
    "A1": "apo",
    "A2": "apo",
    "A3": "apo",
    "H1": "holo",
    "H2": "holo",
    "H3": "holo",
}

MOD_INFO = {
    "m6A":  "DRACH",
    "m1A":  "NNANN",
    "A_I":  "NNANN",
    "m5C":  "NNCNN",
    "hm5C": "NNCNN",
    "m7G":  "NNGNN",
    "G_I":  "NNGNN",
    "psU":  "NNTNN",
}

PROB_CUTOFFS = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]


# =========================
# Helper functions
# =========================

def find_prediction_file(sample, mod, feature_type):
    """
    Expected filename:
    SAMPLE_MOD_FEATURE_prediction.tsv

    Example:
    A1_m6A_DRACH_prediction.tsv
    """
    expected = os.path.join(
        PRED_DIR,
        f"{sample}_{mod}_{feature_type}_prediction.tsv"
    )

    if os.path.exists(expected):
        return expected

    # Fallback: useful if some filenames differ slightly
    pattern = os.path.join(PRED_DIR, f"{sample}_{mod}_*prediction.tsv")
    matches = glob.glob(pattern)

    if len(matches) == 1:
        return matches[0]
    elif len(matches) > 1:
        print(f"WARNING: multiple files found for {sample} {mod}:")
        for m in matches:
            print("  ", m)
        print("Using first one.")
        return matches[0]
    else:
        return None


def read_prediction_file(path):
    """
    Robust reader for TandemMod prediction files.
    Handles files with or without a header.
    """

    with open(path, "r") as f:
        first_line = f.readline().strip()

    first_fields = first_line.split("\t")
    first_lower = [x.lower() for x in first_fields]

    header_keywords = [
        "transcript", "transcript_id", "site", "motif",
        "read", "read_id", "prediction", "label",
        "prob", "probability"
    ]

    has_header = any(
        any(keyword in field for keyword in header_keywords)
        for field in first_lower
    )

    if has_header:
        df = pd.read_csv(path, sep="\t")
    else:
        df = pd.read_csv(path, sep="\t", header=None)

        ncol = df.shape[1]

        if ncol >= 6:
            # Common TandemMod-like format
            base_cols = [
                "transcript_id",
                "site",
                "motif",
                "read_id",
                "prediction",
                "probability"
            ]
            extra_cols = [f"extra_{i}" for i in range(ncol - 6)]
            df.columns = base_cols + extra_cols
        else:
            # Fallback for unexpected format
            df.columns = [f"col_{i}" for i in range(ncol)]

    # Normalize column names
    df.columns = [str(c).strip() for c in df.columns]

    return df


def find_column(df, candidates):
    """
    Find the first matching column from candidate names.
    Case-insensitive.
    """
    col_map = {c.lower(): c for c in df.columns}

    for cand in candidates:
        if cand.lower() in col_map:
            return col_map[cand.lower()]

    # partial match fallback
    for c in df.columns:
        c_lower = c.lower()
        for cand in candidates:
            if cand.lower() in c_lower:
                return c

    return None


def infer_probability_column(df):
    prob_col = find_column(
        df,
        ["probability", "prob", "score", "prediction_probability", "mod_prob"]
    )

    if prob_col is not None:
        return prob_col

    # fallback: choose the last numeric-looking column
    for c in reversed(df.columns):
        values = pd.to_numeric(df[c], errors="coerce")
        if values.notna().mean() > 0.9:
            return c

    return None


def infer_prediction_column(df):
    pred_col = find_column(
        df,
        ["prediction", "predict", "label", "pred_label", "mod_label", "pred"]
    )

    if pred_col is not None:
        return pred_col

    # fallback: often second last column if probability is last
    if df.shape[1] >= 2:
        return df.columns[-2]

    return None


def infer_modified_mask(pred_series, mod):
    """
    Try to identify modified calls.

    Handles labels such as:
    1 / 0
    modified / unmodified
    mod / unmod
    m6A / unmodified
    True / False
    """

    s = pred_series.astype(str).str.strip()

    # numeric labels
    numeric = pd.to_numeric(s, errors="coerce")
    if numeric.notna().mean() > 0.9:
        return numeric == 1

    s_lower = s.str.lower()
    mod_lower = mod.lower()

    modified_terms = {
        "1", "true", "t", "yes", "y",
        "mod", "modified", "m",
        mod_lower,
    }

    # Additional aliases
    if mod == "A_I":
        modified_terms.update(["a_i", "a-to-i", "a2i", "ai"])
    if mod == "G_I":
        modified_terms.update(["g_i", "g-to-i", "g2i", "gi"])
    if mod == "psU":
        modified_terms.update(["psu", "pseudouridine", "pseudo"])

    return s_lower.isin(modified_terms) | s_lower.str.contains(mod_lower, regex=False)


def summarize_one_file(sample, mod, feature_type, path):
    df = read_prediction_file(path)

    pred_col = infer_prediction_column(df)
    prob_col = infer_probability_column(df)

    transcript_col = find_column(df, ["transcript_id", "transcript", "tx", "tx_id"])
    site_col = find_column(df, ["site", "position", "pos", "transcript_pos"])
    motif_col = find_column(df, ["motif", "kmer"])
    read_col = find_column(df, ["read_id", "read", "readname"])

    result = {
        "sample": sample,
        "condition": CONDITION_MAP.get(sample, "unknown"),
        "modification": mod,
        "feature_type": feature_type,
        "file": path,
        "file_size_MB": os.path.getsize(path) / 1024 / 1024,
        "n_rows": len(df),
        "prediction_column": pred_col,
        "probability_column": prob_col,
        "transcript_column": transcript_col,
        "site_column": site_col,
        "motif_column": motif_col,
        "read_column": read_col,
    }

    if transcript_col is not None:
        result["n_unique_transcripts"] = df[transcript_col].nunique()
    else:
        result["n_unique_transcripts"] = np.nan

    if site_col is not None:
        result["n_unique_sites"] = df[site_col].nunique()
    else:
        result["n_unique_sites"] = np.nan

    if transcript_col is not None and site_col is not None:
        result["n_unique_transcript_sites"] = (
            df[transcript_col].astype(str) + ":" + df[site_col].astype(str)
        ).nunique()
    else:
        result["n_unique_transcript_sites"] = np.nan

    if motif_col is not None:
        result["n_unique_motifs"] = df[motif_col].nunique()
    else:
        result["n_unique_motifs"] = np.nan

    if read_col is not None:
        result["n_unique_reads"] = df[read_col].nunique()
    else:
        result["n_unique_reads"] = np.nan

    if prob_col is not None:
        prob = pd.to_numeric(df[prob_col], errors="coerce")
        result["prob_mean"] = prob.mean()
        result["prob_median"] = prob.median()
        result["prob_q25"] = prob.quantile(0.25)
        result["prob_q75"] = prob.quantile(0.75)
        result["prob_q90"] = prob.quantile(0.90)
        result["prob_q95"] = prob.quantile(0.95)
        result["prob_min"] = prob.min()
        result["prob_max"] = prob.max()
    else:
        prob = None
        for x in [
            "prob_mean", "prob_median", "prob_q25", "prob_q75",
            "prob_q90", "prob_q95", "prob_min", "prob_max"
        ]:
            result[x] = np.nan

    if pred_col is not None:
        modified_mask = infer_modified_mask(df[pred_col], mod)
        result["n_modified_all"] = int(modified_mask.sum())
        result["n_unmodified_all"] = int(len(df) - modified_mask.sum())
        result["modified_fraction_all"] = modified_mask.mean()

        # label counts
        label_counts = df[pred_col].astype(str).value_counts()
        label_count_str = "; ".join(
            [f"{label}:{count}" for label, count in label_counts.items()]
        )
        result["prediction_label_counts"] = label_count_str
    else:
        modified_mask = None
        result["n_modified_all"] = np.nan
        result["n_unmodified_all"] = np.nan
        result["modified_fraction_all"] = np.nan
        result["prediction_label_counts"] = "NA"

    # Probability cutoff summaries
    for cutoff in PROB_CUTOFFS:
        tag = str(cutoff).replace(".", "p")

        if prob is not None:
            high_conf_mask = prob >= cutoff
            result[f"n_rows_prob_ge_{tag}"] = int(high_conf_mask.sum())
            result[f"fraction_rows_prob_ge_{tag}"] = high_conf_mask.mean()

            if modified_mask is not None:
                high_conf_modified = high_conf_mask & modified_mask
                result[f"n_modified_prob_ge_{tag}"] = int(high_conf_modified.sum())
                result[f"modified_fraction_total_prob_ge_{tag}"] = high_conf_modified.sum() / len(df)

                if high_conf_mask.sum() > 0:
                    result[f"modified_fraction_within_prob_ge_{tag}"] = (
                        high_conf_modified.sum() / high_conf_mask.sum()
                    )
                else:
                    result[f"modified_fraction_within_prob_ge_{tag}"] = np.nan
            else:
                result[f"n_modified_prob_ge_{tag}"] = np.nan
                result[f"modified_fraction_total_prob_ge_{tag}"] = np.nan
                result[f"modified_fraction_within_prob_ge_{tag}"] = np.nan
        else:
            result[f"n_rows_prob_ge_{tag}"] = np.nan
            result[f"fraction_rows_prob_ge_{tag}"] = np.nan
            result[f"n_modified_prob_ge_{tag}"] = np.nan
            result[f"modified_fraction_total_prob_ge_{tag}"] = np.nan
            result[f"modified_fraction_within_prob_ge_{tag}"] = np.nan

    return result


# =========================
# Main summary
# =========================

all_results = []
missing_files = []

for sample in SAMPLES:
    for mod, feature_type in MOD_INFO.items():
        pred_file = find_prediction_file(sample, mod, feature_type)

        if pred_file is None:
            print(f"Missing file: sample={sample}, mod={mod}, feature={feature_type}")
            missing_files.append({
                "sample": sample,
                "condition": CONDITION_MAP.get(sample, "unknown"),
                "modification": mod,
                "feature_type": feature_type,
                "expected_pattern": f"{sample}_{mod}_{feature_type}_prediction.tsv"
            })
            continue

        print(f"Summarizing: {sample} {mod} {feature_type}")
        one_result = summarize_one_file(sample, mod, feature_type, pred_file)
        all_results.append(one_result)


summary_df = pd.DataFrame(all_results)

summary_out = os.path.join(OUT_DIR, "TandemMod_prediction_sample_modification_summary.tsv")
summary_df.to_csv(summary_out, sep="\t", index=False)

missing_df = pd.DataFrame(missing_files)
missing_out = os.path.join(OUT_DIR, "TandemMod_missing_prediction_files.tsv")
missing_df.to_csv(missing_out, sep="\t", index=False)

print("\nSaved sample/modification summary:")
print(summary_out)

print("\nSaved missing file table:")
print(missing_out)


# =========================
# Condition-level summary
# =========================

numeric_cols = summary_df.select_dtypes(include=[np.number]).columns.tolist()

condition_summary = (
    summary_df
    .groupby(["condition", "modification", "feature_type"], as_index=False)[numeric_cols]
    .agg(["mean", "median", "std"])
)

condition_summary.columns = [
    "_".join([str(x) for x in col if x != ""])
    for col in condition_summary.columns
]

condition_summary_out = os.path.join(
    OUT_DIR,
    "TandemMod_prediction_condition_summary.tsv"
)

condition_summary.to_csv(condition_summary_out, sep="\t", index=False)

print("\nSaved condition-level summary:")
print(condition_summary_out)


# =========================
# Apo vs holo comparison table
# =========================

main_metric = "modified_fraction_total_prob_ge_0p9"

if main_metric in summary_df.columns:
    comparison = (
        summary_df
        .pivot_table(
            index=["modification", "feature_type"],
            columns="condition",
            values=main_metric,
            aggfunc="mean"
        )
        .reset_index()
    )

    if "apo" in comparison.columns and "holo" in comparison.columns:
        comparison["delta_holo_minus_apo"] = comparison["holo"] - comparison["apo"]

    comparison_out = os.path.join(
        OUT_DIR,
        "TandemMod_apo_vs_holo_modified_fraction_prob_ge_0p9.tsv"
    )

    comparison.to_csv(comparison_out, sep="\t", index=False)

    print("\nSaved apo-vs-holo comparison summary:")
    print(comparison_out)


# =========================
# Print compact overview
# =========================

print("\nCompact overview using probability >= 0.9:")
overview_cols = [
    "sample",
    "condition",
    "modification",
    "feature_type",
    "n_rows",
    "n_unique_transcript_sites",
    "prob_median",
    "n_modified_prob_ge_0p8",
    "modified_fraction_total_prob_ge_0p9",
    "modified_fraction_within_prob_ge_0p9",
]

overview_cols = [c for c in overview_cols if c in summary_df.columns]
print(summary_df[overview_cols].to_string(index=False))
