#!/usr/bin/env python3
"""
17_bins_per_cluster.py — Histogram of bins per cluster (plot 02).

Input:
    relational_table_functional_bins_multigenerational_microbiota.txt
    (output of 16_parse_mcl.py)
    Columns: cluster_id  protein_id  bin  generation  sample_type

Output (saved to --outdir):
    02_bins_per_cluster.tiff
        Histogram of how many distinct bins each cluster contains.
        Uses a broken y-axis when the n_bins=1 bar dwarfs the rest,
        so all bars remain visible.

Usage:
    python 17_bins_per_cluster.py \
        --input  functional_annotation/relational_table_functional_bins_multigenerational_microbiota.txt \
        --outdir functional_annotation/mcl_stats
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker


# ---------------------------------------------------------------------------
# Plotting style
# ---------------------------------------------------------------------------

plt.rcParams.update({
    "font.family":    "sans-serif",
    "font.size":      11,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "xtick.labelsize":10,
    "ytick.labelsize":10,
    "savefig.dpi":    300,
    "savefig.format": "tiff",
})


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def save_fig(fig, outdir: str, name: str) -> None:
    path = os.path.join(outdir, name)
    fig.savefig(path, dpi=300, bbox_inches="tight", format="tiff")
    plt.close(fig)
    print(f"  [saved] {path}")


def load_data(input_file: str) -> pd.DataFrame:
    print(f"[INFO] Loading {input_file} ...")
    df = pd.read_csv(input_file, sep="\t", dtype=str)
    required = {"cluster_id", "protein_id", "bin", "generation", "sample_type"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"[ERROR] Missing columns: {missing}")
    df["cluster_id"] = df["cluster_id"].astype(int)
    print(f"  Rows     : {len(df):,}")
    print(f"  Clusters : {df['cluster_id'].nunique():,}")
    print(f"  Bins     : {df['bin'].nunique():,}")
    return df


def build_summary(df: pd.DataFrame) -> pd.DataFrame:
    """One row per cluster with n_bins (distinct bins)."""
    return (df.groupby("cluster_id")["bin"]
              .nunique()
              .reset_index()
              .rename(columns={"bin": "n_bins"}))


# ---------------------------------------------------------------------------
# Plot 02 — bins per cluster histogram with broken y-axis
# ---------------------------------------------------------------------------

def plot_bins_per_cluster(summary: pd.DataFrame, outdir: str) -> None:
    """
    Histogram of bins per cluster using a log-spaced x-axis.
    The x range spans from 1 to the maximum, with log-spaced bins so that
    bars at both ends of the distribution are clearly visible regardless of
    how wide the range is (here 1 to ~900).
    """
    data     = summary["n_bins"]
    max_bins = int(data.max())
    median_val = data.median()

    print(f"  Max bins per cluster : {max_bins}")
    print(f"  Median bins/cluster  : {median_val:.1f}")
    print(f"  Mean bins/cluster    : {data.mean():.1f}")

    # Log-spaced bins: gives visible bar widths across the full 1–max_bins range
    log_bins = np.logspace(0, np.log10(max_bins + 1), 60)

    fig, ax = plt.subplots(figsize=(9, 5))

    ax.hist(data, bins=log_bins, color="#55A868",
            edgecolor="white", linewidth=0.4)
    ax.set_xscale("log")

    ax.set_xlabel("Bins per cluster (log scale)")
    ax.set_ylabel("Number of clusters")
    ax.set_title("Number of bins contributing to each cluster")

    ax.axvline(median_val, color="tomato",  linestyle="--", linewidth=1.5,
               label=f"Median = {median_val:.0f}")
    ax.axvline(data.mean(), color="orange", linestyle=":",  linewidth=1.5,
               label=f"Mean = {data.mean():.1f}")
    ax.legend(fontsize=10)

    # Clean x-tick labels at round numbers
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_major_locator(
        matplotlib.ticker.LogLocator(base=10, numticks=8)
    )

    save_fig(fig, outdir, "02_bins_per_cluster.tiff")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(input_file: str, outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)
    df      = load_data(input_file)
    summary = build_summary(df)
    print("\n[INFO] Generating bins-per-cluster histogram ...")
    plot_bins_per_cluster(summary, outdir)
    print("\n[INFO] Done.")


def parse_args():
    p = argparse.ArgumentParser(
        description="Histogram of bins per cluster with broken y-axis (plot 02).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input",  "-i", required=True,
                   help="Relational table from 16_parse_mcl.py")
    p.add_argument("--outdir", "-o", required=True,
                   help="Output directory for the TIFF figure.")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if not os.path.isfile(args.input):
        sys.exit(f"[ERROR] Input not found: {args.input}")
    main(args.input, args.outdir)

