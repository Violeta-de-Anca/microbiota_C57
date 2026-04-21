#!/usr/bin/env python3
"""
17_heatmap_sharing.py — Cluster-sharing heatmap and cluster composition plot.

Input:
    relational_table_functional_bins_multigenerational_microbiota.txt
    (output of 16_parse_mcl.py)
    Columns: cluster_id  protein_id  bin  generation  sample_type

Outputs (saved to --outdir):
    06_heatmap_clusters_shared.tiff
        Pairwise heatmap: cell(i,j) = % of group i's clusters also in group j.
        Diagonal = 100 % by definition. Groups = generation × sample_type.

    07_cluster_composition.tiff
        Stacked horizontal bar chart of the top N multi-bin clusters.
        Each bar = one cluster, sorted by total bin count (descending).
        Segments show how many bins within that cluster belong to each
        generation × sample_type group.

Usage:
    python 17_heatmap_sharing.py \
        --input  functional_annotation/relational_table_functional_bins_multigenerational_microbiota.txt \
        --outdir functional_annotation/mcl_stats \
        [--top_n 50]
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
import seaborn as sns


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

GROUP_COLORS = {
    "F0_cecum_samples": "#2166AC",
    "F0_last_feces":    "#92C5DE",
    "F1_cecum_samples": "#D6604D",
    "F1_last_feces":    "#FDDBC7",
    "F2_cecum_samples": "#1A9850",
    "F2_last_feces":    "#A6D96A",
}


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
    print(f"  Rows        : {len(df):,}")
    print(f"  Clusters    : {df['cluster_id'].nunique():,}")
    print(f"  Bins        : {df['bin'].nunique():,}")
    print(f"  Generations : {sorted(df['generation'].unique())}")
    print(f"  Sample types: {sorted(df['sample_type'].unique())}")
    return df


# ---------------------------------------------------------------------------
# Plot 06 — pairwise sharing heatmap
# ---------------------------------------------------------------------------

def plot_heatmap_sharing(df: pd.DataFrame, outdir: str) -> None:
    """
    Pairwise heatmap: cell(i, j) = % of group i's clusters that are
    also present in group j. Diagonal = 100 % by definition.
    Groups are defined as generation × sample_type combinations.
    """
    df = df.copy()
    df["group"] = df["generation"] + "_" + df["sample_type"]

    # For each cluster record which groups are present
    cluster_groups = (df.groupby("cluster_id")["group"]
                        .apply(set)
                        .reset_index())

    groups = sorted(df["group"].unique())
    matrix = pd.DataFrame(0.0, index=groups, columns=groups)

    for _, row in cluster_groups.iterrows():
        present = list(row["group"])
        for i in range(len(present)):
            for j in range(len(present)):
                matrix.loc[present[i], present[j]] += 1.0

    # Row-normalise: cell(i,j) = shared(i,j) / own(i) * 100
    # Replace 0 on diagonal with NaN to avoid division by zero
    diag = pd.Series(np.diag(matrix.values), index=matrix.index)
    diag = diag.replace(0, np.nan)
    pct_matrix = matrix.div(diag, axis=0).mul(100).round(1)

    annot_data = pct_matrix.map(lambda x: f"{x:.1f}%" if not np.isnan(x) else "—")

    cmap = LinearSegmentedColormap.from_list(
        "blue_white", ["#f7fbff", "#2171b5"]
    )

    fig, ax = plt.subplots(figsize=(9, 7))

    sns.heatmap(
        pct_matrix,
        ax=ax,
        cmap=cmap,
        annot=annot_data,
        fmt="",
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "% of row-group clusters"},
        vmin=0, vmax=100,
    )

    ax.set_title(
        "Pairwise cluster sharing\n"
        "(generation × sample type | row = % of that group's clusters shared with column)"
    )
    ax.tick_params(axis="x", rotation=45)
    ax.tick_params(axis="y", rotation=0)

    save_fig(fig, outdir, "06_heatmap_clusters_shared.tiff")


# ---------------------------------------------------------------------------
# Plot 07 — cluster composition stacked bar chart
# ---------------------------------------------------------------------------

def plot_cluster_composition(df: pd.DataFrame, outdir: str,
                              top_n: int = 50) -> None:
    """
    Stacked horizontal bar chart of the top N multi-bin clusters (sorted by
    total number of distinct bins, descending).

    Each bar represents one cluster; its total length is the number of distinct
    bins in that cluster. Each coloured segment shows how many of those bins
    come from a given generation × sample_type group.

    Only clusters with >= 2 distinct bins are considered, since single-bin
    clusters carry no cross-group information.
    """
    print(f"[INFO] Building cluster composition chart (top {top_n} multi-bin clusters) ...")

    df = df.copy()
    df["group"] = df["generation"] + "_" + df["sample_type"]

    # Count distinct bins per cluster per group
    counts = (df.groupby(["cluster_id", "group"])["bin"]
                .nunique()
                .reset_index()
                .rename(columns={"bin": "n_bins"}))

    # Total bins per cluster — sort descending
    totals = (counts.groupby("cluster_id")["n_bins"]
                    .sum()
                    .sort_values(ascending=False))

    # Keep only multi-bin clusters
    multi_bin_totals = totals[totals >= 2]
    n_multi = len(multi_bin_totals)
    n_shown = min(top_n, n_multi)
    print(f"  Multi-bin clusters: {n_multi:,}  —  showing top {n_shown}")

    if n_shown == 0:
        print("  [WARNING] No multi-bin clusters found — skipping composition plot.")
        return

    top_ids = multi_bin_totals.head(n_shown).index
    counts_top = counts[counts["cluster_id"].isin(top_ids)]

    # Pivot: rows = clusters (in size order), columns = groups
    pivot = (counts_top
             .pivot_table(index="cluster_id", columns="group",
                          values="n_bins", fill_value=0)
             .reindex(top_ids))

    # Keep columns in GROUP_COLORS order (only those present in data)
    present_groups = [g for g in GROUP_COLORS if g in pivot.columns]
    # Any group in data but not in GROUP_COLORS gets a neutral colour
    extra_groups  = [g for g in pivot.columns if g not in GROUP_COLORS]
    ordered_cols  = present_groups + extra_groups
    pivot = pivot[ordered_cols]

    colors = [GROUP_COLORS.get(g, "#AAAAAA") for g in ordered_cols]

    # Figure height scales with number of bars
    fig_height = max(6, n_shown * 0.32 + 2)
    fig, ax = plt.subplots(figsize=(11, fig_height))

    # Draw stacked bars manually so we can add total labels cleanly
    lefts = np.zeros(n_shown)
    y_pos = np.arange(n_shown)

    for col, color in zip(ordered_cols, colors):
        vals = pivot[col].values.astype(float)
        ax.barh(y_pos, vals, left=lefts,
                color=color, edgecolor="white", linewidth=0.3)
        lefts += vals

    # Y-axis: cluster labels
    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"Cluster {cid}" for cid in pivot.index], fontsize=8)
    ax.invert_yaxis()   # largest at top

    # Total bin count label to the right of each bar
    totals_shown = multi_bin_totals.head(n_shown).values
    x_max = totals_shown.max()
    for i, total in enumerate(totals_shown):
        ax.text(total + x_max * 0.005, i, str(int(total)),
                va="center", ha="left", fontsize=7, color="#444444")

    ax.set_xlabel("Number of bins")
    ax.set_xlim(right=x_max * 1.08)    # leave room for labels
    ax.set_title(
        f"Cluster composition — top {n_shown} multi-bin clusters\n"
        f"(sorted by total bins | coloured by generation × sample type)"
    )

    # Legend
    handles = [mpatches.Patch(color=GROUP_COLORS.get(g, "#AAAAAA"),
                               label=g.replace("_", " "))
               for g in ordered_cols]
    ax.legend(handles=handles, loc="lower right", fontsize=9,
              title="Generation × Sample type", title_fontsize=9,
              framealpha=0.8)

    save_fig(fig, outdir, "07_cluster_composition.tiff")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(input_file: str, outdir: str, top_n: int) -> None:
    os.makedirs(outdir, exist_ok=True)
    df = load_data(input_file)

    print("\n[INFO] Generating heatmap (plot 06) ...")
    plot_heatmap_sharing(df, outdir)

    print("\n[INFO] Generating cluster composition chart (plot 07) ...")
    plot_cluster_composition(df, outdir, top_n=top_n)

    print("\n[INFO] Done.")


def parse_args():
    p = argparse.ArgumentParser(
        description="Cluster-sharing heatmap and cluster composition chart.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input",  "-i", required=True,
                   help="Relational table from 16_parse_mcl.py")
    p.add_argument("--outdir", "-o", required=True,
                   help="Output directory for the TIFF figures.")
    p.add_argument("--top_n", "-n", type=int, default=50,
                   help="Number of top multi-bin clusters to show in plot 07.")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if not os.path.isfile(args.input):
        sys.exit(f"[ERROR] Input not found: {args.input}")
    main(args.input, args.outdir, args.top_n)
