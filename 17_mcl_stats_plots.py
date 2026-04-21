#!/usr/bin/env python3
"""
17_mcl_stats_plots.py — Statistics and visualisations for MCL functional clusters.

Input:
    relational_table_functional_bins_multigenerational_microbiota.txt
    (output of 16_parse_mcl.py)
    Columns: cluster_id  protein_id  bin  generation  sample_type

Outputs (all saved to --outdir):
    Statistics tables (tab-separated):
        cluster_summary.txt         — per-cluster: n_proteins, n_bins, generations, sample_types present
        generation_sharing.txt      — how many clusters are shared across F0/F1/F2 vs unique
        sample_type_sharing.txt     — how many clusters are shared between cecum/feces vs unique
        bins_per_cluster.txt        — distribution of number of bins per cluster

    TIFF figures (300 dpi):
        01_cluster_size_distribution.tiff   — histogram of proteins per cluster
        02_bins_per_cluster.tiff            — histogram of bins per cluster
        03_generation_sharing.tiff          — bar chart: exclusive vs shared across generations
        07_bin_network.tiff                 — network of bins connected by functional repertoire similarity
        04_sample_type_sharing.tiff         — bar chart: exclusive vs shared across sample types
        05_proteins_per_generation.tiff     — boxplot of cluster sizes split by generation
        06_heatmap_clusters_shared.tiff     — heatmap: pairwise cluster sharing between generations x sample types

Usage:
    python 17_mcl_stats_plots.py \
        --input  functional_annotation/relational_table_functional_bins_multigenerational_microbiota.txt \
        --outdir functional_annotation/mcl_stats
"""

import argparse
import os
import sys

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")          # non-interactive backend for cluster use
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
import seaborn as sns
import networkx as nx


# ---------------------------------------------------------------------------
# Plotting style
# ---------------------------------------------------------------------------

plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.size":        11,
    "axes.titlesize":   13,
    "axes.labelsize":   12,
    "xtick.labelsize":  10,
    "ytick.labelsize":  10,
    "figure.dpi":       100,
    "savefig.dpi":      300,
    "savefig.format":   "tiff",
})

PALETTE = {
    "F0": "#4C72B0",
    "F1": "#DD8452",
    "F2": "#55A868",
    "cecum_samples": "#C44E52",
    "last_feces":    "#8172B3",
}


# ---------------------------------------------------------------------------
# Helper — save figure
# ---------------------------------------------------------------------------

def save_fig(fig, outdir: str, name: str) -> None:
    path = os.path.join(outdir, name)
    fig.savefig(path, dpi=300, bbox_inches="tight", format="tiff")
    plt.close(fig)
    print(f"  [saved] {path}")


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

def load_data(input_file: str) -> pd.DataFrame:
    print(f"[INFO] Loading {input_file} ...")
    df = pd.read_csv(input_file, sep="\t", dtype=str)
    required = {"cluster_id", "protein_id", "bin", "generation", "sample_type"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"[ERROR] Missing columns: {missing}")
    df["cluster_id"] = df["cluster_id"].astype(int)
    print(f"  Rows       : {len(df):,}")
    print(f"  Clusters   : {df['cluster_id'].nunique():,}")
    print(f"  Proteins   : {df['protein_id'].nunique():,}")
    print(f"  Bins       : {df['bin'].nunique():,}")
    print(f"  Generations: {sorted(df['generation'].unique())}")
    print(f"  Sample types: {sorted(df['sample_type'].unique())}")
    return df


# ---------------------------------------------------------------------------
# Build per-cluster summary
# ---------------------------------------------------------------------------

def build_cluster_summary(df: pd.DataFrame) -> pd.DataFrame:
    """One row per cluster with key statistics."""
    grp = df.groupby("cluster_id")

    summary = pd.DataFrame({
        "n_proteins":   grp["protein_id"].nunique(),
        "n_bins":       grp["bin"].nunique(),
        "n_generations":grp["generation"].nunique(),
        "generations":  grp["generation"].apply(lambda x: "|".join(sorted(x.unique()))),
        "n_sample_types":grp["sample_type"].nunique(),
        "sample_types": grp["sample_type"].apply(lambda x: "|".join(sorted(x.unique()))),
    }).reset_index()

    return summary


# ---------------------------------------------------------------------------
# Statistics tables
# ---------------------------------------------------------------------------

def stats_generation_sharing(summary: pd.DataFrame) -> pd.DataFrame:
    """Count clusters exclusive to one generation vs shared across 2 or 3."""
    counts = summary["generations"].value_counts().reset_index()
    counts.columns = ["generations_present", "n_clusters"]
    counts["sharing_class"] = counts["generations_present"].apply(
        lambda x: "exclusive" if "|" not in x else f"shared_{x.count('|') + 1}_generations"
    )
    return counts.sort_values("n_clusters", ascending=False)


def stats_sample_type_sharing(summary: pd.DataFrame) -> pd.DataFrame:
    """Count clusters exclusive to one sample type vs shared."""
    counts = summary["sample_types"].value_counts().reset_index()
    counts.columns = ["sample_types_present", "n_clusters"]
    counts["sharing_class"] = counts["sample_types_present"].apply(
        lambda x: "exclusive" if "|" not in x else "shared_both"
    )
    return counts.sort_values("n_clusters", ascending=False)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def plot_cluster_size_distribution(summary: pd.DataFrame, outdir: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))

    data = summary["n_proteins"]
    bins = np.logspace(0, np.log10(data.max() + 1), 50)

    ax.hist(data, bins=bins, color="#4C72B0", edgecolor="white", linewidth=0.4)
    ax.set_xscale("log")
    ax.set_xlabel("Proteins per cluster")
    ax.set_ylabel("Number of clusters")
    ax.set_title("Cluster size distribution")

    # Annotate median and mean
    ax.axvline(data.median(), color="tomato",  linestyle="--", linewidth=1.5,
               label=f"Median = {data.median():.0f}")
    ax.axvline(data.mean(),   color="orange",  linestyle=":",  linewidth=1.5,
               label=f"Mean = {data.mean():.1f}")
    ax.legend(fontsize=10)

    save_fig(fig, outdir, "01_cluster_size_distribution.tiff")


def plot_bins_per_cluster(summary: pd.DataFrame, outdir: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))

    data = summary["n_bins"]
    ax.hist(data, bins=range(1, data.max() + 2), color="#55A868",
            edgecolor="white", linewidth=0.4, align="left")
    ax.set_xlabel("Bins per cluster")
    ax.set_ylabel("Number of clusters")
    ax.set_title("Number of bins contributing to each cluster")
    ax.axvline(data.median(), color="tomato", linestyle="--", linewidth=1.5,
               label=f"Median = {data.median():.0f}")
    ax.legend(fontsize=10)

    save_fig(fig, outdir, "02_bins_per_cluster.tiff")


def plot_generation_sharing(gen_stats: pd.DataFrame, outdir: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))

    labels = gen_stats["generations_present"].tolist()
    values = gen_stats["n_clusters"].tolist()
    colors = [PALETTE.get(l.split("|")[0], "#888888") if "|" not in l
              else "#AAAAAA" for l in labels]

    bars = ax.bar(labels, values, color=colors, edgecolor="white", linewidth=0.5)
    ax.set_xlabel("Generations present in cluster")
    ax.set_ylabel("Number of clusters")
    ax.set_title("Cluster sharing across generations")
    ax.tick_params(axis="x", rotation=30)

    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max(values) * 0.01,
                f"{val:,}", ha="center", va="bottom", fontsize=9)

    save_fig(fig, outdir, "03_generation_sharing.tiff")


def plot_sample_type_sharing(st_stats: pd.DataFrame, outdir: str) -> None:
    fig, ax = plt.subplots(figsize=(6, 5))

    labels = st_stats["sample_types_present"].tolist()
    values = st_stats["n_clusters"].tolist()
    colors = [PALETTE.get(l, "#AAAAAA") for l in labels]

    bars = ax.bar(labels, values, color=colors, edgecolor="white", linewidth=0.5)
    ax.set_xlabel("Sample types present in cluster")
    ax.set_ylabel("Number of clusters")
    ax.set_title("Cluster sharing between sample types")
    ax.tick_params(axis="x", rotation=20)

    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max(values) * 0.01,
                f"{val:,}", ha="center", va="bottom", fontsize=9)

    save_fig(fig, outdir, "04_sample_type_sharing.tiff")


def plot_proteins_per_generation(df: pd.DataFrame, outdir: str) -> None:
    """Boxplot: distribution of cluster sizes split by generation."""

    # For each cluster, count proteins per generation
    data = (df.groupby(["cluster_id", "generation"])["protein_id"]
              .nunique()
              .reset_index()
              .rename(columns={"protein_id": "n_proteins"}))

    generations = sorted(data["generation"].unique())
    colors = [PALETTE.get(g, "#888888") for g in generations]

    fig, ax = plt.subplots(figsize=(7, 5))

    groups = [data.loc[data["generation"] == g, "n_proteins"].values
              for g in generations]

    bp = ax.boxplot(groups, patch_artist=True, labels=generations,
                    medianprops={"color": "black", "linewidth": 1.5},
                    flierprops={"marker": "o", "markersize": 2,
                                "markerfacecolor": "grey", "alpha": 0.4})

    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_yscale("log")
    ax.set_xlabel("Generation")
    ax.set_ylabel("Proteins per cluster (log scale)")
    ax.set_title("Cluster size by generation")

    save_fig(fig, outdir, "05_proteins_per_generation.tiff")


def plot_heatmap_sharing(df: pd.DataFrame, outdir: str) -> None:
    """
    Heatmap: for each pair of (generation x sample_type) groups,
    how many clusters do they share?
    """
    # Create a group label per protein
    df = df.copy()
    df["group"] = df["generation"] + "_" + df["sample_type"]

    # For each cluster, get the set of groups present
    cluster_groups = (df.groupby("cluster_id")["group"]
                        .apply(set)
                        .reset_index())

    groups = sorted(df["group"].unique())
    n = len(groups)
    matrix = pd.DataFrame(0, index=groups, columns=groups)

    for _, row in cluster_groups.iterrows():
        present = list(row["group"])
        for i in range(len(present)):
            for j in range(len(present)):
                matrix.loc[present[i], present[j]] += 1

    fig, ax = plt.subplots(figsize=(9, 7))

    cmap = LinearSegmentedColormap.from_list(
        "blue_white", ["#f7fbff", "#2171b5"]
    )

    sns.heatmap(
        matrix,
        ax=ax,
        cmap=cmap,
        annot=True,
        fmt=",",
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "Shared clusters"},
    )

    ax.set_title("Pairwise cluster sharing\n(generation × sample type)")
    ax.tick_params(axis="x", rotation=45)
    ax.tick_params(axis="y", rotation=0)

    save_fig(fig, outdir, "06_heatmap_clusters_shared.tiff")


# ---------------------------------------------------------------------------
# Network plot — bins connected by functional repertoire similarity
# ---------------------------------------------------------------------------

def compute_bin_similarity(df: pd.DataFrame) -> pd.DataFrame:
    """
    For every pair of bins compute functional repertoire similarity:
        similarity = |shared clusters| / max(|clusters_A|, |clusters_B|)
    This is the same metric used in the paper (Zhu et al. 2015).

    Returns a DataFrame with columns: bin_a, bin_b, similarity, shared_clusters
    """
    print("[INFO] Computing pairwise bin functional similarity ...")

    # For each bin, get the set of cluster IDs it participates in
    bin_clusters = (df.groupby("bin")["cluster_id"]
                      .apply(set)
                      .to_dict())

    bins = list(bin_clusters.keys())
    n    = len(bins)
    print(f"  Bins to compare: {n:,}  ({n*(n-1)//2:,} pairs)")

    records = []
    for i in range(n):
        a      = bins[i]
        set_a  = bin_clusters[a]
        size_a = len(set_a)
        for j in range(i + 1, n):
            b       = bins[j]
            set_b   = bin_clusters[b]
            shared  = len(set_a & set_b)
            if shared == 0:
                continue
            sim = shared / max(size_a, len(set_b))
            records.append((a, b, sim, shared))

    sim_df = pd.DataFrame(records,
                          columns=["bin_a", "bin_b", "similarity", "shared_clusters"])
    print(f"  Non-zero pairs: {len(sim_df):,}")
    return sim_df


def plot_bin_network(df: pd.DataFrame, outdir: str,
                     sim_threshold: float = 0.1) -> None:
    """
    Network where:
      - Nodes  = bins, coloured by generation × sample_type
      - Edges  = functional repertoire similarity >= sim_threshold
      - Node size proportional to number of functional clusters in that bin
      - Edge width / alpha proportional to similarity
    """
    print(f"[INFO] Building bin network (similarity threshold = {sim_threshold}) ...")

    # Bin metadata (generation, sample_type)
    bin_meta = (df[["bin", "generation", "sample_type"]]
                  .drop_duplicates("bin")
                  .set_index("bin"))

    # Bin repertoire size (number of distinct clusters)
    bin_size = df.groupby("bin")["cluster_id"].nunique()

    # Pairwise similarity
    sim_df = compute_bin_similarity(df)
    sim_df = sim_df[sim_df["similarity"] >= sim_threshold]
    print(f"  Edges after threshold filter: {len(sim_df):,}")

    # Build networkx graph
    G = nx.Graph()
    for b in bin_meta.index:
        G.add_node(b)

    for _, row in sim_df.iterrows():
        G.add_edge(row["bin_a"], row["bin_b"],
                   weight=row["similarity"],
                   shared=row["shared_clusters"])

    print(f"  Nodes: {G.number_of_nodes():,}  |  Edges: {G.number_of_edges():,}")

    # Node colour: generation × sample_type
    GROUP_COLORS = {
        "F0_cecum_samples": "#2166AC",
        "F0_last_feces":    "#92C5DE",
        "F1_cecum_samples": "#D6604D",
        "F1_last_feces":    "#FDDBC7",
        "F2_cecum_samples": "#1A9850",
        "F2_last_feces":    "#A6D96A",
    }

    node_colors = []
    node_sizes  = []
    for node in G.nodes():
        if node in bin_meta.index:
            gen = bin_meta.loc[node, "generation"]
            st  = bin_meta.loc[node, "sample_type"]
            node_colors.append(GROUP_COLORS.get(f"{gen}_{st}", "#AAAAAA"))
            # Scale node size: min 20, max 300
            n_clust = bin_size.get(node, 1)
            node_sizes.append(20 + 280 * (n_clust / bin_size.max()))
        else:
            node_colors.append("#AAAAAA")
            node_sizes.append(30)

    # Edge width and alpha proportional to similarity
    edge_weights = [G[u][v]["weight"] for u, v in G.edges()]
    max_w        = max(edge_weights) if edge_weights else 1.0
    edge_widths  = [0.3 + 2.5 * (w / max_w) for w in edge_weights]
    edge_alphas  = [0.15 + 0.65 * (w / max_w) for w in edge_weights]

    # Layout — Kamada-Kawai gives clearer cluster separation than spring
    print("  Computing layout (Kamada-Kawai) ...")
    pos = nx.kamada_kawai_layout(G, weight="weight")

    # Draw
    fig, ax = plt.subplots(figsize=(14, 12))
    ax.set_facecolor("#F8F8F8")
    fig.patch.set_facecolor("#F8F8F8")

    # Draw edges in batches by alpha band for visual clarity
    for (u, v), w, lw, alpha in zip(G.edges(), edge_weights,
                                     edge_widths, edge_alphas):
        ax.plot([pos[u][0], pos[v][0]],
                [pos[u][1], pos[v][1]],
                color="#888888", linewidth=lw, alpha=alpha, zorder=1)

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, ax=ax,
                           node_color=node_colors,
                           node_size=node_sizes,
                           linewidths=0.4,
                           edgecolors="white",
                           alpha=0.9)

    ax.set_title(
        f"Functional repertoire similarity network (threshold ≥ {sim_threshold})\n"
        f"Node size ∝ functional repertoire size | "
        f"Edge width ∝ similarity",
        fontsize=12, pad=12
    )
    ax.axis("off")

    # Legend
    legend_handles = [
        mpatches.Patch(color=color, label=label.replace("_", " "))
        for label, color in GROUP_COLORS.items()
    ]
    ax.legend(handles=legend_handles, loc="lower left",
              fontsize=9, framealpha=0.8, title="Generation × Sample type",
              title_fontsize=9)

    save_fig(fig, outdir, "07_bin_network.tiff")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(input_file: str, outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)

    # Load
    df = load_data(input_file)

    # Per-cluster summary
    print("\n[INFO] Building per-cluster summary ...")
    summary = build_cluster_summary(df)

    # Statistics tables
    print("\n[INFO] Writing statistics tables ...")

    out_summary = os.path.join(outdir, "cluster_summary.txt")
    summary.to_csv(out_summary, sep="\t", index=False)
    print(f"  [saved] {out_summary}")

    gen_stats = stats_generation_sharing(summary)
    out_gen = os.path.join(outdir, "generation_sharing.txt")
    gen_stats.to_csv(out_gen, sep="\t", index=False)
    print(f"  [saved] {out_gen}")

    st_stats = stats_sample_type_sharing(summary)
    out_st = os.path.join(outdir, "sample_type_sharing.txt")
    st_stats.to_csv(out_st, sep="\t", index=False)
    print(f"  [saved] {out_st}")

    bins_dist = (summary["n_bins"].value_counts()
                                  .reset_index()
                                  .rename(columns={"count": "n_clusters"})
                                  .sort_values("n_bins"))
    out_bins = os.path.join(outdir, "bins_per_cluster.txt")
    bins_dist.to_csv(out_bins, sep="\t", index=False)
    print(f"  [saved] {out_bins}")

    # Figures
    print("\n[INFO] Generating figures ...")
    plot_cluster_size_distribution(summary, outdir)
    plot_bins_per_cluster(summary, outdir)
    plot_generation_sharing(gen_stats, outdir)
    plot_sample_type_sharing(st_stats, outdir)
    plot_proteins_per_generation(df, outdir)
    plot_heatmap_sharing(df, outdir)
    plot_bin_network(df, outdir, sim_threshold=0.1)

    print("\n[INFO] All done.")
    print(f"  Output directory: {outdir}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Statistics and TIFF plots for MCL functional clusters.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input",  "-i", required=True,
                   help="Relational table from 16_parse_mcl.py")
    p.add_argument("--outdir", "-o", required=True,
                   help="Output directory for tables and figures.")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if not os.path.isfile(args.input):
        sys.exit(f"[ERROR] Input not found: {args.input}")
    main(args.input, args.outdir)
