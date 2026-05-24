#!/usr/bin/env python3
"""
19_relational_table_abundance_functionality_multigen.py
==================================

Goal
----
Build a wide matrix where:

    rows    = clusters from clusterq_bin_protein_proportion_matrix.tsv
    columns = every (bin, individual) pair — one column for each individual
              that has that bin in their quantification file
    cell    = abundance(bin, individual) * proportion(cluster, bin)

NA semantics
------------
The matrix is wide on (bin x individual), so the column for each bin
exists in every cluster row, even clusters where that bin doesn't belong.
We distinguish two kinds of "no signal":

  * proportion(cluster, bin) == 0  =>  bin NOT in cluster
        Cell written as "NA" (not a real measurement).
  * proportion(cluster, bin) > 0 but abundance(bin, individual) == 0
        =>  bin IS in cluster, but this individual has no abundance
        Cell written as "0" (a real zero).

A cluster row that is all-NA across a bin's columns means the bin doesn't
contribute to that cluster at all. A row with mixed 0s and numbers means
the bin IS in the cluster, and only some individuals carry it.

NO summation across bins. With 938 bins and per-file individual counts
of 17..49, total column count is ~27,138 and the full TSV is ~1 GB
decompressed / ~100 MB gzipped.

Inputs
------
  1) functional_annotatioin/cluster_bin_protein_proportion_matrix.tsv
        rows = clusters; cols = fully qualified bins (F0_cecum_samples_bin.1, ...)
        cells = fraction of that bin's proteins in this cluster.
  2) quantification/bin_abundance_table_<gen>_<cecum|last_feces>_*.tab
        rows = "bin.<NNN>"; cols = per-individual abundance values.
        We prepend "<generation>_<sample_type>_" so bin names line up
        with the proportion-matrix bins. <sample_type> is "cecum_samples"
        when the filename says "cecum" (the cecum files lack the trailing
        "_samples" used in the proteome FASTA), and "last_feces" otherwise.

Outputs (next to the proportion matrix in functional_annotatioin/)
------------------------------------------------------------------
  - cluster_individual_per_bin_weighted_abundance_matrix.tsv.gz
        First column "cluster", then one column per (bin, individual) pair
        with header "<qualified_bin>__<individual>".
  - intermediate per-file blocks (uncompressed TSVs):
        cluster_individual_per_bin_block_<gen>_<sampletype>[_part_NNNN_MMMM].tsv
        cluster_individual_per_bin_block_cluster_index.tsv

Why we split the work into per-file blocks
------------------------------------------
The sandbox shell budget per call is ~45s and the OneDrive mount writes
at ~9 MB/s, so a 1 GB TSV cannot land in a single call. Strategy:
  1) build one block per quantification file (split very large blocks
     by bin-range when they wouldn't fit in 45s),
  2) write the cluster-name column once,
  3) `paste` everything side-by-side and gzip the result.

For each block we ALSO stage in /tmp first (writing/seding directly on
the OneDrive mount is too slow), then do a single cp into the project
folder. /tmp is local (~1 GB/s) so the in-place sed replacing "nan" with
"NA" runs in seconds.

Usage
-----
    python3 build_cluster_individual_matrix.py <idx>                   # full block
    python3 build_cluster_individual_matrix.py <idx> <bstart> <bend>   # bin slice
    python3 build_cluster_individual_matrix.py cluster_index           # cluster col
    python3 build_cluster_individual_matrix.py combine                 # paste+gzip
    python3 build_cluster_individual_matrix.py                         # all in one
"""

import os, re, sys, glob, time, gzip, shutil, subprocess
import pandas as pd
import numpy as np

# Paths
ROOT = "/proj/naiss2024-23-57/C57_female_lineage_microbiota"
PROP_MATRIX_FILE   = f"{ROOT}/functional_annotation/cluster_bin_protein_proportion_matrix.tsv"
QUANT_DIR          = f"{ROOT}/quantification"
OUT_DIR            = f"{ROOT}/quantification"
FINAL_OUT_GZ       = f"{OUT_DIR}/cluster_individual_per_bin_weighted_abundance_matrix.tsv.gz"
CLUSTER_INDEX_FILE = f"{OUT_DIR}/cluster_individual_per_bin_block_cluster_index.tsv"
TMP_DIR            = f"{ROOT}/functional_annotation/mcl/cluster_individual_blocks"
STAGING_TSV        = f"{ROOT}/functional_annotation/mcl/cluster_individual_per_bin_weighted_abundance_matrix.tsv"

FNAME_PATTERN = re.compile(r"^bin_abundance_table_(F\d)_(cecum|last_feces)_.*\.tab$")
SAMPLE_TYPE_MAP = {"cecum": "cecum_samples", "last_feces": "last_feces"}


def list_quant_files():
    return sorted(glob.glob(os.path.join(QUANT_DIR, "bin_abundance_table_*.tab")))


def block_out_path(path):
    fname = os.path.basename(path)
    m = FNAME_PATTERN.match(fname)
    if m is None:
        # Caller passed a path that doesn't match our expected naming;
        # bail out loudly so the bug is obvious instead of silently
        # producing an unindexable None.group(...) error.
        raise ValueError(f"filename does not match expected pattern: {fname}")
    generation, sample_short = m.group(1), m.group(2)
    return os.path.join(
        OUT_DIR,
        f"cluster_individual_per_bin_block_{generation}_{SAMPLE_TYPE_MAP[sample_short]}.tsv",
    )


def build_block(path, prop_df, bin_start=None, bin_end=None):
    """Write one (file, optional bin-slice) block.

    Cell semantics:
        if proportion[c, b] == 0: "NA"  (bin not in cluster)
        else:                     proportion[c, b] * abundance[b, i]

    Pipeline (all in-budget):
      1) compute the dense block (NaN where proportion==0)
      2) np.savetxt to /tmp  (fast, ~1 GB/s)
      3) sed s/nan/NA/g IN /tmp  (fast)
      4) cp /tmp -> OneDrive  (one write at ~9 MB/s)
    """
    t0 = time.time()
    fname = os.path.basename(path)
    m = FNAME_PATTERN.match(fname)
    if not m:
        print(f"  SKIP (bad filename): {fname}", flush=True)
        return
    bin_prefix = f"{m.group(1)}_{SAMPLE_TYPE_MAP[m.group(2)]}_"

    # Load abundance table; coerce to floats.
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip()
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    bins_qualified = [bin_prefix + b for b in df.index]
    indivs = list(df.columns)

    # Drop bins not in proportion matrix (defensive).
    keep = [i for i, b in enumerate(bins_qualified) if b in prop_df.columns]
    if len(keep) != len(bins_qualified):
        bins_qualified = [bins_qualified[i] for i in keep]
        df = df.iloc[keep]

    # Optional bin slice (used to split very large blocks).
    is_partial = (bin_start is not None) or (bin_end is not None)
    if bin_start is None: bin_start = 0
    if bin_end   is None: bin_end   = len(bins_qualified)
    bins_qualified = bins_qualified[bin_start:bin_end]
    df = df.iloc[bin_start:bin_end]

    A = df.to_numpy(dtype=np.float64)                       # (n_bins, n_indivs)
    P = prop_df[bins_qualified].to_numpy(dtype=np.float64)  # (n_clusters, n_bins)

    n_clusters = P.shape[0]
    n_bins   = len(bins_qualified)
    n_indivs = len(indivs)

    # Column names match the array's flatten layout: for each bin in order,
    # all individuals.
    col_names = [f"{b}__{ind}" for b in bins_qualified for ind in indivs]

    final_path = block_out_path(path)
    if is_partial:
        final_path = final_path[:-4] + f"_part_{bin_start:04d}_{bin_end:04d}.tsv"
    print(f"  -> {final_path}", flush=True)
    print(f"     {n_bins} bins x {n_indivs} indivs = {len(col_names)} cols", flush=True)

    # Build dense block bin-by-bin (bounded peak memory).
    # NaN marks "bin not in cluster" -> printed as "nan" -> sed -> "NA".
    out_arr = np.empty((n_clusters, n_bins * n_indivs), dtype=np.float64)
    for j in range(n_bins):
        p_col = P[:, j]
        sub   = p_col[:, None] * A[j:j+1, :]
        sub[p_col == 0.0, :] = np.nan
        out_arr[:, j*n_indivs:(j+1)*n_indivs] = sub

    # Stage in /tmp (fast local writes), then sed.
    # The cp /tmp -> OneDrive happens in a SEPARATE shell call via the
    # `publish` command — that way savetxt+sed and cp each fit in 45s
    # even for the biggest blocks.
    os.makedirs(TMP_DIR, exist_ok=True)
    tmp_path = os.path.join(TMP_DIR, os.path.basename(final_path))
    header_line = "\t".join(col_names) + "\n"
    with open(tmp_path, "wb") as out:
        out.write(header_line.encode())
        np.savetxt(out, out_arr, fmt="%g", delimiter="\t")
    t_writetmp = time.time()
    subprocess.check_call(["sed", "-i", "s/\\bnan\\b/NA/g", tmp_path])
    t_sed = time.time()
    print(f"     STAGED {os.path.getsize(tmp_path)/1e6:.1f} MB to {tmp_path} "
          f"(savetxt={t_writetmp-t0:.1f}s, sed={t_sed-t_writetmp:.1f}s)",
          flush=True)
    print(f"     run `publish` to copy /tmp blocks to {OUT_DIR}", flush=True)


def write_cluster_index():
    prop_df = pd.read_csv(PROP_MATRIX_FILE, sep="\t", index_col=0, usecols=[0])
    with open(CLUSTER_INDEX_FILE, "w") as out:
        out.write("cluster\n")
        for c in prop_df.index:
            out.write(f"{c}\n")
    print(f"  wrote {CLUSTER_INDEX_FILE} ({len(prop_df.index)} clusters)", flush=True)


def cmd_block(idx, bin_start=None, bin_end=None):
    print("Loading proportion matrix ...", flush=True)
    prop_df = pd.read_csv(PROP_MATRIX_FILE, sep="\t", index_col=0)
    files = list_quant_files()
    if idx < 0 or idx >= len(files):
        raise SystemExit(f"index {idx} out of range 0..{len(files)-1}")
    print(f"Building block {idx}: {os.path.basename(files[idx])}", flush=True)
    build_block(files[idx], prop_df, bin_start=bin_start, bin_end=bin_end)


def cmd_combine():
    """Paste cluster-index + every block side-by-side into a /tmp staging
    TSV, gzip it, then copy the .gz to the project folder."""
    all_blocks = sorted(glob.glob(os.path.join(
        OUT_DIR, "cluster_individual_per_bin_block_*.tsv")))
    all_blocks = [p for p in all_blocks if p != CLUSTER_INDEX_FILE]

    # If a block has any "_part_" slices, drop the un-suffixed full file
    # (the slices are authoritative).
    blocks_with_parts = set()
    for p in all_blocks:
        base = os.path.basename(p)
        m = re.match(r"^(cluster_individual_per_bin_block_.+?)_part_\d+_\d+\.tsv$", base)
        if m:
            blocks_with_parts.add(os.path.join(OUT_DIR, m.group(1) + ".tsv"))
    final_blocks = [p for p in all_blocks if p not in blocks_with_parts]

    if not os.path.exists(CLUSTER_INDEX_FILE):
        write_cluster_index()
    parts = [CLUSTER_INDEX_FILE] + final_blocks

    print(f"Pasting {len(parts)} parts into {STAGING_TSV} (in /tmp for speed)",
          flush=True)
    for p in parts:
        print(f"  + {os.path.basename(p)}", flush=True)
    with open(STAGING_TSV, "w") as out:
        subprocess.check_call(["paste"] + parts, stdout=out)
    print(f"  staging TSV: {os.path.getsize(STAGING_TSV)/1e6:.1f} MB", flush=True)

    tmp_gz = STAGING_TSV + ".gz"
    print(f"gzipping -> {tmp_gz}", flush=True)
    with open(STAGING_TSV, "rb") as src, gzip.open(tmp_gz, "wb", compresslevel=1) as dst:
        shutil.copyfileobj(src, dst, length=8 * 1024 * 1024)
    print(f"  gz size: {os.path.getsize(tmp_gz)/1e6:.1f} MB", flush=True)

    print(f"copying {tmp_gz} -> {FINAL_OUT_GZ}", flush=True)
    shutil.copyfile(tmp_gz, FINAL_OUT_GZ)
    print(f"Final file: {os.path.getsize(FINAL_OUT_GZ)/1e6:.1f} MB", flush=True)

    # Clean up staging files so /tmp doesn't fill up across runs.
    for p in (STAGING_TSV, tmp_gz):
        try: os.remove(p)
        except OSError: pass


def cmd_all():
    print("Loading proportion matrix once ...", flush=True)
    prop_df = pd.read_csv(PROP_MATRIX_FILE, sep="\t", index_col=0)
    for i, p in enumerate(list_quant_files()):
        print(f"\n[{i+1}/6] {os.path.basename(p)}", flush=True)
        build_block(p, prop_df)
    write_cluster_index()
    cmd_publish()
    cmd_combine()




def cmd_publish(idx=None):
    """Copy staged block files from /tmp to the project folder.

    With no argument, copies every staged block. With an integer index,
    only copies the block(s) corresponding to that quantification file
    (including any "_part_*" slices).
    """
    if not os.path.isdir(TMP_DIR):
        raise SystemExit(f"no staging dir at {TMP_DIR}")
    if idx is None:
        names = sorted(os.listdir(TMP_DIR))
    else:
        files = list_quant_files()
        if idx < 0 or idx >= len(files):
            raise SystemExit(f"index {idx} out of range")
        # Match the un-suffixed name and any "_part_" slices for it.
        prefix = os.path.basename(block_out_path(files[idx]))[:-4]  # drop ".tsv"
        names = sorted([n for n in os.listdir(TMP_DIR)
                        if n == prefix + ".tsv" or n.startswith(prefix + "_part_")])
    if not names:
        raise SystemExit("nothing to publish")
    for n in names:
        src = os.path.join(TMP_DIR, n)
        dst = os.path.join(OUT_DIR, n)
        t0 = time.time()
        shutil.copyfile(src, dst)
        print(f"  cp {n}: {os.path.getsize(dst)/1e6:.1f} MB in {time.time()-t0:.1f}s",
              flush=True)
        try: os.remove(src)
        except OSError: pass


if __name__ == "__main__":
    arg = sys.argv[1] if len(sys.argv) > 1 else "all"
    if arg == "combine":
        cmd_combine()
    elif arg == "cluster_index":
        write_cluster_index()
    elif arg == "all":
        cmd_all()
    elif arg == "publish":
        # Optional: int index of which file's staged blocks to copy.
        idx = int(sys.argv[2]) if len(sys.argv) > 2 else None
        cmd_publish(idx=idx)
    else:
        try:
            idx = int(arg)
        except ValueError:
            raise SystemExit(f"unknown arg: {arg}")
        bin_start = int(sys.argv[2]) if len(sys.argv) > 2 else None
        bin_end   = int(sys.argv[3]) if len(sys.argv) > 3 else None
        cmd_block(idx, bin_start=bin_start, bin_end=bin_end)
