#!/usr/bin/env python3
"""
Goal
----
Build a wide matrix where:

    rows    = clusters from cluster_bin_protein_proportion_matrix.tsv
              (cluster_1 ... cluster_N)
    columns = every (bin, individual) pair — one column for each individual
              that has that bin in their quantification file
    cell    = abundance(bin, individual) * proportion(cluster, bin)

NO summation across bins.

Inputs
------
  1) functional_annotatioin/cluster_bin_protein_proportion_matrix.tsv
        rows  = clusters, cols = fully qualified bins (F0_cecum_samples_bin.1, ...)
        cells = fraction of that bin's proteins that landed in this cluster.

  2) quantification/bin_abundance_table_<gen>_<cecum|last_feces>_*.tab
        rows  = bins, written as "bin.<NNN>" (no gen / sample-type prefix)
        cols  = per-individual abundance values
        We prepend "<generation>_<sample_type>_" so bin names line up with
        the proportion-matrix columns. <sample_type> is "cecum_samples"
        when the filename says "cecum" (the cecum files lack the trailing
        "_samples" used in the proteome FASTA), and "last_feces" otherwise.

Output (next to the proportion matrix in functional_annotatioin/)
----------------------------------------------------------------
  - cluster_individual_per_bin_weighted_abundance_matrix.tsv.gz
        Wide gzipped TSV. First column "cluster", then one column per
        (bin, individual) pair with header "<qualified_bin>__<individual>"
        (double underscore separator).
  - intermediate per-file blocks (uncompressed TSVs):
        cluster_individual_per_bin_block_<gen>_<sampletype>[_part_NNNN_MMMM].tsv
        cluster_individual_per_bin_block_cluster_index.tsv

Why we split the work into per-file blocks
------------------------------------------
The sandbox shell budget per call is ~45s and writing to the OneDrive
mount runs at ~9 MB/s, so a 1 GB TSV cannot land in a single call. We:
  1) build one block per quantification file (and split very large blocks
     by bin-range when they wouldn't fit in 45s),
  2) write the cluster-name column once as its own file,
  3) `paste` everything side-by-side and gzip the result.

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

# Paths (Linux mount paths — script runs in the sandbox).
ROOT = "/proj/naiss2024-23-57/C57_female_lineage_microbiota"
PROP_MATRIX_FILE   = f"{ROOT}/functional_annotation/cluster_bin_protein_proportion_matrix.tsv"
QUANT_DIR          = f"{ROOT}/quantification"
OUT_DIR            = f"{ROOT}/quantification"
FINAL_OUT_GZ       = f"{OUT_DIR}/cluster_individual_per_bin_weighted_abundance_matrix.tsv.gz"
CLUSTER_INDEX_FILE = f"{OUT_DIR}/relational_table_abundance_functionality_multigen_C57_maternal.txt"
# Pasted-but-uncompressed staging file under /tmp; OneDrive is too slow
# for multi-GB writes within a 45s call.
STAGING_TSV        = "/tmp/cluster_individual_per_bin_weighted_abundance_matrix.tsv"

# Filename pattern for the abundance tables.
FNAME_PATTERN = re.compile(r"^bin_abundance_table_(F\d)_(cecum|last_feces)_.*\.tab$")
# Cecum filenames omit "_samples"; the proportion matrix uses "cecum_samples".
SAMPLE_TYPE_MAP = {"cecum": "cecum_samples", "last_feces": "last_feces"}


def list_quant_files():
    """Discover the abundance tables in deterministic (sorted) order."""
    return sorted(glob.glob(os.path.join(QUANT_DIR, "bin_abundance_table_*.tab")))


def block_out_path(path: str) -> str:
    """Per-file block output path (no bin-slice suffix)."""
    fname = os.path.basename(path)
    m = FNAME_PATTERN.match(fname)
    return os.path.join(
        OUT_DIR,
        f"cluster_individual_per_bin_block_{m.group(1)}_{SAMPLE_TYPE_MAP[m.group(2)]}.tsv",
    )


def build_block(path, prop_df, bin_start=None, bin_end=None):
    """Write one (file, optional bin-slice) block as a header + data TSV.

    Cell: out[c, (b, i)] = proportion[c, b] * abundance[b, i].
    The cluster-name column is NOT in this file — it's written separately
    by write_cluster_index() and pasted in front at the combine step.
    """
    t0 = time.time()
    fname = os.path.basename(path)
    m = FNAME_PATTERN.match(fname)
    if not m:
        print(f"  SKIP (bad filename): {fname}", flush=True)
        return
    bin_prefix = f"{m.group(1)}_{SAMPLE_TYPE_MAP[m.group(2)]}_"

    # Load the abundance table (bins x individuals); coerce to floats.
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip()
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # Qualify bin names so they line up with proportion-matrix columns.
    bins_qualified = [bin_prefix + b for b in df.index]
    indivs = list(df.columns)

    # Drop bins not present in proportion matrix (defensive).
    keep = [i for i, b in enumerate(bins_qualified) if b in prop_df.columns]
    if len(keep) != len(bins_qualified):
        bins_qualified = [bins_qualified[i] for i in keep]
        df = df.iloc[keep]

    # Optional bin slice — used to split very large blocks (F2_cecum) so
    # each piece fits in a 45s shell call.
    is_partial = (bin_start is not None) or (bin_end is not None)
    if bin_start is None: bin_start = 0
    if bin_end   is None: bin_end   = len(bins_qualified)
    bins_qualified = bins_qualified[bin_start:bin_end]
    df = df.iloc[bin_start:bin_end]

    A = df.to_numpy(dtype=np.float64)                       # (n_bins, n_indivs)
    P = prop_df[bins_qualified].to_numpy(dtype=np.float64)  # (n_clusters, n_bins)

    n_clusters = P.shape[0]
    n_bins = len(bins_qualified)
    n_indivs = len(indivs)

    # Output column names match the array's flatten layout: for each bin,
    # all individuals (so column k = bin k//n_indivs, individual k%n_indivs).
    col_names = [f"{b}__{ind}" for b in bins_qualified for ind in indivs]

    out_path = block_out_path(path)
    if is_partial:
        out_path = out_path[:-4] + f"_part_{bin_start:04d}_{bin_end:04d}.tsv"
    print(f"  -> {out_path}", flush=True)
    print(f"     {n_bins} bins x {n_indivs} indivs = {len(col_names)} cols", flush=True)

    # Build dense block bin-by-bin to bound peak memory (avoid 3D tensor).
    out_arr = np.empty((n_clusters, n_bins * n_indivs), dtype=np.float64)
    for j in range(n_bins):
        out_arr[:, j*n_indivs:(j+1)*n_indivs] = P[:, j:j+1] * A[j:j+1, :]

    # numpy.savetxt is C-level; fmt="%g" gives short reps ("0" for zeros).
    header_line = "\t".join(col_names) + "\n"
    with open(out_path, "wb") as out:
        out.write(header_line.encode())
        np.savetxt(out, out_arr, fmt="%g", delimiter="\t")
    print(f"     wrote {os.path.getsize(out_path)/1e6:.1f} MB in {time.time()-t0:.1f}s",
          flush=True)


def write_cluster_index():
    """Write the standalone single-column file holding cluster names."""
    prop_df = pd.read_csv(PROP_MATRIX_FILE, sep="\t", index_col=0, usecols=[0])
    with open(CLUSTER_INDEX_FILE, "w") as out:
        out.write("cluster\n")
        for c in prop_df.index:
            out.write(f"{c}\n")
    print(f"  wrote {CLUSTER_INDEX_FILE} ({len(prop_df.index)} clusters)", flush=True)


def cmd_block(idx, bin_start=None, bin_end=None):
    """Build one block by sorted index (with optional bin slice)."""
    print("Loading proportion matrix ...", flush=True)
    prop_df = pd.read_csv(PROP_MATRIX_FILE, sep="\t", index_col=0)
    files = list_quant_files()
    if idx < 0 or idx >= len(files):
        raise SystemExit(f"index {idx} out of range 0..{len(files)-1}")
    print(f"Building block {idx}: {os.path.basename(files[idx])}", flush=True)
    build_block(files[idx], prop_df, bin_start=bin_start, bin_end=bin_end)


def cmd_combine():
    """Paste cluster-index + every block (or bin-slice) TSV side-by-side
    into a staging TSV in /tmp, then gzip and copy to the project folder.

    Why /tmp first: OneDrive write speed is ~9 MB/s — a 1 GB TSV cannot
    land in a 45s call. /tmp is local (~1 GB/s) so we paste there, gzip
    (compresses ~10x), and copy the smaller .gz to OneDrive in one call.
    """
    all_blocks = sorted(glob.glob(os.path.join(
        OUT_DIR, "cluster_individual_per_bin_block_*.tsv")))
    all_blocks = [p for p in all_blocks if p != CLUSTER_INDEX_FILE]

    # If a block has any "_part_<x>_<y>" slices, drop the un-suffixed full
    # file (it would be missing or incomplete; the slices are authoritative).
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

    # gzip into /tmp, then copy the (much smaller) .gz to the project folder.
    tmp_gz = STAGING_TSV + ".gz"
    print(f"gzipping -> {tmp_gz}", flush=True)
    with open(STAGING_TSV, "rb") as src, gzip.open(tmp_gz, "wb", compresslevel=1) as dst:
        shutil.copyfileobj(src, dst, length=8 * 1024 * 1024)
    print(f"  gz size: {os.path.getsize(tmp_gz)/1e6:.1f} MB", flush=True)

    print(f"copying {tmp_gz} -> {FINAL_OUT_GZ}", flush=True)
    shutil.copyfile(tmp_gz, FINAL_OUT_GZ)
    print(f"Final file: {os.path.getsize(FINAL_OUT_GZ)/1e6:.1f} MB", flush=True)


def cmd_all():
    """Build every block, write cluster-index, and combine.
    NOTE: typically does NOT fit in a single 45s shell call — use the
    per-block CLI mode for long runs.
    """
    print("Loading proportion matrix once ...", flush=True)
    prop_df = pd.read_csv(PROP_MATRIX_FILE, sep="\t", index_col=0)
    for i, p in enumerate(list_quant_files()):
        print(f"\n[{i+1}/6] {os.path.basename(p)}", flush=True)
        build_block(p, prop_df)
    write_cluster_index()
    cmd_combine()


if __name__ == "__main__":
    arg = sys.argv[1] if len(sys.argv) > 1 else "all"
    if arg == "combine":
        cmd_combine()
    elif arg == "cluster_index":
        write_cluster_index()
    elif arg == "all":
        cmd_all()
    else:
        # Single-block mode. Optional bin slice: <idx> <bstart> <bend>.
        try:
            idx = int(arg)
        except ValueError:
            raise SystemExit(f"unknown arg: {arg}")
        bin_start = int(sys.argv[2]) if len(sys.argv) > 2 else None
        bin_end   = int(sys.argv[3]) if len(sys.argv) > 3 else None
        cmd_block(idx, bin_start=bin_start, bin_end=bin_end)
