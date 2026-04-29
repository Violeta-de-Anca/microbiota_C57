#!/usr/bin/env python3
"""
build_cluster_bin_matrix.py
============================

Goal
----
Build a matrix where:
  - Rows    = MCL clusters (one per line in mcl_clusters.txt).
  - Columns = every unique bin across all generations / sample types.
  - Cell    = (number of proteins from that bin found in that cluster)
              / (total number of proteins that bin has in all_proteins_multigen.faa)
            -> a fraction in [0, 1] telling us what proportion of the bin's
               proteome ended up in this cluster.

Inputs
------
  1) mcl_clusters.txt
       - Tab-separated, one cluster per line.
       - Each token is a protein ID, e.g. F1_cecum_samples_bin.100_CPNENCLC_00967

  2) all_proteins_multigen.faa
       - FASTA file with all proteins from all bins / generations.
       - Headers look like:
            >F1_cecum_samples_bin.100_CPNENCLC_00003 Putative NAD(P)H ...
       - The bin identifier is the part of the protein ID *before* the last
         two '_'-separated fields:
              F1_cecum_samples_bin.100_CPNENCLC_00003
              \_____________ bin _____________/  \locus/  \prot#/
         => bin = "F1_cecum_samples_bin.100"

Outputs (saved next to the FASTA, in functional_annotatioin/)
-------------------------------------------------------------
  - bin_total_protein_counts.tsv
       2-column TSV: bin <tab> total_proteins  (one row per unique bin)

  - cluster_bin_protein_proportion_matrix.tsv
       Wide matrix: first column "cluster", then one column per bin (sorted).
       Each cell is the fraction described above (0 if the bin has no
       proteins in that cluster).
"""

import os
from collections import Counter

# ---------------------------------------------------------------------------
# File paths (Linux-side mount paths, since this script runs in the sandbox).
# ---------------------------------------------------------------------------
MCL_FILE   = "/proj/naiss2024-23-57/C57_female_lineage_microbiota/functional_annotation/mcl/mcl_clusters.txt"
FAA_FILE   = "/proj/naiss2024-23-57/C57_female_lineage_microbiota/functional_annotation/all_proteins_multigen.faa"
OUT_DIR    = "/proj/naiss2024-23-57/C57_female_lineage_microbiota/functional_annotation"
OUT_MATRIX = os.path.join(OUT_DIR, "cluster_bin_protein_proportion_matrix.tsv")
OUT_TOTALS = os.path.join(OUT_DIR, "bin_total_protein_counts.tsv")


# ---------------------------------------------------------------------------
# Helper: turn a protein ID into its bin ID.
#
# Protein ID structure:
#   <generation>_<sample_type>_bin.<NNN>_<LOCUSTAG>_<PROTNUM>
# We always strip the last two underscore-separated fields (LOCUSTAG and
# PROTNUM) to recover the bin identifier. This works whether sample_type is
# "cecum_samples" or "last_feces" because both contribute one underscore each
# and the bin token "bin.<NNN>" itself contains a dot, not an underscore.
# ---------------------------------------------------------------------------
def protein_to_bin(protein_id: str) -> str:
    parts = protein_id.split("_")
    # Defensive: if the ID is shorter than expected just return as-is.
    if len(parts) < 3:
        return protein_id
    # Drop last 2 fields (locus tag + protein number) -> rejoin with "_".
    return "_".join(parts[:-2])


# ===========================================================================
# STEP 1 — Count total proteins per bin from all_proteins_multigen.faa
# ===========================================================================
# We scan the FASTA once and look only at header lines (those starting with ">").
# For every header we extract the protein ID (first whitespace-delimited token,
# minus the leading '>'), convert it to its bin, and increment that bin's
# counter. The result is a dict: bin -> total #proteins of that bin.
print("Step 1/3: counting total proteins per bin in all_proteins_multigen.faa ...",
      flush=True)

bin_totals: Counter = Counter()      # bin_id -> total #proteins for that bin
seen_protein_ids = set()              # guard against duplicate FASTA headers
n_headers = 0                         # raw header line count (sanity check)

with open(FAA_FILE, "r") as fh:
    for line in fh:
        # We only care about FASTA header lines.
        if not line.startswith(">"):
            continue
        n_headers += 1

        # Header may look like:  ">F1_cecum_samples_bin.100_CPNENCLC_00003 description..."
        # We take only the first token (the protein ID) and drop the leading ">".
        prot_id = line[1:].split(None, 1)[0]

        # Skip exact duplicates so a protein listed twice in the FASTA is not
        # double-counted. (sort | uniq | wc behaviour from the task description.)
        if prot_id in seen_protein_ids:
            continue
        seen_protein_ids.add(prot_id)

        # Map protein -> bin and increment that bin's total.
        bin_id = protein_to_bin(prot_id)
        bin_totals[bin_id] += 1

print(f"  parsed {n_headers:,} FASTA headers; "
      f"{len(seen_protein_ids):,} unique protein IDs; "
      f"{len(bin_totals):,} unique bins.", flush=True)

# Persist the per-bin totals to disk as a 2-column TSV.
# Useful as a sanity check, and as a standalone lookup file.
with open(OUT_TOTALS, "w") as out:
    out.write("bin\ttotal_proteins\n")
    for b in sorted(bin_totals):
        out.write(f"{b}\t{bin_totals[b]}\n")
print(f"  wrote {OUT_TOTALS}", flush=True)


# ===========================================================================
# STEP 2 — Count, per cluster, how many proteins of each bin are present
# ===========================================================================
# mcl_clusters.txt has one cluster per line; tokens are tab-separated protein
# IDs. For every line we build a Counter mapping bin -> #proteins-of-that-bin
# present in this cluster. We keep one Counter per cluster, in the same order
# as the file (so cluster_counts[0] corresponds to the first line of the file).
print("Step 2/3: counting proteins per bin in each cluster ...", flush=True)

cluster_counts: list = []   # list of Counter objects, one per cluster line
n_clusters = 0

with open(MCL_FILE, "r") as fh:
    for line in fh:
        line = line.rstrip("\n")

        # An empty line is an empty cluster; preserve order with an empty Counter.
        if not line:
            cluster_counts.append(Counter())
            n_clusters += 1
            continue

        # Tokens are tab-separated protein IDs.
        proteins = line.split("\t")

        # For this cluster, count proteins per bin.
        c: Counter = Counter()
        for p in proteins:
            if not p:
                continue
            c[protein_to_bin(p)] += 1

        cluster_counts.append(c)
        n_clusters += 1

print(f"  parsed {n_clusters:,} clusters.", flush=True)

# Sanity check: every bin that shows up in the clusters file should also have
# been seen in the FASTA. If any are missing, something is off (e.g. naming
# mismatch). We just warn rather than fail so the run still produces output.
bins_in_clusters = set()
for c in cluster_counts:
    bins_in_clusters.update(c.keys())

missing = bins_in_clusters - set(bin_totals.keys())
if missing:
    print(f"  WARNING: {len(missing)} bin(s) in clusters not found in FASTA totals "
          f"(first few: {list(missing)[:5]})", flush=True)


# ===========================================================================
# STEP 3 — Build & write the cluster x bin proportion matrix
# ===========================================================================
# Column order: every unique bin from the FASTA, sorted alphabetically (so the
# output is reproducible and bins from the same generation/sample-type cluster
# together visually).
# Row order: cluster_1, cluster_2, ... matching the input file's line order.
# Cell value: (count in this cluster) / (bin's total in FASTA).
# When the bin has zero proteins in the cluster we write a literal "0" to keep
# the file small and easy to parse.
print("Step 3/3: writing cluster x bin proportion matrix ...", flush=True)

all_bins = sorted(bin_totals.keys())
header = "cluster\t" + "\t".join(all_bins) + "\n"

with open(OUT_MATRIX, "w") as out:
    out.write(header)

    # Walk through clusters in their original order and compute one row each.
    for i, c in enumerate(cluster_counts, start=1):
        row_vals = []
        for b in all_bins:
            n = c.get(b, 0)
            if n == 0:
                # Most cells will be 0 — keep it short.
                row_vals.append("0")
            else:
                tot = bin_totals.get(b, 0)
                if tot == 0:
                    # Should not happen given the sanity check above, but be safe.
                    row_vals.append("NA")
                else:
                    # Format with up to 6 significant digits to keep file size sane.
                    row_vals.append(f"{n / tot:.6g}")
        out.write(f"cluster_{i}\t" + "\t".join(row_vals) + "\n")

print(f"  wrote {OUT_MATRIX}", flush=True)
print("Done.", flush=True)
