#!/usr/bin/env python3
"""
15_hssp_chunk.py — Compute HSSP distances for ONE PSI-BLAST chunk file.

Designed to run as a SLURM array job (one job per original query chunk).

Why this is correct without a global merge step:
    Queries were split into 600 non-overlapping chunks before PSI-BLAST.
    Therefore the pair (qseqid, sseqid) can only appear in ONE chunk file.
    The only deduplication needed is within a single chunk (multiple
    PSI-BLAST iterations may produce the same pair more than once).
    After all 600 jobs finish, results are simply concatenated.

HSSP formula (Zhu et al. 2015, Eq. 1):
    L < 11:          HSSP = -99
    11 < L <= 450:   HSSP = Id - 480 * L^(-0.32*(1+exp(-L/1000)))
    L > 450:         HSSP = Id - 19.5
    Threshold >= 10 → functionally similar proteins

Usage (called by 15_submit_hssp_array.sh):
    python 15_hssp_chunk.py --input chunk_0001.txt --output hssp_0001.txt
"""

import pandas as pd
import numpy as np
import argparse
import sys
import os


# ---------------------------------------------------------------------------
# HSSP formula (fully vectorised)
# ---------------------------------------------------------------------------

def hssp_vectorized(L: np.ndarray, Id: np.ndarray) -> np.ndarray:
    L  = np.asarray(L,  dtype=np.float64)
    Id = np.asarray(Id, dtype=np.float64)
    return np.where(
        L < 11,
        -99.0,
        np.where(
            L <= 450,
            Id - 480.0 * np.power(L, -0.32 * (1.0 + np.exp(-L / 1000.0))),
            Id - 19.5,
        ),
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def process_chunk(input_file: str, output_file: str, threshold: float = 10.0) -> None:

    col_names = [
        "qseqid", "sseqid", "pident", "length",
        "mismatch", "gapopen", "qstart", "qend",
        "sstart", "send", "evalue", "bitscore", "nident",
    ]

    # Read only the 4 columns we actually need
    # on_bad_lines='skip' silently drops lines that do not have exactly 13
    # tab-separated fields, e.g. "Search has CONVERGED!" lines still present
    # in the individual (non-merged) chunk files.
    # length is read as float32 (not int32) so that any remaining NA values
    # are represented as NaN rather than raising a ValueError.
    df = pd.read_csv(
        input_file,
        sep="\t",
        header=None,
        names=col_names,
        usecols=["qseqid", "sseqid", "pident", "length"],
        dtype={"qseqid": "string", "sseqid": "string",
               "pident": np.float32, "length": np.float32},
        engine="c",
        on_bad_lines="skip",
    )

    # Drop any rows where pident or length could not be parsed (secondary guard)
    before = len(df)
    df.dropna(subset=["pident", "length"], inplace=True)
    dropped = before - len(df)
    if dropped > 0:
        print(f"[WARN] Dropped {dropped:,} malformed rows in {input_file}",
              file=sys.stderr)

    if df.empty:
        print(f"[WARN] Empty input: {input_file}", file=sys.stderr)
        # Write empty file with header so merge step doesn't fail
        with open(output_file, "w") as fh:
            fh.write("qseqid\tsseqid\tqbin\tsbin\thssp_distance\n")
        return

    # Compute HSSP distance (vectorised, no Python loop)
    df["hssp"] = hssp_vectorized(df["length"].values, df["pident"].values)

    # Keep max HSSP per (query, subject) pair — handles multiple PSI-BLAST iterations
    # sort=False preserves insertion order → faster groupby
    result = (
        df.groupby(["qseqid", "sseqid"], sort=False)["hssp"]
        .max()
        .reset_index()
    )

    # Apply functional similarity threshold
    result = result[result["hssp"] >= threshold]

    # Extract bin name: strip trailing protein number (_NNNNN)
    # e.g. F1_cecum_samples_bin.100_CPNENCLC_00003 → F1_cecum_samples_bin.100_CPNENCLC
    result["qbin"] = result["qseqid"].str.replace(r"_\d+$", "", regex=True)
    result["sbin"] = result["sseqid"].str.replace(r"_\d+$", "", regex=True)

    # Write output (no index, 4-decimal HSSP)
    result[["qseqid", "sseqid", "qbin", "sbin", "hssp"]].to_csv(
        output_file,
        sep="\t",
        index=False,
        header=True,
        float_format="%.4f",
    )

    print(
        f"[OK] {os.path.basename(input_file)} → "
        f"{len(result):,} pairs (HSSP >= {threshold}) written to "
        f"{os.path.basename(output_file)}"
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Compute HSSP distances for one PSI-BLAST chunk file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input",     "-i", required=True,
                   help="PSI-BLAST chunk result file (13-column tab-separated).")
    p.add_argument("--output",    "-o", required=True,
                   help="Output HSSP file for this chunk (tab-separated).")
    p.add_argument("--threshold", "-t", type=float, default=10.0,
                   help="Minimum HSSP distance to retain a pair.")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if not os.path.isfile(args.input):
        sys.exit(f"[ERROR] Input not found: {args.input}")
    process_chunk(args.input, args.output, args.threshold)
