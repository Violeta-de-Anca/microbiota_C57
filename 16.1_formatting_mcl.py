#!/usr/bin/env python3
"""
16_parse_mcl.py — Parse MCL cluster output and add bin/generation/sample metadata.

Takes the human-readable MCL output (one cluster per line, protein IDs
tab-separated) and expands it into a long-format table with one row per
protein, annotated with:

    cluster_id   : integer cluster ID (1-based)
    protein_id   : full protein identifier
    bin          : bin name (protein_id minus trailing _NNNNN)
    generation   : F0, F1, or F2  (parsed from bin name)
    sample_type  : cecum_samples or last_feces (parsed from bin name)

Protein ID format expected:
    {generation}_{sample_type}_bin.{N}_{locus_prefix}_{protein_number}
    e.g.  F1_cecum_samples_bin.100_CPNENCLC_00003

Output columns (tab-separated):
    cluster_id  protein_id  bin  generation  sample_type

This table can be directly joined with your abundance count matrix on the
`bin` column to link functional clusters to per-sample bin abundances.

Usage:
    python 16_parse_mcl.py \
        --input  mcl/mcl_clusters.txt \
        --output mcl/functional_groups.txt
"""

import argparse
import os
import re
import sys


# ---------------------------------------------------------------------------
# Metadata extraction
# ---------------------------------------------------------------------------

# Matches: F0, F1, F2 at the very start of the protein/bin ID
RE_GENERATION = re.compile(r"^(F\d+)_")

# Matches known sample types
SAMPLE_TYPES = ("cecum_samples", "last_feces")


def extract_bin_name(protein_id: str) -> str:
    """Strip trailing _NNNNN protein number from protein ID."""
    idx = protein_id.rfind("_")
    if idx != -1 and protein_id[idx + 1:].isdigit():
        return protein_id[:idx]
    return protein_id


def extract_generation(bin_name: str) -> str:
    """Extract generation label (F0/F1/F2) from bin name."""
    m = RE_GENERATION.match(bin_name)
    return m.group(1) if m else "unknown"


def extract_sample_type(bin_name: str) -> str:
    """Extract sample type (cecum_samples / last_feces) from bin name."""
    for st in SAMPLE_TYPES:
        if st in bin_name:
            return st
    return "unknown"


# ---------------------------------------------------------------------------
# Main parser
# ---------------------------------------------------------------------------

def parse_mcl(input_file: str, output_file: str) -> None:

    n_clusters  = 0
    n_proteins  = 0
    n_singleton = 0

    print(f"[INFO] Parsing MCL output: {input_file}")

    with open(input_file, "r") as fin, open(output_file, "w") as fout:

        # Write header
        fout.write("cluster_id\tprotein_id\tbin\tgeneration\tsample_type\n")

        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue

            proteins = line.split("\t")
            n_clusters += 1
            cluster_id = n_clusters

            if len(proteins) == 1:
                n_singleton += 1

            for protein_id in proteins:
                protein_id = protein_id.strip()
                if not protein_id:
                    continue

                bin_name    = extract_bin_name(protein_id)
                generation  = extract_generation(bin_name)
                sample_type = extract_sample_type(bin_name)

                fout.write(
                    f"{cluster_id}\t{protein_id}\t{bin_name}\t"
                    f"{generation}\t{sample_type}\n"
                )
                n_proteins += 1

    print(f"[INFO] Done.")
    print(f"  Total clusters  : {n_clusters:>10,}")
    print(f"  Singleton clusters (1 protein): {n_singleton:>10,}")
    print(f"  Multi-protein clusters: {n_clusters - n_singleton:>10,}")
    print(f"  Total proteins  : {n_proteins:>10,}")
    print(f"  Output file     : {output_file}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Parse MCL clusters and add bin/generation/sample metadata.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input",  "-i", required=True,
                   help="MCL output file (mcxdump format: one cluster per line, "
                        "proteins tab-separated).")
    p.add_argument("--output", "-o", required=True,
                   help="Output long-format table (tab-separated).")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if not os.path.isfile(args.input):
        sys.exit(f"[ERROR] Input not found: {args.input}")
    parse_mcl(args.input, args.output)
