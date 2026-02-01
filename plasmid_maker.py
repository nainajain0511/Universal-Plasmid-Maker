# =========================================================
# Universal Plasmid Maker
# ORI detection using GC Skew + k-mer enrichment
# =========================================================

import sys
import os
from collections import defaultdict


# ---------------- FASTA READER ----------------
def read_fasta(file_path):
    sequence = []
    if not os.path.exists(file_path):
        sys.exit(f"ERROR: File not found -> {file_path}")

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('>'):
                sequence.append(line.upper())

    if not sequence:
        sys.exit(f"ERROR: Empty FASTA -> {file_path}")

    return "".join(sequence)


# ---------------- ORI DETECTION (GC SKEW + k-MER) ----------------
def find_ori_gc_kmer(sequence, k=8, window_size=5000, step=500):
    """
    ORI detection using:
    1) Cumulative GC skew minimum
    2) k-mer enrichment peak
    """
    n = len(sequence)

    cumulative_skew = 0
    min_skew = float('inf')
    ori_coord = 0

    max_kmer_count = 0
    enriched_coord = 0

    for i in range(0, n - window_size + 1, step):
        window = sequence[i : i + window_size]
        mid_point = i + window_size // 2

        # ---- GC skew ----
        g = window.count('G')
        c = window.count('C')
        skew = (g - c) / (g + c) if (g + c) > 0 else 0
        cumulative_skew += skew

        if cumulative_skew < min_skew:
            min_skew = cumulative_skew
            ori_coord = mid_point

        # ---- k-mer enrichment ----
        counts = defaultdict(int)
        for j in range(len(window) - k + 1):
            kmer = window[j : j + k]
            counts[kmer] += 1

        if counts:
            local_max = max(counts.values())
            if local_max > max_kmer_count:
                max_kmer_count = local_max
                enriched_coord = mid_point

    # Final ORI = average of both signals
    final_ori_coord = (ori_coord + enriched_coord) // 2

    # Extract ORI sequence (~500 bp)
    start = max(0, final_ori_coord - 250)
    end = min(n, final_ori_coord + 250)

    return sequence[start:end]


# ---------------- DESIGN FILE PARSER ----------------
def parse_design(design_file):
    """
    Design file lists ALLOWED restriction sites
    """
    allowed_enzymes = set()

    with open(design_file, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            name, value = map(str.strip, line.split(','))
            if "site" in name.lower():
                allowed_enzymes.add(value)

    return allowed_enzymes


# ---------------- MARKER DICTIONARY ----------------
def load_markers(marker_file):
    marker_dict = {}

    with open(marker_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                marker_dict[parts[0]] = parts[1]

    return marker_dict


# ---------------- RESTRICTION SITES ----------------
RESTRICTION_SITES = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "PstI": "CTGCAG",
    "SphI": "GCATGC",
    "SalI": "GTCGAC",
    "XbaI": "TCTAGA",
    "KpnI": "GGTACC",
    "SacI": "GAGCTC",
    "SmaI": "CCCGGG"
}


# ---------------- REMOVE DISALLOWED SITES ----------------
def remove_disallowed_sites(sequence, allowed_enzymes):
    for enzyme, site in RESTRICTION_SITES.items():
        if enzyme not in allowed_enzymes:
            sequence = sequence.replace(site, "")
    return sequence


# ---------------- FASTA WRITER ----------------
def write_fasta(path, header, sequence):
    with open(path, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i+70] + "\n")


# ---------------- MAIN ----------------
def main():
    organism_fasta = "pUC19.fa"
    design_file = "Design_pUC19.txt"
    marker_file = "markers.tab"
    output_file = "Output.Fa"

    print("Reading organism genome...")
    genome = read_fasta(organism_fasta)

    print("Detecting ORI using GC skew + k-mer enrichment...")
    ori_seq = find_ori_gc_kmer(genome)

    print("Parsing design file...")
    allowed_enzymes = parse_design(design_file)

    print("Loading marker dictionary...")
    markers = load_markers(marker_file)

    print("Removing disallowed restriction sites...")
    plasmid_core = remove_disallowed_sites(genome, allowed_enzymes)

    print("Assembling plasmid...")
    plasmid = ori_seq + plasmid_core

    print("Writing output FASTA...")
    write_fasta(output_file, "Universal_Plasmid | circular", plasmid)

    print("\nSUCCESS")
    print(f"Final plasmid length: {len(plasmid)} bp")

    if "GAATTC" in plasmid:
        print("ERROR: EcoRI site still present")
    else:
        print("EcoRI successfully removed")


if __name__ == "__main__":
    main()
