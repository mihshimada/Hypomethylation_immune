#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script processes 'ALL_cis.tsv' and performs the following steps:

1) Group CpG entries by gene.
2) Split entries into positive- and negative-effect groups based on Fx
   (gene_P for Fx >= 0, gene_N for Fx < 0).
3) Within each gene/sign group, rank CpGs by log10P
   (smaller log10P = stronger association = higher rank),
   using competition ranking (1,2,2,4,...).

Output:
    eQTM_Framingham_cis_1row_gene_rank_NP.txt

Each line has:
    GeneID \t joined CpG info \t #CpGs(in this sign group) \t #CpGs(total for the gene)

The second column is a colon-separated list:
    cg_id*chr*pos*rank : cg_id*chr*pos*rank : ...
"""

from collections import defaultdict, Counter

IN_PATH = "ALL_cis.tsv"
OUT_PATH = "eQTM_Framingham_cis_1row_gene_rank_NP.txt"


def competition_ranks(values):
    """
    Assign competition ranks to a list of numeric values.
    Smaller value -> better rank (1 is the best).
    Ties receive the same rank, and subsequent ranks are skipped
    (i.e. 1, 2, 2, 4, ...).
    """
    cnt = Counter(values)
    cum = 0
    rank_map = {}
    for v in sorted(cnt.keys()):
        rank_map[v] = cum + 1
        cum += cnt[v]
    return [rank_map[v] for v in values]


def main():
    # --------------------------------------------------
    # 1. Read input and collect CpG records per gene.
    #    key: gene (ProbesetID), value: list of items
    #    item = (cg_id, cg_chr, cg_pos, Fx, log10P)
    # --------------------------------------------------
    gene2items = defaultdict(list)

    with open(IN_PATH, "r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n")  # skip header
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")

            # Input columns (example):
            # Rs_ID, ProbesetID, RSq, Fx, T, log10P, ...
            cg_id = cols[0]
            gene_id = cols[1]
            fx = float(cols[3])
            log10p = float(cols[5])

            # Genomic position
            cg_chr = cols[9]    # CpG chromosome
            cg_pos = cols[10]   # CpG position

            gene2items[gene_id].append((cg_id, cg_chr, cg_pos, fx, log10p))

    # --------------------------------------------------
    # 2. For each gene, split CpGs into positive-Fx and negative-Fx
    #    groups and store them in a list for output.
    # --------------------------------------------------
    out_lines = []

    # Header line for the output file
    out_lines.append("GeneID\tCpGs(joined)\t#CpGs(sign)\t#CpGs(total)")

    for gene, items in gene2items.items():
        total_n = len(items)

        pos_items = []  # Fx >= 0
        neg_items = []  # Fx < 0
        for (cg_id, cg_chr, cg_pos, fx, log10p) in items:
            if fx >= 0:
                pos_items.append((cg_id, cg_chr, cg_pos, log10p))
            else:
                neg_items.append((cg_id, cg_chr, cg_pos, log10p))

        # Store gene-level entries; ranking will be added later
        if pos_items:
            out_lines.append((gene + "_P", pos_items, len(pos_items), total_n))
        if neg_items:
            out_lines.append((gene + "_N", neg_items, len(neg_items), total_n))

    # --------------------------------------------------
    # 3. For each stored entry, rank CpGs by log10P and write the result.
    # --------------------------------------------------
    with open(OUT_PATH, "w", encoding="utf-8") as w:
        # write header
        w.write(out_lines[0] + "\n")

        for entry in out_lines[1:]:
            gene_label, items, n_sign, n_total = entry

            # items: list of (cg_id, chr, pos, log10p)
            logps = [x[3] for x in items]
            ranks = competition_ranks(logps)

            joined_parts = []
            for (item, rk) in zip(items, ranks):
                cg_id, cg_chr, cg_pos, _ = item
                # output format: cg_id*chr*pos*rank
                joined_parts.append(f"{cg_id}*{cg_chr}*{cg_pos}*{rk}")

            joined_str = ":".join(joined_parts)
            w.write(f"{gene_label}\t{joined_str}\t{n_sign}\t{n_total}\n")


if __name__ == "__main__":
    main()
