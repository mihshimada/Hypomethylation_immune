#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Input:
    eQTM_Framingham_cis_1row_gene_rank_NP.txt
    (one line per gene-sign: GeneID  CpGs(joined)  #CpGs(sign)  #CpGs(total))
    CpGs(joined) looks like:
        cg00000001*chr1*12345*1:cg00000002*chr1*54321*3:...

Goal:
    Reorganize the data by CpG ID (probe).
    For each CpG, list all gene-level entries in which it appears,
    together with:
        - the gene label (e.g. ENSG00000123456_P),
        - its rank among CpGs for that gene,
        - the number of CpGs in that sign group,
        - the total number of CpGs for the gene.

Output:
    eQTM_Framingham_cis_1row_gene_rank_NP_probe_annotation.txt
    with two columns:
        CpG_ID   gene1*rank/num_in_sign/num_total : gene2*rank/...

"""

from collections import OrderedDict

IN_PATH = "eQTM_Framingham_cis_1row_gene_rank_NP.txt"
OUT_PATH = "eQTM_Framingham_cis_1row_gene_rank_NP_probe_annotation.txt"


def main():
    # CpG_ID -> list of "gene*rank/n_sign/n_total"
    probe_map = OrderedDict()

    with open(IN_PATH, "r", encoding="utf-8") as fin:
        header = fin.readline()  # skip header line

        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t")
            # Expecting at least: GeneID, CpGs(joined), #CpGs(sign), #CpGs(total)
            if len(parts) < 4:
                # skip malformed lines
                continue

            gene_label = parts[0]          # e.g. ENSG00000123456_P
            joined_cpgs = parts[1]         # e.g. cg*chr*pos*rank:cg*...
            num_in_sign = parts[2]         # e.g. "5"
            num_total = parts[3]           # e.g. "10"

            if not joined_cpgs:
                continue

            for cpg_item in joined_cpgs.split(":"):
                # cpg_item: cg_id*chr*pos*rank
                fields = cpg_item.split("*")
                if len(fields) < 4:
                    # skip unexpected format
                    continue
                cpg_id = fields[0]
                rank = fields[3]

                annotation = f"{gene_label}*{rank}/{num_in_sign}/{num_total}"

                if cpg_id not in probe_map:
                    probe_map[cpg_id] = [annotation]
                else:
                    probe_map[cpg_id].append(annotation)

    # write out
    with open(OUT_PATH, "w", encoding="utf-8") as fout:
        for cpg_id, ann_list in probe_map.items():
            joined_ann = ":".join(ann_list)
            fout.write(f"{cpg_id}\t{joined_ann}\n")


if __name__ == "__main__":
    main()
