#!/usr/bin/env python3

import sys
from pathlib import Path

import numpy as np
import pandas as pd


def clonality_from_counts(counts):
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()

    if total <= 0 or len(counts) == 0:
        return {
            "n_clones": 0,
            "top1": np.nan,
            "top10": np.nan,
            "clonality": np.nan,
        }

    p = counts / total
    entropy = -(p * np.log(p)).sum()
    normalized_entropy = entropy / np.log(len(p)) if len(p) > 1 else 0.0

    return {
        "n_clones": int(len(p)),
        "top1": float(p.max()),
        "top10": float(np.sort(p)[::-1][:10].sum()),
        "clonality": float(1 - normalized_entropy),
    }


def filter_productive_locus(df, locus):
    t = df[df["locus"].astype(str) == locus].copy()
    t = t[t["productive"].astype(str).str.lower().isin(["t", "true", "1"])]
    return t


def summarize_clonotypes(df, clonotype_col):
    required_cols = [clonotype_col, "consensus_count"]
    if any(col not in df.columns for col in required_cols):
        return {
            "n_clones": 0,
            "top1": np.nan,
            "top10": np.nan,
            "clonality": np.nan,
        }

    t = df.dropna(subset=[clonotype_col]).copy()
    t["consensus_count"] = pd.to_numeric(t["consensus_count"], errors="coerce").fillna(0)

    counts = t.groupby(clonotype_col)["consensus_count"].sum()
    return clonality_from_counts(counts.values)


def read_sample_list(sample_list_file):
    return [
        line.strip().replace("\r", "")
        for line in sample_list_file.read_text().splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]


def main():
    if len(sys.argv) != 3:
        print("Usage: python 03_trust4_clonality.py <trust4_result_dir> <sample_list.txt>")
        sys.exit(1)

    workdir = Path(sys.argv[1])
    sample_list_file = Path(sys.argv[2])

    airr_suffix = "_TRUST4_airr.tsv"
    out_csv = workdir / "trust4_clonality_TRA_TRB.csv"

    samples = read_sample_list(sample_list_file)
    rows = []

    for sample in samples:
        airr_file = workdir / f"{sample}{airr_suffix}"

        if not airr_file.exists():
            rows.append({"sample": sample, "status": "missing_airr"})
            continue

        try:
            df = pd.read_csv(airr_file, sep="\t")
        except Exception as e:
            rows.append({"sample": sample, "status": f"read_error:{type(e).__name__}"})
            continue

        row = {"sample": sample, "status": "ok"}

        for locus in ["TRA", "TRB"]:
            t = filter_productive_locus(df, locus)

            aa_metrics = summarize_clonotypes(t, "junction_aa")
            nt_metrics = summarize_clonotypes(t, "junction")

            for key, value in aa_metrics.items():
                row[f"{locus}_AA_{key}"] = value

            for key, value in nt_metrics.items():
                row[f"{locus}_NT_{key}"] = value

        rows.append(row)

    out = pd.DataFrame(rows)
    out.to_csv(out_csv, index=False)
    print(f"wrote {out_csv}")


if __name__ == "__main__":
    main()
