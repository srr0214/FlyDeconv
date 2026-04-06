#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import pandas as pd
import json
import os
import psutil

from rnasieve.preprocessing import model_from_raw_counts


def main(args):
    print("==== RNA-Sieve start ====")

    # ---------- memory monitor ----------
    process = psutil.Process(os.getpid())
    peak_mem = 0.0

    def record_mem():
        nonlocal peak_mem
        mem = process.memory_info().rss / 1024 / 1024  # MB
        peak_mem = max(peak_mem, mem)

    record_mem()

    # ---------- 1. 读取输入 ----------
    sc_df = pd.read_csv(args.sc_counts, sep="\t", index_col=0)
    ann_df = pd.read_csv(args.sc_anno, sep="\t", index_col=0)
    bulk_df = pd.read_csv(args.bulk, sep="\t", index_col=0)

    record_mem()

    # ---------- 2. 基因对齐 ----------
    common_genes = sc_df.index.intersection(bulk_df.index)
    if len(common_genes) < 100:
        raise ValueError("Too few common genes (<100), aborting.")

    sc_df = sc_df.loc[common_genes]
    bulk_df = bulk_df.loc[common_genes]

    # ---------- 3. 解析 celltype 列 ----------
    if args.celltype_col not in ann_df.columns:
        celltype_col = ann_df.columns[0]
    else:
        celltype_col = args.celltype_col

    ann_df = ann_df.loc[sc_df.columns]
    cell_types = ann_df[celltype_col].astype(str)

    # ---------- 4. 构造 counts_by_type ----------
    counts_by_type = {}
    for ct in cell_types.unique():
        cells = cell_types[cell_types == ct].index
        counts_by_type[ct] = sc_df[cells].values.astype(np.float32)

    record_mem()

    # ---------- 5. bulk ----------
    psis = bulk_df.values.astype(np.float32)

    # ---------- 6. RNA-Sieve 核心 ----------
    model, cleaned_psis = model_from_raw_counts(
        counts_by_type,
        psis
    )
    record_mem()

    proportions = model.predict(cleaned_psis)
    record_mem()

    # ---------- 7. tensor → numpy ----------
    if hasattr(proportions, "detach"):
        proportions = proportions.detach().cpu().numpy()

    proportions = np.asarray(proportions)

    # ---------- 8. sanity check ----------
    n_sample = psis.shape[1]
    n_celltype = len(counts_by_type)

    if proportions.shape != (n_sample, n_celltype):
        raise ValueError("Shape mismatch in proportions")

    if np.isnan(proportions).any():
        raise ValueError("NA detected in proportions")

    # ---------- 9. 保存结果 ----------
    prop_df = pd.DataFrame(
        proportions,
        index=bulk_df.columns,
        columns=list(counts_by_type.keys())
    )

    out_file = args.out_prefix + "_cell_fraction.tsv"
    prop_df.to_csv(out_file, sep="\t", float_format="%.6f")

    # ---------- 10. 输出 benchmark ----------
    bench_file = args.out_prefix + "_benchmark.json"
    with open(bench_file, "w") as f:
        json.dump(
            {"core_mem_MB": round(peak_mem, 2)},
            f,
            indent=2
        )

    print("Saved result to:", out_file)
    print("Core peak memory (MB):", round(peak_mem, 2))
    print("==== RNA-Sieve finished successfully ====")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run RNA-Sieve using matrix inputs"
    )

    parser.add_argument("--sc_counts", required=True)
    parser.add_argument("--sc_anno", required=True)
    parser.add_argument("--bulk", required=True)
    parser.add_argument("--celltype_col", default="cellType")
    parser.add_argument("--out_prefix", required=True)

    args = parser.parse_args()
    main(args)
