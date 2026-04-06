#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import numpy as np
import pandas as pd
import anndata
import autogenes as ag
import time
import json
import psutil


def main(args):

    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    prefix = args.prefix
    if prefix and not prefix.endswith("_"):
        prefix = prefix + "_"

    # ===============================
    # core memory profiler
    # ===============================
    process = psutil.Process(os.getpid())
    core_mem_peak_mb = 0.0

    def update_mem():
        nonlocal core_mem_peak_mb
        mem_mb = process.memory_info().rss / 1024 / 1024
        if mem_mb > core_mem_peak_mb:
            core_mem_peak_mb = mem_mb

    # ===============================
    # 1. 读取数据（非 core）
    # ===============================
    sc   = pd.read_csv(args.sc_counts, sep="\t", index_col=0)   # gene × cell
    ann  = pd.read_csv(args.sc_anno,   sep="\t", index_col=0)  # cell × meta
    bulk = pd.read_csv(args.bulk,      sep="\t", index_col=0)  # gene × sample

    # ===============================
    # 2. 对齐基因（非 core）
    # ===============================
    common_genes = sc.index.intersection(bulk.index)
    if len(common_genes) == 0:
        raise ValueError("No common genes between scRNA and bulk.")

    sc   = sc.loc[common_genes]
    bulk = bulk.loc[common_genes]

    # ===============================
    # 3. 构造 AnnData（非 core）
    # ===============================
    adata = anndata.AnnData(X=sc.T.values)   # cell × gene
    adata.var_names = sc.index.astype(str)
    adata.obs_names = sc.columns.astype(str)

    adata.obs["annotation"] = ann.loc[
        adata.obs_names, args.celltype_col
    ].astype(str).fillna("Unknown")

    # ===============================
    # ⭐⭐⭐ core algorithm start ⭐⭐⭐
    # ===============================
    t_core_start = time.time()
    update_mem()

    ag.init(
        adata,
        use_highly_variable=args.use_hvg,
        celltype_key="annotation"
    )

    update_mem()

    ag.optimize(
        ngen=args.ngen,
        nfeatures=args.nfeatures,
        seed=args.seed,
        mode="fixed"
    )

    update_mem()

    ag.select(close_to=(1, args.max_genes))

    update_mem()

    coef = ag.deconvolve(
        bulk.T,
        model="nnls"
    )  # sample × celltype

    update_mem()
    t_core_end = time.time()
    # ===============================
    # ⭐⭐⭐ core algorithm end ⭐⭐⭐
    # ===============================

    coef_df = pd.DataFrame(
        coef,
        index=bulk.columns,
        columns=sorted(adata.obs["annotation"].unique())
    )

    coef_norm = coef_df.div(coef_df.sum(axis=1), axis=0)

    out_file = os.path.join(outdir, f"{prefix}AutoGeneS_fraction.tsv")
    coef_norm.to_csv(out_file, sep="\t")

    # ===============================
    # 写 benchmark JSON
    # ===============================
    bench = {
        "core_mem_MB": round(core_mem_peak_mb, 3),
        "core_time_sec": round(t_core_end - t_core_start, 4)
    }

    with open(os.path.join(outdir, "AutoGeneS_benchmark.json"), "w") as f:
        json.dump(bench, f, indent=2)

    print("[AutoGeneS] Finished")
    print("[AutoGeneS] Core peak memory (MB):", bench["core_mem_MB"])
    print("[AutoGeneS] Core time (sec):", bench["core_time_sec"])
    print("[AutoGeneS] Output:", out_file)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--sc_counts", required=True,
                        help="scRNA expression (gene × cell)")
    parser.add_argument("--sc_anno", required=True,
                        help="scRNA annotation (cell × meta)")
    parser.add_argument("--bulk", required=True,
                        help="bulk expression (gene × sample)")
    parser.add_argument("--celltype_col", default="cellType")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--prefix", default="")
    parser.add_argument("--ngen", type=int, default=2000)
    parser.add_argument("--nfeatures", type=int, default=300)
    parser.add_argument("--max_genes", type=int, default=50)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--use_hvg", action="store_true")

    args = parser.parse_args()
    main(args)
