#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import time
import json
import psutil

import numpy as np
import pandas as pd
import scanpy as sc
import torch
import deconv


def main(args):

    print("==== DeconV start ====")

    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    prefix = args.prefix
    if prefix and not prefix.endswith("_"):
        prefix += "_"

    CELL_TYPE_KEY = args.celltype_col
    TOP_N_GENES = 5000

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
    # 1. 读入 bulk（gene × sample）
    # ===============================
    bulk_raw = pd.read_csv(args.bulk, sep="\t", index_col=0)
    bulk_df = bulk_raw.T                      # sample × gene

    # ===============================
    # 2. 读入 scRNA（gene × cell）
    # ===============================
    sc_raw = pd.read_csv(args.sc_counts, sep="\t", index_col=0)
    sc_expr = sc_raw.T                        # cell × gene

    # ===============================
    # 3. 注释
    # ===============================
    ann = pd.read_csv(args.sc_anno, sep="\t", index_col=0)

    common_cells = sc_expr.index.intersection(ann.index)
    sc_expr = sc_expr.loc[common_cells]
    ann = ann.loc[common_cells]

    # ===============================
    # 4. 构建 AnnData
    # ===============================
    adata = sc.AnnData(
        X=sc_expr.values,
        obs=ann,
        var=pd.DataFrame(index=sc_expr.columns)
    )

    if CELL_TYPE_KEY not in adata.obs.columns:
        raise ValueError(f"{CELL_TYPE_KEY} not found in annotation")

    # ===============================
    # 5. 对齐基因
    # ===============================
    common_genes = adata.var_names.intersection(bulk_df.columns)
    if len(common_genes) == 0:
        raise ValueError("No common genes between scRNA and bulk")

    adata = adata[:, common_genes]
    bulk_df = bulk_df[common_genes]

    # ===============================
    # 6. 设备
    # ===============================
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print("Using device:", device)

    # ===============================
    # 7. 初始化模型（non-core）
    # ===============================
    model = deconv.DeconV(
        adata=adata,
        bulk=bulk_df,
        cell_type_key=CELL_TYPE_KEY,
        dropout_type="separate",
        model_type="gamma",
        device=device,
        layer=None,
        top_n_variable_genes=TOP_N_GENES
    )

    # ===============================
    # ⭐⭐⭐ core algorithm start ⭐⭐⭐
    # ===============================
    t_core_start = time.time()
    update_mem()

    model.fit_reference(
        num_epochs=2000,
        lr=0.1,
        lrd=0.999
    )

    update_mem()

    prop = model.deconvolute(
        num_epochs=1000,
        lr=0.1,
        lrd=0.999
    )

    update_mem()
    t_core_end = time.time()
    # ===============================
    # ⭐⭐⭐ core algorithm end ⭐⭐⭐
    # ===============================

    # ===============================
    # 8. 输出比例（统一为 celltype × sample）
    # ===============================
    prop_df = pd.DataFrame(
        prop,
        index=bulk_df.index,      # sample
        columns=model.cell_types  # celltype
    ).T                           # 👉 celltype × sample

    prop_file = os.path.join(outdir, f"{prefix}DeconV_fraction.tsv")
    prop_df.to_csv(prop_file, sep="\t")

    # ===============================
    # 9. 写 benchmark JSON
    # ===============================
    bench = {
        "core_mem_MB": round(core_mem_peak_mb, 3),
        "core_time_sec": round(t_core_end - t_core_start, 4)
    }

    with open(os.path.join(outdir, "DeconV_benchmark.json"), "w") as f:
        json.dump(bench, f, indent=2)

    print("[DeconV] Finished")
    print("[DeconV] Core peak memory (MB):", bench["core_mem_MB"])
    print("[DeconV] Core time (sec):", bench["core_time_sec"])
    print("[DeconV] Output:", prop_file)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--bulk", required=True)
    parser.add_argument("--sc_counts", required=True)
    parser.add_argument("--sc_anno", required=True)
    parser.add_argument("--celltype_col", default="cellType")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--prefix", default="")

    args = parser.parse_args()
    main(args)
