#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import subprocess
import numpy as np
import pandas as pd
import anndata

# ===== 新增：benchmark =====
import json
import time
import psutil


def run(cmd, msg=None):
    """Run shell command with log."""
    if msg:
        print(f"\n[DAISM] {msg}")
    print(">>", " ".join(cmd))
    subprocess.check_call(cmd)


def main(args):

    # =====================================================
    # core memory profiler（Python 内部 RSS）
    # =====================================================
    process = psutil.Process(os.getpid())
    core_mem_peak_mb = 0.0

    def update_mem():
        nonlocal core_mem_peak_mb
        mem_mb = process.memory_info().rss / 1024 / 1024
        if mem_mb > core_mem_peak_mb:
            core_mem_peak_mb = mem_mb

    # =====================================================
    # 0. 路径与环境检查（非 core）
    # =====================================================
    outdir = os.path.abspath(args.outdir)

    try:
        outdir.encode("ascii")
    except UnicodeEncodeError:
        raise RuntimeError(
            f"Output directory contains non-ASCII characters:\n{outdir}\n"
            "Please use an ASCII-only path (no Chinese characters)."
        )

    os.makedirs(outdir, exist_ok=True)

    output_dir = os.path.join(outdir, "output")
    pred_dir   = os.path.join(outdir, "prediction")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(pred_dir, exist_ok=True)

    # =====================================================
    # 1. 读取输入数据（非 core）
    # =====================================================
    print("\n[DAISM] Loading input data...")

    sc   = pd.read_csv(args.sc_counts, sep="\t", index_col=0)   # gene × cell
    ann  = pd.read_csv(args.sc_anno,   sep="\t", index_col=0)  # cell × meta
    bulk = pd.read_csv(args.bulk,      sep="\t", index_col=0)  # gene × sample

    common_genes = sc.index.intersection(bulk.index)
    if len(common_genes) == 0:
        raise ValueError("No common genes between scRNA and bulk.")

    sc   = sc.loc[common_genes]
    bulk = bulk.loc[common_genes]

    # =====================================================
    # 2. 构造 AnnData（非 core）
    # =====================================================
    print("[DAISM] Building AnnData object...")

    adata = anndata.AnnData(X=sc.T.values)
    adata.var_names = sc.index.astype(str)
    adata.obs_names = sc.columns.astype(str)

    if args.celltype_col not in ann.columns:
        raise KeyError(
            f"Column '{args.celltype_col}' not found in sc_anno file."
        )

    adata.obs["cell.type"] = ann.loc[
        adata.obs_names, args.celltype_col
    ].astype(str)

    sc_h5ad = os.path.join(outdir, "sc_for_daism.h5ad")
    adata.write_h5ad(sc_h5ad)

    bulk_file = os.path.join(outdir, "bulk_for_daism.txt")
    bulk.to_csv(bulk_file, sep="\t")

    # =====================================================
    # ⭐⭐⭐ core algorithm starts ⭐⭐⭐
    # =====================================================
    t_core_start = time.time()
    update_mem()

    # ===============================
    # 3. Generic simulation (core)
    # ===============================
    sim_dir = os.path.join(outdir, "generic")
    os.makedirs(sim_dir, exist_ok=True)

    run(
        [
            "daism", "Generic_simulation",
            "-platform", "S",
            "-aug", sc_h5ad,
            "-N", str(args.N),
            "-testexp", bulk_file,
            "-outdir", sim_dir
        ],
        msg="Running Generic_simulation"
    )

    update_mem()

    # ===============================
    # 4. Training (core)
    # ===============================
    print("\n[DAISM] Preparing training data...")

    mixsam = pd.read_csv(
        os.path.join(sim_dir, "output", "Generic_mixsam.txt"),
        sep="\t", index_col=0
    )

    mixsam = mixsam.clip(lower=0)
    mixsam = np.log2(mixsam + 1)

    mixsam_file = os.path.join(outdir, "Generic_mixsam_log.txt")
    mixsam.to_csv(mixsam_file, sep="\t")

    run(
        [
            "daism", "training",
            "-trainexp", mixsam_file,
            "-trainfra", os.path.join(sim_dir, "output", "Generic_mixfra.txt"),
            "-outdir", outdir
        ],
        msg="Training DAISM model"
    )

    update_mem()

    # ===============================
    # 5. Prediction (core)
    # ===============================
    run(
        [
            "daism", "prediction",
            "-testexp", bulk_file,
            "-model", os.path.join(output_dir, "DAISM_model.pkl"),
            "-celltype", os.path.join(output_dir, "DAISM_model_celltypes.txt"),
            "-feature", os.path.join(output_dir, "DAISM_model_feature.txt"),
            "-outdir", pred_dir
        ],
        msg="Running prediction"
    )

    update_mem()
    t_core_end = time.time()
    # =====================================================
    # ⭐⭐⭐ core algorithm ends ⭐⭐⭐
    # =====================================================

    # =====================================================
    # 6. 写 benchmark JSON（给 R 端）
    # =====================================================
    bench = {
        "core_mem_MB": round(core_mem_peak_mb, 3),
        "core_time_sec": round(t_core_end - t_core_start, 4)
    }

    bench_file = os.path.join(outdir, "DAISM_benchmark.json")
    with open(bench_file, "w") as f:
        json.dump(bench, f, indent=2)

    print("\n[DAISM] Pipeline finished successfully!")
    print("[DAISM] Core peak memory (MB):", bench["core_mem_MB"])
    print("[DAISM] Core time (sec):", bench["core_time_sec"])
    print(f"[DAISM] Results saved in:\n{outdir}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Run DAISM using gene×cell matrix instead of h5ad"
    )

    parser.add_argument("--sc_counts", required=True,
                        help="scRNA expression matrix (gene × cell, TSV)")
    parser.add_argument("--sc_anno", required=True,
                        help="scRNA annotation (cell × meta, TSV)")
    parser.add_argument("--bulk", required=True,
                        help="bulk expression (gene × sample, TSV)")
    parser.add_argument("--celltype_col", default="cellType",
                        help="column name of cell type in sc_anno")
    parser.add_argument("--N", type=int, default=20000,
                        help="number of simulated mixtures")
    parser.add_argument("--outdir", required=True,
                        help="output directory (ASCII only!)")

    args = parser.parse_args()
    main(args)
