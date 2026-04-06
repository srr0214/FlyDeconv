#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use("Agg")   # ⭐ 无图形后端

import sys
import argparse
import json
import os
import time
import psutil

import numpy as np
import pandas as pd

# ===== 1. 指定 GLDADec 路径 =====
GLDADEC_SRC = "/project/deconv_py/GLDADec"
if GLDADEC_SRC not in sys.path:
    sys.path.insert(0, GLDADEC_SRC)

from run import pipeline


def main(args):

    print("==== GLDADec start ====")

    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    prefix = args.prefix
    if prefix and not prefix.endswith("_"):
        prefix += "_"

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
    # 1. 读入 bulk（non-core）
    # ===============================
    bulk = pd.read_csv(args.bulk, index_col=0)
    bulk.index = bulk.index.astype(str)

    bulk = bulk[
        ~bulk.index.str.lower().isin(["nan", "none", "na", ""])
    ]
    bulk = bulk[~bulk.index.duplicated()]

    print("Bulk shape after cleaning:", bulk.shape)

    # ===============================
    # 2. 读入 marker gene list（non-core）
    # ===============================
    with open(args.marker_json) as f:
        marker_dict = json.load(f)

    # ===============================
    # 3. 初始化 Pipeline（non-core）
    # ===============================
    pp = pipeline.Pipeline(verbose=True)

    pp.from_predata(
        bulk,
        target_samples=[],
        do_ann=False,
        linear2log=False,
        log2linear=False,
        do_drop=True,
        do_batch_norm=False,
        do_quantile=False,
        remove_noise=False
    )

    # ===============================
    # ⭐⭐⭐ core algorithm start ⭐⭐⭐
    # ===============================
    t_core_start = time.time()
    update_mem()

    # ===== 4. 基因筛选 =====
    pp.gene_selection(
        method="CV",
        outlier=True,
        topn=args.topn
    )
    update_mem()

    # ===== 5. 加 marker =====
    pp.add_marker_genes(
        target_cells=list(marker_dict.keys()),
        add_dic=marker_dict
    )
    update_mem()

    # ===== 6. 准备 LDA =====
    pp.deconv_prep(
        random_sets=[args.seed],
        do_plot=False,
        specific=True,
        prior_norm=True,
        norm_scale=10
    )
    update_mem()

    # ===== 7. 正式去卷积 =====
    pp.deconv(
        n=args.nrep,
        add_topic=0,
        n_iter=args.n_iter,
        alpha=args.alpha,
        eta=args.eta,
        refresh=20,
        initial_conf=1.0,
        seed_conf=1.0,
        other_conf=0.0
    )
    update_mem()

    t_core_end = time.time()
    # ===============================
    # ⭐⭐⭐ core algorithm end ⭐⭐⭐
    # ===============================

    # ===============================
    # 8. 汇总结果（non-core）
    # ===============================
    merge_res = pp.merge_total_res
    deconv_res = sum(merge_res) / len(merge_res)

    # ⚠️ GLDADec 输出通常是 sample × celltype
    # 统一转为 celltype × sample
    if deconv_res.shape[0] == bulk.shape[1]:
        deconv_res = deconv_res.T

    out_file = os.path.join(outdir, f"{prefix}GLDADec_fraction.csv")
    deconv_res.to_csv(out_file)

    # ===============================
    # 9. 写 benchmark JSON
    # ===============================
    bench = {
        "core_mem_MB": round(core_mem_peak_mb, 3),
        "core_time_sec": round(t_core_end - t_core_start, 4)
    }

    with open(os.path.join(outdir, "GLDADec_benchmark.json"), "w") as f:
        json.dump(bench, f, indent=2)

    print("[GLDADec] Finished")
    print("[GLDADec] Core peak memory (MB):", bench["core_mem_MB"])
    print("[GLDADec] Core time (sec):", bench["core_time_sec"])
    print("[GLDADec] Output:", out_file)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--bulk", required=True)
    parser.add_argument("--marker_json", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--prefix", default="")
    parser.add_argument("--topn", type=int, default=150)
    parser.add_argument("--nrep", type=int, default=10)
    parser.add_argument("--n_iter", type=int, default=200)
    parser.add_argument("--alpha", type=float, default=0.1)
    parser.add_argument("--eta", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=42)

    args = parser.parse_args()
    main(args)
