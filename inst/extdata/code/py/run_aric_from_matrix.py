#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import pandas as pd
import time
import json
import psutil

from ARIC import ARIC


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
    mix = pd.read_csv(args.bulk, sep="\t", index_col=0)  # gene × sample
    ref = pd.read_csv(args.ref,  sep="\t", index_col=0)  # gene × celltype

    # ===============================
    # 2. 对齐基因（非 core）
    # ===============================
    common_genes = mix.index.intersection(ref.index)
    if len(common_genes) == 0:
        raise ValueError("No common genes between bulk and reference.")

    mix = mix.loc[common_genes]
    ref = ref.loc[common_genes]

    # ===============================
    # 3. 写 ARIC 临时文件（非 core）
    # ===============================
    mix_file = os.path.join(outdir, f"{prefix}mix_ARIC.csv")
    ref_file = os.path.join(outdir, f"{prefix}ref_ARIC.csv")

    mix.to_csv(mix_file)
    ref.to_csv(ref_file)

    # ===============================
    # ⭐⭐⭐ core algorithm start ⭐⭐⭐
    # ===============================
    t_core_start = time.time()
    update_mem()

    out_file = os.path.join(outdir, f"{prefix}ARIC_fraction.csv")

    ARIC(
        mix_path=mix_file,
        ref_path=ref_file,
        save_path=out_file,
        selected_marker=args.selected_marker,
        scale=args.scale,
        delcol_factor=args.delcol_factor,
        iter_num=args.iter_num,
        confidence=args.confidence,
        unknown=args.unknown,
        is_methylation=False
    )

    update_mem()
    t_core_end = time.time()
    # ===============================
    # ⭐⭐⭐ core algorithm end ⭐⭐⭐
    # ===============================

    # ===============================
    # 写 benchmark JSON
    # ===============================
    bench = {
        "core_mem_MB": round(core_mem_peak_mb, 3),
        "core_time_sec": round(t_core_end - t_core_start, 4)
    }

    with open(os.path.join(outdir, "ARIC_benchmark.json"), "w") as f:
        json.dump(bench, f, indent=2)

    print("[ARIC] Finished")
    print("[ARIC] Core peak memory (MB):", bench["core_mem_MB"])
    print("[ARIC] Core time (sec):", bench["core_time_sec"])
    print("[ARIC] Output:", out_file)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Run ARIC using gene×celltype reference matrix"
    )

    parser.add_argument("--bulk", required=True,
                        help="bulk expression: gene × sample (TSV)")
    parser.add_argument("--ref", required=True,
                        help="reference: gene × celltype (TSV)")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--prefix", default="")
    parser.add_argument("--scale", type=float, default=0.1)
    parser.add_argument("--delcol_factor", type=int, default=10)
    parser.add_argument("--iter_num", type=int, default=10)
    parser.add_argument("--confidence", type=float, default=0.75)
    parser.add_argument("--selected_marker", action="store_true")
    parser.add_argument("--unknown", action="store_true")

    args = parser.parse_args()
    main(args)
