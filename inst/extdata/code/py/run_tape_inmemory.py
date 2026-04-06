#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ======= 禁用所有交互式绘图（防止弹窗阻塞） =======
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
plt.show = lambda *args, **kwargs: None
# ===============================================

import argparse
import pandas as pd
import tempfile
import os
import json
import time
import psutil

from TAPE import Deconvolution


# =========================
# 内存监控工具（核心）
# =========================
_process = psutil.Process(os.getpid())
_core_mem_peak_mb = 0.0

def update_mem():
    """记录 Python 进程 RSS 峰值（MB）"""
    global _core_mem_peak_mb
    mem_mb = _process.memory_info().rss / 1024 / 1024
    if mem_mb > _core_mem_peak_mb:
        _core_mem_peak_mb = mem_mb


def main(args):

    # =========================================================
    # 0. 初始化内存探针（⚠️ 不算 core）
    # =========================================================
    update_mem()

    # ---------- 1. 读入数据（非核心阶段） ----------
    scRef = pd.read_csv(args.scRef, index_col=0, sep=args.sep)
    bulk  = pd.read_csv(args.bulk,  index_col=0, sep=args.sep)

    # scRef: celltype × gene
    # bulk : sample × gene
    common_genes = scRef.columns.intersection(bulk.columns)

    scRef = scRef.loc[:, common_genes]
    bulk  = bulk.loc[:, common_genes]

    print(f"Common genes: {len(common_genes)}")

    update_mem()

    # ---------- 2. 构造 fake GeneLength（非核心） ----------
    fake_len = pd.DataFrame({
        "Gene name": common_genes,
        "Transcript start (bp)": [1] * len(common_genes),
        "Transcript end (bp)": [1000] * len(common_genes)
    })

    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".txt")
    fake_len.to_csv(tmp.name, sep=",", index=False)

    update_mem()

    # =========================================================
    # 3. 核心 TAPE 计算（★★★★★ benchmark 核心 ★★★★★）
    # =========================================================
    t_core_start = time.time()
    update_mem()

    SignatureMatrix, CellFractionPrediction = Deconvolution(
        scRef,
        bulk,
        sep=args.sep,
        datatype=args.datatype,
        genelenfile=tmp.name,
        mode=args.mode,
        adaptive=True,
        variance_threshold=args.variance,
        batch_size=args.batch_size,
        epochs=args.epochs,
        seed=args.seed
    )

    update_mem()
    t_core_end = time.time()

    # ---------- 4. 输出结果（非核心） ----------
    CellFractionPrediction.to_csv(
        args.out_prefix + "_cell_fraction.tsv",
        sep="\t"
    )

    SignatureMatrix.to_csv(
        args.out_prefix + "_signature.tsv",
        sep="\t"
    )

    # =========================================================
    # 5. 输出 benchmark JSON（给 R 读）
    # =========================================================
    bench = {
        "core_mem_MB": round(_core_mem_peak_mb, 3),
        "core_time_sec": round(t_core_end - t_core_start, 4)
    }

    with open(args.out_prefix + "_benchmark.json", "w") as f:
        json.dump(bench, f, indent=2)

    print("==== TAPE finished successfully ====")
    print("Core peak memory (MB):", bench["core_mem_MB"])
    print("Core time (sec):", bench["core_time_sec"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run TAPE deconvolution")

    parser.add_argument("--scRef", required=True,
                        help="scRef matrix (celltype × gene)")
    parser.add_argument("--bulk",  required=True,
                        help="bulk matrix (sample × gene)")
    parser.add_argument("--out_prefix", required=True,
                        help="output prefix")

    parser.add_argument("--sep", default="\t")
    parser.add_argument("--datatype", default="TPM")
    parser.add_argument("--mode", default="overall")
    parser.add_argument("--variance", type=float, default=0.98)
    parser.add_argument("--batch_size", type=int, default=64)
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--seed", type=int, default=1)

    args = parser.parse_args()
    main(args)
