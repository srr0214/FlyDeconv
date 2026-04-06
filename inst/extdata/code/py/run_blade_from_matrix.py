#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import numpy as np
import pandas as pd
import time
import json
import psutil

import BLADE_Deconvolution.BLADE_framework as bf
from BLADE_Deconvolution.BLADE import BLADE


# ======================================================
# 修复后的 BLADE_wrapper（与你当前一致）
# ======================================================
def BLADE_wrapper(X, stdX, Y, Alpha, Alpha0, Kappa0, SY,
                  Init_Fraction, Rep, Crit='E_step', fsel=0.5):

    Ngene, Nsample = Y.shape
    Ncell = X.shape[1]

    Mu0 = X
    logY = np.log(Y + 1e-6)

    SigmaY = np.tile(
        np.std(logY, axis=1)[:, None],
        (1, Nsample)
    ) * SY

    Omega_Init = np.maximum(stdX, 1e-6)
    Beta0 = Alpha0 * Omega_Init

    Nu_Init = np.zeros((Nsample, Ngene, Ncell))
    for i in range(Nsample):
        Nu_Init[i, :, :] = Mu0

    Beta_Init = (
        np.random.gamma(shape=1, size=(Nsample, Ncell))
        + Init_Fraction.T * 10
    )

    obs = BLADE(
        logY, SigmaY, Mu0,
        Alpha, Alpha0, Beta0, Kappa0,
        Nu_Init, Omega_Init, Beta_Init,
        fix_Nu=True,
        fix_Omega=True
    )

    obs.Optimize()
    obs.Fix_par['Nu'] = False
    obs.Fix_par['Omega'] = False
    obs.Optimize()

    setting = {
        'Alpha': Alpha,
        'Alpha0': Alpha0,
        'Kappa0': Kappa0,
        'SigmaY': SY,
        'Rep': Rep
    }

    return obs, setting


bf.BLADE_wrapper = BLADE_wrapper


# ======================================================
# 主流程
# ======================================================
def main(args):

    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

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
    # 1. 读入数据（非 core）
    # ===============================
    bulk = pd.read_csv(args.bulk, sep="\t", index_col=0)          # gene × sample
    X    = pd.read_csv(args.celltype_ref, sep="\t", index_col=0) # gene × celltype
    sc   = pd.read_csv(args.sc_counts, sep="\t", index_col=0)    # gene × cell
    ann  = pd.read_csv(args.sc_anno, sep="\t", index_col=0)      # cell × meta

    common_genes = bulk.index.intersection(X.index)
    if len(common_genes) == 0:
        raise ValueError("No common genes between bulk and reference.")

    bulk = bulk.loc[common_genes]
    X    = X.loc[common_genes]

    # ===============================
    # 2. 计算 stdX（非 core）
    # ===============================
    stdX = pd.DataFrame(index=common_genes, columns=X.columns)

    for ct in X.columns:
        cells = ann.index[ann[args.celltype_col] == ct]
        if len(cells) == 0:
            raise ValueError(f"No cells found for cell type: {ct}")
        stdX[ct] = np.log1p(sc.loc[common_genes, cells]).std(axis=1)

    stdX = stdX.fillna(1e-6)

    # ===============================
    # ⭐⭐⭐ core algorithm start ⭐⭐⭐
    # ===============================
    t_core_start = time.time()
    update_mem()

    Ngene, Nsample = bulk.shape
    Ind_Marker = np.ones(Ngene, dtype=bool)
    Ind_sample = np.ones(Nsample, dtype=bool)

    result = bf.BLADE_framework(
        X=X.values,
        stdX=stdX.values,
        Y=bulk.values,
        Ind_Marker=Ind_Marker,
        Ind_sample=Ind_sample,
        Nrep=3,
        Nrepfinal=10,
        Njob=args.njob
    )

    final = result[0]

    F = pd.DataFrame(
        final.ExpF(final.Beta),
        index=bulk.columns,
        columns=X.columns
    )

    update_mem()
    t_core_end = time.time()
    # ===============================
    # ⭐⭐⭐ core algorithm end ⭐⭐⭐
    # ===============================

    prefix = args.prefix
    if prefix and not prefix.endswith("_"):
        prefix = prefix + "_"

    out_file = os.path.join(outdir, f"{prefix}BLADE_fraction.tsv")
    F.to_csv(out_file, sep="\t")

    # ===============================
    # 写 benchmark JSON
    # ===============================
    bench = {
        "core_mem_MB": round(core_mem_peak_mb, 3),
        "core_time_sec": round(t_core_end - t_core_start, 4)
    }

    with open(os.path.join(outdir, "BLADE_benchmark.json"), "w") as f:
        json.dump(bench, f, indent=2)

    print("[BLADE] Finished successfully")
    print("[BLADE] Core peak memory (MB):", bench["core_mem_MB"])
    print("[BLADE] Core time (sec):", bench["core_time_sec"])
    print("[BLADE] Output:", out_file)


# ======================================================
# CLI
# ======================================================
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Run BLADE using gene×cell / gene×celltype matrices"
    )

    parser.add_argument("--sc_counts", required=True)
    parser.add_argument("--sc_anno", required=True)
    parser.add_argument("--celltype_ref", required=True)
    parser.add_argument("--bulk", required=True)
    parser.add_argument("--celltype_col", default="cellType")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--prefix", default="")
    parser.add_argument("--njob", type=int, default=4)

    args = parser.parse_args()
    main(args)
