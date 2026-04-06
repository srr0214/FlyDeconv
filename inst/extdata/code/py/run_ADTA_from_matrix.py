#!/usr/bin/env python
import argparse
import sys
import os
import numpy as np
import pandas as pd

ADTD_SRC = os.environ.get("ADTD_SRC")
if ADTD_SRC not in sys.path:
    sys.path.insert(0, ADTD_SRC)

from ADTD import ADTD


def main(args):

    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    prefix = args.prefix
    if prefix and not prefix.endswith("_"):
        prefix += "_"

    # ---------- 1. 读入数据 ----------
    bulk = pd.read_csv(args.bulk, sep="\t", index_col=0)   # gene × sample
    ref  = pd.read_csv(args.ref,  sep="\t", index_col=0)   # gene × celltype

    # ---------- 2. 基因对齐 ----------
    genes = bulk.index.intersection(ref.index)
    if len(genes) == 0:
        raise ValueError("No common genes between bulk and reference")

    bulk = bulk.loc[genes]
    ref  = ref.loc[genes]

    X = ref.values.astype(np.double)     # gene × celltype
    Y = bulk.values.astype(np.double)    # gene × sample

    # ---------- 3. gene 权重 gamma ----------
    if args.gamma == "var":
        gamma = np.var(X, axis=1)
        gamma = gamma / np.mean(gamma)
    else:
        gamma = np.ones(X.shape[0])

    # ---------- 4. 运行 ADTD ----------
    C_hat, c_hat, x_hat, Delta = ADTD(
        X, Y, gamma,
        lambda1=args.lambda1,
        lambda2=args.lambda2,
        max_iterations=args.max_iter,
        eps=1e-8,
        quiet=False
    )

    # ---------- 5. 整理输出 ----------
    C_hat[C_hat < 0] = 0

    prop = pd.DataFrame(
        C_hat,
        index=ref.columns,      # cell types
        columns=bulk.columns   # samples
    )

    prop_file = os.path.join(outdir, f"{prefix}ADTD_fraction.csv")
    prop.to_csv(prop_file)

    np.save(os.path.join(outdir, f"{prefix}ADTD_c_hat.npy"), c_hat)
    np.save(os.path.join(outdir, f"{prefix}ADTD_x_hat.npy"), x_hat)
    np.save(os.path.join(outdir, f"{prefix}ADTD_Delta.npy"), Delta)

    print("[ADTD] Finished")
    print("[ADTD] Proportion file:", prop_file)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--bulk", required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--prefix", default="")
    parser.add_argument("--gamma", choices=["var", "uniform"], default="uniform")
    parser.add_argument("--lambda1", type=float, default=1.0)
    parser.add_argument("--lambda2", type=float, default=100.0)
    parser.add_argument("--max_iter", type=int, default=200)

    args = parser.parse_args()
    main(args)
