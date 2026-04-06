#!/usr/bin/python
import sys
import os
import json
import time
import psutil
import random
import numpy as np

# =====================================================
# 1. 可靠地定位 scripts 目录（关键修复）
# =====================================================
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

CANDIDATE_SCRIPT_DIRS = [
    os.path.join(SCRIPT_DIR, "scripts"),
    os.path.join(SCRIPT_DIR, "GEDIT-master", "GEDITv3.0", "scripts")
]

SCRIPT_PATH = None
for d in CANDIDATE_SCRIPT_DIRS:
    if os.path.isdir(d):
        SCRIPT_PATH = d
        break

if SCRIPT_PATH is None:
    raise ImportError(
        "❌ Cannot find GEDIT scripts directory.\n"
        "Tried:\n" + "\n".join(CANDIDATE_SCRIPT_DIRS)
    )

sys.path.insert(0, SCRIPT_PATH)

print("### GEDIT3.py path :", SCRIPT_DIR)
print("### Using scripts  :", SCRIPT_PATH)

# =====================================================
# 2. GEDIT 原始模块
# =====================================================
import getSigGenesModal
import MatrixTools
import HandleInput

# =====================================================
# 3. 主函数
# =====================================================
def main():
    """
    usage default:
    python GEDIT3.py -mix SamplesMat.tsv -ref RefMat.tsv

    user selected parameters:
    python GEDIT3.py -mix SamplesMat.tsv -ref RefMat.tsv -numSigs SigsPerCT -method SigMethod -RS rowscaling
    """

    # =========================
    # core memory profiler
    # =========================
    process = psutil.Process(os.getpid())
    core_mem_peak_mb = 0.0

    def update_mem():
        nonlocal core_mem_peak_mb
        mem_mb = process.memory_info().rss / 1024 / 1024
        if mem_mb > core_mem_peak_mb:
            core_mem_peak_mb = mem_mb

    # =========================
    # 参数解析
    # =========================
    myArgs = HandleInput.checkInputs(sys.argv[1:])
    if myArgs[0] is False:
        print(myArgs[1:])
        return

    rawMix = myArgs[0]
    rawRef = myArgs[1]
    SigsPerCT = myArgs[2]
    SigMethod = myArgs[4]
    RowScaling = myArgs[5]
    MixFName = myArgs[6].split("/")[-1]
    RefFName = myArgs[7].split("/")[-1]
    outFile = myArgs[8]
    SaveFiles = myArgs[9]

    numCTs = len(rawRef[0]) - 1
    TotalSigs = int(SigsPerCT * numCTs)

    SampleNames = rawMix[0]
    CTNames = rawRef[0]

    # =========================
    # 非核心：预处理
    # =========================
    betRef = MatrixTools.remove0s(rawRef)
    normMix, normRef = MatrixTools.qNormMatrices(rawMix, betRef)
    sharedMix, sharedRef = MatrixTools.getSharedRows(normMix, betRef)

    if len(sharedMix) < 1:
        print("error: no gene names match between reference and mixture")
        return

    # =====================================================
    # ⭐⭐⭐ 核心计算开始 ⭐⭐⭐
    # =====================================================
    t_core_start = time.time()
    update_mem()

    SigRef = getSigGenesModal.returnSigMatrix(
        [CTNames] + sharedRef,
        SigsPerCT,
        TotalSigs,
        SigMethod
    )

    update_mem()

    SigMix, SigRef = MatrixTools.getSharedRows(sharedMix, SigRef)

    ScaledRef, ScaledMix = MatrixTools.RescaleRows(
        SigRef[1:], SigMix[1:], RowScaling
    )

    ScaledRef = [CTNames] + ScaledRef
    ScaledMix = [SampleNames] + ScaledMix

    update_mem()

    Predictions = MatrixTools.PerformRegression(ScaledMix, ScaledRef)

    update_mem()
    t_core_end = time.time()
    # =====================================================
    # ⭐⭐⭐ 核心计算结束 ⭐⭐⭐
    # =====================================================

    # =========================
    # 输出路径
    # =========================
    if outFile == "None":
        SaveLocation = "Output/" + MixFName
    elif outFile.endswith("/"):
        SaveLocation = outFile + MixFName
    else:
        SaveLocation = outFile

    # =========================
    # 非核心：写文件
    # =========================
    if SaveFiles != "None":
        MatrixTools.writeMatrix(
            [CTNames] + SigRef,
            SaveLocation + "_SignatureGenes.tsv"
        )

    MatrixTools.writeMatrix(
        Predictions,
        SaveLocation + "_CTPredictions.tsv"
    )

    # =========================
    # 写 benchmark JSON
    # =========================
    bench = {
        "core_mem_MB": round(core_mem_peak_mb, 3),
        "core_time_sec": round(t_core_end - t_core_start, 4)
    }

    with open(SaveLocation + "_benchmark.json", "w") as f:
        json.dump(bench, f, indent=2)

    print("✅ GEDIT finished")
    print("Core peak memory (MB):", bench["core_mem_MB"])
    print("Core time (sec):", bench["core_time_sec"])


if __name__ == "__main__":
    main()
