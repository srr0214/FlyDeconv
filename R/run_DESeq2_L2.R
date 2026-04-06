run_DESeq2_L2 <- function(
    bulk,
    tissue,
    sex = c("mix", "male", "female"),
    bulk_name = NULL,
    out_base = "results"
) {

  ## =========================
  ## 0. 参数处理
  ## =========================
  sex <- match.arg(sex)
  method_name <- "DESeq2_L2"

  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }

  ## =========================
  ## 1. 依赖检查（关键）
  ## =========================
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Method DESeq2_L2 requires package 'DESeq2'.")
  }
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    stop("Method DESeq2_L2 requires package 'peakRAM'.")
  }

  ## ⚠️ 关键：必须 attach，unmix() 才会进入 search path
  suppressPackageStartupMessages({
    library(DESeq2)
  })

  ## =========================
  ## 2. 输出路径 & benchmark 表
  ## =========================
  prop_dir <- file.path(out_base, method_name)
  log_file <- file.path(out_base, "runtime_benchmark.csv")
  dir.create(prop_dir, recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(log_file)) {
    runtime_log <- data.frame(
      method = character(),
      bulk   = character(),
      core_time_sec  = numeric(),
      total_time_sec = numeric(),
      core_mem_MB    = numeric(),
      total_mem_MB   = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    runtime_log <- utils::read.csv(log_file, stringsAsFactors = FALSE)
  }

  ## =========================
  ## 3. 获取 reference（cell-type expression）
  ## =========================
  ref_expr <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "celltype_expr"
  )

  ## =========================
  ## 4. 总时间 & 总内存
  ## =========================
  t_total_start <- Sys.time()
  prop_DESeq2_L2 <- NULL
  core_time <- NA
  core_mem <- NULL

  total_mem <- peakRAM::peakRAM({

    ## ---------- Step 1: 基因交集 ----------
    common_genes <- intersect(
      rownames(bulk),
      rownames(ref_expr)
    )

    if (length(common_genes) == 0) {
      stop("No overlapping genes between bulk and reference.")
    }

    bulk_use <- bulk[common_genes, , drop = FALSE]
    ref_use  <- ref_expr[common_genes, , drop = FALSE]

    ## ---------- Step 2: 核心 DESeq2-L2 ----------
    core_time <- system.time({
      core_mem <- peakRAM::peakRAM({

        mix <- unmix(
          x      = bulk_use,
          pure   = ref_use,
          shift  = 1,      # 与你本地示例一致
          power  = 2,      # L2 loss
          format = "matrix"
        )

      })
    })

    ## ---------- Step 3: 后处理 ----------
    prop_DESeq2_L2 <- t(mix)
    prop_DESeq2_L2[prop_DESeq2_L2 < 0] <- 0
    prop_DESeq2_L2 <- sweep(
      prop_DESeq2_L2,
      2,
      colSums(prop_DESeq2_L2),
      "/"
    )
    prop_DESeq2_L2 <- as.data.frame(prop_DESeq2_L2)
  })

  t_total_end <- Sys.time()

  ## =========================
  ## 5. 写 benchmark
  ## =========================
  runtime_log <- rbind(
    runtime_log,
    data.frame(
      method = method_name,
      bulk   = bulk_name,
      core_time_sec  = unname(core_time["elapsed"]),
      total_time_sec = as.numeric(
        difftime(t_total_end, t_total_start, units = "secs")
      ),
      core_mem_MB    = core_mem$Peak_RAM_Used_MiB[1],
      total_mem_MB   = total_mem$Peak_RAM_Used_MiB[1],
      stringsAsFactors = FALSE
    )
  )

  utils::write.csv(runtime_log, log_file, row.names = FALSE)

  ## =========================
  ## 6. 保存比例矩阵
  ## =========================
  prop_file <- file.path(
    prop_dir,
    paste0("prop_", method_name, "_", bulk_name, ".csv")
  )
  utils::write.csv(prop_DESeq2_L2, prop_file)

  ## =========================
  ## 7. 返回
  ## =========================
  invisible(list(
    prop    = prop_DESeq2_L2,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}
