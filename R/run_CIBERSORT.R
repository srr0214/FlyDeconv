run_CIBERSORT <- function(
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
  method_name <- "CIBERSORT"

  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }

  ## =========================
  ## 1. 依赖检查
  ## =========================
  if (!requireNamespace("CIBERSORT", quietly = TRUE)) {
    stop("Method CIBERSORT requires package 'CIBERSORT'.")
  }
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    stop("Method CIBERSORT requires package 'peakRAM'.")
  }

  ## ⚠️ CIBERSORT 必须 attach
  suppressPackageStartupMessages(library(CIBERSORT))

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
  ## 3. reference（cell-type signature）
  ## =========================
  ref_expr <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "celltype_expr"
  )

  ## =========================
  ## 4. 运行 CIBERSORT（关键修复点）
  ## =========================
  t_total_start <- Sys.time()
  prop <- NULL
  core_time <- NA

  ## 保存并临时修改 mc.cores
  old_mc_cores <- getOption("mc.cores")
  options(mc.cores = 2)   # ⭐⭐⭐ 关键：CIBERSORT 内部硬要求 ≥2

  total_mem <- peakRAM::peakRAM({

    core_time <- system.time({

      res <- cibersort(
        sig_matrix   = as.matrix(ref_expr),
        mixture_file = as.matrix(bulk),
        perm         = 0,        # benchmark：不跑 permutation
        QN           = FALSE     # RNA-seq
      )

    })

    ## =========================
    ## 5. 后处理
    ## =========================
    ## 最后三列是：P-value / Correlation / RMSE
    prop <- t(
      res[, seq_len(ncol(res) - 3), drop = FALSE]
    )

    prop[prop < 0] <- 0
    prop <- sweep(prop, 2, colSums(prop), "/")
    prop <- as.data.frame(prop)
  })

  ## 恢复 mc.cores（非常重要）
  options(mc.cores = old_mc_cores)

  t_total_end <- Sys.time()

  ## =========================
  ## 6. 写 benchmark
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
      core_mem_MB    = total_mem$Peak_RAM_Used_MiB[1],
      total_mem_MB   = total_mem$Peak_RAM_Used_MiB[1],
      stringsAsFactors = FALSE
    )
  )

  utils::write.csv(runtime_log, log_file, row.names = FALSE)

  ## =========================
  ## 7. 保存比例矩阵
  ## =========================
  prop_file <- file.path(
    prop_dir,
    paste0("prop_", method_name, "_", bulk_name, ".csv")
  )
  utils::write.csv(prop, prop_file)

  ## =========================
  ## 8. 返回
  ## =========================
  invisible(list(
    prop    = prop,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}
