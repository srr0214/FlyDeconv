run_LinDeconSeq <- function(
    bulk,
    tissue,
    sex = c("mix", "male", "female"),
    bulk_name = NULL,
    out_base = "results"
) {

  ## =========================
  ## 0. 依赖检查 + attach（关键）
  ## =========================
  if (!requireNamespace("LinDeconSeq", quietly = TRUE)) {
    stop("LinDeconSeq package is required.")
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("LinDeconSeq requires package 'MASS' (for rlm).")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("LinDeconSeq requires package 'dplyr' (for %>%).")
  }
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    stop("Package 'peakRAM' is required.")
  }

  ## ⚠️ 关键：必须 attach
  suppressPackageStartupMessages({
    library(dplyr)      # 提供 %>%
    library(MASS)       # 提供 rlm
  })

  ## =========================
  ## 1. 参数处理
  ## =========================
  sex <- match.arg(sex)
  method_name <- "LinDeconSeq"

  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }

  ## =========================
  ## 2. 输出路径 & benchmark
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
  ## 3. reference
  ## =========================
  ref_expr <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "celltype_expr"
  )

  ## =========================
  ## 4. 主流程
  ## =========================
  t_total_start <- Sys.time()
  prop_LinDeconSeq <- NULL
  core_time <- NA

  total_mem <- peakRAM::peakRAM({

    common_genes <- intersect(
      rownames(bulk),
      rownames(ref_expr)
    )

    bulk_use <- bulk[common_genes, , drop = FALSE]
    signature <- as.data.frame(
      ref_expr[common_genes, , drop = FALSE]
    )

    core_time <- system.time({
      core_mem <- peakRAM::peakRAM({

        ## ⚠️ LinDeconSeq 内部依赖 %>% 和 rlm
        fractions <- LinDeconSeq::deconSeq(
          bulk_use,      # 必须是位置参数
          signature,
          verbose = FALSE
        )

      })
    })

    ## 后处理
    prop_LinDeconSeq <- t(as.matrix(fractions))
    prop_LinDeconSeq[prop_LinDeconSeq < 0] <- 0
    prop_LinDeconSeq <- sweep(
      prop_LinDeconSeq,
      2,
      colSums(prop_LinDeconSeq),
      "/"
    )
    prop_LinDeconSeq <- as.data.frame(prop_LinDeconSeq)
  })

  t_total_end <- Sys.time()

  ## =========================
  ## 5. benchmark
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
      total_mem_MB   = core_mem$Peak_RAM_Used_MiB[1],
      stringsAsFactors = FALSE
    )
  )

  utils::write.csv(runtime_log, log_file, row.names = FALSE)

  ## =========================
  ## 6. 保存结果
  ## =========================
  prop_file <- file.path(
    prop_dir,
    paste0("prop_", method_name, "_", bulk_name, ".csv")
  )
  utils::write.csv(prop_LinDeconSeq, prop_file)

  invisible(list(
    prop    = prop_LinDeconSeq,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}
