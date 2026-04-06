run_MIXTURE <- function(
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
  method_name <- "MIXTURE"

  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }

  ## =========================
  ## 1. 依赖检查（全部隐式依赖）
  ## =========================
  if (!requireNamespace("MIXTURE", quietly = TRUE))
    stop("Method MIXTURE requires package 'MIXTURE'.")
  if (!requireNamespace("BiocParallel", quietly = TRUE))
    stop("Method MIXTURE requires package 'BiocParallel'.")
  if (!requireNamespace("e1071", quietly = TRUE))
    stop("Method MIXTURE requires package 'e1071' (for svm).")
  if (!requireNamespace("abind", quietly = TRUE))
    stop("Method MIXTURE requires package 'abind'.")
  if (!requireNamespace("peakRAM", quietly = TRUE))
    stop("Package 'peakRAM' is required.")

  ## ⚠️ 关键：必须 attach（MIXTURE 内部裸调用）
  suppressPackageStartupMessages({
    library(e1071)   # svm()
    library(abind)   # abind()
  })

  ## =========================
  ## 2. 输出路径 & benchmark
  ## =========================
  prop_dir <- file.path(out_base, method_name)
  log_file <- file.path(out_base, "runtime_benchmark.csv")
  dir.create(prop_dir, recursive = TRUE, showWarnings = FALSE)

  runtime_log <- if (!file.exists(log_file)) {
    data.frame(
      method = character(),
      bulk   = character(),
      core_time_sec  = numeric(),
      total_time_sec = numeric(),
      core_mem_MB    = numeric(),
      total_mem_MB   = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    utils::read.csv(log_file, stringsAsFactors = FALSE)
  }

  ## =========================
  ## 3. reference
  ## =========================
  celltype_expr <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "celltype_expr"
  )

  ## =========================
  ## 4. 修复 BiocParallel（非常关键）
  ## =========================
  old_param <- BiocParallel::bpparam()
  BiocParallel::register(
    BiocParallel::MulticoreParam(workers = 1),
    default = TRUE
  )

  ## =========================
  ## 5. 主流程
  ## =========================
  t_total_start <- Sys.time()
  prop_MIXTURE <- NULL
  core_time <- NA

  total_mem <- peakRAM::peakRAM({

    common_genes <- intersect(
      rownames(bulk),
      rownames(celltype_expr)
    )

    bulk_MIX <- bulk[common_genes, , drop = FALSE]
    ref_MIX  <- celltype_expr[common_genes, , drop = FALSE]

    core_time <- system.time({

      res <- MIXTURE::MIXTURE(
        expressionMatrix = as.matrix(bulk_MIX),
        signatureMatrix  = as.matrix(ref_MIX),
        iter             = 100,
        functionMixture  = MIXTURE::nu.svm.robust.RFE,
        useCores         = 1,
        verbose          = FALSE,
        nullDist         = "SingleSubjectTest"
      )

    })

    prop_MIXTURE <- t(res[["Subjects"]][["MIXprop"]])
    prop_MIXTURE[prop_MIXTURE < 0] <- 0
    prop_MIXTURE <- sweep(
      prop_MIXTURE,
      2,
      colSums(prop_MIXTURE),
      "/"
    )
    prop_MIXTURE <- as.data.frame(prop_MIXTURE)
  })

  t_total_end <- Sys.time()

  ## =========================
  ## 6. 恢复 BiocParallel
  ## =========================
  BiocParallel::register(old_param, default = TRUE)

  ## =========================
  ## 7. benchmark
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
  ## 8. 保存结果
  ## =========================
  prop_file <- file.path(
    prop_dir,
    paste0("prop_", method_name, "_", bulk_name, ".csv")
  )
  utils::write.csv(prop_MIXTURE, prop_file)

  invisible(list(
    prop    = prop_MIXTURE,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}
