run_NNLS <- function(
    bulk,
    tissue,
    sex = c("mix", "male", "female"),
    bulk_name = NULL,
    out_base = "results"
) {

  ## =========================
  ## 0. 依赖检查（关键）
  ## =========================
  if (!requireNamespace("MuSiC", quietly = TRUE)) {
    stop("NNLS requires package 'MuSiC'.")
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("NNLS requires package 'SingleCellExperiment'.")
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("NNLS requires package 'SummarizedExperiment'.")
  }
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    stop("Package 'peakRAM' is required.")
  }

  ## ⚠️ 关键：必须 attach，MuSiC 内部直接调用 counts()
  suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(SummarizedExperiment)
  })

  ## =========================
  ## 1. 参数处理
  ## =========================
  sex <- match.arg(sex)
  method_name <- "NNLS"

  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }

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
  ## 3. 读取 scRNA 数据
  ## =========================
  sc_expr <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "cell_expr"
  )

  sc_ann <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "annotation"
  )

  sc_ann <- sc_ann[colnames(sc_expr), , drop = FALSE]

  if (!"sample" %in% colnames(sc_ann)) {
    stop("NNLS (MuSiC) requires 'sample' column in annotation.")
  }

  ## =========================
  ## 4. 主流程
  ## =========================
  t_total_start <- Sys.time()
  prop_NNLS <- NULL
  core_time <- NA
  core_mem <- NULL

  total_mem <- peakRAM::peakRAM({

    ## ---- 构建 SCE ----
    sce <- SingleCellExperiment(
      assays  = list(counts = as.matrix(sc_expr)),
      colData = sc_ann
    )

    SummarizedExperiment::rowData(sce)$gene.name <- rownames(sce)
    sce$sample.id <- sc_ann$sample

    ## ---- 核心：MuSiC NNLS ----
    core_time <- system.time({
      core_mem <- peakRAM::peakRAM({

        Est <- MuSiC::music_prop(
          bulk.mtx = as.matrix(bulk),
          sc.sce   = sce,
          clusters = "cellType",
          samples  = "sample.id",
          verbose  = FALSE
        )

      })
    })

    ## ---- 后处理 ----
    prop_NNLS <- t(Est$Est.prop.allgene)
    prop_NNLS[prop_NNLS < 0] <- 0
    prop_NNLS <- sweep(prop_NNLS, 2, colSums(prop_NNLS), "/")
    prop_NNLS <- as.data.frame(prop_NNLS)
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
      total_time_sec = as.numeric(difftime(t_total_end, t_total_start, units = "secs")),
      core_mem_MB    = core_mem$Peak_RAM_Used_MiB[1],
      total_mem_MB   = total_mem$Peak_RAM_Used_MiB[1],
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
  utils::write.csv(prop_NNLS, prop_file)

  invisible(list(
    prop    = prop_NNLS,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}
