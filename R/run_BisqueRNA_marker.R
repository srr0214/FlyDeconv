run_BisqueRNA_marker <- function(
  bulk,
  tissue,
  sex = c("mix", "male", "female"),
  bulk_name = NULL,
  out_base = "results"
) {

  ## =========================
  ## 0. 依赖检查
  ## =========================
  if (!requireNamespace("peakRAM", quietly = TRUE))
    stop("Package 'peakRAM' is required.")
  if (!requireNamespace("Biobase", quietly = TRUE))
    stop("Package 'Biobase' is required.")
  if (!requireNamespace("BisqueRNA", quietly = TRUE))
    stop("Package 'BisqueRNA' is required.")

  ## =========================
  ## 1. 参数处理
  ## =========================
  sex <- match.arg(sex)
  method_name <- "BisqueRNA_marker"

  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }

  ## 固定方法参数（benchmark 统一）
  min_gene        <- 5
  max_gene        <- 200
  unique_markers  <- FALSE
  verbose         <- TRUE

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
  ## 3. 获取 marker dataframe
  ## =========================
  marker_long <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "marker_long"
  )

  ## =========================
  ## 3.1 关键修复：BisqueRNA-compatible marker 过滤
  ## =========================
  ## （1）先与 bulk 基因集取交集 —— 这是 BisqueRNA 内部真实逻辑
  marker_long <- marker_long[
    marker_long$gene %in% rownames(bulk),
  ]

  ## （2）重新统计每个 cell type 的 marker 数
  marker_count <- table(marker_long$Simplify.annotations)

  valid_celltypes   <- names(marker_count)[marker_count >= min_gene]
  invalid_celltypes <- names(marker_count)[marker_count <  min_gene]

  if (length(invalid_celltypes) > 0) {
    message(
      "[BisqueRNA_marker] Dropping ",
      length(invalid_celltypes),
      " cell types with insufficient markers after bulk filtering (< ",
      min_gene, "): ",
      paste(invalid_celltypes, collapse = ", ")
    )
  }

  ## （3）只保留合法 cell types
  marker_long <- marker_long[
    marker_long$Simplify.annotations %in% valid_celltypes,
  ]

  ## （4）安全检查
  if (length(unique(marker_long$Simplify.annotations)) < 2) {
    stop(
      "[BisqueRNA_marker] < 2 cell types remain after marker filtering. ",
      "Deconvolution is not meaningful."
    )
  }

  ## =========================
  ## 4. 总时间 & 总内存
  ## =========================
  t_total_start <- Sys.time()

  total_mem <- peakRAM::peakRAM({

    ## -------- 1. 构建 ExpressionSet --------
    bulk.eset <- Biobase::ExpressionSet(
      assayData = as.matrix(bulk)
    )

    ## -------- 2. 核心步骤（BisqueRNA）--------
    core_time <- system.time({
      core_mem <- peakRAM::peakRAM({

        res <- BisqueRNA::MarkerBasedDecomposition(
          bulk.eset = bulk.eset,
          markers   = marker_long,
          ct_col    = "Simplify.annotations",
          gene_col  = "gene",
          min_gene  = min_gene,
          max_gene  = max_gene,
          unique_markers = unique_markers,
          verbose   = verbose
        )

      })
    })

    ## -------- 3. 后处理 --------
    prop_BisqueRNA <- res$bulk.props

    ## min–max scaling
    prop_BisqueRNA <- apply(
      prop_BisqueRNA, 2,
      function(x) {
        if (max(x) == min(x)) return(rep(0, length(x)))
        (x - min(x)) / (max(x) - min(x))
      }
    )

    ## 列归一化
    prop_BisqueRNA <- sweep(
      prop_BisqueRNA,
      2,
      colSums(prop_BisqueRNA),
      FUN = "/"
    )

    prop_BisqueRNA <- as.data.frame(prop_BisqueRNA)
  })

  t_total_end <- Sys.time()

  ## =========================
  ## 5. 时间 & 内存
  ## =========================
  core_time_sec  <- unname(core_time["elapsed"])
  total_time_sec <- as.numeric(
    difftime(t_total_end, t_total_start, units = "secs")
  )

  core_mem_MB  <- core_mem$Peak_RAM_Used_MiB[1]
  total_mem_MB <- total_mem$Peak_RAM_Used_MiB[1]

  ## =========================
  ## 6. 保存结果
  ## =========================
  prop_file <- file.path(
    prop_dir,
    paste0("prop_", method_name, "_", bulk_name, ".csv")
  )

  utils::write.csv(prop_BisqueRNA, prop_file)

  ## =========================
  ## 7. 写 benchmark 表
  ## =========================
  runtime_log <- rbind(
    runtime_log,
    data.frame(
      method = method_name,
      bulk   = bulk_name,
      core_time_sec  = core_time_sec,
      total_time_sec = total_time_sec,
      core_mem_MB    = core_mem_MB,
      total_mem_MB   = total_mem_MB,
      stringsAsFactors = FALSE
    )
  )

  utils::write.csv(runtime_log, log_file, row.names = FALSE)

  ## =========================
  ## 8. 返回
  ## =========================
  invisible(list(
    prop    = prop_BisqueRNA,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex,
    dropped_celltypes = invalid_celltypes
  ))
}
