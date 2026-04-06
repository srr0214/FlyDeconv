#' Run MCPcounter-based deconvolution
#'
#' Run MCPcounter on bulk transcriptomic data using built-in marker genes.
#'
#' @param bulk A gene expression matrix or data frame, with genes in rows
#'   and samples in columns.
#' @param tissue Character. Tissue name.
#' @param sex Character. One of `"mix"`, `"male"`, or `"female"`.
#' @param bulk_name Character. Optional bulk dataset name.
#' @param out_base Character. Output directory.
#'
#' @return An invisible list containing proportion estimates, runtime
#'   information, method name, tissue, and sex.
run_MCPcounter <- function(
    bulk,
    tissue,
    sex = c("mix", "male", "female"),
    bulk_name = NULL,
    out_base = "results"
) {
  ## ===========================
  ## 0. 依赖检查
  ## ===========================
  if (!requireNamespace("MCPcounter", quietly = TRUE)) {
    stop("Package 'MCPcounter' is required but not installed.", call. = FALSE)
  }
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    stop("Package 'peakRAM' is required but not installed.", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.", call. = FALSE)
  }
  
  ## ===========================
  ## 1. 参数处理
  ## ===========================
  sex <- match.arg(sex)
  method_name <- "MCPcounter"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  bulk <- as.matrix(bulk)
  
  if (is.null(rownames(bulk))) {
    stop("bulk must have rownames as gene IDs/symbols.", call. = FALSE)
  }
  
  disable_parallel(1)
  
  ## ===========================
  ## 2. 输出路径 & benchmark 表
  ## ===========================
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
  
  ## ===========================
  ## 3. 获取 internal marker
  ## ===========================
  marker_long <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "long"
  )
  
  if (!is.data.frame(marker_long)) {
    stop("marker_long must be a data.frame.", call. = FALSE)
  }
  
  required_cols <- c("gene", "Simplify.annotations")
  if (!all(required_cols %in% colnames(marker_long))) {
    stop(
      "marker_long must contain columns: ",
      paste(required_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  marker_df <- dplyr::transmute(
    marker_long,
    `HUGO symbols`   = gene,
    `Cell population` = Simplify.annotations
  )
  
  ## 去重并只保留 bulk 中存在的 marker
  marker_df <- unique(marker_df)
  marker_df <- marker_df[marker_df$`HUGO symbols` %in% rownames(bulk), , drop = FALSE]
  
  if (nrow(marker_df) == 0) {
    stop("No marker genes overlap with bulk matrix.", call. = FALSE)
  }
  
  ## ===========================
  ## 4. 总时间 & 总内存
  ## ===========================
  t_total_start <- Sys.time()
  core_time <- NA
  core_mem <- NULL
  prop_MCPcounter <- NULL
  
  total_mem <- peakRAM::peakRAM({
    
    ## -------- core 计算 --------
    core_time <- system.time({
      core_mem <- peakRAM::peakRAM({
        res_MCPcounter <- MCPcounter::MCPcounter.estimate(
          expression   = bulk,
          featuresType = "HUGO_symbols",
          genes        = marker_df
        )
      })
    })
    
    ## -------- 后处理 --------
    res_MCPcounter <- as.matrix(res_MCPcounter)
    res_MCPcounter[res_MCPcounter < 0] <- 0
    
    cs <- colSums(res_MCPcounter, na.rm = TRUE)
    valid_cols <- cs > 0
    
    if (!all(valid_cols)) {
      warning("Some MCPcounter result columns had zero sum and were left unnormalized.", call. = FALSE)
    }
    
    prop_MCPcounter <- res_MCPcounter
    prop_MCPcounter[, valid_cols] <- sweep(
      prop_MCPcounter[, valid_cols, drop = FALSE],
      2,
      cs[valid_cols],
      "/"
    )
    
    prop_MCPcounter <- as.data.frame(prop_MCPcounter)
  })
  
  t_total_end <- Sys.time()
  
  ## ===========================
  ## 5. benchmark 指标
  ## ===========================
  core_time_sec  <- unname(core_time["elapsed"])
  total_time_sec <- as.numeric(
    difftime(t_total_end, t_total_start, units = "secs")
  )
  
  core_mem_MB  <- core_mem$Peak_RAM_Used_MiB[1]
  total_mem_MB <- total_mem$Peak_RAM_Used_MiB[1]
  
  ## ===========================
  ## 6. 保存结果
  ## ===========================
  prop_file <- file.path(
    prop_dir,
    paste0("prop_", method_name, "_", bulk_name, ".csv")
  )
  
  utils::write.csv(prop_MCPcounter, prop_file)
  
  ## ===========================
  ## 7. 写 benchmark 表
  ## ===========================
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
  
  ## ===========================
  ## 8. 返回
  ## ===========================
  invisible(list(
    prop    = prop_MCPcounter,
    runtime = runtime_log[nrow(runtime_log), , drop = FALSE],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}