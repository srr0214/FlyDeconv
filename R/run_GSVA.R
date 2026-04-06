#' Run GSVA-based deconvolution
#'
#' Run GSVA on bulk transcriptomic data using built-in marker gene sets.
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
run_GSVA <- function(
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
  method_name <- "GSVA"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("Method GSVA requires package 'GSVA'.", call. = FALSE)
  }
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    stop("Method GSVA requires package 'peakRAM'.", call. = FALSE)
  }
  
  bulk <- as.matrix(bulk)
  
  if (is.null(rownames(bulk))) {
    stop("bulk must have rownames as gene IDs/symbols.", call. = FALSE)
  }
  
  ## =========================
  ## 1. 输出路径 & benchmark 表
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
  ## 2. 获取 marker list
  ## =========================
  marker_list <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "list"
  )
  
  if (!is.list(marker_list) || length(marker_list) == 0) {
    stop("marker_list is empty or invalid for tissue: ", tissue, call. = FALSE)
  }
  
  ## 只保留在 bulk 中存在的基因
  marker_list <- lapply(marker_list, intersect, rownames(bulk))
  marker_list <- marker_list[lengths(marker_list) > 0]
  
  if (length(marker_list) == 0) {
    stop("No marker genes overlap with bulk matrix for tissue: ", tissue, call. = FALSE)
  }
  
  ## =========================
  ## 3. 总时间 & 总内存
  ## =========================
  t_total_start <- Sys.time()
  core_time <- NA
  core_mem  <- NULL
  prop_gsva <- NULL
  
  total_mem <- peakRAM::peakRAM({
    
    ## -------- 1. 构建 GSVA 参数对象（非 core）--------
    gsvapar <- GSVA::gsvaParam(
      exprData = bulk,
      geneSets = marker_list
    )
    
    ## -------- 2. 核心步骤：GSVA --------
    core_time <- system.time({
      core_mem <- peakRAM::peakRAM({
        gsva_es <- GSVA::gsva(gsvapar)
      })
    })
    
    ## -------- 3. 后处理 --------
    prop_gsva <- apply(
      gsva_es, 1,
      function(x) scale(x, center = TRUE, scale = FALSE)
    )
    
    prop_gsva <- apply(
      prop_gsva, 1,
      function(x) round(x - min(x) + min(abs(x)), 3)
    )
    
    prop_gsva <- sweep(
      prop_gsva, 2,
      colSums(prop_gsva), FUN = "/"
    )
    
    colnames(prop_gsva) <- colnames(bulk)
    prop_gsva <- as.data.frame(prop_gsva)
  })
  
  t_total_end <- Sys.time()
  
  ## =========================
  ## 4. 时间 & 内存统计
  ## =========================
  core_time_sec  <- unname(core_time["elapsed"])
  total_time_sec <- as.numeric(
    difftime(t_total_end, t_total_start, units = "secs")
  )
  
  core_mem_MB  <- core_mem$Peak_RAM_Used_MiB[1]
  total_mem_MB <- total_mem$Peak_RAM_Used_MiB[1]
  
  ## =========================
  ## 5. 保存 prop 结果
  ## =========================
  prop_file <- file.path(
    prop_dir,
    paste0("prop_", method_name, "_", bulk_name, ".csv")
  )
  
  utils::write.csv(prop_gsva, prop_file)
  
  ## =========================
  ## 6. 写入 benchmark 表
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
  
  invisible(list(
    prop    = prop_gsva,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}