#' Run SECRET-based deconvolution
#'
#' Run SECRET on bulk transcriptomic data using built-in cell-type
#' reference expression matrices.
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
run_SECRET <- function(
    bulk,
    tissue,
    sex = c("mix", "male", "female"),
    bulk_name = NULL,
    out_base = "results"
) {
  ## ===========================
  ## 0. 依赖检查
  ## ===========================
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    stop("Package 'peakRAM' is required but not installed.", call. = FALSE)
  }
  if (!requireNamespace("SECRET", quietly = TRUE)) {
    stop("Package 'SECRET' is required but not installed.", call. = FALSE)
  }
  
  ## ===========================
  ## 1. 参数处理
  ## ===========================
  sex <- match.arg(sex)
  method_name <- "SECRET"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  bulk <- as.matrix(bulk)
  
  if (is.null(rownames(bulk))) {
    stop("bulk must have rownames.", call. = FALSE)
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
  ## 3. 获取 internal celltype expression
  ## ===========================
  ref_type <- switch(
    sex,
    mix    = "celltype_all",
    female = "celltype_female",
    male   = "celltype_male"
  )
  
  celltype_expr <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = ref_type
  )
  
  celltype_expr <- as.matrix(celltype_expr)
  
  if (is.null(rownames(celltype_expr))) {
    stop("celltype_expr must have rownames.", call. = FALSE)
  }
  
  ## ===========================
  ## 4. 总时间 & 总内存
  ## ===========================
  t_total_start <- Sys.time()
  core_time <- NA
  core_mem <- NULL
  prop_SECRET <- NULL
  
  total_mem <- peakRAM::peakRAM({
    
    ## -------- 1. 数据准备 --------
    common_genes <- intersect(
      rownames(bulk),
      rownames(celltype_expr)
    )
    
    if (length(common_genes) < 50) {
      stop("Too few common genes between bulk and celltype_expr: ", length(common_genes), call. = FALSE)
    }
    
    bulk_use <- bulk[common_genes, , drop = FALSE]
    ref_use  <- celltype_expr[common_genes, , drop = FALSE]
    
    ## SECRET 必须显式提供 w
    w_vec <- rep(1, length(common_genes))
    names(w_vec) <- common_genes
    
    ## -------- 2. core 计算 --------
    core_time <- system.time({
      core_mem <- peakRAM::peakRAM({
        
        res_SECRET <- SECRET::SECRET(
          bulk        = as.matrix(bulk_use),
          reference   = as.matrix(ref_use),
          withUnknown = TRUE,
          w           = w_vec,
          yNorm       = "none",
          bNorm       = "none"
        )
        
      })
    })
    
    ## -------- 3. 后处理 --------
    if (is.null(res_SECRET) || length(res_SECRET) < 1) {
      stop("SECRET returned invalid result.", call. = FALSE)
    }
    
    prop_SECRET <- t(res_SECRET[[1]])
    prop_SECRET <- as.data.frame(prop_SECRET)
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
  
  utils::write.csv(prop_SECRET, prop_file)
  
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
    prop    = prop_SECRET,
    runtime = runtime_log[nrow(runtime_log), , drop = FALSE],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}