run_SVMDECON <- function(
    bulk,
    tissue,
    sex = c("mix", "male", "female"),
    bulk_name = NULL,
    out_base = "results"
) {
  
  ## ===========================
  ## 0. 依赖检查（不 attach）
  ## ===========================
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    stop("Package 'peakRAM' is required but not installed.")
  }
  if (!requireNamespace("ADAPTS", quietly = TRUE)) {
    stop("Package 'ADAPTS' is required but not installed.")
  }
  
  ## ===========================
  ## 1. 参数处理
  ## ===========================
  sex <- match.arg(sex)
  method_name <- "SVMDECON"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
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
  ## 3. 获取 internal celltype expression（不计 core）
  ## ===========================
  celltype_expr <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "celltype_expr"
  )
  
  ## ===========================
  ## 4. 总时间 & 总内存
  ## ===========================
  t_total_start <- Sys.time()
  
  total_mem <- peakRAM::peakRAM({
    
    ## -------- 1. 数据准备 --------
    common_genes <- intersect(
      rownames(bulk),
      rownames(celltype_expr)
    )
    
    m <- bulk[common_genes, , drop = FALSE]
    B <- celltype_expr[common_genes, , drop = FALSE]
    
    ## 标准化（SVMDECON 要求）
    m_standardized <- scale(m)
    B_standardized <- scale(B)
    
    ## -------- 2. core 计算 --------
    core_time <- system.time({
      core_mem <- peakRAM::peakRAM({
        cellEst_svmdecon <- ADAPTS::estCellPercent.svmdecon(
          refExpr  = B_standardized,
          geneExpr = m_standardized
        )
      })
    })
    
    ## -------- 3. 后处理（归一化为比例）--------
    cellEst_svmdecon_norm <- sweep(
      cellEst_svmdecon,
      2,
      colSums(cellEst_svmdecon),
      FUN = "/"
    )
    
    prop_SVMDECON <- as.data.frame(cellEst_svmdecon_norm)
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
  
  utils::write.csv(prop_SVMDECON, prop_file)
  
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
    prop    = prop_SVMDECON,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}

