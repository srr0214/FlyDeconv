run_deconf <- function(
    bulk,
    tissue,
    sex = c("mix", "male", "female"),
    bulk_name = NULL,
    out_base = "results"
) {
  
  ## =========================
  ## 0. 依赖检查
  ## =========================
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    stop("Package 'peakRAM' is required.")
  }
  if (!requireNamespace("deconf", quietly = TRUE)) {
    stop("Package 'deconf' is required.")
  }
  
  ## =========================
  ## 1. 参数处理
  ## =========================
  sex <- match.arg(sex)
  method_name <- "deconf"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  ## 固定方法参数（benchmark 统一）
  n.iterations    <- 2000
  error.threshold <- 1e-6
  
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
  ## 3. 总时间 & 总内存
  ## =========================
  t_total_start <- Sys.time()
  
  total_mem <- peakRAM::peakRAM({
    
    ## -------- 1. 数据准备（非 core）--------
    I <- as.matrix(bulk)   # 已标准化 bulk（TPM / logCPM）
    
    ## latent component 数：与 cell type 数量对齐（benchmark 约定）
    ann <- get_internal_data(
      tissue = tissue,
      sex    = sex,
      type   = "annotation"
    )
    
    n_ct <- length(unique(ann$cellType))
    
    ## -------- 2. core 计算 --------
    core_time <- system.time({
      core_mem <- peakRAM::peakRAM({
        
        res <- deconf::deconfounding(
          I,
          n.cell.types    = n_ct,
          n.iterations    = n.iterations,
          error.threshold = error.threshold
        )
        
      })
    })
    
    ## -------- 3. 后处理 --------
    prop_deconf <- res[["C"]][["Matrix"]]
    
    rownames(prop_deconf) <- paste0(
      "Latent_", seq_len(nrow(prop_deconf))
    )
    colnames(prop_deconf) <- colnames(bulk)
    
    prop_deconf <- as.data.frame(prop_deconf)
  })
  
  t_total_end <- Sys.time()
  
  ## =========================
  ## 4. benchmark 指标
  ## =========================
  core_time_sec  <- unname(core_time["elapsed"])
  total_time_sec <- as.numeric(
    difftime(t_total_end, t_total_start, units = "secs")
  )
  
  core_mem_MB  <- core_mem$Peak_RAM_Used_MiB[1]
  total_mem_MB <- total_mem$Peak_RAM_Used_MiB[1]
  
  ## =========================
  ## 5. 保存结果
  ## =========================
  prop_file <- file.path(
    prop_dir,
    paste0("prop_", method_name, "_", bulk_name, ".csv")
  )
  
  utils::write.csv(prop_deconf, prop_file)
  
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
  
  ## =========================
  ## 7. 返回
  ## =========================
  invisible(list(
    prop    = prop_deconf,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}



