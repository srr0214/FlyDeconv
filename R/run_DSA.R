
run_DSA <- function(
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
  if (!requireNamespace("DSA", quietly = TRUE)) {
    stop("Package 'DSA' is required but not installed.")
  }
  
  ## ===========================
  ## 1. 参数处理
  ## ===========================
  sex <- match.arg(sex)
  method_name <- "DSA"
  
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
  ## 3. 获取 internal marker list（不计 core）
  ## ===========================
  marker_list <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "marker_list"
  )
  
  ## ===========================
  ## 4. 总时间 & 总内存
  ## ===========================
  t_total_start <- Sys.time()
  
  total_mem <- peakRAM::peakRAM({
    
    ## -------- 1. 预处理 --------
    bulk_log <- log(as.matrix(bulk) + 1)
    
    common_genes <- intersect(
      rownames(bulk_log),
      unique(unlist(marker_list))
    )
    
    bulk_log <- bulk_log[common_genes, , drop = FALSE]
    
    gene_list <- lapply(
      marker_list,
      function(g) intersect(g, rownames(bulk_log))
    )
    
    ## -------- 2. core 计算 --------
    core_time <- system.time({
      core_mem <- peakRAM::peakRAM({
        estimated_weight <- DSA::EstimateWeight(
          mix_ob    = bulk_log,
          gene_list = gene_list,
          method    = "LM"
        )
      })
    })
    
    ## -------- 3. 后处理 --------
    prop_DSA <- estimated_weight$weight
    
    abs_min <- abs(min(prop_DSA, na.rm = TRUE))
    prop_DSA <- ifelse(is.na(prop_DSA), 0, prop_DSA + abs_min)
    prop_DSA <- as.data.frame(prop_DSA)
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
  
  utils::write.csv(prop_DSA, prop_file)
  
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
    prop    = prop_DSA,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}


