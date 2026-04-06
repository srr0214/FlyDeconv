run_qprog <- function(
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
  method_name <- "qprog"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
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
  ## 2. 获取 celltype 表达矩阵
  ## =========================
  celltype_expr <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "celltype_expr"
  )
  
  ## =========================
  ## 3. 总时间 & 总内存
  ## =========================
  t_total_start <- Sys.time()
  
  total_mem <- peakRAM::peakRAM({
    
    ## -------- 1. 基因交集（非 core）--------
    common_genes <- intersect(
      rownames(bulk),
      rownames(celltype_expr)
    )
    
    if (length(common_genes) == 0) {
      stop("No overlapping genes between bulk and celltype reference.")
    }
    
    bulk_use <- bulk[common_genes, , drop = FALSE]
    ref_use  <- celltype_expr[common_genes, , drop = FALSE]
    
    ## -------- 2. 核心步骤（QPROG）--------
    core_time <- system.time({
      core_mem <- peakRAM::peakRAM({
        
        decon <- granulator::deconvolute(
          m         = bulk_use,
          sigMatrix = ref_use,
          methods   = "qprog"
        )
        
      })
    })
    
    ## -------- 3. 后处理（非 core）--------
    prop_qprog <- t(decon$proportions$qprog_sig1)
    
    ## 去负值
    prop_qprog[prop_qprog < 0] <- 0
    
    ## 列归一化
    prop_qprog <- sweep(
      prop_qprog, 2, colSums(prop_qprog), FUN = "/"
    )
    
    prop_qprog <- as.data.frame(prop_qprog)
  })
  
  t_total_end <- Sys.time()
  
  ## =========================
  ## 4. 时间 & 内存
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
  
  utils::write.csv(prop_qprog, prop_file)
  
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
    prop    = prop_qprog,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}
