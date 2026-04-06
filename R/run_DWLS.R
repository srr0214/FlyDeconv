run_DWLS <- function(
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
  method_name <- "DWLS"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  ## =========================
  ## 1. 依赖检查（不 library）
  ## =========================
  if (!requireNamespace("DWLS", quietly = TRUE)) {
    stop("Method DWLS requires package 'DWLS'.")
  }
  if (!requireNamespace("MAST", quietly = TRUE)) {
    stop("Method DWLS requires package 'MAST'.")
  }
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    stop("Method DWLS requires package 'peakRAM'.")
  }
  
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
    runtime_log <- read.csv(log_file, stringsAsFactors = FALSE)
  }
  
  ## =========================
  ## 3. 获取 internal reference（celltype 表达矩阵）
  ## =========================
  Signature <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "celltype_expr"
  )
  
  ## =========================
  ## 4. 总体计时 & 内存
  ## =========================
  t_total_start <- Sys.time()
  
  total_mem <- peakRAM::peakRAM({
    
    Signature  <- as.matrix(Signature)
    bulkMatrix <- as.matrix(bulk)
    
    sample_names <- colnames(bulkMatrix)
    result_list  <- list()
    
    ## ---------- 核心 DWLS（逐 sample） ----------
    core_time <- system.time({
      core_mem <- peakRAM::peakRAM({
        
        for (sample in sample_names) {
          
          bulk_vec <- bulkMatrix[, sample]
          names(bulk_vec) <- rownames(bulkMatrix)
          
          trimmed <- DWLS::trimData(Signature, bulk_vec)
          S <- trimmed$sig
          B <- trimmed$bulk
          
          if (nrow(S) < 50) {
            warning(
              sprintf(
                "DWLS: sample %s skipped (only %d overlapping genes).",
                sample, nrow(S)
              )
            )
            next
          }
          
          res <- DWLS::solveDampenedWLS(S, B)
          result_list[[sample]] <- res
        }
      })
    })
    
    ## ---------- 后处理 ----------
    if (length(result_list) == 0) {
      stop("DWLS failed: no samples produced valid results.")
    }
    
    prop_DWLS <- do.call(cbind, result_list)
    
    ## 非负 + 归一化
    prop_DWLS[prop_DWLS < 0] <- 0
    prop_DWLS <- sweep(prop_DWLS, 2, colSums(prop_DWLS), "/")
    prop_DWLS <- as.data.frame(prop_DWLS)
  })
  
  t_total_end <- Sys.time()
  
  ## =========================
  ## 5. 时间 & 内存统计
  ## =========================
  core_time_sec  <- unname(core_time["elapsed"])
  total_time_sec <- as.numeric(
    difftime(t_total_end, t_total_start, units = "secs")
  )
  
  core_mem_MB  <- core_mem$Peak_RAM_Used_MiB[1]
  total_mem_MB <- total_mem$Peak_RAM_Used_MiB[1]
  
  ## =========================
  ## 6. 保存比例矩阵
  ## =========================
  prop_file <- file.path(
    prop_dir,
    paste0("prop_", method_name, "_", bulk_name, ".csv")
  )
  write.csv(prop_DWLS, prop_file)
  
  ## =========================
  ## 7. 写入 benchmark 表
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
  
  write.csv(runtime_log, log_file, row.names = FALSE)
  
  invisible(list(
    prop    = prop_DWLS,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}


