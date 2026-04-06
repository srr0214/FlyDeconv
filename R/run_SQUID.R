run_SQUID <- function(
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
  if (!requireNamespace("SQUID", quietly = TRUE)) {
    stop("Package 'SQUID' is required but not installed.")
  }
  if (!requireNamespace("DWLS", quietly = TRUE)) {
    stop("Package 'DWLS' is required by SQUID but not installed.")
  }
  if (!requireNamespace("magrittr", quietly = TRUE)) {
    stop("Package 'magrittr' is required by SQUID but not installed.")
  }
  
  ## ===========================
  ## 1. 参数处理
  ## ===========================
  sex <- match.arg(sex)
  method_name <- "SQUID"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  
  help_r <- file.path("code", "R", "helper_functions.R")
  
  if (!file.exists(help_r))
    stop("❌ SQUID helper_functions.R not found at: ", help_r)
  source(help_r, local = .GlobalEnv)


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
  ## 3. 获取 internal 单细胞数据（不计 core）
  ## ===========================
  cell_expr <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "cell_expr"
  )
  
  ann <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "annotation"
  )
  
  ## ===========================
  ## 4. 总时间 & 总内存
  ## ===========================
  t_total_start <- Sys.time()
  
  total_mem <- peakRAM::peakRAM({
    
    ## -------- 1. 基因交集 --------
    common_gene <- intersect(
      rownames(bulk),
      rownames(cell_expr)
    )
    
    bulk_use <- bulk[common_gene, , drop = FALSE]
    expr_use <- cell_expr[common_gene, , drop = FALSE]
    
    ## -------- 2. 构建先验比例矩阵 P --------
    cell_frac <- table(ann$cellType)
    cell_frac <- cell_frac / sum(cell_frac)
    
    bulk_samples <- colnames(bulk_use)
    
    P <- matrix(
      rep(cell_frac, length(bulk_samples)),
      nrow = length(cell_frac),
      ncol = length(bulk_samples)
    )
    rownames(P) <- names(cell_frac)
    colnames(P) <- bulk_samples
    
    ## -------- 3. core 计算（⚠️ SQUID 强依赖 attach） --------
    core_time <- system.time({
      core_mem <- peakRAM::peakRAM({
        
        ## SQUID 内部裸用 trimData() 和 %>%
        suppressPackageStartupMessages({
          library(DWLS)
          library(magrittr)
        })
        
        RESULTS <- SQUID::SQUID(
          B           = bulk_use,
          scC         = expr_use,
          scMeta      = ann,
          pB          = NULL,
          P           = P,
          LeaveOneOut = FALSE
        )
      })
    })
    
    ## -------- 4. 后处理（base R） --------
    df <- RESULTS[, c("cellType", "sample.id", "observed_fraction")]
    
    prop_SQUID <- reshape(
      df,
      timevar   = "sample.id",
      idvar     = "cellType",
      direction = "wide"
    )
    
    rownames(prop_SQUID) <- prop_SQUID$cellType
    prop_SQUID$cellType <- NULL
    
    colnames(prop_SQUID) <- sub(
      "^observed_fraction\\.",
      "",
      colnames(prop_SQUID)
    )
    
    prop_SQUID <- prop_SQUID[order(rownames(prop_SQUID)), , drop = FALSE]
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
  
  utils::write.csv(prop_SQUID, prop_file)
  
  ## ===========================
  ## 7. 写入 benchmark 表
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
    prop    = prop_SQUID,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}
