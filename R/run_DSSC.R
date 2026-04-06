run_DSSC <- function(
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
  method_name <- "DSSC"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  ## =========================
  ## 1. 依赖检查（不 library）
  ## =========================
  pkgs_needed <- c("scuttle", "dplyr", "limma", "peakRAM")
  for (p in pkgs_needed) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Method DSSC requires package: ", p)
    }
  }
  
  ## =========================
  ## 2. source DSSC 核心代码
  ## =========================
  


  dssc_r <- file.path("code", "R", "DSSC.R")

  if (!file.exists(dssc_r)) {
    stop("❌ DSSC DSSC.R not found at: ", dssc_r)
  }

  source(dssc_r, local = .GlobalEnv)  

  help_r <- file.path("code", "R", "function_help.R")
  
  if (!file.exists(help_r))
    stop("❌ DSSC function_help.R not found at: ", help_r)
    
  source(help_r, local = .GlobalEnv)  

  
  ## =========================
  ## 3. 输出路径 & benchmark 表
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
  ## 4. 获取 reference（cell-type expression）
  ## =========================
  ref_expr <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "celltype_expr"
  )
  
  ## =========================
  ## 5. 总体计时 & 内存
  ## =========================
  t_total_start <- Sys.time()
  
  total_mem <- peakRAM::peakRAM({
    
    ## ---------- Step 1: gene 交集 ----------
    common_gene <- intersect(
      rownames(ref_expr),
      rownames(bulk)
    )
    
    if (length(common_gene) < 100) {
      stop("DSSC: too few overlapping genes (", length(common_gene), ").")
    }
    
    T <- bulk[common_gene, , drop = FALSE]
    C <- ref_expr[common_gene, , drop = FALSE]
    
    ## ---------- Step 2: log1p ----------
    T_log <- log1p(T)
    C_log <- log1p(C)
    
    ## ---------- Step 3: 去除零方差基因 ----------
    keep <- apply(T_log, 1, stats::sd) > 0 &
      apply(C_log, 1, stats::sd) > 0
    
    T_use <- T_log[keep, , drop = FALSE]
    C_use <- C_log[keep, , drop = FALSE]
    
    ## ---------- Step 4: 构造正则矩阵 ----------
    Ss_mat <- diag(ncol(T_use))
    Sg_mat <- diag(nrow(T_use))
    
    ## ---------- Step 5: 核心 DSSC ----------
    core_time <- system.time({
      core_mem <- peakRAM::peakRAM({
        
        result <- DSSC(
          data_bulk = T_use,
          data_ref  = C_use,
          lambda1 = 0,
          lambda2 = 0,
          lambdaC = 1000,
          Ss = Ss_mat,
          Sg = Sg_mat
        )
        
      })
    })
    
    ## ---------- Step 6: 后处理 ----------
    prop_DSSC <- result[["p"]]
    
    rownames(prop_DSSC) <- colnames(C)
    prop_DSSC[prop_DSSC < 0] <- 0
    prop_DSSC <- sweep(prop_DSSC, 2, colSums(prop_DSSC), "/")
    prop_DSSC <- as.data.frame(prop_DSSC)
  })
  
  t_total_end <- Sys.time()
  
  ## =========================
  ## 6. 时间 & 内存统计
  ## =========================
  core_time_sec  <- unname(core_time["elapsed"])
  total_time_sec <- as.numeric(
    difftime(t_total_end, t_total_start, units = "secs")
  )
  
  core_mem_MB  <- core_mem$Peak_RAM_Used_MiB[1]
  total_mem_MB <- total_mem$Peak_RAM_Used_MiB[1]
  
  ## =========================
  ## 7. 保存比例矩阵
  ## =========================
  prop_file <- file.path(
    prop_dir,
    paste0("prop_", method_name, "_", bulk_name, ".csv")
  )
  write.csv(prop_DSSC, prop_file)
  
  ## =========================
  ## 8. 写 benchmark
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
    prop    = prop_DSSC,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}


