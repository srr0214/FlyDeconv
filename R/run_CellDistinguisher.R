run_CellDistinguisher <- function(
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
  method_name <- "CellDistinguisher"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  ## =========================
  ## 1. 依赖检查（不 attach）
  ## =========================
  if (!requireNamespace("CellDistinguisher", quietly = TRUE)) {
    stop("Method CellDistinguisher requires package 'CellDistinguisher'.")
  }
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    stop("Method CellDistinguisher requires package 'peakRAM'.")
  }
  
  ## =========================
  ## 2. 输出路径 & benchmark
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
  ## 3. 总体时间 & 内存
  ## =========================
  t_total_start <- Sys.time()
  
  total_mem <- peakRAM::peakRAM({
    
    ## ---------- Step 1: 输入准备（非 core） ----------
    exprLinear <- as.matrix(bulk)
    
    ## ---------- Step 2: distinguisher 识别（非 core） ----------
    distinguishers <- CellDistinguisher::gecd_CellDistinguisher(
      exprLinear = exprLinear,
      genesymb   = rownames(exprLinear),
      numCellClasses = 4,
      minDistinguisherAlternatives = 1,
      maxDistinguisherAlternatives = 50
    )
    
    ## ---------- Step 3: 核心去卷积 ----------
    core_time <- system.time({
      core_mem <- peakRAM::peakRAM({
        
        deconv <- CellDistinguisher::gecd_DeconvolutionByDistinguishers(
          exprLinear,
          distinguishers$bestDistinguishers,
          nonNegativeOnly = TRUE,
          convexSolution = TRUE
        )
        
      })
    })
    
    prop <- as.data.frame(deconv$sampleCompositions)
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
  ## 5. 保存比例矩阵
  ## =========================
  prop_file <- file.path(
    prop_dir,
    paste0("prop_", method_name, "_", bulk_name, ".csv")
  )
  write.csv(prop, prop_file)
  
  ## =========================
  ## 6. 写 benchmark
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
  
  ## =========================
  ## 7. 返回
  ## =========================
  invisible(list(
    prop    = prop,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}



