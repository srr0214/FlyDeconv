run_SCDC <- function(
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
  method_name <- "SCDC"
  
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
  ## 2. 获取 internal scRNA 数据
  ## =========================
  sc_expr <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "cell_expr"
  )
  
  sc_ann <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "annotation"
  )
  
  ## 确保 annotation 与表达矩阵对齐
  sc_ann <- sc_ann[colnames(sc_expr), , drop = FALSE]
  
  ## =========================
  ## 3. 总时间 & 总内存
  ## =========================
  t_total_start <- Sys.time()
  
  total_mem <- peakRAM::peakRAM({
    
    ## -------- 1. bulk ExpressionSet（非 core）--------
    bulk.eset <- Biobase::ExpressionSet(
      assayData = as.matrix(bulk)
    )
    
    ## -------- 2. scRNA ExpressionSet（非 core）--------
    sc.eset <- Biobase::ExpressionSet(
      assayData = as.matrix(sc_expr),
      phenoData = Biobase::AnnotatedDataFrame(sc_ann)
    )
    
    ## -------- 3. 参数准备（非 core）--------
    ct.varname <- "cellType"
    sample     <- "sample"
    ct.sub     <- unique(sc_ann[[ct.varname]])
    
    ## -------- 4. 核心步骤：SCDC --------
    core_time <- system.time({
      core_mem <- peakRAM::peakRAM({
        
        res <- SCDC::SCDC_prop(
          bulk.eset = bulk.eset,
          sc.eset   = sc.eset,
          ct.varname = ct.varname,
          sample     = sample,
          ct.sub     = ct.sub,
          weight.basis     = TRUE,
          Transform_bisque = TRUE
        )
        
      })
    })
    
    ## -------- 5. 后处理（非 core）--------
    prop_SCDC <- res$prop.est.mvw
    prop_SCDC <- t(prop_SCDC)
    prop_SCDC <- as.data.frame(prop_SCDC)
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
  
  utils::write.csv(prop_SCDC, prop_file)
  
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
    prop    = prop_SCDC,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}


