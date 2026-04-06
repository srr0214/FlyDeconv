run_BisqueRNA_refer <- function(
    bulk,
    tissue,
    sex = c("mix", "male", "female"),
    bulk_name = NULL,
    out_base = "results"
) {
  
  ## =========================
  ## 0. 依赖检查
  ## =========================
  if (!requireNamespace("peakRAM", quietly = TRUE))
    stop("Package 'peakRAM' is required.")
  if (!requireNamespace("Biobase", quietly = TRUE))
    stop("Package 'Biobase' is required.")
  if (!requireNamespace("BisqueRNA", quietly = TRUE))
    stop("Package 'BisqueRNA' is required.")
  
  ## =========================
  ## 1. 参数处理
  ## =========================
  sex <- match.arg(sex)
  method_name <- "BisqueRNA_refer"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  ## 固定方法参数（benchmark 统一）
  verbose <- FALSE
  pseudo_scale <- 10
  
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
  ## 3. 获取 internal scRNA 数据（非 core）
  ## =========================
  sc_counts <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "cell_expr"
  )
  
  sc_anno <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "annotation"
  )
  
  stopifnot(all(colnames(sc_counts) == rownames(sc_anno)))
  
  ## =========================
  ## 4. 总时间 & 总内存
  ## =========================
  t_total_start <- Sys.time()
  
  total_mem <- peakRAM::peakRAM({
    
    ## -------- 1. gene 对齐（非 core）--------
    common_genes <- intersect(
      rownames(bulk),
      rownames(sc_counts)
    )
    
    bulk_use <- bulk[common_genes, , drop = FALSE]
    sc_use   <- sc_counts[common_genes, , drop = FALSE]
    
    ## -------- 2. pseudo-count（BisqueRNA 要求）--------
    pseudo_count <- round(bulk_use * pseudo_scale)
    pseudo_count[pseudo_count < 0] <- 0
    
    bulk.eset <- Biobase::ExpressionSet(
      assayData = as.matrix(pseudo_count)
    )
    
    sc.eset <- Biobase::ExpressionSet(
      assayData = as.matrix(sc_use),
      phenoData = Biobase::AnnotatedDataFrame(sc_anno)
    )
    
    ## -------- 3. 核心步骤（BisqueRNA reference）--------
    core_time <- system.time({
      core_mem <- peakRAM::peakRAM({
        
        res <- BisqueRNA::ReferenceBasedDecomposition(
          bulk.eset       = bulk.eset,
          sc.eset         = sc.eset,
          markers         = NULL,
          cell.types      = "cellType",
          subject.names   = "sample",
          use.overlap     = FALSE,
          verbose         = verbose
        )
        
      })
    })
    
    ## -------- 4. 后处理（与你原逻辑一致）--------
    prop_BisqueRNA <- res$bulk.props
    
    ## min–max scaling
    prop_BisqueRNA <- apply(
      prop_BisqueRNA, 2,
      function(x) {
        if (max(x) == min(x)) return(rep(0, length(x)))
        (x - min(x)) / (max(x) - min(x))
      }
    )
    
    ## 列归一化
    prop_BisqueRNA <- sweep(
      prop_BisqueRNA,
      2,
      colSums(prop_BisqueRNA),
      FUN = "/"
    )
    
    prop_BisqueRNA <- as.data.frame(prop_BisqueRNA)
  })
  
  t_total_end <- Sys.time()
  
  ## =========================
  ## 5. 时间 & 内存
  ## =========================
  core_time_sec  <- unname(core_time["elapsed"])
  total_time_sec <- as.numeric(
    difftime(t_total_end, t_total_start, units = "secs")
  )
  
  core_mem_MB  <- core_mem$Peak_RAM_Used_MiB[1]
  total_mem_MB <- total_mem$Peak_RAM_Used_MiB[1]
  
  ## =========================
  ## 6. 保存结果
  ## =========================
  prop_file <- file.path(
    prop_dir,
    paste0("prop_", method_name, "_", bulk_name, ".csv")
  )
  
  utils::write.csv(prop_BisqueRNA, prop_file)
  
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
  
  utils::write.csv(runtime_log, log_file, row.names = FALSE)
  
  ## =========================
  ## 8. 返回
  ## =========================
  invisible(list(
    prop    = prop_BisqueRNA,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}

