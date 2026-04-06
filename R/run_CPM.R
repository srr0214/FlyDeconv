run_CPM <- function(
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
  method_name <- "CPM"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  ## =========================
  ## 1. 依赖检查（不 attach）
  ## =========================
  if (!requireNamespace("cpm", quietly = TRUE)) {
    stop("Method CPM requires package 'cpm'.")
  }
  if (!requireNamespace("scBio", quietly = TRUE)) {
    stop("Method CPM requires package 'scBio'.")
  }
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    stop("Method CPM requires package 'peakRAM'.")
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
  ## 3. 获取 single-cell 数据
  ## =========================
  sc_expr <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "cell_expr"
  )
  
  ann <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "annotation"
  )
  
  ## =========================
  ## 4. 总体流程
  ## =========================
  t_total_start <- Sys.time()
  
  total_mem <- peakRAM::peakRAM({
    
    ## ---------- 数据准备（非 core） ----------
    SCData   <- log2(sc_expr + 1)
    BulkData <- log2(bulk + 1)
    
    common_cells <- intersect(
      colnames(SCData),
      rownames(ann)
    )
    
    SCData   <- SCData[, common_cells, drop = FALSE]
    SCLabels <- ann[common_cells, "cellType"]
    
    genes <- intersect(
      rownames(SCData),
      rownames(BulkData)
    )
    
    SCData   <- SCData[genes, , drop = FALSE]
    BulkData <- BulkData[genes, , drop = FALSE]
    
    ## ---------- 构建 cellSpace（非 core） ----------
    hvg <- head(
      order(apply(SCData, 1, var), decreasing = TRUE),
      2000
    )
    
    pca <- stats::prcomp(
      t(SCData[hvg, ]),
      center = TRUE,
      scale. = TRUE
    )
    
    cellSpace <- pca$x[, 1:2]
    
    ## =========================
    ## 核心步骤：CPM
    ## =========================
    core_time <- system.time({
      core_mem <- peakRAM::peakRAM({
        
        res <- scBio::CPM(
          SCData,
          SCLabels,
          BulkData,
          cellSpace,
          neighborhoodSize = 5,
          modelSize = 30,
          quantifyTypes = TRUE,
          no_cores = 4
        )
        
      })
    })
    
    ## ---------- 后处理 ----------
    prop <- t(res$cellTypePredictions)
    prop[prop < 0] <- 0
    prop <- sweep(prop, 2, colSums(prop), "/")
    prop <- as.data.frame(prop)
  })
  
  t_total_end <- Sys.time()
  
  ## =========================
  ## 5. 写 benchmark
  ## =========================
  runtime_log <- rbind(
    runtime_log,
    data.frame(
      method = method_name,
      bulk   = bulk_name,
      core_time_sec  = unname(core_time["elapsed"]),
      total_time_sec = as.numeric(
        difftime(t_total_end, t_total_start, units = "secs")
      ),
      core_mem_MB    = core_mem$Peak_RAM_Used_MiB[1],
      total_mem_MB   = total_mem$Peak_RAM_Used_MiB[1],
      stringsAsFactors = FALSE
    )
  )
  
  write.csv(runtime_log, log_file, row.names = FALSE)
  
  ## =========================
  ## 6. 保存比例矩阵
  ## =========================
  prop_file <- file.path(
    prop_dir,
    paste0("prop_", method_name, "_", bulk_name, ".csv")
  )
  write.csv(prop, prop_file)
  
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