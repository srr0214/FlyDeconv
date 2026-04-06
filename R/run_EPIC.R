run_EPIC <- function(
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
  method_name <- "EPIC"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  ## =========================
  ## 1. 依赖检查（不 library）
  ## =========================
  if (!requireNamespace("EPIC", quietly = TRUE)) {
    stop("Method EPIC requires package 'EPIC'.")
  }
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    stop("Method EPIC requires package 'peakRAM'.")
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
  ## 3. 获取 reference（celltype expression）
  ## =========================
  ref_expr <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "celltype_expr"
  )
  
  ## =========================
  ## 4. 总体时间 & 内存
  ## =========================
  t_total_start <- Sys.time()
  
  total_mem <- peakRAM::peakRAM({
    
    ## ---------- Step 1: 基因交集（非 core） ----------
    common_genes <- intersect(
      rownames(bulk),
      rownames(ref_expr)
    )
    
    Y <- bulk[common_genes, , drop = FALSE]
    refProfiles <- ref_expr[common_genes, , drop = FALSE]
    
    ## EPIC 要求的 reference list
    refList <- list(
      refProfiles = refProfiles,
      sigGenes    = rownames(refProfiles)
    )
    
    ## ---------- Step 2: 核心 EPIC ----------
    core_time <- system.time({
      core_mem <- peakRAM::peakRAM({
        res <- EPIC::EPIC(
          bulk = Y,
          reference = refList,
          scaleExprs = TRUE,
          withOtherCells = TRUE,
          constrainedSum = TRUE
        )
      })
    })
    
    ## ---------- Step 3: 后处理 ----------
    prop_EPIC <- t(res$cellFractions)
    prop_EPIC[prop_EPIC < 0] <- 0
    prop_EPIC <- sweep(
      prop_EPIC, 2, colSums(prop_EPIC), "/"
    )
    prop_EPIC <- as.data.frame(prop_EPIC)
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
  ## 6. 保存比例矩阵
  ## =========================
  prop_file <- file.path(
    prop_dir,
    paste0("prop_", method_name, "_", bulk_name, ".csv")
  )
  write.csv(prop_EPIC, prop_file)
  
  ## =========================
  ## 7. 写 benchmark
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
    prop    = prop_EPIC,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}



