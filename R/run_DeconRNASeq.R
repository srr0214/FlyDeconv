run_DeconRNASeq <- function(
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
  method_name <- "DeconRNASeq"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  ## =========================
  ## 1. 依赖加载（⚠️ 特例：必须 attach）
  ## =========================
  suppressPackageStartupMessages({
    library(DeconRNASeq)
    library(peakRAM)
  })
  
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
  ## 3. reference
  ## =========================
  ref_expr <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "celltype_expr"
  )
  
  ## =========================
  ## 4. 总体计时
  ## =========================
  t_total_start <- Sys.time()
  
  total_mem <- peakRAM({
    
    common_genes <- intersect(rownames(bulk), rownames(ref_expr))
    
    bulk_use <- as.data.frame(bulk[common_genes, , drop = FALSE])
    ref_use  <- as.data.frame(ref_expr[common_genes, , drop = FALSE])
    
    ## ---- core ----
    core_time <- system.time({
      core_mem <- peakRAM({
        res <- DeconRNASeq(
          datasets   = bulk_use,
          signatures = ref_use,
          checksig   = FALSE,
          known.prop = FALSE,
          use.scale  = TRUE,
          fig        = FALSE
        )
      })
    })
    
    ## ---- post ----
    prop <- t(res$out.all)
    prop[prop < 0] <- 0
    prop <- sweep(prop, 2, colSums(prop), "/")
    prop <- as.data.frame(prop)
  })
  
  t_total_end <- Sys.time()
  
  ## =========================
  ## 5. benchmark
  ## =========================
  runtime_log <- rbind(
    runtime_log,
    data.frame(
      method = method_name,
      bulk   = bulk_name,
      core_time_sec  = unname(core_time["elapsed"]),
      total_time_sec = as.numeric(difftime(t_total_end, t_total_start, units = "secs")),
      core_mem_MB    = core_mem$Peak_RAM_Used_MiB[1],
      total_mem_MB   = total_mem$Peak_RAM_Used_MiB[1]
    )
  )
  
  write.csv(runtime_log, log_file, row.names = FALSE)
  
  write.csv(
    prop,
    file.path(prop_dir, paste0("prop_", method_name, "_", bulk_name, ".csv"))
  )

  invisible(list(
    prop    = prop,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name
  ))
}

