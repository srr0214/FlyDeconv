run_RNASieve <- function(
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
  method_name <- "RNA-Sieve"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  ## =========================
  ## 1. 固定的内部配置（不对外暴露）
  ## =========================
  ## 依赖检查（非常加分）
  if (!requireNamespace("processx", quietly = TRUE)) {
    stop("RNA-Sieve requires package: processx")
  }
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    stop("RNA-Sieve requires package: peakRAM")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("RNA-Sieve requires package: jsonlite")
  }
  
  ## =========================
  ## 2. internal data
  ## =========================
  cell_expr <- get_internal_data(tissue, sex, "cell_expr")
  ann       <- get_internal_data(tissue, sex, "annotation")
  
  ## =========================
  ## 3. 输出目录
  ## =========================
  prop_dir <- file.path(out_base, method_name)
  mid_dir  <- file.path(out_base, method_name, "python_mid")
  
  dir.create(prop_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(mid_dir,  recursive = TRUE, showWarnings = FALSE)
  
  ## =========================
  ## 4. 临时输入文件
  ## =========================
  sc_file   <- tempfile(fileext = ".tsv")
  ann_file  <- tempfile(fileext = ".tsv")
  bulk_file <- tempfile(fileext = ".tsv")
  
  write.table(cell_expr, sc_file, sep = "\t", quote = FALSE, col.names = NA)
  write.table(ann,       ann_file, sep = "\t", quote = FALSE, col.names = NA)
  write.table(bulk,      bulk_file, sep = "\t", quote = FALSE, col.names = NA)
  
  out_prefix <- file.path(mid_dir, bulk_name)
  
  ## =========================
  ## 5. 运行 + benchmark
  ## =========================
  total_peak <- peakRAM::peakRAM({
    
    total_time <- system.time({
      
      core_time <- system.time({
        
        res <- processx::run(
          "conda",
          args = c(
            "run", "-n", "cpdb",
            "python", "/project/code/py/run_rnasieve.py",
            "--sc_counts", sc_file,
            "--sc_anno", ann_file,
            "--bulk", bulk_file,
            "--celltype_col", "cellType",
            "--out_prefix", out_prefix
          ),
          echo = TRUE,
          error_on_status = FALSE
        )
        
        if (res$status != 0) {
          stop("RNA-Sieve failed:\n", paste(res$stderr, collapse = "\n"))
        }
        
      })
      
      ## =========================
      ## 6. 读取结果
      ## =========================
      prop_file <- paste0(out_prefix, "_cell_fraction.tsv")
      
      prop_RNASieve <- read.table(
        prop_file,
        header = TRUE,
        sep = "\t",
        row.names = 1
      )
      
      prop_RNASieve <- t(prop_RNASieve)
      prop_RNASieve <- sweep(prop_RNASieve, 2, colSums(prop_RNASieve), "/")
      prop_RNASieve <- as.data.frame(prop_RNASieve)
      
      prop_outfile <- file.path(
        prop_dir,
        paste0("prop_", method_name, "_", bulk_name, ".csv")
      )
      write.csv(prop_RNASieve, prop_outfile)
      
    })
    
  })
  
  total_mem_MB <- max(total_peak$Peak_RAM_Used_MiB)
  
  ## =========================
  ## 7. Python core memory
  ## =========================
  bench_json <- paste0(out_prefix, "_benchmark.json")
  core_mem_MB <- if (file.exists(bench_json)) {
    jsonlite::fromJSON(bench_json)$core_mem_MB
  } else {
    NA_real_
  }
  
  ## =========================
  ## 8. benchmark 记录
  ## =========================
  bench_row <- data.frame(
    method         = method_name,
    bulk           = bulk_name,
    core_time_sec  = unname(core_time["elapsed"]),
    total_time_sec = unname(total_time["elapsed"]),
    core_mem_MB    = core_mem_MB,
    total_mem_MB   = total_mem_MB,
    stringsAsFactors = FALSE
  )
  
  bench_file <- file.path(out_base, "runtime_benchmark.csv")
  if (!file.exists(bench_file)) {
    write.csv(bench_row, bench_file, row.names = FALSE)
  } else {
    write.table(
      bench_row, bench_file,
      sep = ",", row.names = FALSE,
      col.names = FALSE, append = TRUE
    )
  }
  
  invisible(list(
    prop    = prop_RNASieve,
    runtime = bench_row,
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}


