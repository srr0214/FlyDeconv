run_TAPE <- function(
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
  method_name <- "TAPE"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  ## =========================
  ## 1. 固定内部配置
  ## =========================
  conda_env     <- "cpdb"
  python_script <- "/project/code/py/run_tape_inmemory.py"
  
  if (!requireNamespace("processx", quietly = TRUE)) {
    stop("TAPE requires package: processx")
  }
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    stop("TAPE requires package: peakRAM")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("TAPE requires package: jsonlite")
  }
  
  ## =========================
  ## 2. internal data
  ## =========================
  celltype_expr <- get_internal_data(tissue, sex, "celltype_expr")
  
  ## =========================
  ## 3. 输出目录
  ## =========================
  prop_dir <- normalizePath(
    file.path(out_base, method_name),
    winslash = "/",
    mustWork = FALSE
  )
  mid_dir <- normalizePath(
    file.path(out_base, method_name, "python_mid"),
    winslash = "/",
    mustWork = FALSE
  )
  
  dir.create(prop_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(mid_dir,  recursive = TRUE, showWarnings = FALSE)
  
  ## =========================
  ## 4. 临时文件
  ## =========================
  scRef_file <- tempfile(fileext = ".tsv")
  bulk_file  <- tempfile(fileext = ".tsv")
  
  write.table(t(celltype_expr), scRef_file,
              sep = "\t", quote = FALSE, col.names = NA)
  write.table(t(bulk), bulk_file,
              sep = "\t", quote = FALSE, col.names = NA)
  
  out_prefix <- normalizePath(
    file.path(mid_dir, bulk_name),
    winslash = "/",
    mustWork = FALSE
  )
  
  ## =========================
  ## 5. 运行 + benchmark
  ## =========================
  total_peak <- peakRAM::peakRAM({
    
    total_time <- system.time({
      
      core_time <- system.time({
        
        res <- processx::run(
          "conda",
          args = c(
            "run", "-n", conda_env,
            "python", python_script,
            "--scRef", scRef_file,
            "--bulk",  bulk_file,
            "--out_prefix", out_prefix
          ),
          echo = TRUE,
          error_on_status = FALSE
        )
        
        if (res$status != 0) {
          stop(
            "TAPE failed!\n",
            paste(res$stderr, collapse = "\n")
          )
        }
        
      })
      
      ## =========================
      ## 6. 读取比例
      ## =========================
      prop_file <- paste0(out_prefix, "_cell_fraction.tsv")
      
      prop_TAPE <- read.table(
        prop_file,
        header = TRUE,
        sep = "\t",
        row.names = 1
      )
      
      prop_TAPE <- t(prop_TAPE)
      prop_TAPE <- sweep(prop_TAPE, 2, colSums(prop_TAPE), "/")
      prop_TAPE <- as.data.frame(prop_TAPE)
      
      write.csv(
        prop_TAPE,
        file.path(prop_dir, paste0("prop_", method_name, "_", bulk_name, ".csv"))
      )
    })
  })
  
  ## =========================
  ## 7. 读取 Python core memory
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
    total_mem_MB   = max(total_peak$Peak_RAM_Used_MiB),
    stringsAsFactors = FALSE
  )
  
  bench_file <- file.path(out_base, "runtime_benchmark.csv")
  if (!file.exists(bench_file)) {
    write.csv(bench_row, bench_file, row.names = FALSE)
  } else {
    write.table(
      bench_row,
      bench_file,
      sep = ",",
      row.names = FALSE,
      col.names = FALSE,
      append = TRUE
    )
  }
  
  invisible(list(
    prop    = prop_TAPE,
    runtime = bench_row,
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}

