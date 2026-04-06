run_GLDADec <- function(
    bulk,
    tissue,
    sex = c("mix", "male", "female"),
    bulk_name = NULL,
    out_base = "results"
) {
  
  ## =========================
  ## 0. 参数
  ## =========================
  sex <- match.arg(sex)
  method_name <- "GLDADec"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  conda_env     <- "cpdb"
  python_script <- "/project/code/py/run_gldadec_from_matrix.py"
  
  ## =========================
  ## 1. 依赖
  ## =========================
  if (!requireNamespace("processx", quietly = TRUE))
    stop("need processx")
  if (!requireNamespace("peakRAM", quietly = TRUE))
    stop("need peakRAM")
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("need jsonlite")
  
  ## =========================
  ## 2. internal data
  ## =========================
  marker_list <- get_internal_data(tissue, sex, "marker_list")
  
  ## =========================
  ## 3. 输出目录
  ## =========================
  prop_dir <- file.path(out_base, method_name)
  mid_dir  <- file.path(out_base, method_name, "python_mid", bulk_name)
  
  dir.create(prop_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(mid_dir,  recursive = TRUE, showWarnings = FALSE)
  
  ## =========================
  ## 4. 临时文件
  ## =========================
  bulk_file   <- tempfile(fileext = ".csv")
  marker_file <- tempfile(fileext = ".json")
  
  write.csv(bulk, bulk_file)
  jsonlite::write_json(marker_list, marker_file, auto_unbox = TRUE)
  
  ## =========================
  ## 5. total time + total mem
  ## =========================
  total_peak <- peakRAM::peakRAM({
    
    total_time <- system.time({
      
      res <- processx::run(
        "conda",
        args = c(
          "run", "-n", conda_env,
          "python", python_script,
          "--bulk", bulk_file,
          "--marker_json", marker_file,
          "--prefix", paste0(bulk_name, "_"),
          "--outdir", mid_dir
        ),
        echo = TRUE,
        error_on_status = FALSE
      )
      
      if (res$status != 0) {
        stop("GLDADec failed:\n", paste(res$stderr, collapse = "\n"))
      }
    })
  })
  
  ## =========================
  ## 6. 读取比例（celltype × sample）
  ## =========================
  prop_file <- file.path(mid_dir, paste0(bulk_name, "_GLDADec_fraction.csv"))
  
  prop_GLDADec <- read.table(
    prop_file,
    header = TRUE,
    sep = ",",
    row.names = 1,
    check.names = FALSE
  )
  
  prop_GLDADec <- as.data.frame(prop_GLDADec)
  
  write.csv(
    prop_GLDADec,
    file.path(prop_dir, paste0("prop_GLDADec_", bulk_name, ".csv"))
  )
  
  ## =========================
  ## 7. core benchmark
  ## =========================
  bench <- jsonlite::fromJSON(
    file.path(mid_dir, "GLDADec_benchmark.json")
  )
  
  bench_row <- data.frame(
    method = method_name,
    bulk   = bulk_name,
    core_time_sec  = bench$core_time_sec,
    total_time_sec = unname(total_time["elapsed"]),
    core_mem_MB    = bench$core_mem_MB,
    total_mem_MB   = max(total_peak$Peak_RAM_Used_MiB),
    stringsAsFactors = FALSE
  )
  
  bench_file <- file.path(out_base, "runtime_benchmark.csv")
  
  if (!file.exists(bench_file)) {
    write.csv(bench_row, bench_file, row.names = FALSE)
  } else {
    write.table(
      bench_row, bench_file,
      sep = ",",
      row.names = FALSE,
      col.names = FALSE,
      append = TRUE
    )
  }
  
  invisible(list(
    prop    = prop_GLDADec,
    runtime = bench_row,
    method  = method_name
  ))
}


