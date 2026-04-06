run_DAISM <- function(
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
  method_name <- "DAISM"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  ## =========================
  ## 1. 固定内部配置
  ## =========================
  conda_env     <- "cpdb"
  python_script <- "/project/code/py/run_daism_from_matrix.py"
  
  ## 依赖检查
  if (!requireNamespace("processx", quietly = TRUE)) {
    stop("DAISM requires package: processx")
  }
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    stop("DAISM requires package: peakRAM")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("DAISM requires package: jsonlite")
  }
  
  ## =========================
  ## 2. internal data
  ## =========================
  ## DAISM 使用：
  ##   sc_expr : gene × cell
  ##   ann     : cell × meta（含 cellType）
  sc_expr <- get_internal_data(tissue, sex, "cell_expr")
  ann     <- get_internal_data(tissue, sex, "annotation")
  
  ## =========================
  ## 3. 输出目录
  ## =========================
  prop_dir <- normalizePath(
    file.path(out_base, method_name),
    winslash = "/", mustWork = FALSE
  )
  
  mid_dir <- normalizePath(
    file.path(out_base, method_name, "python_mid", bulk_name),
    winslash = "/", mustWork = FALSE
  )
  
  dir.create(prop_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(mid_dir,  recursive = TRUE, showWarnings = FALSE)
  
  ## =========================
  ## 4. 临时输入文件
  ## =========================
  sc_file   <- tempfile(fileext = ".tsv")
  ann_file  <- tempfile(fileext = ".tsv")
  bulk_file <- tempfile(fileext = ".tsv")
  
  write.table(sc_expr, sc_file, sep = "\t", quote = FALSE, col.names = NA)
  write.table(ann,     ann_file, sep = "\t", quote = FALSE, col.names = NA)
  write.table(bulk,    bulk_file, sep = "\t", quote = FALSE, col.names = NA)
  
  ## =========================
  ## 5. 运行 + benchmark
  ## =========================
  total_peak <- peakRAM::peakRAM({
    
    total_time <- system.time({
      
      core_time <- system.time({
        
        message("🧬 [DAISM] running...")
        message("    outdir : ", mid_dir)
        
        res <- processx::run(
          "conda",
          args = c(
            "run", "-n", conda_env,
            "python", python_script,
            "--sc_counts", sc_file,
            "--sc_anno",   ann_file,
            "--bulk",      bulk_file,
            "--celltype_col", "cellType",
            "--outdir", mid_dir
          ),
          echo = TRUE,
          error_on_status = FALSE
        )
        
        ## ========= 艾玛：失败直接掀桌 =========
        if (res$status != 0) {
          stop(
            "😱 DAISM failed!\n",
            "---- STDOUT ----\n",
            paste(res$stdout, collapse = "\n"),
            "\n---- STDERR ----\n",
            paste(res$stderr, collapse = "\n")
          )
        }
        
      })  # core_time
      
      ## =========================
      ## 6. 读取结果
      ## =========================
      prop_file <- file.path(
        mid_dir,
        "prediction",
        "output",
        "DAISM_result.txt"
      )
      
      if (!file.exists(prop_file)) {
        stop("😱 DAISM finished but output not found:\n", prop_file)
      }
      
      prop_DAISM <- read.table(
        prop_file,
        header = TRUE,
        sep = "\t",
        row.names = 1,
        check.names = FALSE
      )
      
      ## sample × celltype
      prop_DAISM <- sweep(prop_DAISM, 2, colSums(prop_DAISM), "/")
      prop_DAISM <- as.data.frame(prop_DAISM)
      
      write.csv(
        prop_DAISM,
        file.path(
          prop_dir,
          paste0("prop_", method_name, "_", bulk_name, ".csv")
        )
      )
      
    })  # total_time
    
  })  # peakRAM
  
  ## =========================
  ## 7. Python core benchmark
  ## =========================
  bench_json <- file.path(mid_dir, "DAISM_benchmark.json")
  
  core_mem_MB  <- NA_real_
  core_time_py <- NA_real_
  
  if (file.exists(bench_json)) {
    bench <- jsonlite::fromJSON(bench_json)
    core_mem_MB  <- bench$core_mem_MB
    core_time_py <- bench$core_time_sec
  }
  
  ## =========================
  ## 8. runtime_benchmark 记录
  ## =========================
  bench_row <- data.frame(
    method         = method_name,
    bulk           = bulk_name,
    core_time_sec  = core_time_py,
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
  
  message("✅ [DAISM] finished successfully")
  
  invisible(list(
    prop    = prop_DAISM,
    runtime = bench_row,
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}
