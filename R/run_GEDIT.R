run_GEDIT <- function(
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
  method_name <- "GEDIT"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  ## =========================
  ## 1. 固定内部配置
  ## =========================
  conda_env <- "cpdb"

  python_script <- "/project/code/py/GEDIT3.py"
  
  ## =========================
  ## 2. 依赖检查
  ## =========================
  if (!requireNamespace("processx", quietly = TRUE)) {
    stop("GEDIT requires package: processx")
  }
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    stop("GEDIT requires package: peakRAM")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("GEDIT requires package: jsonlite")
  }
  
  ## =========================
  ## 3. internal data
  ## =========================
  ## GEDIT 使用 gene × celltype reference
  celltype_expr <- get_internal_data(tissue, sex, "celltype_expr")
  
  ## =========================
  ## 4. 输出目录（绝对路径）
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
  ## 5. 临时输入文件
  ## =========================
  bulk_file  <- tempfile(fileext = ".tsv")
  scRef_file <- tempfile(fileext = ".tsv")
  
  write.table(
    bulk,
    bulk_file,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  
  write.table(
    celltype_expr,
    scRef_file,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  
  ## GEDIT 输出前缀（必须是绝对路径）
  out_prefix <- normalizePath(
    file.path(prop_dir, bulk_name),
    winslash = "/",
    mustWork = FALSE
  )
  
  ## =========================
  ## 6. 运行 + benchmark
  ## =========================
  total_peak <- peakRAM::peakRAM({
    
    total_time <- system.time({
      
      core_time <- system.time({
        
        message("🧬 [GEDIT] running...")
        message("    workdir : ", python_script)
        message("    out     : ", out_prefix)
        
        res <- processx::run(
          "conda",
          args = c(
            "run", "-n", conda_env,
            "python", python_script,
            "-mix", bulk_file,
            "-ref", scRef_file,
            "-out", out_prefix,
            "-NumSigs", "150",
            "-MinSigs", "100",
            "-SigMethod", "fsDiff",
            "-RowScaling", "0.2"
          ),
          echo = TRUE,
          error_on_status = FALSE
        )
        
        ## ---------- 错误处理（艾玛级） ----------
        if (res$status != 0) {
          stop(
            "😱 GEDIT failed!\n",
            "---- STDOUT ----\n",
            paste(res$stdout, collapse = "\n"),
            "\n---- STDERR ----\n",
            paste(res$stderr, collapse = "\n")
          )
        }
      })
      
      ## =========================
      ## 7. 读取结果
      ## =========================
      prop_file <- paste0(out_prefix, "_CTPredictions.tsv")
      
      if (!file.exists(prop_file)) {
        stop(
          "😱 GEDIT finished but output not found:\n",
          prop_file
        )
      }
      
      prop_GEDIT <- read.table(
        prop_file,
        header = TRUE,
        sep = "\t",
        row.names = 1
      )
      
      ## sample × celltype
      prop_GEDIT <- t(prop_GEDIT)
      prop_GEDIT <- sweep(prop_GEDIT, 2, colSums(prop_GEDIT), "/")
      prop_GEDIT <- as.data.frame(prop_GEDIT)
      
      write.csv(
        prop_GEDIT,
        file.path(
          prop_dir,
          paste0("prop_", method_name, "_", bulk_name, ".csv")
        )
      )
    })
  })
  
  ## =========================
  ## 8. 读取 Python core memory
  ## =========================
  bench_json <- paste0(out_prefix, "_benchmark.json")
  
  core_mem_MB <- if (file.exists(bench_json)) {
    jsonlite::fromJSON(bench_json)$core_mem_MB
  } else {
    NA_real_
  }
  
  ## =========================
  ## 9. benchmark 记录
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
  
  message("✅ [GEDIT] finished successfully")
  
  invisible(list(
    prop    = prop_GEDIT,
    runtime = bench_row,
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}



