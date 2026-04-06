run_DeconV <- function(
  bulk,
  tissue,
  sex = c("mix", "male", "female"),
  bulk_name = NULL,
  out_base = "results"
) {

  sex <- match.arg(sex)
  method_name <- "DeconV"

  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }

  if (!requireNamespace("processx", quietly = TRUE))
    stop("need processx")
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("need jsonlite")

  ## =====================================================
  ## 固定配置（与你服务器容器对齐）
  ## =====================================================
  PYTHON <- "/opt/conda/envs/py310/bin/python"
  python_script <- "/project/code/py/run_deconv_from_matrix.py"

  if (!file.exists(PYTHON)) stop("Python not found: ", PYTHON)
  if (!file.exists(python_script)) stop("Python script not found")

  ## =====================================================
  ## internal data
  ## =====================================================
  sc_expr <- get_internal_data(tissue, sex, "cell_expr")
  ann     <- get_internal_data(tissue, sex, "annotation")

  ## =====================================================
  ## 输出目录
  ## =====================================================
  prop_dir <- file.path(out_base, method_name)
  mid_dir  <- file.path(out_base, method_name, "python_mid", bulk_name)

  dir.create(prop_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(mid_dir,  recursive = TRUE, showWarnings = FALSE)

  ## =====================================================
  ## 文件（sc / ann 缓存；bulk 每次写）
  ## =====================================================
  sc_file   <- file.path(mid_dir, "sc_expr.tsv")
  ann_file  <- file.path(mid_dir, "annotation.tsv")
  bulk_file <- file.path(mid_dir, "bulk.tsv")

  if (!file.exists(sc_file)) {
    write.table(sc_expr, sc_file, sep = "\t",
                quote = FALSE, col.names = NA)
    write.table(ann, ann_file, sep = "\t",
                quote = FALSE, col.names = NA)
  }

  write.table(bulk, bulk_file, sep = "\t",
              quote = FALSE, col.names = NA)

  ## =====================================================
  ## stdout / stderr 日志
  ## =====================================================
  stdout_file <- file.path(mid_dir, "deconv_stdout.log")
  stderr_file <- file.path(mid_dir, "deconv_stderr.log")

  ## =====================================================
  ## 运行 DeconV
  ## =====================================================
  total_time <- system.time({

    res <- processx::run(
      command = PYTHON,
      args = c(
        python_script,
        "--bulk", bulk_file,
        "--sc_counts", sc_file,
        "--sc_anno", ann_file,
        "--celltype_col", "cellType",
        "--prefix", paste0(bulk_name, "_"),
        "--outdir", mid_dir
      ),
      echo = FALSE,
      stdout = stdout_file,
      stderr = stderr_file,
      error_on_status = FALSE,
      timeout = 24 * 3600
    )

    if (res$status != 0) {
      stop(
        "DeconV failed.\n",
        "See logs:\n",
        "  ", stdout_file, "\n",
        "  ", stderr_file
      )
    }
  })

  ## =====================================================
  ## 读取结果
  ## =====================================================
  prop_file <- file.path(mid_dir,
                         paste0(bulk_name, "_DeconV_fraction.tsv"))

  if (!file.exists(prop_file)) {
    stop("DeconV output not found. Check logs.")
  }

  prop_DeconV <- read.table(
    prop_file,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE
  )

  write.csv(
    prop_DeconV,
    file.path(prop_dir,
              paste0("prop_DeconV_", bulk_name, ".csv"))
  )

  ## =====================================================
  ## Python benchmark
  ## =====================================================
  bench_json <- file.path(mid_dir, "DeconV_benchmark.json")
  bench <- if (file.exists(bench_json))
    jsonlite::fromJSON(bench_json) else list()

  bench_row <- data.frame(
    method = method_name,
    bulk   = bulk_name,
    core_time_sec  =
      if (!is.null(bench$core_time_sec)) bench$core_time_sec else NA_real_,
    total_time_sec = unname(total_time["elapsed"]),
    core_mem_MB    =
      if (!is.null(bench$core_mem_MB)) bench$core_mem_MB else NA_real_,
    total_mem_MB   = NA_real_,
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
    prop    = prop_DeconV,
    runtime = bench_row,
    method  = method_name
  ))
}
