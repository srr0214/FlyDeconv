run_AutoGeneS <- function(
    bulk,
    tissue,
    sex = c("mix", "male", "female"),
    bulk_name = NULL,
    out_base = "results"
) {
  
  sex <- match.arg(sex)
  method_name <- "AutoGeneS"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  conda_env     <- "cpdb"
  python_script <- "/project/code/py/run_autogenes_from_matrix.py"
  
  if (!requireNamespace("processx", quietly = TRUE)) stop("need processx")
  if (!requireNamespace("peakRAM", quietly = TRUE))  stop("need peakRAM")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("need jsonlite")
  
  ## internal data
  sc_expr <- get_internal_data(tissue, sex, "cell_expr")
  ann     <- get_internal_data(tissue, sex, "annotation")
  
  ## output dirs
  prop_dir <- file.path(out_base, method_name)
  mid_dir  <- file.path(out_base, method_name, "python_mid", bulk_name)
  
  dir.create(prop_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(mid_dir,  recursive = TRUE, showWarnings = FALSE)
  
  ## temp files
  sc_file   <- tempfile(fileext = ".tsv")
  ann_file  <- tempfile(fileext = ".tsv")
  bulk_file <- tempfile(fileext = ".tsv")
  
  write.table(sc_expr, sc_file, sep="\t", quote=FALSE, col.names=NA)
  write.table(ann,     ann_file, sep="\t", quote=FALSE, col.names=NA)
  write.table(bulk,    bulk_file, sep="\t", quote=FALSE, col.names=NA)
  
  total_peak <- peakRAM::peakRAM({
    
    total_time <- system.time({
      
      res <- processx::run(
        "conda",
        args = c(
          "run", "-n", conda_env,
          "python", python_script,
          "--sc_counts", sc_file,
          "--sc_anno", ann_file,
          "--bulk", bulk_file,
          "--celltype_col", "cellType",
          "--prefix", paste0(bulk_name, "_"),
          "--outdir", mid_dir
        ),
        echo = TRUE,
        error_on_status = FALSE
      )
      
      if (res$status != 0) {
        stop("AutoGeneS failed:\n",
             paste(res$stderr, collapse = "\n"))
      }
    })
  })
  
  ## read result
  prop_file <- file.path(mid_dir, paste0(bulk_name, "_AutoGeneS_fraction.tsv"))
  prop_AutoGeneS <- read.table(prop_file, header=TRUE, sep="\t", row.names=1)
  prop_AutoGeneS <- t(prop_AutoGeneS)
  prop_AutoGeneS <- as.data.frame(prop_AutoGeneS)
  
  write.csv(
    prop_AutoGeneS,
    file.path(prop_dir, paste0("prop_AutoGeneS_", bulk_name, ".csv"))
  )
  
  ## python benchmark
  bench <- jsonlite::fromJSON(file.path(mid_dir, "AutoGeneS_benchmark.json"))
  
  bench_row <- data.frame(
    method = method_name,
    bulk = bulk_name,
    core_time_sec  = bench$core_time_sec,
    total_time_sec = unname(total_time["elapsed"]),
    core_mem_MB    = bench$core_mem_MB,
    total_mem_MB   = max(total_peak$Peak_RAM_Used_MiB),
    stringsAsFactors = FALSE
  )
  
  bench_file <- file.path(out_base, "runtime_benchmark.csv")
  if (!file.exists(bench_file)) {
    write.csv(bench_row, bench_file, row.names=FALSE)
  } else {
    write.table(
      bench_row, bench_file,
      sep=",", row.names=FALSE,
      col.names=FALSE, append=TRUE
    )
  }
  
  invisible(list(
    prop = prop_AutoGeneS,
    runtime = bench_row,
    method = method_name
  ))
}

