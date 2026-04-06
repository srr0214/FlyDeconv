run_ARIC <- function(
    bulk,
    tissue,
    sex = c("mix", "male", "female"),
    bulk_name = NULL,
    out_base = "results"
) {
  
  sex <- match.arg(sex)
  method_name <- "ARIC"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  conda_env     <- "cpdb"
  python_script <- "/project/code/py/run_aric_from_matrix.py"
  
  if (!requireNamespace("processx", quietly = TRUE)) stop("need processx")
  if (!requireNamespace("peakRAM", quietly = TRUE))  stop("need peakRAM")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("need jsonlite")
  
  ## internal data
  ct_expr <- get_internal_data(tissue, sex, "celltype_expr")
  
  ## output dirs
  prop_dir <- file.path(out_base, method_name)
  mid_dir  <- file.path(out_base, method_name, "python_mid", bulk_name)
  
  dir.create(prop_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(mid_dir,  recursive = TRUE, showWarnings = FALSE)
  
  ## temp files
  bulk_file <- tempfile(fileext = ".tsv")
  ref_file  <- tempfile(fileext = ".tsv")
  
  write.table(bulk,    bulk_file, sep="\t", quote=FALSE, col.names=NA)
  write.table(ct_expr, ref_file,  sep="\t", quote=FALSE, col.names=NA)
  
  total_peak <- peakRAM::peakRAM({
    
    total_time <- system.time({
      
      res <- processx::run(
        "conda",
        args = c(
          "run", "-n", conda_env,
          "python", python_script,
          "--bulk", bulk_file,
          "--ref", ref_file,
          "--prefix", paste0(bulk_name, "_"),
          "--outdir", mid_dir
        ),
        echo = TRUE,
        error_on_status = FALSE
      )
      
      if (res$status != 0) {
        stop("ARIC failed:\n",
             paste(res$stderr, collapse = "\n"))
      }
    })
  })
  
  ## read result
  prop_file <- file.path(mid_dir, paste0(bulk_name, "_ARIC_fraction.csv"))
  prop_ARIC <- read.csv(prop_file, row.names = 1, check.names = FALSE)
  prop_ARIC <- as.data.frame(prop_ARIC)
  
  write.csv(
    prop_ARIC,
    file.path(prop_dir, paste0("prop_ARIC_", bulk_name, ".csv"))
  )
  
  ## python benchmark
  bench <- jsonlite::fromJSON(file.path(mid_dir, "ARIC_benchmark.json"))
  
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
    prop = prop_ARIC,
    runtime = bench_row,
    method = method_name
  ))
}
