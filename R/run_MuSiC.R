run_MuSiC <- function(
    bulk,
    tissue,
    sex = c("mix", "male", "female"),
    bulk_name = NULL,
    out_base = "results"
) {

  sex <- match.arg(sex)
  method_name <- "MuSiC"

  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }

  ## =========================
  ## 1. 依赖检查 + attach
  ## =========================
  if (!requireNamespace("MuSiC", quietly = TRUE))
    stop("MuSiC package is required.")
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE))
    stop("MuSiC requires SingleCellExperiment.")
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
    stop("MuSiC requires SummarizedExperiment.")
  if (!requireNamespace("peakRAM", quietly = TRUE))
    stop("peakRAM is required.")

  suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(SummarizedExperiment)
  })

  ## =========================
  ## 2. internal scRNA data
  ## =========================
  sc_expr <- get_internal_data(tissue, sex, "cell_expr")
  sc_ann  <- get_internal_data(tissue, sex, "annotation")

  sc_ann <- sc_ann[colnames(sc_expr), , drop = FALSE]

  if (!"sample" %in% colnames(sc_ann)) {
    stop("MuSiC requires 'sample' column in annotation.")
  }

  ## =========================
  ## 3. total time + total mem
  ## =========================
  t_total_start <- Sys.time()

  core_time <- NA
  core_mem  <- NA
  prop_MuSiC <- NULL

  total_mem <- peakRAM::peakRAM({

    ## ---------- build SCE ----------
    sce <- SingleCellExperiment(
      assays  = list(counts = as.matrix(sc_expr)),
      colData = sc_ann
    )
    rowData(sce)$gene.name <- rownames(sce)
    sce$sample.id <- sc_ann$sample

    ## ---------- core: MuSiC ----------
    core_mem <- peakRAM::peakRAM({
      core_time <- system.time({
        Est <- MuSiC::music_prop(
          bulk.mtx = as.matrix(bulk),
          sc.sce   = sce,
          clusters = "cellType",
          samples  = "sample.id",
          verbose  = FALSE
        )
      })
    })

    ## ---------- post ----------
    prop_MuSiC <- t(Est$Est.prop.weighted)
    prop_MuSiC[prop_MuSiC < 0] <- 0
    prop_MuSiC <- sweep(prop_MuSiC, 2, colSums(prop_MuSiC), "/")
    prop_MuSiC <- as.data.frame(prop_MuSiC)
  })

  t_total_end <- Sys.time()

  ## =========================
  ## 4. extract metrics
  ## =========================
  core_time_sec  <- unname(core_time["elapsed"])
  total_time_sec <- as.numeric(difftime(t_total_end, t_total_start, units = "secs"))

  core_mem_MB  <- core_mem$Peak_RAM_Used_MiB[1]
  total_mem_MB <- total_mem$Peak_RAM_Used_MiB[1]

  ## =========================
  ## 5. save outputs
  ## =========================
  prop_dir <- file.path(out_base, method_name)
  dir.create(prop_dir, recursive = TRUE, showWarnings = FALSE)

  write.csv(
    prop_MuSiC,
    file.path(prop_dir, paste0("prop_", method_name, "_", bulk_name, ".csv"))
  )

  log_file <- file.path(out_base, "runtime_benchmark.csv")

  runtime_log <- if (file.exists(log_file)) {
    read.csv(log_file, stringsAsFactors = FALSE)
  } else {
    data.frame(
      method = character(),
      bulk   = character(),
      core_time_sec  = numeric(),
      total_time_sec = numeric(),
      core_mem_MB    = numeric(),
      total_mem_MB   = numeric(),
      stringsAsFactors = FALSE
    )
  }

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
    prop    = prop_MuSiC,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}
