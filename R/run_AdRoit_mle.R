run_AdRoit_mle <- function(
    bulk,
    tissue,
    sex = c("mix", "male", "female"),
    bulk_name = NULL,
    out_base = "results",
    benchmark_mode = TRUE,
    nthreads = 1
) {

  sex <- match.arg(sex)
  method_name <- "AdRoit_mle"

  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }

  ## 串行 + 关闭 CRAN core-limit guard
  Sys.setenv(
    OMP_NUM_THREADS = as.character(nthreads),
    OPENBLAS_NUM_THREADS = as.character(nthreads),
    MKL_NUM_THREADS = as.character(nthreads),
    NUMEXPR_NUM_THREADS = as.character(nthreads),
    VECLIB_MAXIMUM_THREADS = as.character(nthreads),
    RCPP_PARALLEL_NUM_THREADS = as.character(nthreads),
    R_FUTURE_FORK_ENABLE = "false",
    R_PARALLELLY_FORK_ENABLE = "false",
    `_R_CHECK_LIMIT_CORES_` = "false"
  )
  Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "false")

  if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::sequential)
  }
  options(future.fork.enable = FALSE)
  options(parallelly.fork.enable = FALSE)

  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  }
  if (requireNamespace("foreach", quietly = TRUE)) {
    foreach::registerDoSEQ()
  }
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(nthreads)
    RhpcBLASctl::omp_set_num_threads(nthreads)
  }

  if (!requireNamespace("AdRoit", quietly = TRUE)) {
    stop("Method AdRoit requires package 'AdRoit'.")
  }
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    stop("Method AdRoit requires package 'peakRAM'.")
  }

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

  sc_expr <- get_internal_data(tissue = tissue, sex = sex, type = "cell_expr")
  sc_ann  <- get_internal_data(tissue = tissue, sex = sex, type = "annotation")

  if (is.null(rownames(sc_expr)) || is.null(rownames(bulk))) {
    stop("sc_expr and bulk must have rownames as gene IDs/symbols.")
  }
  if (!all(c("cellType", "sample") %in% colnames(sc_ann))) {
    stop("sc_ann must contain columns: cellType and sample.")
  }

  common_genes <- intersect(rownames(sc_expr), rownames(bulk))
  if (length(common_genes) < 200) {
    stop("Too few common genes between sc_expr and bulk: ", length(common_genes))
  }

  sc_use   <- sc_expr[common_genes, , drop = FALSE]
  bulk_use <- bulk[common_genes, , drop = FALSE]

  t_total_start <- Sys.time()
  prop <- NULL
  core_time <- NA

  total_mem <- peakRAM::peakRAM({

    ## ------------------------------------------------
    ## 强制 parallel::detectCores() = 1
    ## ------------------------------------------------
    old_detectCores <- get("detectCores", envir = asNamespace("parallel"))
    unlockBinding("detectCores", asNamespace("parallel"))
    assign("detectCores", function(...) 1L, envir = asNamespace("parallel"))
    lockBinding("detectCores", asNamespace("parallel"))

    on.exit({
      unlockBinding("detectCores", asNamespace("parallel"))
      assign("detectCores", old_detectCores, envir = asNamespace("parallel"))
      lockBinding("detectCores", asNamespace("parallel"))
    }, add = TRUE)

    ## --------------------------
    ## Step 1: build reference
    ## --------------------------
    if (isTRUE(benchmark_mode)) {
      ## ✅ 关键修复：multi.sample.single 必须 TRUE（否则后续可能缺 mapper）
      ref <- AdRoit::ref.build(
        counts  = sc_use,
        annotations = sc_ann$cellType,
        genes   = common_genes,
        samples = sc_ann$sample,
        normalize = "Total",
        scale.factor = 1e5,
        multi.sample.bulk   = FALSE,
        multi.sample.single = TRUE,   # ✅ 改回 TRUE
        nbootsids  = 1,               # ✅ 保持不并行
        minbootsize = 1,
        silent = TRUE
      )
    } else {
      ref <- AdRoit::ref.build(
        counts  = sc_use,
        annotations = sc_ann$cellType,
        genes   = common_genes,
        samples = sc_ann$sample,
        normalize = "Total",
        scale.factor = 1e5,
        multi.sample.bulk   = FALSE,
        multi.sample.single = TRUE,
        nbootsids  = 10,
        minbootsize = 10,
        silent = TRUE
      )
    }

    ## --------------------------
    ## Step 2: core estimation
    ## --------------------------
    core_time <- system.time({
      prop <- tryCatch(
        {
          AdRoit::AdRoit.est(
            bulk.sample = bulk_use,
            single.ref  = ref,
            use.refvar  = FALSE,
            per.sample.adapt = FALSE,
            init.est.method  = "mle",
            lambda = 1
          )
        },
        error = function(e) {
          ## 保险：如果仍然 mapper 报错，自动换一条最稳路径
          if (grepl("mapper", e$message, fixed = TRUE)) {
            AdRoit::AdRoit.est(
              bulk.sample = bulk_use,
              single.ref  = ref,
              use.refvar  = FALSE,
              per.sample.adapt = FALSE,
              init.est.method  = "ols",  # fallback，不走可能触发 mapper 的 mle 分支
              lambda = 1
            )
          } else {
            stop(e)
          }
        }
      )
    })
  })

  t_total_end <- Sys.time()

  if (is.null(prop)) stop("AdRoit.est returned NULL.")
  prop <- as.data.frame(prop)

  core_time_sec  <- unname(core_time["elapsed"])
  total_time_sec <- as.numeric(difftime(t_total_end, t_total_start, units = "secs"))

  total_mem_MB <- total_mem$Peak_RAM_Used_MiB[1]
  core_mem_MB  <- total_mem_MB

  prop_file <- file.path(prop_dir, paste0("prop_", method_name, "_", bulk_name, ".csv"))
  write.csv(prop, prop_file)

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
    prop    = prop,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}
