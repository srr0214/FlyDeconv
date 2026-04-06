#' Run AdRoit with NNLS initialization
#'
#' Run AdRoit-based deconvolution on bulk transcriptomic data using
#' built-in single-cell reference data, with NNLS as the initial estimator.
#'
#' @param bulk A gene expression matrix or data frame, with genes in rows
#'   and samples in columns.
#' @param tissue Character. Tissue name.
#' @param sex Character. One of `"mix"`, `"male"`, or `"female"`.
#' @param bulk_name Character. Optional bulk dataset name.
#' @param out_base Character. Output directory.
#' @param benchmark_mode Logical. Whether to use benchmark-safe settings.
#' @param nthreads Integer. Number of threads to use. Defaults to 1.
#'
#' @return An invisible list containing proportion estimates, runtime
#'   information, method name, tissue, and sex.
run_AdRoit_nnls <- function(
    bulk,
    tissue,
    sex = c("mix", "male", "female"),
    bulk_name = NULL,
    out_base = "results",
    benchmark_mode = TRUE,
    nthreads = 1
) {
  ## =========================
  ## 0. 参数处理
  ## =========================
  sex <- match.arg(sex)
  method_name <- "AdRoit_nnls"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  bulk <- as.matrix(bulk)
  
  if (is.null(rownames(bulk))) {
    stop("bulk must have rownames.", call. = FALSE)
  }
  
  ## =========================
  ## 0.1 强制串行
  ## =========================
  disable_parallel(nthreads)
  
  ## =========================
  ## 1. 依赖检查
  ## =========================
  if (!requireNamespace("AdRoit", quietly = TRUE)) {
    stop("Method AdRoit requires package 'AdRoit'.", call. = FALSE)
  }
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    stop("Method AdRoit requires package 'peakRAM'.", call. = FALSE)
  }
  
  ## =========================
  ## 2. 输出路径 & benchmark
  ## =========================
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
    runtime_log <- utils::read.csv(log_file, stringsAsFactors = FALSE)
  }
  
  ## =========================
  ## 3. 获取 single-cell 数据
  ## =========================
  sc_expr <- get_internal_data(tissue = tissue, sex = sex, type = "cell_expr")
  sc_ann  <- get_internal_data(tissue = tissue, sex = sex, type = "annotation")
  
  sc_expr <- as.matrix(sc_expr)
  
  if (is.null(rownames(sc_expr))) {
    stop("sc_expr must have rownames.", call. = FALSE)
  }
  if (!all(c("cellType", "sample") %in% colnames(sc_ann))) {
    stop("sc_ann must contain columns: cellType and sample.", call. = FALSE)
  }
  
  ## =========================
  ## 4. 基因交集
  ## =========================
  common_genes <- intersect(rownames(sc_expr), rownames(bulk))
  if (length(common_genes) < 200) {
    stop("Too few common genes: ", length(common_genes), call. = FALSE)
  }
  
  sc_use   <- sc_expr[common_genes, , drop = FALSE]
  bulk_use <- bulk[common_genes, , drop = FALSE]
  
  ## =========================
  ## 5. 计时 + 内存
  ## =========================
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
      ref <- AdRoit::ref.build(
        counts  = sc_use,
        annotations = sc_ann$cellType,
        genes   = common_genes,
        samples = sc_ann$sample,
        normalize = "Total",
        scale.factor = 1e5,
        multi.sample.bulk   = FALSE,
        multi.sample.single = TRUE,
        nbootsids  = 1,
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
            init.est.method  = "nnls",
            lambda = 1
          )
        },
        error = function(e) {
          if (grepl("mapper", e$message, fixed = TRUE)) {
            AdRoit::AdRoit.est(
              bulk.sample = bulk_use,
              single.ref  = ref,
              use.refvar  = FALSE,
              per.sample.adapt = FALSE,
              init.est.method  = "ols",
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
  
  if (is.null(prop)) {
    stop("AdRoit.est returned NULL.", call. = FALSE)
  }
  prop <- as.data.frame(prop)
  
  ## =========================
  ## 6. 时间 & 内存
  ## =========================
  core_time_sec  <- unname(core_time["elapsed"])
  total_time_sec <- as.numeric(difftime(t_total_end, t_total_start, units = "secs"))
  
  total_mem_MB <- total_mem$Peak_RAM_Used_MiB[1]
  core_mem_MB  <- total_mem_MB
  
  ## =========================
  ## 7. 保存结果
  ## =========================
  prop_file <- file.path(
    prop_dir,
    paste0("prop_", method_name, "_", bulk_name, ".csv")
  )
  utils::write.csv(prop, prop_file)
  
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
  utils::write.csv(runtime_log, log_file, row.names = FALSE)
  
  invisible(list(
    prop    = prop,
    runtime = runtime_log[nrow(runtime_log), , drop = FALSE],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}