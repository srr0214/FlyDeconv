run_CAM3 <- function(
    bulk,
    tissue,
    sex = c("mix", "male", "female"),
    bulk_name = NULL,
    out_base = "results",
    K = NULL,
    dim.rdc = 20,
    thres.low = 0.30,
    thres.high = 0.95,
    top_var_genes = 3000
) {

  sex <- match.arg(sex)
  method_name <- "CAM3"

  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }

  if (!requireNamespace("CAM3", quietly = TRUE))
    stop("Method CAM3 requires package 'CAM3'.")
  if (!requireNamespace("peakRAM", quietly = TRUE))
    stop("Method CAM3 requires package 'peakRAM'.")

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

  t_total_start <- Sys.time()
  prop <- NULL
  core_time <- NA

  total_mem <- peakRAM::peakRAM({

    ## ---------- 1) 基本 QC ----------
    X <- as.matrix(bulk)
    X <- X[rowSums(X) > 0, , drop = FALSE]
    v <- apply(X, 1, stats::var)
    X <- X[v > 0, , drop = FALSE]

    if (!is.null(top_var_genes) && nrow(X) > top_var_genes) {
      ord <- order(v[v > 0], decreasing = TRUE)
      X <- X[names(v[v > 0])[ord[seq_len(top_var_genes)]], , drop = FALSE]
    }

    ## ---------- 2) 自动 K（安全策略） ----------
    ns <- ncol(X)
    if (is.null(K)) {
      # CAM3 对 K 很敏感：限制在 [2, 8]
      K_use <- min(8, max(2, ns - 1))
    } else {
      K_use <- K
    }

    dim_rdc_use <- max(dim.rdc, K_use + 2)

    message(">>> [CAM3] bulk = ", bulk_name,
            " | K = ", K_use,
            " | dim.rdc = ", dim_rdc_use)

    ## ---------- 3) Core: CAM3 ----------
    core_time <- system.time({
      rCAM3 <- CAM3::CAM3Run(
        data       = X,
        K          = K_use,
        dim.rdc    = dim_rdc_use,
        thres.low  = thres.low,
        thres.high = thres.high
      )
    })

    if (is.null(rCAM3) ||
        is.null(rCAM3@PrepResult) ||
        is.null(rCAM3@PrepResult@W)) {
      stop("CAM3 failed: no valid W returned.")
    }

    ## ---------- 4) Extract proportion ----------
    prop <- t(rCAM3@PrepResult@W)

    rownames(prop) <- paste0("Latent_", seq_len(nrow(prop)))
    colnames(prop) <- colnames(X)

    prop <- as.data.frame(prop)
  })

  t_total_end <- Sys.time()

  runtime_log <- rbind(
    runtime_log,
    data.frame(
      method = method_name,
      bulk   = bulk_name,
      core_time_sec  = unname(core_time["elapsed"]),
      total_time_sec = as.numeric(
        difftime(t_total_end, t_total_start, units = "secs")
      ),
      core_mem_MB    = total_mem$Peak_RAM_Used_MiB[1],
      total_mem_MB   = total_mem$Peak_RAM_Used_MiB[1],
      stringsAsFactors = FALSE
    )
  )

  write.csv(runtime_log, log_file, row.names = FALSE)

  prop_file <- file.path(
    prop_dir,
    paste0("latent_", method_name, "_", bulk_name, ".csv")
  )
  write.csv(prop, prop_file)

  invisible(list(
    prop    = prop,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex,
    K       = K_use
  ))
}
