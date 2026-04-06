run_NITUMID <- function(
    bulk,
    tissue,
    sex = c("mix", "male", "female"),
    bulk_name = NULL,
    out_base = "results"
) {

  ## =========================
  ## 0. 依赖检查
  ## =========================
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    stop("Package 'peakRAM' is required.")
  }
  if (!requireNamespace("NITUMID", quietly = TRUE)) {
    stop("Package 'NITUMID' is required.")
  }

  ## =========================
  ## 1. 参数处理
  ## =========================
  sex <- match.arg(sex)
  method_name <- "NITUMID"

  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }

  ## 固定方法参数（benchmark 统一）
  max.itr <- 2000
  tol     <- 1e-5

  ## =========================
  ## 2. 输出路径 & benchmark 表
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
  ## 3. 一些稳健预处理工具
  ## =========================
  .sanitize_Y <- function(Y) {
    Y <- as.matrix(Y)
    storage.mode(Y) <- "double"

    ## 去掉 NA/Inf
    Y[!is.finite(Y)] <- NA_real_
    if (anyNA(Y)) Y[is.na(Y)] <- 0

    ## 去掉全 0 基因
    keep1 <- rowSums(Y) > 0
    Y <- Y[keep1, , drop = FALSE]

    ## 去掉 0 方差基因（关键：避免 scale 产生 NA）
    v <- apply(Y, 1, stats::var)
    v[!is.finite(v)] <- 0
    keep2 <- v > 0
    Y <- Y[keep2, , drop = FALSE]

    if (nrow(Y) < 10) {
      stop("NITUMID: too few genes left after filtering (need > 10).")
    }
    if (ncol(Y) < 2) {
      stop("NITUMID: need at least 2 bulk samples for stable scaling.")
    }
    Y
  }

  .build_A <- function(ref_mean) {
    ## A: genes x celltypes
    A <- apply(ref_mean, 2, function(x) {
      x <- as.numeric(x)
      x[!is.finite(x)] <- NA_real_
      if (all(is.na(x))) return(rep(0, length(x)))

      q_high <- stats::quantile(x, 0.75, na.rm = TRUE, names = FALSE, type = 7)
      q_low  <- stats::quantile(x, 0.25, na.rm = TRUE, names = FALSE, type = 7)

      ## 防御：当 q_high/q_low 都是 NA 或者相等时，避免产生全 NA / 奇怪逻辑
      if (!is.finite(q_high) || !is.finite(q_low)) return(rep(0, length(x)))

      ifelse(
        x >= q_high,  1,
        ifelse(x <= q_low, -1, 0)
      )
    })

    A <- as.matrix(A)
    storage.mode(A) <- "double"
    A[!is.finite(A)] <- 0
    A
  }

  ## =========================
  ## 4. 总时间 & 总内存
  ## =========================
  t_total_start <- Sys.time()

  prop_NITUMID <- NULL
  core_time <- NA
  core_mem <- NULL

  total_mem <- peakRAM::peakRAM({

    ## -------- 1. 获取 celltype 表达矩阵（非 core）--------
    celltype_expr <- get_internal_data(
      tissue = tissue,
      sex    = sex,
      type   = "celltype_expr"
    )

    ## -------- 2. 基因交集（非 core）--------
    common_genes <- intersect(
      rownames(celltype_expr),
      rownames(bulk)
    )

    if (length(common_genes) == 0) {
      stop("No overlapping genes between bulk and celltype reference.")
    }

    Y_raw <- bulk[common_genes, , drop = FALSE]
    ref_mean <- celltype_expr[common_genes, , drop = FALSE]

    ## -------- 3. 清洗 Y（关键）--------
    Y <- .sanitize_Y(Y_raw)

    ## 让 ref_mean 与 Y 同步（否则 A/Y 行不一致会埋雷）
    ref_mean <- ref_mean[rownames(Y), , drop = FALSE]

    ## -------- 4. 构建指导矩阵 A（非 core）--------
    A <- .build_A(ref_mean)

    ## 再过滤掉 A 全 0 的行（关键：避免后面 row_mean / 内部索引出空）
    keepA <- rowSums(A != 0) > 0
    Y <- Y[keepA, , drop = FALSE]
    A <- A[keepA, , drop = FALSE]

    if (nrow(Y) < 10) {
      stop("NITUMID: too few genes left after filtering by A (need > 10).")
    }

    row_index <- seq_len(nrow(A))
    row_mean  <- rowMeans(A != 0)
    num_cell  <- ncol(A)

    ## -------- 5. 核心步骤（NITUMID）--------
    core_time <- system.time({
      core_mem <- peakRAM::peakRAM({

        res <- NITUMID::NITUMID(
          Y = Y,                # 已清洗
          A = A,                # 已清洗
          row_index = row_index,
          if.bulk   = TRUE,
          row_mean  = row_mean,
          num_cell  = num_cell,
          max.itr   = max.itr,
          tol       = tol
        )

      })
    })

    ## -------- 6. 后处理（非 core）--------
    ## res$consistency_table 可能全 NA 或长度为 0：必须兜底
    ct <- res$consistency_table
    if (is.null(ct) || length(ct) == 0) {
      stop("NITUMID returned empty consistency_table.")
    }
    ct2 <- ct
    ct2[!is.finite(ct2)] <- -Inf
    best_i <- which.max(ct2)
    if (!is.finite(ct2[best_i])) {
      stop("NITUMID consistency_table is all NA/Inf; cannot select best solution.")
    }

    ## result 结构也要兜底
    if (is.null(res$result) || length(res$result) < best_i || is.null(res$result[[best_i]]$H)) {
      stop("NITUMID returned no valid H in result.")
    }

    prop_NITUMID <- res$result[[best_i]]$H
    prop_NITUMID <- as.matrix(prop_NITUMID)
    storage.mode(prop_NITUMID) <- "double"
    prop_NITUMID[!is.finite(prop_NITUMID)] <- 0
    prop_NITUMID[prop_NITUMID < 0] <- 0

    ## 列归一化（防 0 列）
    cs <- colSums(prop_NITUMID)
    cs[cs == 0] <- 1
    prop_NITUMID <- sweep(prop_NITUMID, 2, cs, FUN = "/")

    ## 给行列名：行=celltypes, 列=bulk samples（如果维度匹配）
    if (nrow(prop_NITUMID) == ncol(ref_mean)) {
      rownames(prop_NITUMID) <- colnames(ref_mean)
    } else if (is.null(rownames(prop_NITUMID))) {
      rownames(prop_NITUMID) <- paste0("Latent_", seq_len(nrow(prop_NITUMID)))
    }
    colnames(prop_NITUMID) <- colnames(bulk)

    prop_NITUMID <- as.data.frame(prop_NITUMID)
  })

  t_total_end <- Sys.time()

  ## =========================
  ## 5. 时间 & 内存
  ## =========================
  core_time_sec  <- unname(core_time["elapsed"])
  total_time_sec <- as.numeric(difftime(t_total_end, t_total_start, units = "secs"))

  core_mem_MB  <- if (!is.null(core_mem)) core_mem$Peak_RAM_Used_MiB[1] else NA_real_
  total_mem_MB <- total_mem$Peak_RAM_Used_MiB[1]

  ## =========================
  ## 6. 保存 prop 结果
  ## =========================
  prop_file <- file.path(
    prop_dir,
    paste0("prop_", method_name, "_", bulk_name, ".csv")
  )
  utils::write.csv(prop_NITUMID, prop_file)

  ## =========================
  ## 7. 写入 benchmark 表
  ## =========================
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

  ## =========================
  ## 8. 返回
  ## =========================
  invisible(list(
    prop    = prop_NITUMID,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}
