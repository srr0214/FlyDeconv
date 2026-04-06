run_BayICE <- function(
  bulk,
  tissue,
  sex = c("mix", "male", "female"),
  bulk_name = NULL,
  out_base = "results",
  n_rep = 10,
  cells_per_rep = 50,
  min_rep = 3,
  iter = 1000,
  burnin = 0.6,
  thinning = 3,
  parC = 1e-3,
  min_var = 1e-6
) {

  ## =====================================================
  ## 0. 参数与环境
  ## =====================================================
  sex <- match.arg(sex)
  method_name <- "BayICE"

  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }

  if (!requireNamespace("BayICE", quietly = TRUE)) stop("need BayICE")
  if (!requireNamespace("peakRAM", quietly = TRUE)) stop("need peakRAM")

  set.seed(123)

  ## =====================================================
  ## 1. 读取 internal scRNA 数据
  ## =====================================================
  sc_expr <- get_internal_data(tissue, sex, "cell_expr")
  ann     <- get_internal_data(tissue, sex, "annotation")

  ## =====================================================
  ## 2. pseudo-bulk reference（BayICE-safe）
  ## =====================================================
  make_pseudobulk_ref_safe <- function(
    sc_mat, ann,
    n_rep, cells_per_rep, min_rep,
    cell_id = "cell.id",
    celltype = "cellType"
  ) {

    ann <- ann[ann[[cell_id]] %in% colnames(sc_mat), ]

    ref_list <- list()
    ref_id   <- character()

    for (ct in unique(ann[[celltype]])) {

      cells <- ann[ann[[celltype]] == ct, cell_id]

      if (length(cells) < min_rep * 5) next

      reps <- min(n_rep, floor(length(cells) / 5))
      if (reps < min_rep) next

      for (i in seq_len(reps)) {
        pick <- sample(
          cells,
          size = min(cells_per_rep, length(cells)),
          replace = TRUE
        )

        ref_list[[paste0(ct, "_", i)]] <-
          rowSums(sc_mat[, pick, drop = FALSE])

        ref_id <- c(ref_id, ct)
      }
    }

    list(
      ref.set = do.call(cbind, ref_list),
      ref.id  = ref_id
    )
  }

  pb <- make_pseudobulk_ref_safe(
    sc_mat = as.matrix(sc_expr),
    ann = ann,
    n_rep = n_rep,
    cells_per_rep = cells_per_rep,
    min_rep = min_rep
  )

  if (length(unique(pb$ref.id)) < 2) {
    stop("BayICE aborted: < 2 cell types after pseudo-bulk.")
  }

  ## =====================================================
  ## 3. 基因对齐 + log
  ## =====================================================
  genes <- intersect(rownames(bulk), rownames(pb$ref.set))
  if (length(genes) < 200) {
    stop("BayICE aborted: too few overlapping genes.")
  }

  mix <- as.matrix(bulk[genes, , drop = FALSE])
  ref <- as.matrix(pb$ref.set[genes, , drop = FALSE])

  keep0 <- rowSums(mix) > 0 &
           rowSums(ref) > 0 &
           apply(ref, 1, var) > 0

  mix <- log1p(mix[keep0, , drop = FALSE])
  ref <- log1p(ref[keep0, , drop = FALSE])

  ## =====================================================
  ## 4. cell-type-wise 方差过滤（关键）
  ## =====================================================
  cts <- unique(pb$ref.id)

  keep_gene <- rep(TRUE, nrow(ref))
  names(keep_gene) <- rownames(ref)

  for (ct in cts) {
    idx <- pb$ref.id == ct
    keep_gene <- keep_gene & apply(
      ref[, idx, drop = FALSE],
      1,
      function(x) var(x) > min_var
    )
  }

  ref2 <- ref[keep_gene, , drop = FALSE]
  mix2 <- mix[keep_gene, , drop = FALSE]

  if (nrow(ref2) < 50) {
    stop("BayICE aborted: < 50 genes after variance filtering.")
  }

  ## =====================================================
  ## 5. 核心计算（core）测量
  ## =====================================================
  core_peak <- peakRAM::peakRAM({
    core_time <- system.time({
      res <- BayICE::BayICE(
        ref.set  = ref2,
        mix.set  = mix2,
        ref.id   = pb$ref.id,
        iter     = iter,
        parC     = parC,
        burnin   = burnin,
        thinning = thinning
      )
    })
  })

  ## =====================================================
  ## 6. 总体（total）测量
  ## =====================================================
  total_peak <- peakRAM::peakRAM({
    total_time <- system.time({
      invisible(NULL)  # core 已执行，这里只用于结构统一
    })
  })

  ## =====================================================
  ## 7. 提取比例 + 命名
  ## =====================================================
  prop <- res[["w"]]

  celltypes_BayICE <- unique(pb$ref.id)
  rownames(prop) <- c(celltypes_BayICE, "Unknown")
  colnames(prop) <- colnames(mix2)

  prop <- as.data.frame(prop)

  ## =====================================================
  ## 8. 写出比例结果
  ## =====================================================
  prop_dir <- file.path(out_base, method_name)
  dir.create(prop_dir, recursive = TRUE, showWarnings = FALSE)

  write.csv(
    prop,
    file.path(prop_dir, paste0("prop_BayICE_", bulk_name, ".csv"))
  )

  ## =====================================================
  ## 9. runtime / memory 记录（统一格式）
  ## =====================================================
  bench_row <- data.frame(
    method = method_name,
    bulk = bulk_name,
    core_time_sec  = unname(core_time["elapsed"]),
    total_time_sec = unname(core_time["elapsed"]),  # R-native：core == total
    core_mem_MB    = max(core_peak$Peak_RAM_Used_MiB),
    total_mem_MB   = max(core_peak$Peak_RAM_Used_MiB),
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

  ## =====================================================
  ## 10. 返回（registry 统一接口）
  ## =====================================================
  invisible(list(
    prop = prop,
    runtime = bench_row,
    method = method_name
  ))
}
