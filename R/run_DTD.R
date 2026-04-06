run_DTD <- function(
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
  method_name <- "DTD"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  ## =========================
  ## 1. 依赖检查（不 attach）
  ## =========================
  if (!requireNamespace("DTD", quietly = TRUE)) {
    stop("Method DTD requires package 'DTD'.")
  }
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    stop("Method DTD requires package 'peakRAM'.")
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
    runtime_log <- read.csv(log_file, stringsAsFactors = FALSE)
  }
  
  ## =========================
  ## 3. reference 数据
  ## =========================
  ref_expr <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "celltype_expr"
  )
  
  sc_expr <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "cell_expr"
  )
  
  ann <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "annotation"
  )
  
  ## =========================
  ## 4. 总体流程
  ## =========================
  t_total_start <- Sys.time()
  
  total_mem <- peakRAM::peakRAM({
    
    ## ---------- 基因交集（非 core） ----------
    common_gene <- intersect(
      rownames(sc_expr),
      rownames(bulk)
    )
    
    sc_expr2  <- sc_expr[common_gene, , drop = FALSE]
    bulk_use  <- bulk[common_gene, , drop = FALSE]
    ref_use   <- ref_expr[common_gene, , drop = FALSE]
    
    pheno <- ann[colnames(sc_expr2), "cellType"]
    names(pheno) <- colnames(sc_expr2)
    
    ## ---------- pseudo-bulk 构建（非 core） ----------
    train_data <- DTD::mix_samples(
      expr.data       = sc_expr2,
      pheno           = pheno,
      included.in.X   = unique(pheno),
      n.samples       = 300,
      n.per.mixture   = 100
    )
    
    train_list <- list(
      mixtures   = train_data$mixtures,
      quantities = train_data$quantities
    )
    
    tweak0 <- rep(1, nrow(ref_use))
    names(tweak0) <- rownames(ref_use)
    
    ## =========================
    ## 核心步骤：DTD 训练 + 预测
    ## =========================
    core_time <- system.time({
      core_mem <- peakRAM::peakRAM({
        
        trained <- DTD::train_deconvolution_model(
          tweak              = tweak0,
          X.matrix           = ref_use,
          train.data.list    = train_list,
          estimate.c.type    = "non_negative",
          n_iterations       = 200,
          stepsize           = 0.001
        )
        
        prop <- DTD::estimate_c(
          new.data = bulk_use,
          DTD.model = trained
        )
        
      })
    })
    
    ## ---------- 后处理 ----------
    prop[prop < 0] <- 0
    prop <- sweep(prop, 2, colSums(prop), "/")
    prop <- as.data.frame(prop)
  })
  
  t_total_end <- Sys.time()
  
  ## =========================
  ## 5. 写 benchmark
  ## =========================
  runtime_log <- rbind(
    runtime_log,
    data.frame(
      method = method_name,
      bulk   = bulk_name,
      core_time_sec  = unname(core_time["elapsed"]),
      total_time_sec = as.numeric(
        difftime(t_total_end, t_total_start, units = "secs")
      ),
      core_mem_MB    = core_mem$Peak_RAM_Used_MiB[1],
      total_mem_MB   = total_mem$Peak_RAM_Used_MiB[1],
      stringsAsFactors = FALSE
    )
  )
  
  write.csv(runtime_log, log_file, row.names = FALSE)
  
  ## =========================
  ## 6. 保存比例矩阵
  ## =========================
  prop_file <- file.path(
    prop_dir,
    paste0("prop_", method_name, "_", bulk_name, ".csv")
  )
  write.csv(prop, prop_file)
  
  ## =========================
  ## 7. 返回
  ## =========================
  invisible(list(
    prop    = prop,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}
