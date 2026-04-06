run_MOMF_KL <- function(
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
  method_name <- "MOMF_KL"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  ## ---- 固定算法参数（不暴露给总控）----
  seed     <- 123
  rho      <- 2
  num_iter <- 5000
  set.seed(seed)
  
  ## =========================
  ## 1. 输出路径 & benchmark
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
  ## 2. 获取 internal reference
  ## =========================
  celltype_expr <- get_internal_data(
    tissue = tissue,
    sex    = sex,
    type   = "celltype_expr"
  )
  
  ## =========================
  ## 3. 总时间 & 总内存
  ## =========================
  t_total_start <- Sys.time()
  
  total_mem <- peakRAM::peakRAM({
    
    ## ---------- Step 1: gene 交集 ----------
    common_gene <- intersect(
      rownames(bulk),
      rownames(celltype_expr)
    )
    
    bulk_MOMF <- bulk[common_gene, , drop = FALSE]
    ref_MOMF  <- celltype_expr[common_gene, , drop = FALSE]
    
    ## ---------- Step 2: 构造 list ----------
    DataX <- list(as.matrix(bulk_MOMF))
    DataW <- list(as.matrix(ref_MOMF))
    
    ## ---------- Step 3: 初始化 H（固定策略） ----------
    K  <- ncol(ref_MOMF)
    nS <- ncol(bulk_MOMF)
    
    H_init <- matrix(
      stats::runif(K * nS),
      nrow = K,
      ncol = nS
    )
    H_init <- apply(H_init, 2, function(x) x / sum(x))
    DataH  <- list(H_init)
    
    ## ---------- Step 4: 核心 MOMF（KL） ----------
    core_time <- system.time({
      core_mem <- peakRAM::peakRAM({
        
        res <- MOMF::momf.fit(
          DataX = DataX,
          DataW = DataW,
          DataH = DataH,
          DataU = NULL,
          method = "KL",
          rho = rho,
          num_iter = num_iter
        )
        
      })
    })
    
    ## ---------- Step 5: 后处理 ----------
    prop_MOMF <- res[["H"]]
    rownames(prop_MOMF) <- colnames(ref_MOMF)
    colnames(prop_MOMF) <- colnames(bulk_MOMF)
    
    prop_MOMF[prop_MOMF < 0] <- 0
    prop_MOMF <- sweep(prop_MOMF, 2, colSums(prop_MOMF), "/")
    prop_MOMF <- as.data.frame(prop_MOMF)
  })
  
  t_total_end <- Sys.time()
  
  ## =========================
  ## 4. 时间 & 内存
  ## =========================
  core_time_sec  <- unname(core_time["elapsed"])
  total_time_sec <- as.numeric(
    difftime(t_total_end, t_total_start, units = "secs")
  )
  
  core_mem_MB  <- core_mem$Peak_RAM_Used_MiB[1]
  total_mem_MB <- total_mem$Peak_RAM_Used_MiB[1]
  
  ## =========================
  ## 5. 保存比例矩阵
  ## =========================
  prop_file <- file.path(
    prop_dir,
    paste0("prop_", method_name, "_", bulk_name, ".csv")
  )
  utils::write.csv(prop_MOMF, prop_file)
  
  ## =========================
  ## 6. 写 benchmark
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

  invisible(list(
    prop    = prop_MOMF,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}
