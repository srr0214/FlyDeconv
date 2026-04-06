run_CAMTHC <- function(
    bulk,
    tissue,
    sex = c("mix", "male", "female"),
    bulk_name = NULL,
    out_base = "results"
) {
  
  sex <- match.arg(sex)
  method_name <- "CAMTHC"
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  if (!requireNamespace("CAMTHC", quietly = TRUE)) {
    stop("Method CAMTHC requires package 'CAMTHC'.")
  }
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    stop("Method CAMTHC requires package 'peakRAM'.")
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
  
  t_total_start <- Sys.time()
  
  total_mem <- peakRAM::peakRAM({
    
    core_time <- system.time({
      core_mem <- peakRAM::peakRAM({
        
        ## ✅ 关键修正：不用 X=
        rCAM <- CAMTHC::CAM(
          as.matrix(bulk),
          K = 2
        )
        
      })
    })
    
    prop <- t(rCAM@PrepResult@centers)
    colnames(prop) <- colnames(bulk)
    rownames(prop) <- paste0("Latent_", seq_len(nrow(prop)))
    prop <- as.data.frame(prop)
  })
  
  t_total_end <- Sys.time()
  
  runtime_log <- rbind(
    runtime_log,
    data.frame(
      method = method_name,
      bulk   = bulk_name,
      core_time_sec  = unname(core_time["elapsed"]),
      total_time_sec = as.numeric(difftime(
        t_total_end, t_total_start, units = "secs"
      )),
      core_mem_MB    = core_mem$Peak_RAM_Used_MiB[1],
      total_mem_MB   = total_mem$Peak_RAM_Used_MiB[1],
      stringsAsFactors = FALSE
    )
  )
  
  write.csv(runtime_log, log_file, row.names = FALSE)
  
  prop_file <- file.path(
    prop_dir,
    paste0("prop_", method_name, "_", bulk_name, ".csv")
  )
  write.csv(prop, prop_file)
  
  invisible(list(
    prop    = prop,
    runtime = runtime_log[nrow(runtime_log), ],
    method  = method_name,
    tissue  = tissue,
    sex     = sex
  ))
}



