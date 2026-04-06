#' Run deconvolution analysis with FlyDeconv
#'
#' Unified interface for running one or multiple deconvolution methods
#' on bulk transcriptomic data.
#'
#' @param bulk A gene expression matrix or data frame, with genes in rows
#'   and samples in columns.
#' @param tissue Character. Tissue name.
#' @param sex Character. One of `"mix"`, `"male"`, or `"female"`.
#' @param method Character. Method name(s) or `"all"`.
#' @param bulk_name Character. Optional bulk dataset name.
#' @param out_base Character. Output directory.
#' @param show_result Logical. Whether to automatically print results.
#' @param n_show Integer. Number of rows to display for cell proportions.
#'
#' @return A named list of method-specific results.
#' @export
flydeconv <- function(
    bulk,
    tissue,
    sex = c("mix", "male", "female"),
    method = "all",
    bulk_name = NULL,
    out_base = "results",
    show_result = TRUE,
    n_show = 6
) {
  sex <- match.arg(sex)
  
  if (is.null(bulk_name)) {
    bulk_name <- deparse(substitute(bulk))
  }
  
  disable_parallel()
  
  registry <- get_method_registry()
  
  if (identical(method, "all")) {
    methods_to_run <- names(registry)
  } else {
    unknown_methods <- setdiff(method, names(registry))
    if (length(unknown_methods) > 0) {
      stop(
        "Unknown method(s): ",
        paste(unknown_methods, collapse = ", "),
        call. = FALSE
      )
    }
    methods_to_run <- method
  }
  
  message("Running methods: ", paste(methods_to_run, collapse = ", "))
  
  res <- list()
  failed <- character()
  
  for (m in methods_to_run) {
    message("=== [", m, "] ===")
    
    runner <- tryCatch(
      get_method_runner(m),
      error = function(e) {
        message("❌ ", m, " runner not available: ", e$message)
        failed <<- c(failed, m)
        NULL
      }
    )
    
    if (is.null(runner)) {
      res[[m]] <- NULL
      next
    }
    
    out <- tryCatch(
      runner(
        bulk = bulk,
        tissue = tissue,
        sex = sex,
        bulk_name = bulk_name,
        out_base = out_base
      ),
      error = function(e) {
        message("❌ ", m, " failed: ", e$message)
        failed <<- c(failed, m)
        NULL
      }
    )
    
    res[[m]] <- out
    
    ## =========================
    ## 自动展示结果
    ## =========================
    if (isTRUE(show_result) && !is.null(out)) {
      
      cat("\n")
      cat("========================================\n")
      cat("Method:", m, "\n")
      cat("========================================\n")
      
      ## 1. 展示细胞比例
      if (!is.null(out$prop)) {
        cat("\n[Cell proportion result]\n")
        
        prop_show <- out$prop
        if (is.data.frame(prop_show) || is.matrix(prop_show)) {
          nr <- min(n_show, nrow(prop_show))
          print(prop_show[seq_len(nr), , drop = FALSE])
        } else {
          print(prop_show)
        }
      }
      
      ## 2. 展示时间 / 内存
      if (!is.null(out$runtime)) {
        cat("\n[Runtime summary]\n")
        print(out$runtime)
      }
      
      cat("\n")
    }
  }
  
  message("✔ Finished: ", length(methods_to_run) - length(failed))
  if (length(failed) > 0) {
    message("✖ Failed: ", paste(failed, collapse = ", "))
  }
  
  return(res)
}