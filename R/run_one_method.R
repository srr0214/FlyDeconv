#' Run a single deconvolution method
#'
#' Convenience wrapper for running exactly one method through
#' the unified FlyDeconv interface.
#'
#' @param bulk A gene expression matrix or data frame, with genes in rows
#'   and samples in columns.
#' @param tissue Character. Tissue name.
#' @param sex Character. One of `"mix"`, `"male"`, or `"female"`.
#' @param method Character. A single method name. Cannot be `"all"`.
#' @param bulk_name Character. Optional bulk dataset name.
#' @param out_base Character. Output directory.
#'
#' @return Method-specific result returned by [flydeconv()].
#' @export
run_one_method <- function(
    bulk,
    tissue,
    sex = c("mix", "male", "female"),
    method,
    bulk_name = NULL,
    out_base = "results"
) {
  sex <- match.arg(sex)
  
  # 必须是单个方法
  if (missing(method) || length(method) != 1L) {
    stop("'method' must be a single method name.", call. = FALSE)
  }
  
  # 禁止 all
  if (identical(method, "all")) {
    stop("'run_one_method()' does not accept method = 'all'.", call. = FALSE)
  }
  
  # 检查方法是否合法
  registry <- get_method_registry()
  if (!method %in% names(registry)) {
    stop("Unknown method: ", method, call. = FALSE)
  }
  
  # 调用统一入口（关键）
  flydeconv(
    bulk = bulk,
    tissue = tissue,
    sex = sex,
    method = method,
    bulk_name = bulk_name,
    out_base = out_base
  )
}