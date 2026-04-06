#' Disable parallel execution
#'
#' Force single-threaded execution for benchmarking and reproducibility.
#'
#' @param nthreads Integer. Number of threads to use. Defaults to 1.
#'
#' @return Invisibly returns TRUE.
#' @keywords internal
disable_parallel <- function(nthreads = 1L) {
  nthreads <- as.integer(nthreads)
  
  options(
    mc.cores = nthreads,
    parallelly.maxWorkers.localhost = nthreads,
    parallelly.maxWorkers = nthreads
  )
  
  Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
  Sys.setenv(
    `_R_CHECK_LIMIT_CORES_` = "false",
    OMP_NUM_THREADS = as.character(nthreads),
    OPENBLAS_NUM_THREADS = as.character(nthreads),
    MKL_NUM_THREADS = as.character(nthreads),
    NUMEXPR_NUM_THREADS = as.character(nthreads),
    VECLIB_MAXIMUM_THREADS = as.character(nthreads),
    RCPP_PARALLEL_NUM_THREADS = as.character(nthreads),
    R_FUTURE_FORK_ENABLE = "false",
    R_PARALLELLY_FORK_ENABLE = "false"
  )
  
  if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::sequential)
  }
  
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
  
  invisible(TRUE)
}