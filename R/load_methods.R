load_deconv_methods <- function(
    path = "deconv_function",
    internal_code_dir = "code"
) {

  ## =====================================================
  ## 0. Disable all parallel backends (CRITICAL for Docker)
  ## =====================================================
  if (file.exists("disable_parallel.R")) {
    source("disable_parallel.R", local = .GlobalEnv)
  } else {
    warning("⚠ disable_parallel.R not found; methods may run in parallel")
  }

  ## =====================================================
  ## 1. Load shared internal data utilities (once)
  ## =====================================================
  internal_file <- "get_internal_data.R"

  if (!file.exists(internal_file)) {
    stop("❌ get_internal_data.R not found at: ", internal_file)
  }

  source(internal_file, local = .GlobalEnv)

  ## =====================================================
  ## 2. Load all run_* methods
  ## =====================================================
  files <- list.files(
    path,
    pattern = "^run_.*\\.R$",
    full.names = TRUE
  )

  for (f in files) {
    source(f, local = .GlobalEnv)
  }

  invisible(files)
}
