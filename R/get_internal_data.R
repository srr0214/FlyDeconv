#' Get internal reference data for FlyDeconv
#'
#' Retrieve built-in reference objects for a given tissue, sex, and data type.
#'
#' @param tissue Character. Tissue name.
#' @param sex Character. One of `"mix"`, `"male"`, or `"female"`.
#' @param type Character. Data type.
#'
#' @return An internal reference object.
#' @export
get_internal_data <- function(
    tissue,
    sex = c("mix", "male", "female"),
    type = c(
      "cell_expr", "annotation",
      "celltype_all", "celltype_female", "celltype_male",
      "list", "long"
    )
) {
  sex <- match.arg(sex)
  type <- match.arg(type)
  
  tissue_key <- tolower(gsub(" ", "_", tissue))
  
  ## =========================
  ## 轻量数据：从 data/ 读取
  ## =========================
  data("all_tissues_scRNA", package = "FlyDeconv", envir = environment())
  
  if (!exists("all_tissues_scRNA", inherits = FALSE)) {
    stop("Internal dataset 'all_tissues_scRNA' is not available.", call. = FALSE)
  }
  
  if (!tissue_key %in% names(all_tissues_scRNA)) {
    stop("Unsupported tissue: ", tissue, call. = FALSE)
  }
  
  tissue_data <- all_tissues_scRNA[[tissue_key]]
  
  ## =========================
  ## 大数据：从 inst/extdata 读取
  ## =========================
  if (type %in% c("cell_expr", "annotation")) {
    big_file <- system.file("extdata", "all_tissues_scRNA_big.rda", package = "FlyDeconv")
    
    if (big_file == "" || !file.exists(big_file)) {
      stop("Large internal dataset file not found: all_tissues_scRNA_big.rda", call. = FALSE)
    }
    
    big_env <- new.env()
    load(big_file, envir = big_env)
    
    if (!exists("all_tissues_scRNA_big", envir = big_env, inherits = FALSE)) {
      stop("Object 'all_tissues_scRNA_big' not found in external dataset file.", call. = FALSE)
    }
    
    big_data <- big_env$all_tissues_scRNA_big
    
    if (!tissue_key %in% names(big_data)) {
      stop("Unsupported tissue in big dataset: ", tissue, call. = FALSE)
    }
    
    tissue_big <- big_data[[tissue_key]]
    
    if (type == "cell_expr") {
      obj_name <- switch(
        sex,
        mix = "cell_all",
        female = "cell_female",
        male = "cell_male"
      )
      
      if (!obj_name %in% names(tissue_big)) {
        stop("Cell expression data '", obj_name, "' is not available for tissue '", tissue, "'.", call. = FALSE)
      }
      
      return(tissue_big[[obj_name]])
    }
    
    if (type == "annotation") {
      obj_name <- switch(
        sex,
        mix = "ann_all",
        female = "ann_female",
        male = "ann_male"
      )
      
      if (!obj_name %in% names(tissue_big)) {
        stop("Annotation data '", obj_name, "' is not available for tissue '", tissue, "'.", call. = FALSE)
      }
      
      return(tissue_big[[obj_name]])
    }
  }
  
  ## =========================
  ## 轻量对象直接返回
  ## =========================
  if (!type %in% names(tissue_data)) {
    stop(
      "Data type '", type, "' is not available for tissue '", tissue, "'.",
      call. = FALSE
    )
  }
  
  tissue_data[[type]]
}