# Internal registry of deconvolution methods
.deconv_registry <- list(
  
  AdRoit_mle = list(
    language = "R",
    runner   = "run_AdRoit_mle",
    requires_ref = TRUE
  ),
  
  AdRoit_nnls = list(
    language = "R",
    runner   = "run_AdRoit_nnls",
    requires_ref = TRUE
  ),
  
  ADTD = list(
    language = "Python",
    runner   = "run_ADTD",
    requires_ref = TRUE
  ),
  
  ARIC = list(
    language = "R",
    runner   = "run_ARIC",
    requires_ref = FALSE
  ),
  
  AutoGeneS = list(
    language = "Python",
    runner   = "run_AutoGeneS",
    requires_ref = TRUE
  ),
  
  BayesPrism = list(
    language = "Python",
    runner   = "run_BayesPrism",
    requires_ref = TRUE
  ),
  
  BayICE = list(
    language = "Python",
    runner   = "run_BayICE",
    requires_ref = TRUE
  ),
  
  BisqueRNA_marker = list(
    language = "Python",
    runner   = "run_BisqueRNA_marker",
    requires_ref = TRUE
  ),
  
  BisqueRNA_refer = list(
    language = "Python",
    runner   = "run_BisqueRNA_refer",
    requires_ref = TRUE
  ),
  
  BLADE = list(
    language = "Python",
    runner   = "run_BLADE",
    requires_ref = TRUE
  ),
  
  CAM3 = list(
    language = "Python",
    runner   = "run_CAM3",
    requires_ref = TRUE
  ),
  
  CAMTHC = list(
    language = "Python",
    runner   = "run_CAMTHC",
    requires_ref = TRUE
  ),
  
  CDSeq = list(
    language = "Python",
    runner   = "run_CDSeq",
    requires_ref = TRUE
  ),
  
  CellDistinguisher = list(
    language = "Python",
    runner   = "run_CellDistinguisher",
    requires_ref = TRUE
  ),
  
  CIBERSORT = list(
    language = "Python",
    runner   = "run_CIBERSORT",
    requires_ref = TRUE
  ),
  
  CPM = list(
    language = "Python",
    runner   = "run_CPM",
    requires_ref = TRUE
  ),
  
  DAISM = list(
    language = "Python",
    runner   = "run_DAISM",
    requires_ref = TRUE
  ),
  
  DCQ = list(
    language = "Python",
    runner   = "run_DCQ",
    requires_ref = TRUE
  ),
  
  deconf = list(
    language = "Python",
    runner   = "run_deconf",
    requires_ref = TRUE
  ),
  
  DeconRNASeq = list(
    language = "Python",
    runner   = "run_DeconRNASeq",
    requires_ref = TRUE
  ),
  
  DeconV = list(
    language = "Python",
    runner   = "run_DeconV",
    requires_ref = TRUE
  ),
  
  DESeq2_L1 = list(
    language = "Python",
    runner   = "run_DESeq2_L1",
    requires_ref = TRUE
  ),
  
  DESeq2_L2 = list(
    language = "Python",
    runner   = "run_DESeq2_L2",
    requires_ref = TRUE
  ),
  
  DSA = list(
    language = "Python",
    runner   = "run_DSA",
    requires_ref = TRUE
  ),
  
  DSSC = list(
    language = "Python",
    runner   = "run_DSSC",
    requires_ref = TRUE
  ),
  
  dtangle = list(
    language = "Python",
    runner   = "run_dtangle",
    requires_ref = TRUE
  ),
  
  DTD = list(
    language = "Python",
    runner   = "run_DTD",
    requires_ref = TRUE
  ),
  
  DWLS = list(
    language = "Python",
    runner   = "run_DWLS",
    requires_ref = TRUE
  ),
  
  EPIC = list(
    language = "Python",
    runner   = "run_EPIC",
    requires_ref = TRUE
  ),
  
  FARDEEP = list(
    language = "Python",
    runner   = "run_FARDEEP",
    requires_ref = TRUE
  ),
  
  GEDIT = list(
    language = "Python",
    runner   = "run_GEDIT",
    requires_ref = TRUE
  ),
  
  GLDADec = list(
    language = "Python",
    runner   = "run_GLDADec",
    requires_ref = TRUE
  ),
  
  GSVA = list(
    language = "Python",
    runner   = "run_GSVA",
    requires_ref = TRUE
  ),
  
  LinDeconSeq = list(
    language = "Python",
    runner   = "run_LinDeconSeq",
    requires_ref = TRUE
  ),
  
  MCPcounter = list(
    language = "Python",
    runner   = "run_MCPcounter",
    requires_ref = TRUE
  ),
  
  MIXTURE = list(
    language = "Python",
    runner   = "run_MIXTURE",
    requires_ref = TRUE
  ),
  
  MOMF_IS = list(
    language = "Python",
    runner   = "run_MOMF_IS",
    requires_ref = TRUE
  ),
  
  MOMF_KL = list(
    language = "Python",
    runner   = "run_MOMF_KL",
    requires_ref = TRUE
  ),
  
  MuSiC = list(
    language = "Python",
    runner   = "run_MuSiC",
    requires_ref = TRUE
  ),
  
  NITUMID = list(
    language = "Python",
    runner   = "run_NITUMID",
    requires_ref = TRUE
  ),
  
  NNLS = list(
    language = "Python",
    runner   = "run_NNLS",
    requires_ref = TRUE
  ),
  
  ols = list(
    language = "Python",
    runner   = "run_ols",
    requires_ref = TRUE
  ),
  
  PREDE = list(
    language = "Python",
    runner   = "run_PREDE",
    requires_ref = TRUE
  ),
  
  ProportionsInAdmixture = list(
    language = "Python",
    runner   = "run_ProportionsInAdmixture",
    requires_ref = TRUE
  ),
  
  qprog = list(
    language = "Python",
    runner   = "run_qprog",
    requires_ref = TRUE
  ),
  
  rls = list(
    language = "Python",
    runner   = "run_rls",
    requires_ref = TRUE
  ),
  
  RNASieve = list(
    language = "Python",
    runner   = "run_RNASieve",
    requires_ref = TRUE
  ),
  
  SCDC = list(
    language = "Python",
    runner   = "run_SCDC",
    requires_ref = TRUE
  ),
  
  SECRET = list(
    language = "Python",
    runner   = "run_SECRET",
    requires_ref = TRUE
  ),
  
  SQUID = list(
    language = "Python",
    runner   = "run_SQUID",
    requires_ref = TRUE
  ),
  
  SVMDECON = list(
    language = "Python",
    runner   = "run_SVMDECON",
    requires_ref = TRUE
  ),
  
  TAPE = list(
    language = "Python",
    runner   = "run_TAPE",
    requires_ref = TRUE
  )
)

get_method_registry <- function() {
  .deconv_registry
}

get_method_info <- function(method) {
  reg <- get_method_registry()
  
  if (!method %in% names(reg)) {
    stop("Unsupported method: ", method, call. = FALSE)
  }
  
  reg[[method]]
}

get_method_runner <- function(method) {
  info <- get_method_info(method)
  runner_name <- info$runner
  
  if (!exists(runner_name, mode = "function", inherits = TRUE)) {
    stop("Runner function not found: ", runner_name, call. = FALSE)
  }
  
  get(runner_name, mode = "function", inherits = TRUE)
}

list_supported_methods <- function() {
  names(get_method_registry())
}