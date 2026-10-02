#' Train Final Multi-Environment BGLR Model with GxE Interaction
#'
#' @description
#' This function trains a final multi-environment genomic prediction model
#' using all available phenotypic, environmental, and marker data.
#'
#' The model supports GBLUP and nonlinear genomic kernels, including Gaussian,
#' Laplacian, Polynomial, ANOVA, Tanh and Bessel kernels. The model includes
#' fixed environmental effects, main genomic effects, and Genotype by
#' Environment (GxE) interaction effects.
#'
#' For each genomic kernel, the observation-level genomic kernel is obtained
#' by expanding the genotype kernel according to the genotype incidence
#' matrix. The GxE kernel is then calculated as the Hadamard product between
#' the observation-level genomic kernel and the environmental relationship
#' matrix.
#'
#' @param SNPs A numeric matrix or data frame of SNP genotypes, with genotypes
#'   in rows and markers in columns. Row names must correspond to genotype IDs.
#' @param y A numeric vector of phenotypic values.
#' @param IDs A character vector indicating the genotype identity for each
#'   observation in y.
#' @param env A character vector indicating the environment for each
#'   observation in y.
#' @param EZ Optional incidence matrix for fixed environmental effects.
#'   If NULL, it is automatically generated from env.
#' @param model A character vector indicating the genomic kernel or kernels
#'   to fit. Options are "gblup", "gaussian", "laplacian", "polynomial",
#'   "anova", "tanh" and "bessel". If more than one model is provided,
#'   they are fitted as combined kernels.
#' @param gaussian_sigma Sigma parameter for the Gaussian/RBF kernel.
#'   Default is 0.001.
#' @param laplacian_sigma Sigma parameter for the Laplacian kernel.
#'   Default is 0.01.
#' @param polynomial_degree Degree parameter for the Polynomial kernel.
#'   Default is 2.
#' @param polynomial_scale Scale parameter for the Polynomial kernel.
#'   Default is 2.
#' @param polynomial_offset Offset parameter for the Polynomial kernel.
#'   Default is 2.
#' @param anova_sigma Sigma parameter for the ANOVA kernel.
#'   Default is 0.1.
#' @param anova_degree Degree parameter for the ANOVA kernel.
#'   Default is 2.
#' @param tanh_scale Scale parameter for the hyperbolic tangent kernel.
#'   Default is 1.
#' @param tanh_offset Offset parameter for the hyperbolic tangent kernel.
#'   Default is 1.
#' @param bessel_sigma Sigma parameter for the Bessel kernel.
#'   Default is 0.1.
#' @param bessel_order Order parameter for the Bessel kernel.
#'   Default is 1.
#' @param bessel_degree Degree parameter for the Bessel kernel.
#'   Default is 2.
#' @param ploidy Integer. Ploidy level used in the VanRaden genomic
#'   relationship matrix. Default is 2.
#' @param maf Numeric. Minor allele frequency threshold used by AGHmatrix
#'   for GBLUP. Default is 0.05.
#' @param nIter Total number of iterations for the BGLR Gibbs sampler.
#'   Default is 10000.
#' @param burnIn Number of burn-in iterations. Default is 4000.
#' @param thin Thinning interval. Default is 10.
#' @param save_model Logical value indicating whether to save the fitted
#'   model as an RDS file. Default is TRUE.
#' @param file_name Character string specifying the name of the RDS file.
#'   Default is "final_env_gxe_bglr_model.rds".
#'
#' @return A list containing the fitted BGLR model, fitted values, residuals,
#' environmental incidence matrix, genotype incidence matrix, genomic kernels,
#' GxE kernels, eigen decompositions, model information, hyperparameters,
#' and training data.
#'
#' @export

train_final_multige <- function(
    SNPs,
    y,
    IDs,
    env,
    EZ = NULL,
    model = "gblup",

    gaussian_sigma = 0.001,

    laplacian_sigma = 0.01,

    polynomial_degree = 2,
    polynomial_scale = 2,
    polynomial_offset = 2,

    anova_sigma = 0.1,
    anova_degree = 2,

    tanh_scale = 1,
    tanh_offset = 1,

    bessel_sigma = 0.1,
    bessel_order = 1,
    bessel_degree = 2,

    ploidy = 2,
    maf = 0.05,

    nIter = 10000,
    burnIn = 4000,
    thin = 10,

    save_model = TRUE,
    file_name = "final_env_gxe_bglr_model.rds"
) {

  # ============================================================
  # CHECK PACKAGES
  # ============================================================

  if (!requireNamespace("BGLR", quietly = TRUE)) {
    stop("Package 'BGLR' is required.")
  }

  if (!requireNamespace("kernlab", quietly = TRUE)) {
    stop("Package 'kernlab' is required.")
  }

  if (!requireNamespace("AGHmatrix", quietly = TRUE)) {
    stop("Package 'AGHmatrix' is required.")
  }

  # ============================================================
  # DATA
  # ============================================================

  SNPs <- as.matrix(SNPs)
  storage.mode(SNPs) <- "numeric"

  y <- as.numeric(y)
  IDs <- as.character(IDs)
  env <- as.character(env)

  if (length(y) != length(IDs) ||
      length(y) != length(env)) {

    stop(
      "The length of y, IDs, and env must be the same."
    )
  }

  if (any(is.na(y))) {

    stop(
      "This function trains the final model using all available data. ",
      "Remove or impute missing values in y before fitting."
    )
  }

  if (is.null(rownames(SNPs))) {

    stop(
      "SNPs must have row names corresponding to genotype IDs."
    )
  }

  if (!all(unique(IDs) %in% rownames(SNPs))) {

    stop(
      "Some genotype IDs are not present in rownames(SNPs)."
    )
  }

  # ============================================================
  # VALID MODELS
  # ============================================================

  valid_models <- c(
    "gblup",
    "gaussian",
    "laplacian",
    "polynomial",
    "anova",
    "tanh",
    "bessel"
  )

  if (!all(model %in% valid_models)) {

    stop(
      "Invalid model. Use one or more of: ",
      paste(valid_models, collapse = ", ")
    )
  }

  # ============================================================
  # BASIC INFORMATION
  # ============================================================

  n <- length(y)

  uIDs <- unique(IDs)
  uenv <- unique(env)

  # ============================================================
  # ENVIRONMENTAL INCIDENCE MATRIX
  # ============================================================

  if (is.null(EZ)) {

    EZ <- model.matrix(
      ~ factor(env) - 1
    )

    colnames(EZ) <- paste0(
      "Env_",
      levels(factor(env))
    )
  }

  EZ <- as.matrix(EZ)

  if (nrow(EZ) != n) {

    stop(
      "EZ must have the same number of rows as the length of y."
    )
  }

  # ============================================================
  # GENOTYPE INCIDENCE MATRIX
  # ============================================================

  IDs_factor <- factor(
    IDs,
    levels = rownames(SNPs)
  )

  GZ <- model.matrix(
    ~ IDs_factor - 1
  )

  GZ <- as.matrix(GZ)

  colnames(GZ) <- rownames(SNPs)

  # ============================================================
  # OBSERVATION NAMES
  # ============================================================

  obs_names <- paste0(
    IDs,
    "_",
    env,
    "_",
    seq_along(y)
  )

  rownames(EZ) <- obs_names
  rownames(GZ) <- obs_names

  # ============================================================
  # ENVIRONMENTAL RELATIONSHIP MATRIX
  # ============================================================

  E <- EZ %*% t(EZ)

  rownames(E) <- obs_names
  colnames(E) <- obs_names

  # ============================================================
  # MAKE POSITIVE SEMI-DEFINITE
  # ============================================================

  make_psd <- function(K) {

    K <- as.matrix(K)

    K <- (K + t(K)) / 2

    eig <- eigen(
      K,
      symmetric = TRUE
    )

    eig$values[eig$values < 0] <- 0

    K_psd <- eig$vectors %*%
      diag(
        eig$values,
        nrow = length(eig$values)
      ) %*%
      t(eig$vectors)

    K_psd <- (K_psd + t(K_psd)) / 2

    return(K_psd)
  }

  # ============================================================
  # NORMALIZE KERNEL
  # ============================================================

  normalize_kernel <- function(K) {

    K <- as.matrix(K)

    K <- make_psd(K)

    tr <- sum(diag(K))

    if (is.finite(tr) && tr > 0) {

      K <- K / tr * nrow(K)

    }

    return(K)
  }

  # ============================================================
  # MAKE GENOMIC KERNEL
  # ============================================================

  make_genotype_kernel <- function(model_i) {

    # ----------------------------------------------------------
    # GBLUP
    # ----------------------------------------------------------

    if (model_i == "gblup") {

      K <- AGHmatrix::Gmatrix(
        SNPmatrix = SNPs,
        method = "VanRaden",
        ploidy = ploidy,
        maf = maf
      )

      K <- normalize_kernel(K)

      rownames(K) <- rownames(SNPs)
      colnames(K) <- rownames(SNPs)

      return(
        list(
          K = K,
          method = "gblup"
        )
      )
    }

    # ----------------------------------------------------------
    # GAUSSIAN
    # ----------------------------------------------------------

    if (model_i == "gaussian") {

      K <- kernlab::kernelMatrix(
        kernlab::rbfdot(
          sigma = gaussian_sigma
        ),
        SNPs
      )

      K <- normalize_kernel(K)

      rownames(K) <- rownames(SNPs)
      colnames(K) <- rownames(SNPs)

      return(
        list(
          K = K,
          method = "gaussian"
        )
      )
    }

    # ----------------------------------------------------------
    # LAPLACIAN
    # ----------------------------------------------------------

    if (model_i == "laplacian") {

      K <- kernlab::kernelMatrix(
        kernlab::laplacedot(
          sigma = laplacian_sigma
        ),
        SNPs
      )

      K <- normalize_kernel(K)

      rownames(K) <- rownames(SNPs)
      colnames(K) <- rownames(SNPs)

      return(
        list(
          K = K,
          method = "laplacian"
        )
      )
    }

    # ----------------------------------------------------------
    # POLYNOMIAL
    # ----------------------------------------------------------

    if (model_i == "polynomial") {

      K <- kernlab::kernelMatrix(
        kernlab::polydot(
          degree = polynomial_degree,
          scale = polynomial_scale,
          offset = polynomial_offset
        ),
        SNPs
      )

      K <- normalize_kernel(K)

      rownames(K) <- rownames(SNPs)
      colnames(K) <- rownames(SNPs)

      return(
        list(
          K = K,
          method = "polynomial"
        )
      )
    }

    # ----------------------------------------------------------
    # ANOVA
    # ----------------------------------------------------------

    if (model_i == "anova") {

      K <- kernlab::kernelMatrix(
        kernlab::anovadot(
          sigma = anova_sigma,
          degree = anova_degree
        ),
        SNPs
      )

      K <- normalize_kernel(K)

      rownames(K) <- rownames(SNPs)
      colnames(K) <- rownames(SNPs)

      return(
        list(
          K = K,
          method = "anova"
        )
      )
    }

    # ----------------------------------------------------------
    # TANH
    # ----------------------------------------------------------

    if (model_i == "tanh") {

      K <- kernlab::kernelMatrix(
        kernlab::tanhdot(
          scale = tanh_scale,
          offset = tanh_offset
        ),
        SNPs
      )

      K <- normalize_kernel(K)

      rownames(K) <- rownames(SNPs)
      colnames(K) <- rownames(SNPs)

      return(
        list(
          K = K,
          method = "tanh"
        )
      )
    }

    # ----------------------------------------------------------
    # BESSEL
    # ----------------------------------------------------------

    if (model_i == "bessel") {

      K <- kernlab::kernelMatrix(
        kernlab::besseldot(
          sigma = bessel_sigma,
          order = bessel_order,
          degree = bessel_degree
        ),
        SNPs
      )

      K <- normalize_kernel(K)

      rownames(K) <- rownames(SNPs)
      colnames(K) <- rownames(SNPs)

      return(
        list(
          K = K,
          method = "bessel"
        )
      )
    }
  }

  # ============================================================
  # EXPAND GENOTYPE KERNEL TO OBSERVATIONS
  # ============================================================

  expand_kernel <- function(K_genotype) {

    K_obs <- GZ %*%
      K_genotype %*%
      t(GZ)

    K_obs <- normalize_kernel(K_obs)

    rownames(K_obs) <- obs_names
    colnames(K_obs) <- obs_names

    return(K_obs)
  }

  # ============================================================
  # GxE KERNEL
  # ============================================================

  make_gxe_kernel <- function(
    K_obs,
    E) {

    K_gxe <- K_obs * E

    K_gxe <- normalize_kernel(
      K_gxe
    )

    rownames(K_gxe) <- obs_names
    colnames(K_gxe) <- obs_names

    return(K_gxe)
  }

  # ============================================================
  # EIGEN DECOMPOSITION
  # ============================================================

  decompose_kernel <- function(K) {

    eig <- eigen(
      K,
      symmetric = TRUE
    )

    list(
      values = pmax(
        eig$values,
        0
      ),
      vectors = eig$vectors
    )
  }

  # ============================================================
  # INFORMATION
  # ============================================================

  cat(
    "\n============================================\n"
  )

  cat(
    "FINAL MULTI-ENVIRONMENT BGLR MODEL\n"
  )

  cat(
    "============================================\n"
  )

  cat(
    "Models:",
    paste(
      model,
      collapse = " + "
    ),
    "\n"
  )

  cat(
    "Number of observations:",
    n,
    "\n"
  )

  cat(
    "Number of genotypes:",
    length(uIDs),
    "\n"
  )

  cat(
    "Number of environments:",
    length(uenv),
    "\n"
  )

  cat(
    "============================================\n"
  )

  # ============================================================
  # COMPUTE GENOMIC KERNELS
  # ============================================================

  genotype_kernel_outputs <- lapply(
    model,
    make_genotype_kernel
  )

  names(
    genotype_kernel_outputs
  ) <- model

  genotype_kernels <- lapply(
    genotype_kernel_outputs,
    function(x) x$K
  )

  # ============================================================
  # EXPAND TO OBSERVATION LEVEL
  # ============================================================

  observation_kernels <- lapply(
    genotype_kernels,
    expand_kernel
  )

  # ============================================================
  # COMPUTE GxE KERNELS
  # ============================================================

  gxe_kernels <- lapply(
    observation_kernels,
    make_gxe_kernel,
    E = E
  )

  # ============================================================
  # DECOMPOSE GENOMIC KERNELS
  # ============================================================

  genomic_decompositions <- lapply(
    observation_kernels,
    decompose_kernel
  )

  # ============================================================
  # DECOMPOSE GxE KERNELS
  # ============================================================

  gxe_decompositions <- lapply(
    gxe_kernels,
    decompose_kernel
  )

  # ============================================================
  # BGLR ETA
  # ============================================================

  ETA <- list(

    Environment = list(
      X = EZ,
      model = "FIXED"
    )
  )

  # ============================================================
  # ADD GENOMIC AND GxE COMPONENTS
  # ============================================================

  for (i in seq_along(model)) {

    # ----------------------------------------------------------
    # MAIN GENOMIC EFFECT
    # ----------------------------------------------------------

    ETA[
      [paste0(
        model[i],
        "_G"
      )]
    ] <- list(

      V = genomic_decompositions[[i]]$vectors,

      d = genomic_decompositions[[i]]$values,

      model = "RKHS"
    )

    # ----------------------------------------------------------
    # GxE EFFECT
    # ----------------------------------------------------------

    ETA[
      [paste0(
        model[i],
        "_GxE"
      )]
    ] <- list(

      V = gxe_decompositions[[i]]$vectors,

      d = gxe_decompositions[[i]]$values,

      model = "RKHS"
    )
  }

  # ============================================================
  # FIT FINAL BGLR MODEL
  # ============================================================

  cat(
    "\nFitting BGLR model...\n"
  )

  fit <- BGLR::BGLR(

    y = y,

    ETA = ETA,

    nIter = nIter,

    burnIn = burnIn,

    thin = thin,

    verbose = FALSE
  )

  # ============================================================
  # FITTED VALUES
  # ============================================================

  fitted_values <- fit$yHat

  residuals <- y -
    fitted_values

  # ============================================================
  # MODEL OBJECT
  # ============================================================

  model_object <- list(

    # ----------------------------------------------------------
    # BGLR
    # ----------------------------------------------------------

    fit = fit,

    yHat = fitted_values,

    residuals = residuals,

    # ----------------------------------------------------------
    # MODEL INFORMATION
    # ----------------------------------------------------------

    model = model,

    ETA_components = model,

    ETA = ETA,

    # ----------------------------------------------------------
    # ENVIRONMENT
    # ----------------------------------------------------------

    EZ = EZ,

    E = E,

    # ----------------------------------------------------------
    # GENOTYPE
    # ----------------------------------------------------------

    GZ = GZ,

    # ----------------------------------------------------------
    # TRAINING DATA
    # ----------------------------------------------------------

    IDs_train = IDs,

    env_train = env,

    SNPs_train = SNPs,

    y_train = y,

    # ----------------------------------------------------------
    # GENOMIC KERNELS
    # ----------------------------------------------------------

    genotype_kernels = genotype_kernels,

    observation_kernels = observation_kernels,

    gxe_kernels = gxe_kernels,

    # ----------------------------------------------------------
    # EIGEN DECOMPOSITIONS
    # ----------------------------------------------------------

    genomic_decompositions =
      genomic_decompositions,

    gxe_decompositions =
      gxe_decompositions,

    # ----------------------------------------------------------
    # HYPERPARAMETERS
    # ----------------------------------------------------------

    hyperparameters = list(

      gaussian_sigma =
        gaussian_sigma,

      laplacian_sigma =
        laplacian_sigma,

      polynomial_degree =
        polynomial_degree,

      polynomial_scale =
        polynomial_scale,

      polynomial_offset =
        polynomial_offset,

      anova_sigma =
        anova_sigma,

      anova_degree =
        anova_degree,

      tanh_scale =
        tanh_scale,

      tanh_offset =
        tanh_offset,

      bessel_sigma =
        bessel_sigma,

      bessel_order =
        bessel_order,

      bessel_degree =
        bessel_degree,

      ploidy =
        ploidy,

      maf =
        maf,

      nIter =
        nIter,

      burnIn =
        burnIn,

      thin =
        thin
    ),

    # ----------------------------------------------------------
    # CALL
    # ----------------------------------------------------------

    call = match.call()
  )

  # ============================================================
  # CLASS
  # ============================================================

  class(model_object) <-
    "train_final_env_gxe"

  # ============================================================
  # SAVE
  # ============================================================

  if (save_model) {

    saveRDS(
      model_object,
      file = file_name
    )

    cat(
      "\nModel saved as:",
      file_name,
      "\n"
    )
  }

  # ============================================================
  # RETURN
  # ============================================================

  return(
    model_object
  )
}
