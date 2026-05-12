#' Train Final Multi-Environment BGLR Model with GxE Interaction
#'
#' @description
#' This function trains a final multi-environment genomic prediction model using
#' all available phenotypic, environmental, and marker data. It supports GBLUP,
#' PCA, nonlinear kernels, kernel PCA, and combinations of kernels using the RKHS
#' framework implemented in the BGLR package. The model includes fixed
#' environmental effects, main genomic effects, and genotype-by-environment
#' interaction effects.
#'
#' @param SNPs A numeric matrix or data frame of SNP genotypes, with genotypes in rows
#'   and markers in columns. Row names must correspond to genotype IDs.
#' @param y A numeric vector of phenotypic values.
#' @param IDs A character vector indicating the genotype identity for each observation in y.
#' @param env A character vector indicating the environment for each observation in y.
#' @param EZ Optional incidence matrix for fixed environmental effects. If NULL,
#'   it is automatically generated from env.
#' @param model A character vector indicating the model or models to fit. Options are
#'   "gblup", "pca", "gaussian", "laplacian", "polynomial", "bessel",
#'   "pca_gaussian", "pca_laplacian", "pca_polynomial", and "pca_bessel".
#'   If more than one model is provided, they are fitted as combined kernels.
#' @param gaussian_sigma Sigma parameter for the Gaussian kernel. Default is 0.001.
#' @param laplacian_sigma Sigma parameter for the Laplacian kernel. Default is 0.01.
#' @param polynomial_degree Degree parameter for the polynomial kernel. Default is 2.
#' @param polynomial_scale Scale parameter for the polynomial kernel. Default is 2.
#' @param polynomial_offset Offset parameter for the polynomial kernel. Default is 2.
#' @param bessel_sigma Sigma parameter for the Bessel kernel. Default is 0.1.
#' @param bessel_order Order parameter for the Bessel kernel. Default is 1.
#' @param bessel_degree Degree parameter for the Bessel kernel. Default is 2.
#' @param var_threshold Minimum proportion of variance explained required for a principal
#'   component to be retained in PCA or kernel PCA. Default is 0.01.
#' @param nPC Optional number of principal components to retain. If NULL, the number
#'   of PCs is selected using var_threshold.
#' @param scale_pca Logical value indicating whether SNPs should be scaled before PCA.
#'   Default is TRUE.
#' @param ploidy Integer. Ploidy level used in the VanRaden genomic relationship matrix.
#'   Default is 2.
#' @param maf Numeric. Minor allele frequency threshold used by AGHmatrix for GBLUP.
#'   Default is 0.05.
#' @param nIter Total number of iterations for the BGLR Gibbs sampler. Default is 10000.
#' @param burnIn Number of burn-in iterations. Default is 4000.
#' @param thin Thinning interval. Default is 10.
#' @param save_model Logical value indicating whether to save the fitted model as an RDS file.
#'   Default is TRUE.
#' @param file_name Character string specifying the name of the RDS file.
#'   Default is "final_env_gxe_bglr_model.rds".
#'
#' @return A list containing the fitted BGLR model, fitted values, residuals,
#' environmental incidence matrix, genotype incidence matrix, genomic kernels,
#' GxE kernels, eigen decompositions, model information, hyperparameters, and
#' training data.
#'
#' @examples
#' \dontrun{
#' fit <- train_final_multige(
#'   SNPs = SNPs,
#'   y = y,
#'   IDs = IDs,
#'   env = env,
#'   model = c("laplacian", "bessel"),
#'   laplacian_sigma = 0.01,
#'   bessel_sigma = 0.1,
#'   bessel_order = 1,
#'   bessel_degree = 2,
#'   file_name = "final_env_laplacian_bessel_gxe.rds"
#' )
#' }
#'
#' @export

train_final_multige <- function(SNPs,
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
                                bessel_sigma = 0.1,
                                bessel_order = 1,
                                bessel_degree = 2,
                                var_threshold = 0.01,
                                nPC = NULL,
                                scale_pca = TRUE,
                                ploidy = 2,
                                maf = 0.05,
                                nIter = 10000,
                                burnIn = 4000,
                                thin = 10,
                                save_model = TRUE,
                                file_name = "final_env_gxe_bglr_model.rds") {

  if (!requireNamespace("BGLR", quietly = TRUE)) {
    stop("Package 'BGLR' is required.")
  }

  if (!requireNamespace("kernlab", quietly = TRUE)) {
    stop("Package 'kernlab' is required.")
  }

  if (!requireNamespace("AGHmatrix", quietly = TRUE)) {
    stop("Package 'AGHmatrix' is required.")
  }

  SNPs <- as.matrix(SNPs)
  storage.mode(SNPs) <- "numeric"

  y <- as.numeric(y)
  IDs <- as.character(IDs)
  env <- as.character(env)

  if (length(y) != length(IDs) || length(y) != length(env)) {
    stop("The length of y, IDs, and env must be the same.")
  }

  if (any(is.na(y))) {
    stop("This function trains the final model using all available data. Remove or impute missing values in y before fitting.")
  }

  if (is.null(rownames(SNPs))) {
    stop("SNPs must have row names corresponding to genotype IDs.")
  }

  if (!all(unique(IDs) %in% rownames(SNPs))) {
    stop("Some genotype IDs are not present in rownames(SNPs).")
  }

  valid_models <- c(
    "gblup",
    "pca",
    "gaussian",
    "laplacian",
    "polynomial",
    "bessel",
    "pca_gaussian",
    "pca_laplacian",
    "pca_polynomial",
    "pca_bessel"
  )

  if (!all(model %in% valid_models)) {
    stop(
      "Invalid model. Use one or more of: ",
      paste(valid_models, collapse = ", ")
    )
  }

  n <- length(y)

  if (is.null(EZ)) {
    EZ <- model.matrix(~ factor(env) - 1)
    colnames(EZ) <- paste0("Env_", levels(factor(env)))
  }

  EZ <- as.matrix(EZ)

  if (nrow(EZ) != n) {
    stop("EZ must have the same number of rows as the length of y.")
  }

  IDs_factor <- factor(IDs, levels = rownames(SNPs))

  GZ <- model.matrix(~ IDs_factor - 1)
  GZ <- as.matrix(GZ)

  colnames(GZ) <- rownames(SNPs)

  obs_names <- paste0(IDs, "_", env, "_", seq_along(y))

  rownames(EZ) <- obs_names
  rownames(GZ) <- obs_names

  E <- EZ %*% t(EZ)

  rownames(E) <- obs_names
  colnames(E) <- obs_names

  make_psd <- function(K) {

    K <- as.matrix(K)
    K <- (K + t(K)) / 2

    eig <- eigen(K, symmetric = TRUE)
    eig$values[eig$values < 0] <- 0

    K_psd <- eig$vectors %*%
      diag(eig$values, nrow = length(eig$values)) %*%
      t(eig$vectors)

    K_psd <- (K_psd + t(K_psd)) / 2

    return(K_psd)
  }

  normalize_kernel <- function(K) {

    K <- as.matrix(K)
    K <- make_psd(K)

    tr <- sum(diag(K))

    if (is.finite(tr) && tr > 0) {
      K <- K / tr * nrow(K)
    }

    return(K)
  }

  select_npc <- function(eig_values, var_threshold, nPC = NULL) {

    eig_values <- eig_values[eig_values > 0]

    if (length(eig_values) == 0) {
      stop("No positive eigenvalues were found.")
    }

    if (!is.null(nPC)) {

      nPC_selected <- min(nPC, length(eig_values))

    } else {

      var_explained <- eig_values / sum(eig_values)
      nPC_selected <- sum(var_explained > var_threshold)
    }

    if (nPC_selected < 1) {
      nPC_selected <- 1
    }

    return(nPC_selected)
  }

  make_pca_kernel <- function(SNPs, var_threshold, nPC, scale_pca) {

    pca_obj <- stats::prcomp(
      SNPs,
      center = TRUE,
      scale. = scale_pca
    )

    eig_values <- pca_obj$sdev^2
    nPC_selected <- select_npc(eig_values, var_threshold, nPC)

    scores <- pca_obj$x[, seq_len(nPC_selected), drop = FALSE]

    K <- tcrossprod(scores) / ncol(scores)
    K <- normalize_kernel(K)

    rownames(K) <- rownames(SNPs)
    colnames(K) <- rownames(SNPs)

    list(
      K = K,
      pca_object = pca_obj,
      nPC = nPC_selected
    )
  }

  make_kpca_kernel <- function(SNPs,
                               kernel_name,
                               kpar,
                               var_threshold,
                               nPC) {

    kpca_temp <- kernlab::kpca(
      x = SNPs,
      kernel = kernel_name,
      kpar = kpar,
      features = 0
    )

    eig_values <- kernlab::eig(kpca_temp)
    nPC_selected <- select_npc(eig_values, var_threshold, nPC)

    kpca_obj <- kernlab::kpca(
      x = SNPs,
      kernel = kernel_name,
      kpar = kpar,
      features = nPC_selected
    )

    embedding <- as.matrix(kpca_obj@rotated)

    K <- tcrossprod(embedding) / ncol(embedding)
    K <- normalize_kernel(K)

    rownames(K) <- rownames(SNPs)
    colnames(K) <- rownames(SNPs)

    list(
      K = K,
      kpca_object = kpca_obj,
      nPC = nPC_selected
    )
  }

  make_genotype_kernel <- function(model_i) {

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

      return(list(
        K = K,
        method = "gblup",
        object = NULL,
        nPC = NA
      ))
    }

    if (model_i == "gaussian") {

      K <- kernlab::kernelMatrix(
        kernlab::rbfdot(sigma = gaussian_sigma),
        SNPs
      )

      K <- normalize_kernel(K)

      rownames(K) <- rownames(SNPs)
      colnames(K) <- rownames(SNPs)

      return(list(
        K = K,
        method = "gaussian",
        object = NULL,
        nPC = NA
      ))
    }

    if (model_i == "laplacian") {

      K <- kernlab::kernelMatrix(
        kernlab::laplacedot(sigma = laplacian_sigma),
        SNPs
      )

      K <- normalize_kernel(K)

      rownames(K) <- rownames(SNPs)
      colnames(K) <- rownames(SNPs)

      return(list(
        K = K,
        method = "laplacian",
        object = NULL,
        nPC = NA
      ))
    }

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

      return(list(
        K = K,
        method = "polynomial",
        object = NULL,
        nPC = NA
      ))
    }

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

      return(list(
        K = K,
        method = "bessel",
        object = NULL,
        nPC = NA
      ))
    }

    if (model_i == "pca") {

      out <- make_pca_kernel(
        SNPs = SNPs,
        var_threshold = var_threshold,
        nPC = nPC,
        scale_pca = scale_pca
      )

      return(list(
        K = out$K,
        method = "pca",
        object = out$pca_object,
        nPC = out$nPC
      ))
    }

    if (model_i == "pca_gaussian") {

      out <- make_kpca_kernel(
        SNPs = SNPs,
        kernel_name = "rbfdot",
        kpar = list(sigma = gaussian_sigma),
        var_threshold = var_threshold,
        nPC = nPC
      )

      return(list(
        K = out$K,
        method = "pca_gaussian",
        object = out$kpca_object,
        nPC = out$nPC
      ))
    }

    if (model_i == "pca_laplacian") {

      out <- make_kpca_kernel(
        SNPs = SNPs,
        kernel_name = "laplacedot",
        kpar = list(sigma = laplacian_sigma),
        var_threshold = var_threshold,
        nPC = nPC
      )

      return(list(
        K = out$K,
        method = "pca_laplacian",
        object = out$kpca_object,
        nPC = out$nPC
      ))
    }

    if (model_i == "pca_polynomial") {

      out <- make_kpca_kernel(
        SNPs = SNPs,
        kernel_name = "polydot",
        kpar = list(
          degree = polynomial_degree,
          scale = polynomial_scale,
          offset = polynomial_offset
        ),
        var_threshold = var_threshold,
        nPC = nPC
      )

      return(list(
        K = out$K,
        method = "pca_polynomial",
        object = out$kpca_object,
        nPC = out$nPC
      ))
    }

    if (model_i == "pca_bessel") {

      out <- make_kpca_kernel(
        SNPs = SNPs,
        kernel_name = "besseldot",
        kpar = list(
          sigma = bessel_sigma,
          order = bessel_order,
          degree = bessel_degree
        ),
        var_threshold = var_threshold,
        nPC = nPC
      )

      return(list(
        K = out$K,
        method = "pca_bessel",
        object = out$kpca_object,
        nPC = out$nPC
      ))
    }
  }

  expand_kernel <- function(K_genotype) {

    K_obs <- GZ %*% K_genotype %*% t(GZ)

    K_obs <- normalize_kernel(K_obs)

    rownames(K_obs) <- obs_names
    colnames(K_obs) <- obs_names

    return(K_obs)
  }

  make_gxe_kernel <- function(K_obs, E) {

    K_gxe <- K_obs * E

    K_gxe <- normalize_kernel(K_gxe)

    rownames(K_gxe) <- obs_names
    colnames(K_gxe) <- obs_names

    return(K_gxe)
  }

  decompose_kernel <- function(K) {

    eig <- eigen(K, symmetric = TRUE)

    list(
      values = pmax(eig$values, 0),
      vectors = eig$vectors
    )
  }

  cat("\nTraining final multi-environment BGLR model with GxE interaction\n")
  cat("Model:", paste(model, collapse = " + "), "\n")

  genotype_kernel_outputs <- lapply(model, make_genotype_kernel)
  names(genotype_kernel_outputs) <- model

  genotype_kernels <- lapply(genotype_kernel_outputs, function(x) x$K)

  observation_kernels <- lapply(genotype_kernels, expand_kernel)

  gxe_kernels <- lapply(observation_kernels, make_gxe_kernel, E = E)

  genomic_decompositions <- lapply(observation_kernels, decompose_kernel)

  gxe_decompositions <- lapply(gxe_kernels, decompose_kernel)

  ETA <- list(
    Environment = list(
      X = EZ,
      model = "FIXED"
    )
  )

  for (i in seq_along(model)) {

    ETA[[paste0(model[i], "_G")]] <- list(
      V = genomic_decompositions[[i]]$vectors,
      d = genomic_decompositions[[i]]$values,
      model = "RKHS"
    )

    ETA[[paste0(model[i], "_GxE")]] <- list(
      V = gxe_decompositions[[i]]$vectors,
      d = gxe_decompositions[[i]]$values,
      model = "RKHS"
    )
  }

  fit <- BGLR::BGLR(
    y = y,
    ETA = ETA,
    nIter = nIter,
    burnIn = burnIn,
    thin = thin,
    verbose = FALSE
  )

  fitted_values <- fit$yHat
  residuals <- y - fitted_values

  model_object <- list(
    fit = fit,
    yHat = fitted_values,
    residuals = residuals,
    model = model,
    ETA_components = model,
    ETA = ETA,
    EZ = EZ,
    E = E,
    GZ = GZ,
    IDs_train = IDs,
    env_train = env,
    SNPs_train = SNPs,
    y_train = y,
    genotype_kernels = genotype_kernels,
    observation_kernels = observation_kernels,
    gxe_kernels = gxe_kernels,
    genomic_decompositions = genomic_decompositions,
    gxe_decompositions = gxe_decompositions,
    pca_objects = lapply(genotype_kernel_outputs, function(x) x$object),
    nPC = sapply(genotype_kernel_outputs, function(x) x$nPC),
    hyperparameters = list(
      gaussian_sigma = gaussian_sigma,
      laplacian_sigma = laplacian_sigma,
      polynomial_degree = polynomial_degree,
      polynomial_scale = polynomial_scale,
      polynomial_offset = polynomial_offset,
      bessel_sigma = bessel_sigma,
      bessel_order = bessel_order,
      bessel_degree = bessel_degree,
      var_threshold = var_threshold,
      nPC = nPC,
      scale_pca = scale_pca,
      ploidy = ploidy,
      maf = maf,
      nIter = nIter,
      burnIn = burnIn,
      thin = thin
    ),
    call = match.call()
  )

  class(model_object) <- "train_final_env_gxe"

  if (save_model) {
    saveRDS(model_object, file = file_name)
    cat("\nModel saved as:", file_name, "\n")
  }

  return(model_object)
}
