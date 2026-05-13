#' Train Final Single Environment Model for Genomic Prediction
#'
#' @description
#' This function trains a final genomic prediction model using all available
#' phenotypic and marker data. It supports GBLUP, PCA, nonlinear kernels,
#' kernel PCA, and combinations of kernels using the RKHS framework implemented
#' in the BGLR package.
#'
#' @param SNPs A numeric matrix or data frame of SNP genotypes, with individuals in rows and markers in columns.
#' @param y A numeric vector of phenotypic values corresponding to the individuals.
#' @param model A character vector indicating the model or models to fit. Options are
#' "gblup", "pca", "gaussian", "laplacian", "polynomial", "bessel",
#' "pca_gaussian", "pca_laplacian", "pca_polynomial", and "pca_bessel".
#' If more than one model is provided, they are fitted as combined kernels in the BGLR ETA.
#' @param gaussian_sigma Sigma parameter for the Gaussian kernel. Default is 0.001.
#' @param laplacian_sigma Sigma parameter for the Laplacian kernel. Default is 0.01.
#' @param polynomial_degree Degree parameter for the polynomial kernel. Default is 2.
#' @param polynomial_scale Scale parameter for the polynomial kernel. Default is 2.
#' @param polynomial_offset Offset parameter for the polynomial kernel. Default is 2.
#' @param bessel_sigma Sigma parameter for the Bessel kernel. Default is 0.1.
#' @param bessel_order Order parameter for the Bessel kernel. Default is 1.
#' @param bessel_degree Degree parameter for the Bessel kernel. Default is 2.
#' @param var_threshold Minimum proportion of variance explained required for a principal component to be retained in PCA or kernel PCA. Default is 0.01.
#' @param nPC Optional number of principal components to retain. If NULL, the number of PCs is selected using var_threshold.
#' @param scale_pca Logical value indicating whether SNPs should be scaled before PCA. Default is TRUE.
#' @param nIter Total number of iterations for the BGLR model. Default is 10000.
#' @param burnIn Number of burn-in iterations for the BGLR model. Default is 4000.
#' @param thin Thinning interval for the BGLR model. Default is 10.
#' @param save_model Logical value indicating whether to save the fitted model as an RDS file. Default is TRUE.
#' @param file_name Character string specifying the name of the RDS file. Default is "final_bglr_model.rds".
#'
#' @return A list containing the fitted BGLR model, fitted values, residuals,
#' kernel matrices, model information, selected PCs, hyperparameters, and training data.
#'
#' @examples
#' \dontrun{
#' fit <- train_final_single(
#'   SNPs = SNPs,
#'   y = y,
#'   model = c("laplacian", "bessel"),
#'   laplacian_sigma = 0.01,
#'   bessel_sigma = 0.1,
#'   bessel_order = 1,
#'   bessel_degree = 2,
#'   file_name = "final_laplacian_bessel_gxe.rds"
#' )
#' }
#'
#' @export

train_final_single <- function(SNPs,
                             y,
                             model = "gaussian",
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
                             nIter = 10000,
                             burnIn = 4000,
                             thin = 10,
                             save_model = TRUE,
                             file_name = "final_bglr_model.rds") {

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

  if (nrow(SNPs) != length(y)) {
    stop("Number of rows in SNPs must match length of y.")
  }

  if (any(is.na(y))) {
    stop("This function trains the final model using all available data. Remove or impute missing values in y before fitting.")
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

  make_psd <- function(K) {
    K <- as.matrix(K)
    K <- (K + t(K)) / 2

    eig <- eigen(K, symmetric = TRUE)
    eig$values[eig$values < 0] <- 0

    K_psd <- eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
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

    list(
      K = K,
      kpca_object = kpca_obj,
      nPC = nPC_selected
    )
  }

  make_kernel <- function(model_i) {

    if (model_i == "gblup") {
      K <- AGHmatrix::Gmatrix(
        SNPmatrix = SNPs,
        method = "VanRaden",
        ploidy = 2
      )

      return(list(
        K = normalize_kernel(K),
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

      return(list(
        K = normalize_kernel(K),
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

      return(list(
        K = normalize_kernel(K),
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

      return(list(
        K = normalize_kernel(K),
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

      return(list(
        K = normalize_kernel(K),
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

  cat("\nTraining final BGLR model using all available data\n")
  cat("Model:", paste(model, collapse = " + "), "\n")

  kernel_outputs <- lapply(model, make_kernel)
  names(kernel_outputs) <- model

  ETA <- lapply(seq_along(kernel_outputs), function(i) {
    list(
      K = kernel_outputs[[i]]$K,
      model = "RKHS"
    )
  })

  names(ETA) <- model

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
    ETA = ETA,
    kernels = lapply(kernel_outputs, function(x) x$K),
    pca_objects = lapply(kernel_outputs, function(x) x$object),
    nPC = sapply(kernel_outputs, function(x) x$nPC),
    SNPs_train = SNPs,
    y_train = y,
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
      nIter = nIter,
      burnIn = burnIn,
      thin = thin
    ),
    call = match.call()
  )

  class(model_object) <- "train_final_single"

  if (save_model) {
    saveRDS(model_object, file = file_name)
    cat("\nModel saved as:", file_name, "\n")
  }

  return(model_object)
}
