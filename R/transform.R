#' Transform SNP Matrix for Genomic Prediction
#'
#' @description
#' This function transforms a SNP marker matrix using one or more genomic
#' relationship, nonlinear kernel, PCA, or kernel PCA approaches. It does not
#' require phenotypic values and does not fit a BGLR model.
#'
#' @param SNPs A numeric matrix or data frame of SNP genotypes, with individuals in rows and markers in columns.
#' @param model Character vector indicating the transformation method or methods. Options are
#' "gblup", "pca", "gaussian", "laplacian", "polynomial", "bessel",
#' "pca_gaussian", "pca_laplacian", "pca_polynomial", and "pca_bessel".
#' If more than one model is provided, the function returns a separate matrix for each model.
#' @param gaussian_sigma Sigma parameter for the Gaussian kernel. Default is 0.001.
#' @param laplacian_sigma Sigma parameter for the Laplacian kernel. Default is 0.01.
#' @param polynomial_degree Degree parameter for the polynomial kernel. Default is 2.
#' @param polynomial_scale Scale parameter for the polynomial kernel. Default is 2.
#' @param polynomial_offset Offset parameter for the polynomial kernel. Default is 2.
#' @param bessel_sigma Sigma parameter for the Bessel kernel. Default is 0.1.
#' @param bessel_order Order parameter for the Bessel kernel. Default is 1.
#' @param bessel_degree Degree parameter for the Bessel kernel. Default is 2.
#' @param var_threshold Minimum proportion of variance explained required for retaining PCs. Default is 0.01.
#' @param nPC Optional number of principal components to retain. If NULL, PCs are selected using var_threshold.
#' @param scale_pca Logical value indicating whether SNPs should be scaled before PCA. Default is TRUE.
#' @param save_transform Logical value indicating whether to save the transformed object as an RDS file. Default is FALSE.
#' @param file_name Character string specifying the name of the RDS file. Default is "snp_transform.rds".
#'
#' @return A list containing the transformed SNP matrices or kernel matrices, selected model(s),
#' PCA/kernel PCA objects when applicable, selected PCs, hyperparameters, and original SNP matrix.
#'
#' @examples
#' \dontrun{
#' obj <- transform_snps(
#'   SNPs = SNPs,
#'   model = c("laplacian", "bessel"),
#'   laplacian_sigma = 0.01,
#'   bessel_sigma = 0.1,
#'   bessel_order = 1,
#'   bessel_degree = 2
#' )
#'
#' K_laplacian <- obj$K$laplacian
#' K_bessel <- obj$K$bessel
#' }
#'
#' @export

transform_snps <- function(SNPs,
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
                           save_transform = FALSE,
                           file_name = "snp_transform.rds") {

  if (!requireNamespace("kernlab", quietly = TRUE)) {
    stop("Package 'kernlab' is required.")
  }

  if (!requireNamespace("AGHmatrix", quietly = TRUE)) {
    stop("Package 'AGHmatrix' is required.")
  }

  SNPs <- as.matrix(SNPs)
  storage.mode(SNPs) <- "numeric"

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

  make_pca_transform <- function(SNPs, var_threshold, nPC, scale_pca) {
    pca_obj <- stats::prcomp(
      SNPs,
      center = TRUE,
      scale. = scale_pca
    )

    eig_values <- pca_obj$sdev^2
    nPC_selected <- select_npc(
      eig_values = eig_values,
      var_threshold = var_threshold,
      nPC = nPC
    )

    scores <- pca_obj$x[, seq_len(nPC_selected), drop = FALSE]

    K <- tcrossprod(scores) / ncol(scores)
    K <- normalize_kernel(K)

    return(list(
      transformed = scores,
      K = K,
      object = pca_obj,
      nPC = nPC_selected
    ))
  }

  make_kpca_transform <- function(SNPs,
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

    nPC_selected <- select_npc(
      eig_values = eig_values,
      var_threshold = var_threshold,
      nPC = nPC
    )

    kpca_obj <- kernlab::kpca(
      x = SNPs,
      kernel = kernel_name,
      kpar = kpar,
      features = nPC_selected
    )

    embedding <- as.matrix(kpca_obj@rotated)

    K <- tcrossprod(embedding) / ncol(embedding)
    K <- normalize_kernel(K)

    return(list(
      transformed = embedding,
      K = K,
      object = kpca_obj,
      nPC = nPC_selected
    ))
  }

  make_transform <- function(model_i) {

    if (model_i == "gblup") {
      K <- AGHmatrix::Gmatrix(
        SNPmatrix = SNPs,
        method = "VanRaden",
        ploidy = 2
      )

      K <- normalize_kernel(K)

      return(list(
        transformed = K,
        K = K,
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

      return(list(
        transformed = K,
        K = K,
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

      return(list(
        transformed = K,
        K = K,
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

      return(list(
        transformed = K,
        K = K,
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

      return(list(
        transformed = K,
        K = K,
        object = NULL,
        nPC = NA
      ))
    }

    if (model_i == "pca") {
      return(
        make_pca_transform(
          SNPs = SNPs,
          var_threshold = var_threshold,
          nPC = nPC,
          scale_pca = scale_pca
        )
      )
    }

    if (model_i == "pca_gaussian") {
      return(
        make_kpca_transform(
          SNPs = SNPs,
          kernel_name = "rbfdot",
          kpar = list(sigma = gaussian_sigma),
          var_threshold = var_threshold,
          nPC = nPC
        )
      )
    }

    if (model_i == "pca_laplacian") {
      return(
        make_kpca_transform(
          SNPs = SNPs,
          kernel_name = "laplacedot",
          kpar = list(sigma = laplacian_sigma),
          var_threshold = var_threshold,
          nPC = nPC
        )
      )
    }

    if (model_i == "pca_polynomial") {
      return(
        make_kpca_transform(
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
      )
    }

    if (model_i == "pca_bessel") {
      return(
        make_kpca_transform(
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
      )
    }
  }

  transform_outputs <- lapply(model, make_transform)
  names(transform_outputs) <- model

  transform_object <- list(
    transformed = lapply(transform_outputs, function(x) x$transformed),
    K = lapply(transform_outputs, function(x) x$K),
    model = model,
    objects = lapply(transform_outputs, function(x) x$object),
    nPC = sapply(transform_outputs, function(x) x$nPC),
    SNPs_original = SNPs,
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
      scale_pca = scale_pca
    ),
    call = match.call()
  )

  class(transform_object) <- "snp_transform"

  if (save_transform) {
    saveRDS(transform_object, file = file_name)
    cat("\nTransformed SNP object saved as:", file_name, "\n")
  }

  return(transform_object)
}
