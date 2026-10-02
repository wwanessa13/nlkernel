#' Transform SNP Matrix for Genomic Prediction
#'
#' @description
#' This function transforms a SNP marker matrix using genomic relationship
#' and nonlinear kernel approaches. It does not require phenotypic values
#' and does not fit a BGLR model.
#'
#' @param SNPs A numeric matrix or data frame of SNP genotypes, with
#' individuals in rows and markers in columns.
#' @param model Character vector indicating the transformation method or
#' methods. Options are "gblup", "gaussian", "laplacian", "polynomial",
#' "anova", "tanh", and "bessel".
#' If more than one model is provided, the function returns a separate
#' matrix for each model.
#' @param gaussian_sigma Sigma parameter for the Gaussian kernel.
#' Default is 0.001.
#' @param laplacian_sigma Sigma parameter for the Laplacian kernel.
#' Default is 0.01.
#' @param polynomial_degree Degree parameter for the polynomial kernel.
#' Default is 2.
#' @param polynomial_scale Scale parameter for the polynomial kernel.
#' Default is 2.
#' @param polynomial_offset Offset parameter for the polynomial kernel.
#' Default is 2.
#' @param anova_sigma Sigma parameter for the ANOVA kernel.
#' Default is 0.01.
#' @param anova_degree Degree parameter for the ANOVA kernel.
#' Default is 2.
#' @param tanh_scale Scale parameter for the hyperbolic tangent kernel.
#' Default is 0.01.
#' @param tanh_offset Offset parameter for the hyperbolic tangent kernel.
#' Default is 1.
#' @param bessel_sigma Sigma parameter for the Bessel kernel.
#' Default is 0.1.
#' @param bessel_order Order parameter for the Bessel kernel.
#' Default is 1.
#' @param bessel_degree Degree parameter for the Bessel kernel.
#' Default is 2.
#' @param save_transform Logical value indicating whether to save the
#' transformed object as an RDS file. Default is FALSE.
#' @param file_name Character string specifying the name of the RDS file.
#' Default is "snp_transform.rds".
#'
#' @return A list containing the transformed SNP matrices or kernel matrices,
#' selected model(s), kernel objects, hyperparameters, and original SNP matrix.
#'
#' @examples
#' \dontrun{
#' obj <- transform_snps(
#'   SNPs = SNPs,
#'   model = c("anova", "tanh", "gaussian")
#' )
#'
#' K_anova <- obj$K$anova
#' K_tanh <- obj$K$tanh
#' K_gaussian <- obj$K$gaussian
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
                           anova_sigma = 0.01,
                           anova_degree = 2,
                           tanh_scale = 0.01,
                           tanh_offset = 1,
                           bessel_sigma = 0.1,
                           bessel_order = 1,
                           bessel_degree = 2,
                           ploidy = 2,
                           save_transform = FALSE,
                           file_name = "snp_transform.rds") {

  SNPs <- as.matrix(SNPs)
  storage.mode(SNPs) <- "numeric"

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

  make_psd <- function(K) {

    K <- as.matrix(K)
    K <- (K + t(K)) / 2

    eig <- eigen(K, symmetric = TRUE)

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

  normalize_kernel <- function(K) {

    K <- as.matrix(K)

    K <- make_psd(K)

    tr <- sum(diag(K))

    if (is.finite(tr) && tr > 0) {

      K <- K / tr * nrow(K)

    }

    return(K)
  }

  make_transform <- function(model_i) {

    # ============================================================
    # GBLUP
    # ============================================================

    if (model_i == "gblup") {

      K <- AGHmatrix::Gmatrix(
        SNPmatrix = SNPs,
        method = "VanRaden",
        ploidy = ploidy
      )

      K <- normalize_kernel(K)

      return(list(
        transformed = K,
        K = K,
        object = NULL
      ))
    }

    # ============================================================
    # GAUSSIAN
    # ============================================================

    if (model_i == "gaussian") {

      K <- kernelMatrix(
        rbfdot(
          sigma = gaussian_sigma
        ),
        SNPs
      )

      K <- normalize_kernel(K)

      return(list(
        transformed = K,
        K = K,
        object = NULL
      ))
    }

    # ============================================================
    # LAPLACIAN
    # ============================================================

    if (model_i == "laplacian") {

      K <- kernelMatrix(
        laplacedot(
          sigma = laplacian_sigma
        ),
        SNPs
      )

      K <- normalize_kernel(K)

      return(list(
        transformed = K,
        K = K,
        object = NULL
      ))
    }

    # ============================================================
    # POLYNOMIAL
    # ============================================================

    if (model_i == "polynomial") {

      K <- kernelMatrix(
        polydot(
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
        object = NULL
      ))
    }

    # ============================================================
    # ANOVA
    # ============================================================

    if (model_i == "anova") {

      K <- kernelMatrix(
        anovadot(
          sigma = anova_sigma,
          degree = anova_degree
        ),
        SNPs
      )

      K <- normalize_kernel(K)

      return(list(
        transformed = K,
        K = K,
        object = NULL
      ))
    }

    # ============================================================
    # TANH / SIGMOID
    # ============================================================

    if (model_i == "tanh") {

      K <- kernelMatrix(
        tanhdot(
          scale = tanh_scale,
          offset = tanh_offset
        ),
        SNPs
      )

      K <- normalize_kernel(K)

      return(list(
        transformed = K,
        K = K,
        object = NULL
      ))
    }

    # ============================================================
    # BESSEL
    # ============================================================

    if (model_i == "bessel") {

      K <- kernelMatrix(
        besseldot(
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
        object = NULL
      ))
    }

    stop("Unknown model: ", model_i)
  }

  # ==============================================================
  # TRANSFORM KERNELS
  # ==============================================================

  transform_outputs <- lapply(
    model,
    make_transform
  )

  names(transform_outputs) <- model

  # ==============================================================
  # OUTPUT OBJECT
  # ==============================================================

  transform_object <- list(

    transformed = lapply(
      transform_outputs,
      function(x) x$transformed
    ),

    K = lapply(
      transform_outputs,
      function(x) x$K
    ),

    model = model,

    objects = lapply(
      transform_outputs,
      function(x) x$object
    ),

    SNPs_original = SNPs,

    hyperparameters = list(

      gaussian_sigma = gaussian_sigma,

      laplacian_sigma = laplacian_sigma,

      polynomial_degree = polynomial_degree,

      polynomial_scale = polynomial_scale,

      polynomial_offset = polynomial_offset,

      anova_sigma = anova_sigma,

      anova_degree = anova_degree,

      tanh_scale = tanh_scale,

      tanh_offset = tanh_offset,

      bessel_sigma = bessel_sigma,

      bessel_order = bessel_order,

      bessel_degree = bessel_degree

    ),

    call = match.call()

  )

  class(transform_object) <- "snp_transform"

  # ==============================================================
  # SAVE
  # ==============================================================

  if (save_transform) {

    saveRDS(
      transform_object,
      file = file_name
    )

    cat(
      "\nTransformed SNP object saved as:",
      file_name,
      "\n"
    )
  }

  return(transform_object)
}
