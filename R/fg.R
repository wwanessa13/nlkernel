#' Train Final Single Environment Model for Genomic Prediction
#'
#' @description
#' This function trains a final genomic prediction model using all available
#' phenotypic and marker data. It supports GBLUP and nonlinear kernels,
#' including Gaussian/RBF, Laplacian, Polynomial, Bessel, ANOVA, and
#' hyperbolic tangent (tanh) kernels, using the RKHS framework implemented
#' in the BGLR package.
#'
#' If more than one model is provided, they are fitted simultaneously as
#' combined kernels in the BGLR ETA.
#'
#' @param SNPs A numeric matrix or data frame of SNP genotypes, with individuals
#'   in rows and markers in columns.
#' @param y A numeric vector of phenotypic values corresponding to the individuals.
#' @param model A character vector indicating the model or models to fit.
#'   Options are "gblup", "gaussian", "laplacian", "polynomial", "bessel",
#'   "anova", and "tanh".
#'   If more than one model is provided, they are fitted as combined kernels
#'   in the BGLR ETA.
#' @param gaussian_sigma Sigma parameter for the Gaussian/RBF kernel.
#'   Default is 0.001.
#' @param laplacian_sigma Sigma parameter for the Laplacian kernel.
#'   Default is 0.01.
#' @param polynomial_degree Degree parameter for the polynomial kernel.
#'   Default is 2.
#' @param polynomial_scale Scale parameter for the polynomial kernel.
#'   Default is 2.
#' @param polynomial_offset Offset parameter for the polynomial kernel.
#'   Default is 2.
#' @param bessel_sigma Sigma parameter for the Bessel kernel.
#'   Default is 0.1.
#' @param bessel_order Order parameter for the Bessel kernel.
#'   Default is 1.
#' @param bessel_degree Degree parameter for the Bessel kernel.
#'   Default is 2.
#' @param anova_sigma Sigma parameter for the ANOVA kernel.
#'   Default is 0.1.
#' @param anova_degree Degree parameter for the ANOVA kernel.
#'   Default is 2.
#' @param tanh_scale Scale parameter for the hyperbolic tangent kernel.
#'   Default is 0.01.
#' @param tanh_offset Offset parameter for the hyperbolic tangent kernel.
#'   Default is 1.
#' @param nIter Total number of iterations for the BGLR model.
#'   Default is 10000.
#' @param burnIn Number of burn-in iterations for the BGLR model.
#'   Default is 4000.
#' @param thin Thinning interval for the BGLR model.
#'   Default is 10.
#' @param save_model Logical value indicating whether to save the fitted model
#'   as an RDS file. Default is TRUE.
#' @param file_name Character string specifying the name of the RDS file.
#'   Default is "final_bglr_model.rds".
#'
#' @return A list containing the fitted BGLR model, fitted values, residuals,
#' kernel matrices, model information, hyperparameters, and training data.
#'
#' @examples
#' \dontrun{
#' fit <- train_final_single(
#'   SNPs = SNPs,
#'   y = y,
#'   model = c("laplacian", "anova", "tanh"),
#'   laplacian_sigma = 0.01,
#'   anova_sigma = 0.1,
#'   anova_degree = 2,
#'   tanh_scale = 0.01,
#'   tanh_offset = 1,
#'   file_name = "final_laplacian_anova_tanh.rds"
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
                               anova_sigma = 0.1,
                               anova_degree = 2,
                               tanh_scale = 0.01,
                               tanh_offset = 1,
                               nIter = 10000,
                               burnIn = 4000,
                               thin = 10,
                               save_model = TRUE,
                               file_name = "final_bglr_model.rds") {

  # ============================================================
  # VERIFICAR PACOTES
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
  # PREPARAR DADOS
  # ============================================================

  SNPs <- as.matrix(SNPs)
  storage.mode(SNPs) <- "numeric"

  y <- as.numeric(y)

  if (nrow(SNPs) != length(y)) {
    stop("Number of rows in SNPs must match length of y.")
  }

  if (any(is.na(y))) {
    stop(
      "This function trains the final model using all available data. ",
      "Remove or impute missing values in y before fitting."
    )
  }

  # ============================================================
  # MODELOS DISPONÍVEIS
  # ============================================================

  valid_models <- c(
    "gblup",
    "gaussian",
    "laplacian",
    "polynomial",
    "bessel",
    "anova",
    "tanh"
  )

  if (!all(model %in% valid_models)) {

    stop(
      "Invalid model. Use one or more of: ",
      paste(valid_models, collapse = ", ")
    )
  }

  # ============================================================
  # FUNÇÃO PARA GARANTIR MATRIZ PSD
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
  # NORMALIZAR KERNEL
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
  # CONSTRUIR KERNEL
  # ============================================================

  make_kernel <- function(model_i) {

    # ----------------------------------------------------------
    # GBLUP
    # ----------------------------------------------------------

    if (model_i == "gblup") {

      K <- AGHmatrix::Gmatrix(
        SNPmatrix = SNPs,
        method = "VanRaden",
        ploidy = 2
      )

      K <- normalize_kernel(K)

      return(
        list(
          K = K,
          method = "gblup"
        )
      )
    }

    # ----------------------------------------------------------
    # GAUSSIAN / RBF
    # ----------------------------------------------------------

    if (model_i == "gaussian") {

      K <- kernlab::kernelMatrix(
        kernlab::rbfdot(
          sigma = gaussian_sigma
        ),
        SNPs
      )

      K <- normalize_kernel(K)

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

      return(
        list(
          K = K,
          method = "polynomial"
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

      return(
        list(
          K = K,
          method = "bessel"
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

      return(
        list(
          K = K,
          method = "anova"
        )
      )
    }

    # ----------------------------------------------------------
    # TANH / SIGMOID
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

      return(
        list(
          K = K,
          method = "tanh"
        )
      )
    }
  }

  # ============================================================
  # TREINAMENTO
  # ============================================================

  cat(
    "\nTraining final single-environment BGLR model ",
    "using all available data\n"
  )

  cat(
    "Model:",
    paste(model, collapse = " + "),
    "\n"
  )

  # ============================================================
  # CONSTRUIR KERNELS
  # ============================================================

  kernel_outputs <- lapply(
    model,
    make_kernel
  )

  names(kernel_outputs) <- model

  kernels <- lapply(
    kernel_outputs,
    function(x) x$K
  )

  # ============================================================
  # CONSTRUIR ETA
  # ============================================================

  ETA <- lapply(
    seq_along(kernel_outputs),
    function(i) {

      list(
        K = kernel_outputs[[i]]$K,
        model = "RKHS"
      )
    }
  )

  names(ETA) <- model

  # ============================================================
  # BGLR
  # ============================================================

  fit <- BGLR::BGLR(
    y = y,
    ETA = ETA,
    nIter = nIter,
    burnIn = burnIn,
    thin = thin,
    verbose = FALSE
  )

  # ============================================================
  # VALORES AJUSTADOS
  # ============================================================

  fitted_values <- fit$yHat

  residuals <- y - fitted_values

  # ============================================================
  # OBJETO FINAL
  # ============================================================

  model_object <- list(

    fit = fit,

    yHat = fitted_values,

    residuals = residuals,

    model = model,

    ETA = ETA,

    kernels = kernels,

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

      anova_sigma = anova_sigma,

      anova_degree = anova_degree,

      tanh_scale = tanh_scale,

      tanh_offset = tanh_offset,

      nIter = nIter,

      burnIn = burnIn,

      thin = thin
    ),

    call = match.call()
  )

  class(model_object) <- "train_final_single"

  # ============================================================
  # SALVAR MODELO
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

  return(model_object)
}
