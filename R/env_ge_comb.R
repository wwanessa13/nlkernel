#' Combined Kernel Models for Multi-Environment Genomic Prediction with GxE Interaction
#'
#' This function fits combined kernel models for multi-environment genomic
#' prediction using the RKHS framework implemented in the BGLR package. It
#' evaluates pairwise combinations of nonlinear kernels and GBLUP, accounting
#' for fixed environmental effects, main genomic effects, and Genotype by
#' Environment (GxE) interaction effects. For each genomic kernel, the GxE
#' kernel is computed as the Hadamard product between the observation-level
#' genomic kernel and the environmental relationship matrix.
#'
#' @param SNPs A numeric matrix of SNP genotypes, with genotypes in rows and markers in columns.
#'   Row names must correspond to genotype IDs.
#' @param y A numeric vector of phenotypic values.
#' @param IDs A vector of genotype IDs corresponding to each phenotypic observation.
#' @param env A vector of environment labels corresponding to each phenotypic observation.
#' @param EZ Optional matrix of fixed environmental effects. If \code{NULL}, it is created from \code{env}.
#' @param CV Cross-validation scheme. One of \code{"CV1"}, \code{"CV2"}, or \code{"CV0"}.
#' @param polynomial_degree Degree parameter for the polynomialnomial kernel. Default is 2.
#' @param polynomial_scale Scale parameter for the polynomialnomial kernel. Default is 2.
#' @param polynomial_offset Offset parameter for the polynomialnomial kernel. Default is 2.
#' @param laplacian_sigma Sigma parameter for the Laplacian kernel. Default is 0.01.
#' @param bessel_sigma Sigma parameter for the Bessel kernel. Default is 0.1.
#' @param bessel_order Order parameter for the Bessel kernel. Default is 1.
#' @param bessel_degree Degree parameter for the Bessel kernel. Default is 2.
#' @param gaussian_sigma Sigma parameter for the Gaussian/RBF kernel. Default is 0.001.
#' @param ploidy Ploidy level used in \code{AGHmatrix::Gmatrix}. Default is 2.
#' @param nIter Total number of iterations for the BGLR model. Default is 10000.
#' @param burnIn Number of burn-in iterations for the BGLR model. Default is 4000.
#' @param thin Thinning interval for the BGLR model. Default is 10.
#' @param seed Integer value used to set the random seed for reproducibility in
#'   CV1 and CV2. Different seed values generate different random partitions of
#'   the dataset into cross-validation folds. Default is 123.
#' @param save_xlsx A logical value indicating whether to save results in an Excel file. Default is \code{TRUE}.
#' @param file_name Character string specifying the name of the Excel file. If \code{NULL}, a default name is used.
#'
#' @examples
#' \dontrun{
#' results <- env_ge_comb(
#'   SNPs = X,
#'   y = phen$yield,
#'   IDs = phen$genotype,
#'   env = phen$Env,
#'   CV = "CV2"
#' )
#' }
#'
#' @export

env_ge_combinations <- function(SNPs, y, IDs, env,
                        EZ = NULL,
                        CV = c("CV1", "CV2", "CV0"),
                        polynomial_degree = 2,
                        polynomial_scale = 2,
                        polynomial_offset = 2,
                        laplacian_sigma = 0.01,
                        bessel_sigma = 0.1,
                        bessel_order = 1,
                        bessel_degree = 2,
                        gaussian_sigma = 0.001,
                        ploidy = 2,
                        nIter = 10000,
                        burnIn = 4000,
                        thin = 10,
                        seed = 123,
                        save_xlsx = TRUE,
                        file_name = NULL) {

  CV <- match.arg(CV)

  SNPs <- as.matrix(SNPs)
  storage.mode(SNPs) <- "numeric"

  y <- as.numeric(y)
  IDs <- as.character(IDs)
  env <- as.character(env)

  if (length(y) != length(IDs) || length(y) != length(env)) {
    stop("The length of y, IDs, and env must be the same.")
  }

  if (is.null(rownames(SNPs))) {
    stop("SNPs must have row names corresponding to genotype IDs.")
  }

  n <- length(y)
  uIDs <- unique(IDs)
  uenv <- unique(env)

  if (!all(uIDs %in% rownames(SNPs))) {
    stop("Some genotype IDs are not present in rownames(SNPs).")
  }

  if (is.null(EZ)) {
    EZ <- model.matrix(~ factor(env) - 1)
    colnames(EZ) <- paste0("Env_", levels(factor(env)))
  }

  EZ <- as.matrix(EZ)

  if (nrow(EZ) != n) {
    stop("EZ must have the same number of rows as the length of y.")
  }

  Y <- data.frame(
    ID = IDs,
    Env = env,
    y = y
  )

  if (CV == "CV1") {
    set.seed(seed)
    n_folds <- 5
    fold_id <- rep(1:n_folds, length.out = length(uIDs))
    fold_id <- sample(fold_id)

    names(fold_id) <- uIDs

    Y$Fold <- fold_id[Y$ID]
  }

  if (CV == "CV2") {
    set.seed(seed)
    n_folds <- 5
    Y$Fold <- NA

    for (id in uIDs) {

      idx <- which(Y$ID == id)
      ni <- length(idx)

      Y$Fold[idx] <- sample(
        1:n_folds,
        size = ni,
        replace = ni > n_folds
      )
    }
  }

  if (CV == "CV0") {

    n_folds_env <- length(uenv)

    fold_env <- sample(
      1:n_folds_env,
      size = n_folds_env
    )

    names(fold_env) <- uenv

    Y$Fold <- fold_env[Y$Env]
  }

  folds_run <- sort(unique(Y$Fold))

  IDs_factor <- factor(IDs, levels = rownames(SNPs))

  GZ <- as.matrix(model.matrix(~ IDs_factor - 1))

  colnames(GZ) <- rownames(SNPs)

  obs_names <- paste0(IDs, "_", env, "_", seq_along(y))

  cat("Computing environmental relationship matrix\n")

  E <- EZ %*% t(EZ)

  rownames(E) <- obs_names
  colnames(E) <- obs_names

  kernels <- list(
    polynomial = function(SNPs) {
      kernlab::kernelMatrix(
        kernlab::polydot(
          degree = polynomial_degree,
          scale = polynomial_scale,
          offset = polynomial_offset
        ),
        SNPs
      )
    },

    laplacian = function(SNPs) {
      kernlab::kernelMatrix(
        kernlab::laplacedot(sigma = laplacian_sigma),
        SNPs
      )
    },

    bessel = function(SNPs) {
      kernlab::kernelMatrix(
        kernlab::besseldot(
          sigma = bessel_sigma,
          order = bessel_order,
          degree = bessel_degree
        ),
        SNPs
      )
    },

    gaussian = function(SNPs) {
      kernlab::kernelMatrix(
        kernlab::rbfdot(sigma = gaussian_sigma),
        SNPs
      )
    },

    GBLUP = function(SNPs) {
      AGHmatrix::Gmatrix(
        SNPmatrix = SNPs,
        method = "VanRaden",
        ploidy = ploidy,
        maf = 0.05
      )
    }
  )

  comb_list <- list(
    c("polynomial", "gaussian"),
    c("polynomial", "bessel"),
    c("polynomial", "laplacian"),
    c("gaussian", "bessel"),
    c("gaussian", "laplacian"),
    c("bessel", "laplacian"),
    c("polynomial", "GBLUP"),
    c("gaussian", "GBLUP"),
    c("bessel", "GBLUP"),
    c("laplacian", "GBLUP")
  )

  KDec_by_kernel <- list()

  for (kname in names(kernels)) {

    cat("Computing kernel:", kname, "\n")

    Gn <- kernels[[kname]](SNPs)
    Gn <- as.matrix(Gn)

    rownames(Gn) <- rownames(SNPs)
    colnames(Gn) <- rownames(SNPs)

    cat("Expanding", kname, "kernel to observation level\n")

    G <- GZ %*% Gn %*% t(GZ)

    rownames(G) <- obs_names
    colnames(G) <- obs_names

    cat("Computing", kname, "GxE interaction kernel\n")

    GxE <- G * E

    rownames(GxE) <- obs_names
    colnames(GxE) <- obs_names

    cat("Eigen decomposition of", kname, "G\n")

    GDec <- eigen(G, symmetric = TRUE)

    GDec$values <- pmax(GDec$values, 0)
    rownames(GDec$vectors) <- rownames(G)

    cat("Eigen decomposition of", kname, "GxE\n")

    GxEDec <- eigen(GxE, symmetric = TRUE)

    GxEDec$values <- pmax(GxEDec$values, 0)
    rownames(GxEDec$vectors) <- rownames(GxE)

    KDec_by_kernel[[kname]] <- list(
      G_values = GDec$values,
      G_vectors = GDec$vectors,
      GxE_values = GxEDec$values,
      GxE_vectors = GxEDec$vectors
    )
  }

  list_metrics <- list()

  for (comb in comb_list) {

    kname <- paste(comb, collapse = "_")

    cat("\nRunning combination:", kname, "\n")

    ETA_kernels <- unlist(
      lapply(comb, function(k) {
        list(
          list(
            V = KDec_by_kernel[[k]]$G_vectors,
            d = KDec_by_kernel[[k]]$G_values,
            model = "RKHS"
          ),
          list(
            V = KDec_by_kernel[[k]]$GxE_vectors,
            d = KDec_by_kernel[[k]]$GxE_values,
            model = "RKHS"
          )
        )
      }),
      recursive = FALSE
    )

    ETA_i <- c(
      list(
        list(
          X = EZ,
          model = "FIXED"
        )
      ),
      ETA_kernels
    )

    for (fold in folds_run) {

      cat("  Processing fold", fold, "\n")

      testing <- which(Y$Fold == fold)

      yNA <- y
      yNA[testing] <- NA

      fm <- BGLR::BGLR(
        y = yNA,
        ETA = ETA_i,
        nIter = nIter,
        burnIn = burnIn,
        thin = thin,
        verbose = FALSE
      )

      yHat <- fm$yHat

      for (a in uenv) {

        idx_env <- which(env == a)
        join <- intersect(idx_env, testing)

        if (length(join) > 1) {

          if (sd(yHat[join], na.rm = TRUE) > 0 &&
              sd(y[join], na.rm = TRUE) > 0) {

            cor_val <- cor(
              yHat[join],
              y[join],
              use = "complete.obs"
            )

          } else {

            cor_val <- NA
          }

          list_metrics[[length(list_metrics) + 1]] <-
            data.frame(
              Model = "Combined_Kernels_GxE",
              CV = CV,
              Kernel = kname,
              Fold = fold,
              Environment = a,
              Predictive_Capacity = cor_val
            )
        }
      }
    }
  }

  df_raw <- do.call(rbind, list_metrics)

  df_metrics <- aggregate(
    Predictive_Capacity ~ Model + CV + Kernel + Environment,
    data = df_raw,
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  names(df_metrics)[names(df_metrics) == "Predictive_Capacity"] <- "pred"

  df_metrics <- df_metrics[order(-df_metrics$pred), ]

  if (save_xlsx) {

    if (is.null(file_name)) {
      file_name <- paste0("combined_kernels_gxe_", CV, ".xlsx")
    }

    writexl::write_xlsx(df_metrics, file_name)
  }

  return(df_metrics)
}
