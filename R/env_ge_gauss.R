#' Multi-Environment Genomic Prediction via Gaussian Kernel with GxE Interaction
#'
#' This function performs genomic prediction across multiple environments by
#' accounting for fixed environmental effects, main genomic effects modeled with
#' a Gaussian kernel, and Genotype by Environment (GxE) interaction. The GxE
#' kernel is computed as the Hadamard product between the observation-level
#' genomic kernel and the environmental relationship matrix.
#'
#' @param SNPs A numeric matrix of SNP genotypes (individuals in rows, markers in columns).
#'   Must have \code{rownames} corresponding to the genotype IDs.
#' @param y A numeric vector of phenotypic values.
#' @param IDs A character vector indicating the genotype identity for each observation in \code{y}.
#' @param env A character vector indicating the environment for each observation in \code{y}.
#' @param EZ An incidence matrix for fixed environmental effects. If \code{NULL},
#'   it is automatically generated from the \code{env} vector.
#' @param CV A character string specifying the cross-validation scheme:
#' "CV1": Prediction of unobserved genotypes in observed environments.
#' "CV2": Prediction of genotypes observed in only a subset of environments.
#' "CV0": Prediction of observed genotypes in completely unobserved environments.
#' @param sg A numeric vector of sigma values for the Gaussian kernel.
#'   Default is \code{c(0.001, 0.01, 0.05, 0.1)}.
#' @param nIter Total number of iterations for the BGLR Gibbs sampler. Default is 10000.
#' @param burnIn Number of burn-in iterations to be discarded. Default is 4000.
#' @param thin Thinning interval for the MCMC chain. Default is 10.
#' @param save_xlsx Logical. If \code{TRUE}, saves the predictive capacity results to an Excel file. Default is \code{TRUE}.
#' @param file_name Character string for the Excel file name. If \code{NULL}, a name
#'   is automatically generated as "gblup_CV(1, 2 or 0).xlsx". Default is \code{NULL}.
#'
#' @return A dataframe containing the predictive capacity (mean Pearson correlation)
#'   for each Gaussian sigma value and environment, accounting for the GxE
#'   interaction model.
#'
#' @details
#' The model implemented is:
#' \deqn{y = Xb + Zg + Zi + e}
#' where \eqn{Xb} represents fixed environmental effects, \eqn{Zg} represents the
#' main genomic effect modeled with the Gaussian kernel, and \eqn{Zi} represents
#' the GxE interaction effect. The GxE kernel is computed as the Hadamard product
#' between the observation-level Gaussian genomic kernel (G) and the environmental
#' relationship matrix (E).
#'
#' @examples
#' \dontrun{
#' results <- env_ge_gaussian(
#'   SNPs = X,
#'   y = phen$yield,
#'   IDs = phen$genotype,
#'   env = phen$Env,
#'   CV = "CV2"
#' )
#' }
#'
#' @export

env_ge_gaussian <- function(SNPs, y, IDs, env,
                             EZ = NULL,
                             CV = c("CV1", "CV2", "CV0"),
                             sg = c(0.001, 0.01, 0.05, 0.1),
                             nIter = 10000,
                             burnIn = 4000,
                             thin = 10,
                             save_xlsx = TRUE,
                             file_name = NULL) {

  CV <- match.arg(CV)
  set.seed(1)

  SNPs <- as.matrix(SNPs)
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

    n_folds <- 5

    fold_id <- rep(1:n_folds, length.out = length(uIDs))
    fold_id <- sample(fold_id)

    names(fold_id) <- uIDs

    Y$Fold <- fold_id[Y$ID]
  }

  if (CV == "CV2") {

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

  GDec_list <- vector("list", length(sg))

  for (i in seq_along(sg)) {

    sigma_i <- sg[i]

    cat("Computing Gaussian kernel for sigma =", sigma_i, "\n")

    Gn <- kernlab::kernelMatrix(
      kernlab::rbfdot(sigma = sigma_i),
      SNPs
    )

    Gn <- as.matrix(Gn)

    rownames(Gn) <- rownames(SNPs)
    colnames(Gn) <- rownames(SNPs)

    cat("Expanding Gaussian genomic kernel to observation level\n")

    G <- GZ %*% Gn %*% t(GZ)

    rownames(G) <- obs_names
    colnames(G) <- obs_names

    cat("Computing Gaussian GxE interaction kernel\n")

    GxE <- G * E

    rownames(GxE) <- obs_names
    colnames(GxE) <- obs_names

    cat("Eigen decomposition of Gaussian G\n")

    GDec <- eigen(G, symmetric = TRUE)

    GDec$values <- pmax(GDec$values, 0)
    rownames(GDec$vectors) <- rownames(G)

    cat("Eigen decomposition of Gaussian GxE\n")

    GxEDec <- eigen(GxE, symmetric = TRUE)

    GxEDec$values <- pmax(GxEDec$values, 0)
    rownames(GxEDec$vectors) <- rownames(GxE)

    GDec_list[[i]] <- list(
      sigma = sigma_i,
      G_values = GDec$values,
      G_vectors = GDec$vectors,
      GxE_values = GxEDec$values,
      GxE_vectors = GxEDec$vectors
    )
  }

  names(GDec_list) <- paste0("sigma_", sg)

  list_metrics <- list()

  for (i in seq_along(GDec_list)) {

    GDec_i <- GDec_list[[i]]

    cat("\nRunning Gaussian for sigma =", GDec_i$sigma, "\n")

    ETA_i <- list(
      list(
        X = EZ,
        model = "FIXED"
      ),
      list(
        V = GDec_i$G_vectors,
        d = GDec_i$G_values,
        model = "RKHS"
      ),
      list(
        V = GDec_i$GxE_vectors,
        d = GDec_i$GxE_values,
        model = "RKHS"
      )
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
              Model = "Gaussian_GxE",
              CV = CV,
              Sigma = GDec_i$sigma,
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
    Predictive_Capacity ~ Model + CV + Sigma + Environment,
    data = df_raw,
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  names(df_metrics)[names(df_metrics) == "Predictive_Capacity"] <- "pred"

  if (save_xlsx) {

    if (is.null(file_name)) {
      file_name <- paste0("gaussian_gxe_", CV, ".xlsx")
    }

    writexl::write_xlsx(df_metrics, file_name)
  }

  return(df_metrics)
}
