#' Multi-Environment Genomic Prediction via GBLUP
#'
#' This function performs genomic prediction across multiple environments using the
#' Genomic Best Linear Unbiased Prediction (GBLUP) method. It constructs a genomic
#' relationship matrix based on the VanRaden (2008) method and implements
#' a main effects model through the BGLR Gibbs sampler.
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
#' @param ploidy Integer. The ploidy level of the species. Default is 2 (diploid).
#' @param maf Numeric. Minor Allele Frequency threshold to filter SNPs. Default is 0.05.
#' @param nIter Total number of iterations for the BGLR Gibbs sampler. Default is 10000.
#' @param burnIn Number of burn-in iterations to be discarded. Default is 4000.
#' @param thin Thinning interval for the MCMC chain. Default is 10.
#' @param save_xlsx Logical. If \code{TRUE}, saves the predictive capacity results to an Excel file. Default is \code{TRUE}.
#' @param file_name Character string for the Excel file name. If \code{NULL}, a name
#'   is automatically generated as "gblup_CV(1, 2 or 0).xlsx". Default is \code{NULL}.
#'
#' @return A dataframe containing the predictive capacity (mean Pearson correlation)
#'   for each environment under the specified cross-validation scheme.
#'
#' @details
#' The function uses the \code{AGHmatrix} package to compute the VanRaden GRM and
#' \code{BGLR} for the Bayesian MCMC estimation. The model includes fixed environmental
#' effects and random genomic effects.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' results <- env_g_gblup(
#'   SNPs = X,
#'   y = data$yield,
#'   IDs = data$id,
#'   env = data$location,
#'   CV = "CV1"
#' )
#' }
#'
#' @export

env_g_gblup <- function(SNPs, y, IDs, env,
                        EZ = NULL,
                        CV = c("CV1", "CV2", "CV0"),
                        ploidy = 2,
                        maf = 0.05,
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

  if (!all(unique(IDs) %in% rownames(SNPs))) {
    stop("Some genotype IDs are not present in rownames(SNPs).")
  }

  n <- length(y)
  uIDs <- unique(IDs)
  uenv <- unique(env)

  if (is.null(EZ)) {
    EZ <- model.matrix(~ factor(env) - 1)
    colnames(EZ) <- paste0("Env_", uenv)
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
    fold_env <- sample(1:n_folds_env, size = n_folds_env)

    names(fold_env) <- uenv
    Y$Fold <- fold_env[Y$Env]
  }

  folds_run <- sort(unique(Y$Fold))

  IDs_factor <- factor(IDs, levels = rownames(SNPs))
  GZ <- model.matrix(~ IDs_factor - 1)

  cat("Computing VanRaden genomic relationship matrix\n")

  Gn <- AGHmatrix::Gmatrix(
    SNPs,
    method = "VanRaden",
    ploidy = ploidy,
    maf = maf
  )

  Gn <- as.matrix(Gn)

  G <- GZ %*% Gn %*% t(GZ)

  GDec <- eigen(G, symmetric = TRUE)

  values <- pmax(GDec$values, 0)

  GDec_list <- list(
    GBLUP = list(
      values = values,
      vectors = GDec$vectors
    )
  )

  list_metrics <- list()

  for (i in seq_along(GDec_list)) {

    GDec_i <- GDec_list[[i]]

    cat("\nRunning GBLUP\n")

    ETA_i <- list(
      list(X = EZ, model = "FIXED"),
      list(
        V = GDec_i$vectors,
        d = GDec_i$values,
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
              Model = "GBLUP",
              CV = CV,
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
    Predictive_Capacity ~ Model + CV + Environment,
    data = df_raw,
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  names(df_metrics)[names(df_metrics) == "Predictive_Capacity"] <- "pred"

  if (save_xlsx) {

    if (is.null(file_name)) {
      file_name <- paste0("gblup_", CV, ".xlsx")
    }

    writexl::write_xlsx(df_metrics, file_name)
  }

  return(df_metrics)
}

