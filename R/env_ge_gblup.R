#' Multi-Environment Genomic Prediction via GBLUP with GxE Interaction
#'
#' This function performs genomic prediction across multiple environments by
#' accounting for main effects (Genotype and Environment) and the Genotype by
#' Environment (GxE) interaction. It uses a VanRaden-based genomic relationship
#' matrix and an environmental kernel to model the interaction via a Hadamard product.
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
#' @param ploidy Integer. The ploidy level of the species. Default is 2.
#' @param nIter Total number of iterations for the BGLR Gibbs sampler. Default is 10000.
#' @param burnIn Number of burn-in iterations to be discarded. Default is 4000.
#' @param thin Thinning interval for the MCMC chain. Default is 10.
#' @param save_xlsx Logical. If \code{TRUE}, saves the predictive capacity results to an Excel file. Default is \code{TRUE}.
#' @param file_name Character string for the Excel file name. If \code{NULL}, a name
#'   is automatically generated as "gblup_CV(1, 2 or 0).xlsx". Default is \code{NULL}.
#'
#' @return A dataframe containing the predictive capacity (mean Pearson correlation)
#'   for each environment, accounting for the GxE interaction model.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' results <- env_ge_gblup(
#'   SNPs = X,
#'   y = phen$yield,
#'   IDs = phen$genotype,
#'   env = phen$Env,
#'   CV = "CV2"
#' )
#' }
#'
#' @export

env_ge_gblup <- function(SNPs, y, IDs, env,
                        EZ = NULL,
                        CV = c("CV1", "CV2", "CV0"),
                        ploidy = 2,
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

  cat("Computing VanRaden genomic relationship matrix\n")

  Gn <- AGHmatrix::Gmatrix(
    SNPmatrix = SNPs,
    method = "VanRaden",
    ploidy = ploidy,
  )

  Gn <- as.matrix(Gn)

  IDs_factor <- factor(IDs, levels = rownames(Gn))

  GZ <- as.matrix(model.matrix(~ IDs_factor - 1))

  colnames(GZ) <- rownames(Gn)

  cat("Expanding genomic relationship matrix to observation level\n")

  G <- GZ %*% Gn %*% t(GZ)

  obs_names <- paste0(IDs, "_", env, "_", seq_along(y))

  rownames(G) <- obs_names
  colnames(G) <- obs_names

  cat("Computing environmental relationship matrix\n")

  E <- EZ %*% t(EZ)

  rownames(E) <- obs_names
  colnames(E) <- obs_names

  cat("Computing GxE interaction kernel\n")

  GxE <- G * E

  rownames(GxE) <- obs_names
  colnames(GxE) <- obs_names

  cat("Eigen decomposition of G\n")

  GDec <- eigen(G, symmetric = TRUE)

  GDec$values <- pmax(GDec$values, 0)
  rownames(GDec$vectors) <- rownames(G)

  cat("Eigen decomposition of GxE\n")

  GxEDec <- eigen(GxE, symmetric = TRUE)

  GxEDec$values <- pmax(GxEDec$values, 0)
  rownames(GxEDec$vectors) <- rownames(GxE)

  ETA <- list(
    list(
      X = EZ,
      model = "FIXED"
    ),
    list(
      V = GDec$vectors,
      d = GDec$values,
      model = "RKHS"
    ),
    list(
      V = GxEDec$vectors,
      d = GxEDec$values,
      model = "RKHS"
    )
  )

  list_metrics <- list()

  cat("\nRunning GBLUP + GxE\n")

  for (fold in folds_run) {

    cat("  Processing fold", fold, "\n")

    testing <- which(Y$Fold == fold)

    yNA <- y
    yNA[testing] <- NA

    fm <- BGLR::BGLR(
      y = yNA,
      ETA = ETA,
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
            Model = "GBLUP_GxE",
            CV = CV,
            Fold = fold,
            Environment = a,
            Predictive_Capacity = cor_val
          )
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
      file_name <- paste0("gblup_gxe_", CV, ".xlsx")
    }

    writexl::write_xlsx(df_metrics, file_name)
  }

  return(df_metrics)
}
