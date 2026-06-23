#' Multi-Environment Genomic Prediction via Polynomial Kernel
#'
#' This function performs genomic prediction across multiple environments by
#' integrating genotype and environmental information using a Polynomial Kernel under a main effects framework.
#' It supports three cross-validation schemes (CV1, CV2, CV0) and performs a
#' grid search to optimize kernel hyperparameters.
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
#' @param degree A numeric vector of degree values for the Polynomial kernel. Default is 2 and 3.
#' @param scale A numeric vector of scale values for the Polynomial kernel. Default is 0.5, 1 and 2.
#' @param offset A numeric vector of offset values for the Polynomial kernel. Default is 0, 1 and 2.
#' @param nIter Total number of iterations for the BGLR Gibbs sampler. Default is 10000.
#' @param burnIn Number of burn-in iterations to be discarded. Default is 4000.
#' @param thin Thinning interval for the MCMC chain. Default is 10.
#' @param seed Integer value used to set the random seed for reproducibility in
#'   CV1 and CV2. Different seed values generate different random partitions of
#'   the dataset into cross-validation folds. Default is 123.
#' @param save_xlsx Logical. If TRUE, saves the results to an Excel file.Default is TRUE.
#' @param file_name Character string for the Excel file name. If NULL, a name is
#'   automatically generated based on the CV scheme. Default is NULL.
#'
#' @return A dataframe containing the predictive capacity (Pearson correlation)
#'   averaged across folds for each combination of hyperparameters and environment.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' results <- env_g_polynomial(SNPs = X, y = data$yield, IDs = data$id,
#'                             env = data$location, CV = "CV1")
#' }
#'
#' @export

env_g_polynomial <- function(SNPs, y, IDs, env,
                             EZ = NULL,
                             CV = c("CV1", "CV2", "CV0"),
                             degree = c(2, 3),
                             scale = c(0.5, 1, 2),
                             offset = c(0, 1, 2),
                             nIter = 10000,
                             burnIn = 4000,
                             thin = 10,
                             seed = 123,
                             save_xlsx = TRUE,
                             file_name = NULL) {

  CV <- match.arg(CV)

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
    fold_env <- sample(1:n_folds_env, size = n_folds_env)

    names(fold_env) <- uenv
    Y$Fold <- fold_env[Y$Env]
  }

  folds_run <- sort(unique(Y$Fold))

  IDs_factor <- factor(IDs, levels = rownames(SNPs))
  GZ <- model.matrix(~ IDs_factor - 1)

  # Grid de hiperparâmetros
  grid <- expand.grid(
    degree = degree,
    scale = scale,
    offset = offset
  )

  GDec_list <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {

    deg <- grid$degree[i]
    sc <- grid$scale[i]
    off <- grid$offset[i]

    cat(
      "Computing Polynomial kernel for degree =", deg,
      "scale =", sc,
      "offset =", off, "\n"
    )

    Gn <- kernlab::kernelMatrix(
      kernlab::polydot(
        degree = deg,
        scale = sc,
        offset = off
      ),
      SNPs
    )

    Gn <- as.matrix(Gn)

    G <- GZ %*% Gn %*% t(GZ)

    GDec <- eigen(G, symmetric = TRUE)

    values <- pmax(GDec$values, 0)

    GDec_list[[i]] <- list(
      degree = deg,
      scale = sc,
      offset = off,
      values = values,
      vectors = GDec$vectors
    )
  }

  names(GDec_list) <- apply(grid, 1, function(x) {
    paste0("degree_", x[1], "_scale_", x[2], "_offset_", x[3])
  })

  list_metrics <- list()

  for (i in seq_along(GDec_list)) {

    GDec_i <- GDec_list[[i]]

    cat(
      "\nRunning degree =", GDec_i$degree,
      "scale =", GDec_i$scale,
      "offset =", GDec_i$offset, "\n"
    )

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
              Model = "Polynomial",
              CV = CV,
              Degree = GDec_i$degree,
              Scale = GDec_i$scale,
              Offset = GDec_i$offset,
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
    Predictive_Capacity ~ Model + CV + Degree + Scale + Offset + Environment,
    data = df_raw,
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  names(df_metrics)[names(df_metrics) == "Predictive_Capacity"] <- "pred"

  if (save_xlsx) {

    if (is.null(file_name)) {
      file_name <- paste0("Polynomial_", CV, ".xlsx")
    }

    writexl::write_xlsx(df_metrics, file_name)
  }

  return(df_metrics)
}
