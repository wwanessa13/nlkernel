#' Multi-Environment Genomic Prediction via Polynomial Kernel with GxE Interaction
#'
#' This function performs genomic prediction across multiple environments by
#' accounting for fixed environmental effects, main genomic effects modeled with
#' a Polynomial kernel, and Genotype by Environment (GxE) interaction. The GxE kernel
#' is computed as the Hadamard product between the observation-level genomic
#' kernel and the environmental relationship matrix.
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
#' @param save_xlsx Logical. If \code{TRUE}, saves the predictive capacity results to an Excel file. Default is \code{TRUE}.
#' @param file_name Character string for the Excel file name. If \code{NULL}, a name
#'   is automatically generated as "gblup_CV(1, 2 or 0).xlsx". Default is \code{NULL}.
#'
#' @return A dataframe containing the predictive capacity (mean Pearson correlation)
#'   for each combination of Polynomial kernel hyperparameters and environment,
#'   accounting for the GxE interaction model.
#'
#' @examples
#' \dontrun{
#' results <- env_ge_polynomial(
#'   SNPs = X,
#'   y = phen$yield,
#'   IDs = phen$genotype,
#'   env = phen$Env,
#'   CV = "CV2"
#' )
#' }
#'
#' @export

env_ge_polynomial <- function(SNPs, y, IDs, env,
                          EZ = NULL,
                          CV = c("CV1", "CV2", "CV0"),
                          degree = c(2, 3),
                          scale = c(0.5, 1, 2),
                          offset = c(0, 1, 2),
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

  grid <- expand.grid(
    degree = degree,
    offset = offset,
    scale = scale
  )

  GDec_list <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {

    deg <- grid$degree[i]
    or <- grid$offset[i]
    sig <- grid$scale[i]

    cat(
      "Computing Polynomial kernel for degree =", deg,
      "| offset =", or,
      "| scale =", sig, "\n"
    )

    K_poly <- kernlab::kernelMatrix(
      kernlab::polydot(
        degree = deg,
        offset = or,
        scale = sig
      ),
      SNPs
    )

    K_poly <- as.matrix(K_poly)

    rownames(K_poly) <- rownames(SNPs)
    colnames(K_poly) <- rownames(SNPs)

    cat("Expanding Polynomial genomic kernel to observation level\n")

    G <- GZ %*% K_poly %*% t(GZ)

    rownames(G) <- obs_names
    colnames(G) <- obs_names

    cat("Computing Polynomial GxE interaction kernel\n")

    GxE <- G * E

    rownames(GxE) <- obs_names
    colnames(GxE) <- obs_names

    cat("Eigen decomposition of Polynomial G\n")

    GDec <- eigen(G, symmetric = TRUE)

    GDec$values <- pmax(GDec$values, 0)
    rownames(GDec$vectors) <- rownames(G)

    cat("Eigen decomposition of Polynomial GxE\n")

    GxEDec <- eigen(GxE, symmetric = TRUE)

    GxEDec$values <- pmax(GxEDec$values, 0)
    rownames(GxEDec$vectors) <- rownames(GxE)

    GDec_list[[i]] <- list(
      degree = deg,
      offset = or,
      scale = sig,
      G_values = GDec$values,
      G_vectors = GDec$vectors,
      GxE_values = GxEDec$values,
      GxE_vectors = GxEDec$vectors
    )
  }

  names(GDec_list) <- apply(grid, 1, function(x) {
    paste0("degree_", x[1], "_offset_", x[2], "_scale_", x[3])
  })

  list_metrics <- list()

  for (i in seq_along(GDec_list)) {

    GDec_i <- GDec_list[[i]]

    cat(
      "\nRunning Polynomial for degree =", GDec_i$degree,
      "| offset =", GDec_i$offset,
      "| scale =", GDec_i$scale, "\n"
    )

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
              Model = "Polynomial_GxE",
              CV = CV,
              Degree = GDec_i$degree,
              Offset = GDec_i$offset,
              Scale = GDec_i$scale,
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
    Predictive_Capacity ~ Model + CV + Degree + Offset + Scale + Environment,
    data = df_raw,
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  names(df_metrics)[names(df_metrics) == "Predictive_Capacity"] <- "pred"

  if (save_xlsx) {

    if (is.null(file_name)) {
      file_name <- paste0("polynomial_gxe_", CV, ".xlsx")
    }

    writexl::write_xlsx(df_metrics, file_name)
  }

  return(df_metrics)
}
