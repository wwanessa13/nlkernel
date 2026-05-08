#' Kernel PCA with Bessel Kernel for Multi-Environment Genomic Prediction
#'
#' This function fits kernel PCA models using a Bessel kernel for
#' multi-environment genomic prediction. It evaluates different combinations of
#' sigma, order, and degree parameters. Kernel principal components are selected
#' according to a variance-explained threshold, and a genomic kernel is
#' constructed from the selected component scores. The kernel is expanded to the
#' observation level and combined with fixed environmental effects. Predictive
#' capacity is evaluated by environment using CV1, CV2, or CV0 cross-validation
#' schemes within the RKHS framework implemented in BGLR.
#'
#' @param SNPs A numeric matrix of SNP genotypes, with genotypes in rows and markers in columns.
#' Row names must correspond to genotype IDs.
#' @param y A numeric vector of phenotypic values.
#' @param IDs A vector of genotype IDs corresponding to each phenotypic observation.
#' @param env A vector of environment labels corresponding to each phenotypic observation.
#' @param EZ Optional matrix of fixed environmental effects. If NULL, it is created from env.
#' @param CV Cross-validation scheme. One of "CV1", "CV2", or "CV0".
#' @param sg A numeric vector of sigma values for the Bessel kernel. Default is c(0.1, 0.5, 1).
#' @param ord A numeric vector of order values for the Bessel kernel. Default is c(0, 1, 2).
#' @param dg A numeric vector of degree values for the Bessel kernel. Default is c(2, 3).
#' @param var_threshold Minimum proportion of variance explained required for a principal component to be retained. Default is 0.01.
#' @param nIter Total number of iterations for the BGLR model. Default is 10000.
#' @param burnIn Number of burn-in iterations for the BGLR model. Default is 4000.
#' @param thin Thinning interval for the BGLR model. Default is 10.
#' @param save_xlsx A logical value indicating whether to save results in an Excel file. Default is TRUE.
#' @param file_name Character string specifying the name of the Excel file. If NULL, a default name is used.
#'
#' @return A data frame with the mean predictive capacity by model, CV scheme,
#' Bessel kernel parameters, number of selected PCs, variance threshold, and environment.
#'
#' @export

env_g_pca_bessel <- function(SNPs, y, IDs, env,
                              EZ = NULL,
                              CV = c("CV1", "CV2", "CV0"),
                              sg = c(0.1, 0.5, 1),
                              ord = c(0, 1, 2),
                              dg = c(2, 3),
                              var_threshold = 0.01,
                              nIter = 10000,
                              burnIn = 4000,
                              thin = 10,
                              save_xlsx = TRUE,
                              file_name = NULL) {

  CV <- match.arg(CV)
  set.seed(1)

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

  grid <- expand.grid(
    degree = dg,
    order = ord,
    sigma = sg
  )

  GDec_list <- list()
  counter <- 1

  for (i in seq_len(nrow(grid))) {

    sig <- grid$sigma[i]
    or <- grid$order[i]
    deg <- grid$degree[i]

    cat(
      "Computing KPCA Bessel kernel for sigma =", sig,
      "| order =", or,
      "| degree =", deg, "\n"
    )

    kpca_temp <- kernlab::kpca(
      x = SNPs,
      kernel = "besseldot",
      kpar = list(
        sigma = sig,
        order = or,
        degree = deg
      ),
      features = 0
    )

    eig_vals <- kernlab::eig(kpca_temp)
    var_explained <- eig_vals / sum(eig_vals)

    nPC <- sum(var_explained > var_threshold)

    if (nPC < 2) {
      warning(
        paste(
          "Sigma", sig,
          "Order", or,
          "Degree", deg,
          "selected fewer than 2 PCs. Skipping."
        )
      )
      next
    }

    cat("Number of PCs selected:", nPC, "\n")

    kpca_model <- kernlab::kpca(
      x = SNPs,
      kernel = "besseldot",
      kpar = list(
        sigma = sig,
        order = or,
        degree = deg
      ),
      features = nPC
    )

    embedding <- kpca_model@rotated
    embedding <- as.matrix(embedding)

    Gn <- tcrossprod(embedding) / ncol(embedding)

    rownames(Gn) <- rownames(SNPs)
    colnames(Gn) <- rownames(SNPs)

    G <- GZ %*% Gn %*% t(GZ)

    GDec <- eigen(G, symmetric = TRUE)

    values <- pmax(GDec$values, 0)

    GDec_list[[counter]] <- list(
      sigma = sig,
      order = or,
      degree = deg,
      nPC = nPC,
      var_threshold = var_threshold,
      values = values,
      vectors = GDec$vectors
    )

    counter <- counter + 1
  }

  if (length(GDec_list) == 0) {
    stop("No KPCA Bessel model was fitted. All parameter combinations selected fewer than 2 PCs.")
  }

  names(GDec_list) <- vapply(GDec_list, function(x) {
    paste0(
      "degree_", x$degree,
      "_order_", x$order,
      "_sigma_", x$sigma
    )
  }, character(1))

  list_metrics <- list()

  for (i in seq_along(GDec_list)) {

    GDec_i <- GDec_list[[i]]

    cat(
      "\nRunning KPCA Bessel model with sigma =", GDec_i$sigma,
      "| order =", GDec_i$order,
      "| degree =", GDec_i$degree,
      "| nPC =", GDec_i$nPC, "\n"
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
              Model = "KPCA_Bessel",
              CV = CV,
              Sigma = GDec_i$sigma,
              Order = GDec_i$order,
              Degree = GDec_i$degree,
              nPC = GDec_i$nPC,
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
    Predictive_Capacity ~ Model + CV + Sigma + Order + Degree +
      nPC + Environment,
    data = df_raw,
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  names(df_metrics)[names(df_metrics) == "Predictive_Capacity"] <- "pred"

  df_metrics <- df_metrics[order(-df_metrics$pred), ]

  if (save_xlsx) {

    if (is.null(file_name)) {
      file_name <- paste0("PCA_Bessel_", CV, ".xlsx")
    }

    writexl::write_xlsx(df_metrics, file_name)
  }

  return(df_metrics)
}
