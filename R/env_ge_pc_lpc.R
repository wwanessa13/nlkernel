#' Kernel PCA with Laplacian Kernel for Multi-Environment Genomic Prediction with GxE Interaction
#'
#' This function fits kernel PCA models using a Laplacian kernel for
#' multi-environment genomic prediction. It evaluates different sigma values.
#' Kernel principal components are selected according to a variance-explained
#' threshold, and a genomic kernel is constructed from the selected component
#' scores. The kernel is expanded to the observation level and combined with
#' fixed environmental effects and a Genotype by Environment (GxE) interaction
#' kernel.
#'
#' @param SNPs A numeric matrix of SNP genotypes, with genotypes in rows and markers in columns.
#'   Row names must correspond to genotype IDs.
#' @param y A numeric vector of phenotypic values.
#' @param IDs A vector of genotype IDs corresponding to each phenotypic observation.
#' @param env A vector of environment labels corresponding to each phenotypic observation.
#' @param EZ Optional matrix of fixed environmental effects. If \code{NULL}, it is created from \code{env}.
#' @param CV Cross-validation scheme. One of \code{"CV1"}, \code{"CV2"}, or \code{"CV0"}.
#' @param sg A numeric vector of sigma values for the Laplacian kernel.
#'   Default is \code{c(0.001, 0.01, 0.05, 0.1)}.
#' @param var_threshold Minimum proportion of variance explained required for a principal component to be retained. Default is 0.01.
#' @param nIter Total number of iterations for the BGLR model. Default is 10000.
#' @param burnIn Number of burn-in iterations for the BGLR model. Default is 4000.
#' @param thin Thinning interval for the BGLR model. Default is 10.
#' @param save_xlsx A logical value indicating whether to save results in an Excel file. Default is \code{TRUE}.
#' @param file_name Character string specifying the name of the Excel file. If \code{NULL}, a default name is used.
#'
#' @return A data frame with the mean predictive capacity by model, CV scheme,
#' sigma value, number of selected PCs, and environment, accounting for the
#' GxE interaction model.
#'
#' @details
#' The model implemented is:
#' \deqn{y = Xb + Zg + Zi + e}
#' where \eqn{Xb} represents fixed environmental effects, \eqn{Zg} represents the
#' main genomic effect obtained from the KPCA Laplacian kernel, and \eqn{Zi}
#' represents the GxE interaction effect. The GxE kernel is computed as the
#' Hadamard product between the observation-level KPCA Laplacian genomic kernel
#' and the environmental relationship matrix.
#'
#' @examples
#' \dontrun{
#' results <- env_ge_pca_laplacian(
#'   SNPs = X,
#'   y = phen$yield,
#'   IDs = phen$genotype,
#'   env = phen$Env,
#'   CV = "CV2"
#' )
#' }
#'
#' @export

env_ge_pca_laplian <- function(SNPs, y, IDs, env,
                               EZ = NULL,
                               CV = c("CV1", "CV2", "CV0"),
                               sg = c(0.001, 0.01, 0.05, 0.1),
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

  GDec_list <- list()
  counter <- 1

  for (sigma_i in sg) {

    cat(
      "\nComputing KPCA Laplacian kernel for sigma =",
      sigma_i, "\n"
    )

    kpca_temp <- kernlab::kpca(
      x = SNPs,
      kernel = "laplacedot",
      kpar = list(sigma = sigma_i),
      features = 0
    )

    eig_vals <- kernlab::eig(kpca_temp)
    var_explained <- eig_vals / sum(eig_vals)

    nPC <- sum(var_explained > var_threshold)

    if (nPC < 2) {
      warning(
        paste(
          "Sigma", sigma_i,
          "selected fewer than 2 PCs. Skipping this sigma."
        )
      )
      next
    }

    cat("Number of PCs selected:", nPC, "\n")

    kpca_model <- kernlab::kpca(
      x = SNPs,
      kernel = "laplacedot",
      kpar = list(sigma = sigma_i),
      features = nPC
    )

    embedding <- kpca_model@rotated
    embedding <- as.matrix(embedding)

    Gn <- tcrossprod(embedding) / ncol(embedding)

    rownames(Gn) <- rownames(SNPs)
    colnames(Gn) <- rownames(SNPs)

    cat("Expanding KPCA Laplacian genomic kernel to observation level\n")

    G <- GZ %*% Gn %*% t(GZ)

    rownames(G) <- obs_names
    colnames(G) <- obs_names

    cat("Computing KPCA Laplacian GxE interaction kernel\n")

    GxE <- G * E

    rownames(GxE) <- obs_names
    colnames(GxE) <- obs_names

    cat("Eigen decomposition of KPCA Laplacian G\n")

    GDec <- eigen(G, symmetric = TRUE)

    GDec$values <- pmax(GDec$values, 0)
    rownames(GDec$vectors) <- rownames(G)

    cat("Eigen decomposition of KPCA Laplacian GxE\n")

    GxEDec <- eigen(GxE, symmetric = TRUE)

    GxEDec$values <- pmax(GxEDec$values, 0)
    rownames(GxEDec$vectors) <- rownames(GxE)

    GDec_list[[counter]] <- list(
      sigma = sigma_i,
      nPC = nPC,
      G_values = GDec$values,
      G_vectors = GDec$vectors,
      GxE_values = GxEDec$values,
      GxE_vectors = GxEDec$vectors
    )

    counter <- counter + 1
  }

  if (length(GDec_list) == 0) {
    stop("No KPCA Laplacian model was fitted. All sigma values selected fewer than 2 PCs.")
  }

  names(GDec_list) <- vapply(GDec_list, function(x) {
    paste0("sigma_", x$sigma)
  }, character(1))

  list_metrics <- list()

  for (i in seq_along(GDec_list)) {

    GDec_i <- GDec_list[[i]]

    cat(
      "\nRunning KPCA Laplacian + GxE model with sigma =",
      GDec_i$sigma,
      "| nPC =", GDec_i$nPC, "\n"
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
              Model = "KPCA_Laplacian_GxE",
              CV = CV,
              Sigma = GDec_i$sigma,
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
    Predictive_Capacity ~ Model + CV + Sigma + nPC + Environment,
    data = df_raw,
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  names(df_metrics)[names(df_metrics) == "Predictive_Capacity"] <- "pred"

  df_metrics <- df_metrics[order(-df_metrics$pred), ]

  if (save_xlsx) {

    if (is.null(file_name)) {
      file_name <- paste0("kpca_laplacian_gxe_", CV, ".xlsx")
    }

    writexl::write_xlsx(df_metrics, file_name)
  }

  return(df_metrics)
}
