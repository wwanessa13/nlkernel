#' Combined Kernel Models for Multi-Environment Genomic Prediction
#'
#' This function fits combined kernel models for multi-environment genomic
#' prediction using the RKHS framework implemented in the BGLR package.
#' It evaluates pairwise combinations of nonlinear kernels, GBLUP, and
#' the hyperbolic tangent (tanh) kernel.
#'
#' @param SNPs A numeric matrix of SNP genotypes, with genotypes in rows
#'   and markers in columns. Row names must correspond to genotype IDs.
#' @param y A numeric vector of phenotypic values.
#' @param IDs A vector of genotype IDs corresponding to each phenotypic observation.
#' @param env A vector of environment labels corresponding to each phenotypic observation.
#' @param EZ Optional matrix of fixed environmental effects. If NULL, it is
#'   created from env.
#' @param CV A character string specifying the cross-validation scheme:
#'   "CV1", "CV2", or "CV0".
#' @param polynomial_degree Degree parameter for the polynomial kernel.
#' @param polynomial_scale Scale parameter for the polynomial kernel.
#' @param polynomial_offset Offset parameter for the polynomial kernel.
#' @param laplacian_sigma Sigma parameter for the Laplacian kernel.
#' @param bessel_sigma Sigma parameter for the Bessel kernel.
#' @param bessel_order Order parameter for the Bessel kernel.
#' @param bessel_degree Degree parameter for the Bessel kernel.
#' @param anova_sigma Sigma parameter for the ANOVA kernel.
#' @param anova_degree Degree parameter for the ANOVA kernel.
#' @param gaussian_sigma Sigma parameter for the Gaussian/RBF kernel.
#' @param tanh_scale Scale parameter for the hyperbolic tangent kernel.
#' @param tanh_offset Offset parameter for the hyperbolic tangent kernel.
#' @param ploidy Ploidy level used in AGHmatrix::Gmatrix.
#' @param nIter Total number of iterations for the BGLR model.
#' @param burnIn Number of burn-in iterations for the BGLR model.
#' @param thin Thinning interval for the MCMC chain.
#' @param seed Integer value used to set the random seed for reproducibility.
#' @param save_xlsx Logical. If TRUE, saves results to an Excel file.
#' @param file_name Character string for the Excel file name.
#'
#' @return A data frame with the mean predictive capacity by kernel
#'   combination, CV scheme, and environment.
#'
#' @export

env_g_combinations <- function(
    SNPs, y, IDs, env,
    EZ = NULL,
    CV = c("CV1", "CV2", "CV0"),

    polynomial_degree = 2,
    polynomial_scale = 2,
    polynomial_offset = 2,

    laplacian_sigma = 0.01,

    bessel_sigma = 0.1,
    bessel_order = 1,
    bessel_degree = 2,

    anova_sigma = 0.1,
    anova_degree = 2,

    gaussian_sigma = 0.001,

    tanh_scale = 1,
    tanh_offset = 1,

    ploidy = 2,

    nIter = 10000,
    burnIn = 4000,
    thin = 10,

    seed = 123,

    save_xlsx = TRUE,
    file_name = NULL) {


  # ============================================================
  # PACOTES
  # ============================================================

  library(kernlab)
  library(BGLR)
  library(AGHmatrix)
  library(writexl)


  # ============================================================
  # CROSS-VALIDATION
  # ============================================================

  CV <- match.arg(CV)


  # ============================================================
  # DADOS
  # ============================================================

  SNPs <- as.matrix(SNPs)
  storage.mode(SNPs) <- "numeric"

  y <- as.numeric(y)
  IDs <- as.character(IDs)
  env <- as.character(env)


  # ============================================================
  # CHECAGEM DOS DADOS
  # ============================================================

  if (
    length(y) != length(IDs) ||
    length(y) != length(env)
  ) {

    stop(
      "The length of y, IDs, and env must be the same."
    )
  }


  if (is.null(rownames(SNPs))) {

    stop(
      "SNPs must have row names corresponding to genotype IDs."
    )
  }


  if (!all(unique(IDs) %in% rownames(SNPs))) {

    stop(
      "Some genotype IDs are not present in rownames(SNPs)."
    )
  }


  n <- length(y)

  uIDs <- unique(IDs)
  uenv <- unique(env)


  # ============================================================
  # EFEITO FIXO DE AMBIENTE
  # ============================================================

  if (is.null(EZ)) {

    EZ <- model.matrix(
      ~ factor(env) - 1
    )

    colnames(EZ) <- paste0(
      "Env_",
      uenv
    )
  }

  EZ <- as.matrix(EZ)


  if (nrow(EZ) != n) {

    stop(
      "EZ must have the same number of rows as the length of y."
    )
  }


  # ============================================================
  # DATAFRAME AUXILIAR
  # ============================================================

  Y <- data.frame(
    ID = IDs,
    Env = env,
    y = y
  )


  # ============================================================
  # CV1
  # ============================================================

  if (CV == "CV1") {

    set.seed(seed)

    n_folds <- 5

    fold_id <- rep(
      1:n_folds,
      length.out = length(uIDs)
    )

    fold_id <- sample(fold_id)

    names(fold_id) <- uIDs

    Y$Fold <- fold_id[Y$ID]
  }


  # ============================================================
  # CV2
  # ============================================================

  if (CV == "CV2") {

    set.seed(seed)

    n_folds <- 5

    Y$Fold <- NA

    for (id in uIDs) {

      idx <- which(
        Y$ID == id
      )

      ni <- length(idx)

      Y$Fold[idx] <- sample(
        1:n_folds,
        size = ni,
        replace = ni > n_folds
      )
    }
  }


  # ============================================================
  # CV0
  # ============================================================

  if (CV == "CV0") {

    set.seed(seed)

    n_folds_env <- length(uenv)

    fold_env <- sample(
      1:n_folds_env,
      size = n_folds_env
    )

    names(fold_env) <- uenv

    Y$Fold <- fold_env[Y$Env]
  }


  folds_run <- sort(
    unique(Y$Fold)
  )


  # ============================================================
  # MATRIZ DE INCIDÊNCIA GENOTÍPICA
  # ============================================================

  IDs_factor <- factor(
    IDs,
    levels = rownames(SNPs)
  )

  GZ <- model.matrix(
    ~ IDs_factor - 1
  )


  # ============================================================
  # KERNELS
  # ============================================================

  kernels <- list(

    # ----------------------------------------------------------
    # POLYNOMIAL
    # ----------------------------------------------------------

    polynomial = function(SNPs) {

      kernelMatrix(

        polydot(
          degree = polynomial_degree,
          scale = polynomial_scale,
          offset = polynomial_offset
        ),

        SNPs
      )
    },


    # ----------------------------------------------------------
    # LAPLACIAN
    # ----------------------------------------------------------

    laplacian = function(SNPs) {

      kernelMatrix(

        laplacedot(
          sigma = laplacian_sigma
        ),

        SNPs
      )
    },


    # ----------------------------------------------------------
    # BESSEL
    # ----------------------------------------------------------

    bessel = function(SNPs) {

      kernelMatrix(

        besseldot(
          sigma = bessel_sigma,
          order = bessel_order,
          degree = bessel_degree
        ),

        SNPs
      )
    },


    # ----------------------------------------------------------
    # ANOVA
    # ----------------------------------------------------------

    anova = function(SNPs) {

      kernelMatrix(

        anovadot(
          sigma = anova_sigma,
          degree = anova_degree
        ),

        SNPs
      )
    },


    # ----------------------------------------------------------
    # GAUSSIAN
    # ----------------------------------------------------------

    gaussian = function(SNPs) {

      kernelMatrix(

        rbfdot(
          sigma = gaussian_sigma
        ),

        SNPs
      )
    },


    # ----------------------------------------------------------
    # TANH
    # ----------------------------------------------------------

    tanh = function(SNPs) {

      kernelMatrix(

        tanhdot(
          scale = tanh_scale,
          offset = tanh_offset
        ),

        SNPs
      )
    },


    # ----------------------------------------------------------
    # GBLUP
    # ----------------------------------------------------------

    GBLUP = function(SNPs) {

      AGHmatrix::Gmatrix(
        SNPs,
        method = "VanRaden",
        ploidy = ploidy,
        maf = 0.05
      )
    }
  )


  # ============================================================
  # COMBINAÇÕES DOS KERNELS
  # ============================================================

  comb_list <- list(

    c("polynomial", "gaussian"),
    c("polynomial", "bessel"),
    c("polynomial", "laplacian"),
    c("polynomial", "anova"),
    c("polynomial", "tanh"),
    c("polynomial", "GBLUP"),

    c("gaussian", "bessel"),
    c("gaussian", "laplacian"),
    c("gaussian", "anova"),
    c("gaussian", "tanh"),
    c("gaussian", "GBLUP"),

    c("bessel", "laplacian"),
    c("bessel", "anova"),
    c("bessel", "tanh"),
    c("bessel", "GBLUP"),

    c("laplacian", "anova"),
    c("laplacian", "tanh"),
    c("laplacian", "GBLUP"),

    c("anova", "tanh"),
    c("anova", "GBLUP"),

    c("tanh", "GBLUP")
  )


  # ============================================================
  # KERNELS EXPANDIDOS E DECOMPOSTOS
  # ============================================================

  KDec_by_kernel <- list()


  for (kname in names(kernels)) {

    cat(
      "\nComputing kernel:",
      kname,
      "\n"
    )


    Gn <- kernels[[kname]](SNPs)

    Gn <- as.matrix(Gn)


    rownames(Gn) <- rownames(SNPs)
    colnames(Gn) <- rownames(SNPs)


    # ----------------------------------------------------------
    # EXPANDIR PARA NÍVEL DAS OBSERVAÇÕES
    # ----------------------------------------------------------

    G <- GZ %*%
      Gn %*%
      t(GZ)


    # ----------------------------------------------------------
    # DECOMPOSIÇÃO ESPECTRAL
    # ----------------------------------------------------------

    GDec <- eigen(
      G,
      symmetric = TRUE
    )


    KDec_by_kernel[[kname]] <- list(

      values = pmax(
        GDec$values,
        0
      ),

      vectors = GDec$vectors
    )
  }


  # ============================================================
  # RESULTADOS
  # ============================================================

  list_metrics <- list()


  # ============================================================
  # LOOP SOBRE AS COMBINAÇÕES
  # ============================================================

  for (comb in comb_list) {

    kname <- paste(
      comb,
      collapse = "_"
    )


    cat(
      "\n============================================\n"
    )

    cat(
      "Running combination:",
      kname,
      "\n"
    )

    cat(
      "============================================\n"
    )


    # ==========================================================
    # ETA DOS KERNELS
    # ==========================================================

    ETA_kernels <- lapply(
      comb,
      function(k) {

        list(

          V = KDec_by_kernel[[k]]$vectors,

          d = KDec_by_kernel[[k]]$values,

          model = "RKHS"
        )
      }
    )


    # ==========================================================
    # EFEITO AMBIENTAL + KERNELS
    # ==========================================================

    ETA_i <- c(

      list(
        list(
          X = EZ,
          model = "FIXED"
        )
      ),

      ETA_kernels
    )


    # ==========================================================
    # CROSS-VALIDATION
    # ==========================================================

    for (fold in folds_run) {

      cat(
        "  Processing fold",
        fold,
        "\n"
      )


      # --------------------------------------------------------
      # TESTE
      # --------------------------------------------------------

      testing <- which(
        Y$Fold == fold
      )


      # --------------------------------------------------------
      # OCULTAR FENÓTIPOS
      # --------------------------------------------------------

      yNA <- y

      yNA[testing] <- NA


      # --------------------------------------------------------
      # BGLR
      # --------------------------------------------------------

      fm <- BGLR::BGLR(

        y = yNA,

        ETA = ETA_i,

        nIter = nIter,

        burnIn = burnIn,

        thin = thin,

        verbose = FALSE
      )


      # --------------------------------------------------------
      # PREDIÇÕES
      # --------------------------------------------------------

      yHat <- fm$yHat


      # ========================================================
      # CAPACIDADE PREDITIVA POR AMBIENTE
      # ========================================================

      for (a in uenv) {

        idx_env <- which(
          env == a
        )

        join <- intersect(
          idx_env,
          testing
        )


        if (length(join) > 1) {

          if (
            sd(
              yHat[join],
              na.rm = TRUE
            ) > 0 &&
            sd(
              y[join],
              na.rm = TRUE
            ) > 0
          ) {

            cor_val <- cor(

              yHat[join],

              y[join],

              use = "complete.obs"
            )

          } else {

            cor_val <- NA
          }


          list_metrics[[length(list_metrics) + 1]] <- data.frame(
            Model = "Combined_Kernels",
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


  # ============================================================
  # COMBINAR RESULTADOS
  # ============================================================

  df_raw <- do.call(
    rbind,
    list_metrics
  )


  # ============================================================
  # MÉDIA ENTRE FOLDS
  # ============================================================

  df_metrics <- aggregate(

    Predictive_Capacity ~

      Model +
      CV +
      Kernel +
      Environment,

    data = df_raw,

    FUN = function(x) {

      if (all(is.na(x))) {

        return(NA_real_)

      } else {

        return(
          mean(
            x,
            na.rm = TRUE
          )
        )
      }
    }
  )


  # ============================================================
  # RENOMEAR
  # ============================================================

  names(df_metrics)[
    names(df_metrics) ==
      "Predictive_Capacity"
  ] <- "pred"


  # ============================================================
  # ORDENAR
  # ============================================================

  df_metrics <- df_metrics[
    order(
      -df_metrics$pred
    ),
  ]


  # ============================================================
  # SALVAR
  # ============================================================

  if (save_xlsx) {

    if (is.null(file_name)) {

      file_name <- paste0(
        "Combined_Kernels_",
        CV,
        ".xlsx"
      )
    }


    writexl::write_xlsx(
      df_metrics,
      file_name
    )
  }


  # ============================================================
  # RETORNAR
  # ============================================================

  return(df_metrics)
}
