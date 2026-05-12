library(shiny)
library(readr)
library(dplyr)
library(DT)
library(kernlab)
library(MASS)

predict_saved_bglr <- function(model_object,
                               SNPs_new,
                               IDs_new = NULL,
                               env_new = NULL) {

  if (!requireNamespace("kernlab", quietly = TRUE)) {
    stop("Package 'kernlab' is required.")
  }

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required.")
  }

  SNPs_train <- model_object$SNPs_train

  SNPs_new <- as.matrix(SNPs_new)
  storage.mode(SNPs_new) <- "numeric"

  if (ncol(SNPs_new) != ncol(SNPs_train)) {
    stop("The number of markers in SNPs_new must match SNPs_train.")
  }

  model_class <- class(model_object)[1]

  fit <- model_object$fit
  mu <- fit$mu
  hp <- model_object$hyperparameters

  ETA_components <- model_object$ETA_components

  if (is.null(ETA_components)) {
    ETA_components <- model_object$model
  }

  make_cross_kernel <- function(component) {

    if (component == "gaussian") {

      return(as.matrix(
        kernlab::kernelMatrix(
          kernlab::rbfdot(sigma = hp$gaussian_sigma),
          SNPs_new,
          SNPs_train
        )
      ))
    }

    if (component == "laplacian") {

      return(as.matrix(
        kernlab::kernelMatrix(
          kernlab::laplacedot(sigma = hp$laplacian_sigma),
          SNPs_new,
          SNPs_train
        )
      ))
    }

    if (component == "polynomial") {

      return(as.matrix(
        kernlab::kernelMatrix(
          kernlab::polydot(
            degree = hp$polynomial_degree,
            scale = hp$polynomial_scale,
            offset = hp$polynomial_offset
          ),
          SNPs_new,
          SNPs_train
        )
      ))
    }

    if (component == "bessel") {

      return(as.matrix(
        kernlab::kernelMatrix(
          kernlab::besseldot(
            sigma = hp$bessel_sigma,
            order = hp$bessel_order,
            degree = hp$bessel_degree
          ),
          SNPs_new,
          SNPs_train
        )
      ))
    }

    if (component == "gblup") {

      M_train <- SNPs_train
      M_new <- SNPs_new

      p <- colMeans(M_train, na.rm = TRUE) / 2

      Z_train <- sweep(M_train, 2, 2 * p, "-")
      Z_new <- sweep(M_new, 2, 2 * p, "-")

      denom <- 2 * sum(p * (1 - p), na.rm = TRUE)

      return(as.matrix(
        Z_new %*% t(Z_train) / denom
      ))
    }

    if (component == "pca") {

      pca_obj <- model_object$pca_objects[[component]]

      if (is.null(pca_obj)) {
        pca_obj <- model_object$pca_objects[[which(ETA_components == component)]]
      }

      nPC <- model_object$nPC[[component]]

      if (is.null(nPC)) {
        nPC <- model_object$nPC[[which(ETA_components == component)]]
      }

      scores_train <- pca_obj$x[, seq_len(nPC), drop = FALSE]

      scores_new <- predict(
        pca_obj,
        newdata = SNPs_new
      )[, seq_len(nPC), drop = FALSE]

      return(as.matrix(
        scores_new %*% t(scores_train) / ncol(scores_train)
      ))
    }

    if (component %in% c("pca_gaussian", "pca_laplacian", "pca_polynomial", "pca_bessel")) {

      kpca_obj <- model_object$pca_objects[[component]]

      if (is.null(kpca_obj)) {
        kpca_obj <- model_object$pca_objects[[which(ETA_components == component)]]
      }

      scores_train <- as.matrix(kpca_obj@rotated)

      scores_new <- as.matrix(
        predict(kpca_obj, SNPs_new)
      )

      return(as.matrix(
        scores_new %*% t(scores_train) / ncol(scores_train)
      ))
    }

    stop("Unknown model component: ", component)
  }

  predict_single <- function() {

    y_pred <- rep(mu, nrow(SNPs_new))

    for (i in seq_along(ETA_components)) {

      component <- ETA_components[i]

      K_train <- model_object$kernels[[component]]

      if (is.null(K_train)) {
        K_train <- model_object$kernels[[i]]
      }

      K_new <- make_cross_kernel(component)

      u_train <- fit$ETA[[component]]$u

      if (is.null(u_train)) {
        u_train <- fit$ETA[[i]]$u
      }

      alpha <- MASS::ginv(K_train) %*% u_train

      y_pred <- y_pred + as.vector(K_new %*% alpha)
    }

    return(y_pred)
  }

  predict_env_g <- function() {

    if (is.null(IDs_new) || is.null(env_new)) {
      stop("For multi-environment models, IDs_new and env_new must be provided.")
    }

    IDs_new <- as.character(IDs_new)
    env_new <- as.character(env_new)

    if (length(IDs_new) != nrow(SNPs_new)) {
      stop("IDs_new must have the same length as the number of rows in SNPs_new.")
    }

    if (length(env_new) != nrow(SNPs_new)) {
      stop("env_new must have the same length as the number of rows in SNPs_new.")
    }

    env_train_levels <- sub("^Env_", "", colnames(model_object$EZ))

    if (!all(env_new %in% env_train_levels)) {
      stop(
        "Some environments in env_new were not present in the training data. ",
        "Available environments are: ",
        paste(env_train_levels, collapse = ", ")
      )
    }

    EZ_new <- model.matrix(
      ~ factor(env_new, levels = env_train_levels) - 1
    )

    colnames(EZ_new) <- paste0("Env_", env_train_levels)

    beta_env <- fit$ETA$Environment$b

    y_pred <- as.vector(mu + EZ_new %*% beta_env)

    for (i in seq_along(ETA_components)) {

      component <- ETA_components[i]

      K_obs_train <- model_object$observation_kernels[[component]]

      if (is.null(K_obs_train)) {
        K_obs_train <- model_object$observation_kernels[[i]]
      }

      K_new_gen <- make_cross_kernel(component)

      GZ_train <- model_object$GZ

      K_new_obs <- K_new_gen %*% t(GZ_train)

      u_train <- fit$ETA[[component]]$u

      if (is.null(u_train)) {
        u_train <- fit$ETA[[i + 1]]$u
      }

      alpha <- MASS::ginv(K_obs_train) %*% u_train

      y_pred <- y_pred + as.vector(K_new_obs %*% alpha)
    }

    return(y_pred)
  }

  predict_env_gxe <- function() {

    if (is.null(IDs_new) || is.null(env_new)) {
      stop("For multi-environment GxE models, IDs_new and env_new must be provided.")
    }

    IDs_new <- as.character(IDs_new)
    env_new <- as.character(env_new)

    if (length(IDs_new) != nrow(SNPs_new)) {
      stop("IDs_new must have the same length as the number of rows in SNPs_new.")
    }

    if (length(env_new) != nrow(SNPs_new)) {
      stop("env_new must have the same length as the number of rows in SNPs_new.")
    }

    env_train <- model_object$env_train
    env_train_levels <- sub("^Env_", "", colnames(model_object$EZ))

    if (!all(env_new %in% env_train_levels)) {
      stop(
        "Some environments in env_new were not present in the training data. ",
        "Available environments are: ",
        paste(env_train_levels, collapse = ", ")
      )
    }

    EZ_new <- model.matrix(
      ~ factor(env_new, levels = env_train_levels) - 1
    )

    colnames(EZ_new) <- paste0("Env_", env_train_levels)

    beta_env <- fit$ETA$Environment$b

    y_pred <- as.vector(mu + EZ_new %*% beta_env)

    E_new_train <- outer(
      env_new,
      env_train,
      FUN = "=="
    )

    E_new_train <- matrix(
      as.numeric(E_new_train),
      nrow = length(env_new),
      ncol = length(env_train)
    )

    for (i in seq_along(ETA_components)) {

      component <- ETA_components[i]

      K_obs_train <- model_object$observation_kernels[[component]]

      if (is.null(K_obs_train)) {
        K_obs_train <- model_object$observation_kernels[[i]]
      }

      K_gxe_train <- model_object$gxe_kernels[[component]]

      if (is.null(K_gxe_train)) {
        K_gxe_train <- model_object$gxe_kernels[[i]]
      }

      K_new_gen <- make_cross_kernel(component)

      GZ_train <- model_object$GZ

      K_new_obs <- K_new_gen %*% t(GZ_train)

      K_new_gxe <- K_new_obs * E_new_train

      u_g <- fit$ETA[[paste0(component, "_G")]]$u
      u_gxe <- fit$ETA[[paste0(component, "_GxE")]]$u

      alpha_g <- MASS::ginv(K_obs_train) %*% u_g
      alpha_gxe <- MASS::ginv(K_gxe_train) %*% u_gxe

      y_pred <- y_pred +
        as.vector(K_new_obs %*% alpha_g) +
        as.vector(K_new_gxe %*% alpha_gxe)
    }

    return(y_pred)
  }

  if (model_class == "train_final_single") {
    return(predict_single())
  }

  if (model_class == "train_final_env") {
    return(predict_env_g())
  }

  if (model_class == "train_final_env_gxe") {
    return(predict_env_gxe())
  }

  stop(
    "Unknown model object class. Expected one of: ",
    "train_final_single, train_final_env, or train_final_env_gxe."
  )
}


ui <- fluidPage(

  titlePanel("Genomic Prediction Application"),

  sidebarLayout(

    sidebarPanel(

      fileInput(
        inputId = "model_file",
        label = "Upload trained model (.rds)",
        accept = ".rds"
      ),

      fileInput(
        inputId = "snp_file",
        label = "Upload marker data for new individuals (.csv)",
        accept = c(".csv")
      ),

      actionButton(
        inputId = "predict_btn",
        label = "Predict genomic values",
        class = "btn-primary"
      ),

      br(),
      br(),

      downloadButton(
        outputId = "download_predictions",
        label = "Download predictions"
      ),

      br(),
      br(),

      tags$hr(),

      h5("Input format"),

      p("For single-environment models, the CSV file must contain ID and marker columns."),

      tags$code("ID,SNP1,SNP2,SNP3"),

      br(),
      br(),

      p("For multi-environment models, the CSV file must contain ID, Env, and marker columns."),

      tags$code("ID,Env,SNP1,SNP2,SNP3")
    ),

    mainPanel(

      h4("Model information"),
      verbatimTextOutput("model_info"),

      br(),

      h4("Prediction results"),
      DTOutput("pred_table")
    )
  )
)


server <- function(input, output, session) {

  final_model <- reactive({

    req(input$model_file)

    readRDS(input$model_file$datapath)
  })

  output$model_info <- renderPrint({

    req(final_model())

    model <- final_model()

    cat("Model class:", class(model)[1], "\n")

    if (!is.null(model$model)) {
      cat("Model:", paste(model$model, collapse = " + "), "\n")
    }

    if (!is.null(model$ETA_components)) {
      cat("ETA components:", paste(model$ETA_components, collapse = " + "), "\n")
    }

    if (!is.null(model$env_train)) {
      cat("Training environments:", paste(unique(model$env_train), collapse = ", "), "\n")
    }

    cat("\nHyperparameters:\n")
    print(model$hyperparameters)
  })

  new_snps_data <- reactive({

    req(input$snp_file)
    req(final_model())

    dat <- readr::read_csv(
      input$snp_file$datapath,
      show_col_types = FALSE
    )

    model <- final_model()
    model_class <- class(model)[1]

    if (model_class == "train_final_single") {

      if (!("ID" %in% colnames(dat))) {
        stop("For single-environment models, the file must contain a column named ID.")
      }

      IDs <- dat$ID

      SNPs_new <- dat |>
        dplyr::select(-ID)

      SNPs_new <- as.matrix(SNPs_new)
      storage.mode(SNPs_new) <- "numeric"

      return(list(
        IDs = IDs,
        env = NULL,
        SNPs_new = SNPs_new
      ))
    }

    if (model_class %in% c("train_final_env", "train_final_env_gxe")) {

      if (!all(c("ID", "Env") %in% colnames(dat))) {
        stop("For multi-environment models, the file must contain columns named ID and Env.")
      }

      IDs <- dat$ID
      env <- dat$Env

      SNPs_new <- dat |>
        dplyr::select(-ID, -Env)

      SNPs_new <- as.matrix(SNPs_new)
      storage.mode(SNPs_new) <- "numeric"

      return(list(
        IDs = IDs,
        env = env,
        SNPs_new = SNPs_new
      ))
    }

    stop(
      "Unknown model object class. Expected one of: ",
      "train_final_single, train_final_env, or train_final_env_gxe."
    )
  })

  predictions <- eventReactive(input$predict_btn, {

    req(final_model())
    req(new_snps_data())

    model <- final_model()
    new_data <- new_snps_data()

    y_pred <- predict_saved_bglr(
      model_object = model,
      SNPs_new = new_data$SNPs_new,
      IDs_new = new_data$IDs,
      env_new = new_data$env
    )

    if (is.null(new_data$env)) {

      data.frame(
        ID = new_data$IDs,
        Predicted = y_pred
      )

    } else {

      data.frame(
        ID = new_data$IDs,
        Env = new_data$env,
        Predicted = y_pred
      )
    }
  })

  output$pred_table <- renderDT({

    req(predictions())

    datatable(
      predictions(),
      options = list(
        pageLength = 10,
        scrollX = TRUE
      )
    )
  })

  output$download_predictions <- downloadHandler(

    filename = function() {
      "genomic_predictions.csv"
    },

    content = function(file) {
      readr::write_csv(predictions(), file)
    }
  )
}


shinyApp(ui = ui, server = server)
