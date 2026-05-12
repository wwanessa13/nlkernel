#' Launch the Genomic Prediction Shiny Application
#'
#' @description
#' This function launches the Shiny application included in the \code{nlkernel}
#' package. The application provides a user-friendly interface for applying
#' trained genomic prediction models to new individuals. Users can upload a
#' previously saved model and a marker matrix for new genotypes to obtain
#' predicted genomic values without refitting the model.
#'
#' @details
#' The application supports prediction from models trained using the final-model
#' functions available in the package, including single-environment models,
#' multi-environment models with genomic main effects, and multi-environment
#' models including genotype-by-environment interaction. For multi-environment
#' models, the input file must include both genotype identifiers and environment
#' labels.
#'
#' Please note that the model is based on a transformed marker matrix, and
#' therefore the same preprocessing steps must be applied to any new marker data
#' before prediction.
#'
#' @return
#' This function launches the Shiny application in the user's default web browser.
#' It does not return an R object.
#'
#' @examples
#' \dontrun{
#' run_app()
#' }
#'
#' @export
run_app <- function() {
  app_dir <- system.file("shiny", "app", package = "nlkernel")

  if (app_dir == "") {
    stop(
      "Could not find the Shiny app. Try reinstalling the package.",
      call. = FALSE
    )
  }

  shiny::runApp(app_dir, display.mode = "normal")
}
