#' Launch Shiny App for orthoVisualizer
#'
#' A function that launches the Shiny app for orthoVisualizer.
#' The purpose of this app is only to illustrate how a Shiny
#' app works. The code has been placed in \code{./inst/shiny-scripts}.
#'
#' @return No return value but open up a Shiny page.
#'
#' @examples
#' \dontrun{
#' orthoVisualizer::runorthoVisualizer()
#' }
#'
#' @references
#' Grolemund, G. (2015). Learn Shiny - Video Tutorials.
#' \href{https://shiny.rstudio.com/tutorial/}{Link}
#'
#' @importFrom shiny runApp
#' @import shinyalert
#' @export
runorthoVisualizer <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "orthoVisualizer")
  actionShiny <- shiny::runApp(appDir, display.mode = "normal")
  return(actionShiny)
}
# [END]
