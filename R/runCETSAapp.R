#' @title Run the app to use IMPRINTS.CETSA.app package
#'
#' @description
#'  \code{runCETSAapp} works exactly the same as runExample from \code{\link{shiny}} package.
#'
#' @export
runCETSAapp <- function() {

 appDir <- system.file("shiny-examples", "myapp", package = "IMPRINTS.CETSA.app")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `IMPRINTS.CETSA.app`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
