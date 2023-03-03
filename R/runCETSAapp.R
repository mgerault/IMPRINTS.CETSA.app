#' @title Run the app to use mineCETSAapp package
#'
#' @description
#'  \code{runCETSAapp} works exactly the same as runExample from \code{\link{shiny}} package.
#'
#' @export
runCETSAapp <- function() {

 appDir <- system.file("shiny-examples", "myapp", package = "mineCETSAapp")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mineCETSAapp`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
