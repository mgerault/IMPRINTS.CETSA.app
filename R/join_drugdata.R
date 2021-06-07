#' join_drugdata
#'
#' Function to join a list of data.frames. It is the same as the join_all function
#' from the plyr package but allow to use the full_join function from the dplyr package
#' and so to add suffixes to duplicated names not selecting by the argument by.
#'
#' @param dfs A list of data frames
#' @param by A character vector of variables to join by.
#'           If NULL, the default, full_join() will perform a natural join, using all variables in common across the data frames.
#'           A message lists the variables so that you can check they're correct; suppress the message by supplying
#'           by explicitly.
#'
#'           To join by different variables, use a named vector. For example, by = c("a" = "b") will match x$a to y$b.
#'           To join by multiple variables, use a vector with length > 1. For example, by = c("a", "b") will match
#'           x$a to y$a and x$b to y$b. Use a named vector to match different variables in x and y.
#'           For example, by = c("a" = "b", "c" = "d") will match x$a to y$b and x$c to y$d.
#'           x being the data frame named 'joined' and y the data frame i from dfs.
#'
#'
#' @return The joined dataframe
#'
#' @seealso \code{\link{dplyr}}
#'
#' @export
#'


join_drugdata <- function (dfs, by = NULL){
  if (length(dfs) == 1)
    return(dfs[[1]])

  joined <- dfs[[1]]
  for (i in 2:length(dfs)) {
    joined <- dplyr::full_join(joined, dfs[[i]], by = by)
  }
  return(joined)
}
