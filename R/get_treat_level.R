#' get_treat_level
#'
#' Function to get the treatment level of the CETSA data.
#'
#' @param data Dataset after ms_2D_caldiff
#'
#' @return A character vector which contains all the treatment level from the data
#'
#' @export
#'
#'
get_treat_level <- function(data){
  lev <- names(data)
  lev <- str_subset(lev, "^\\d{2}")
  lev <- str_split(lev, "_")
  lev <- lapply(lev, function(x) x[length(x)]) #take the last element after a '_', so in theory the condition
  lev <- unique(unlist(lev))

  return(lev)
}
