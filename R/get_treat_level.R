#' get_treat_level
#'
#' Function to get the treatment level of the CETSA data.
#'
#' @param data An IMPRINTS-CETSA dataset; after imprints_normalization or imprints_caldiff for example.
#'
#' @return A character vector which contains all the treatment level from the data
#'
#' @export
#'
#'
get_treat_level <- function(data){
  lev <- names(data)
  lev <- grep("^\\d{2}", lev, value = TRUE)
  lev <- strsplit(lev, "_")
  lev <- lapply(lev, function(x) x[length(x)]) #take the last element after a '_', so in theory the treatment
  lev <- unique(unlist(lev))

  return(lev)
}
