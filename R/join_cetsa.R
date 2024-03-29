#' join_cetsa
#'
#' Function to join a list of data frame from IMPRINTS-CETSA experiment and renamed the column.
#' (Allow to avoid the adding of ".x" or ".y" at the end of the column names)
#'
#' @param list_data A list of data frames. Usually the data frames are the output from the imprints_caldiff function.
#' @param new_names The new character element you want to add at the end of the columns name.
#'
#' @return The joined data frame (by 'id' and 'description')
#'
#' @export
#'

join_cetsa <- function(list_data, new_names = c("1h", "6h")){
  if(length(grep("_|/", new_names))){
    stop("Please enter a valid suffix. The character '_' and '/' are not allowed.")
  }
  if(inherits(list_data, "data.frame")){
    list_data <- list(list_data)
  }
  n <- length(list_data)

  if(length(new_names) != n){
    stop("You need to provide as many new_names as you have data.frame in your list !")
  }

  list_data <- mapply(function(x,y){ colnames(x)[!(colnames(x) %in% c("id", "description"))] <-
                                       paste0(colnames(x)[!(colnames(x) %in% c("id", "description"))], y); x},
                      list_data, new_names, SIMPLIFY = FALSE)

  is_OX <- unlist(lapply(list_data,
                         function(z){
                           z <- z$description
                           z <- grepl("OX=\\d{1,}", z);
                           any(z)
                         })
  )
  if(any(is_OX)){   # if in description OX is precised and for other not, remove it
    list_data[is_OX] <- lapply(list_data[is_OX],
                                 function(z){
                                   z$description <- gsub("OX=\\d{1,} ", "", z$description);
                                   z
                              })
  }
  df <- plyr::join_all(list_data, by = c("id", "description"), type = "full")

  return(df)
}


