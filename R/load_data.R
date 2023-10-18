#' load_data
#'
#' Function to load the data saved in the drug_data folder. Will import in the drug_data format.
#'
#'
#'
#' @return The list drug_data
#'
#' @export

load_data <- function() {
  # Read all the files into a list
  files <- list.files("drug_data", full.names = TRUE)
  if(length(files)){
    names(files) <- str_remove_all(list.files("drug_data", full.names = FALSE), "\\.xlsx")
    files <- as.list(files)
    data <- lapply(files, import_list)

    message("Data loaded with success !")
    return(data)
  }
  else{
    message("No folder named 'drug_data' was found; returned NULL")
    return(NULL)
  }
}
