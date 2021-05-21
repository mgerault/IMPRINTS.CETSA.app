#' loadData
#'
#' Function to load the data saved in the drug_data folder. Will import in the drug_data format.
#'
#'
#'
#' @return The list drug_data
#'
#' @export

loadData <- function() {
  # Read all the files into a list
  files <- list.files("drug_data", full.names = TRUE)
  names(files) <- str_remove_all(list.files("drug_data", full.names = FALSE), "\\.xlsx")
  files <- as.list(files)
  data <- lapply(files, import_list)

  # Concatenate all data together into one data.frame
  data$treat_level <- lapply(data$treat_level, function(x) {x <- append(unname(unlist(as.list(x))),
                                                                        names(x), after = 0);x})
  message("Data loaded with success !")
  return(data)
}
