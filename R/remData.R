#' remData
#'
#' Function to remove a drug from the data in the folder drug_data.
#'
#' @param data A list of data (the data saved in the folder drug_data)
#' @param old_drug The name of the drug you want to remove (has to be in the drug_data folder)
#'
#'
#' @export
#'

remData <- function(data, old_drug){
  if(!is.null(old_drug)){
    for (i in 1:length(data)){
      f <- file.path(getwd(), "drug_data", paste0(names(data)[i], ".xlsx"))
      OUT <- openxlsx::loadWorkbook(f)

      openxlsx::removeWorksheet(OUT, old_drug)

      openxlsx::saveWorkbook(OUT, f, overwrite = TRUE)
    }
  }
}

