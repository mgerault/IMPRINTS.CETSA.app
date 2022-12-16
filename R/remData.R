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
    f <- file.path(getwd(), "drug_data")
    if(!file.exists(f)){
      path <- getwd()
      dir.create(file.path(path, "drug_data"))

      # Write the file to the local system
      export_list(data,
                  file.path(getwd(), "drug_data", paste0("%s", ".xlsx"))
      )
    }
    for (i in 1:length(data)){
      f <- file.path(getwd(), "drug_data", paste0(names(data)[i], ".xlsx"))
      OUT <- openxlsx::loadWorkbook(f)

      openxlsx::removeWorksheet(OUT, old_drug)

      openxlsx::saveWorkbook(OUT, f, overwrite = TRUE)
    }
  }
}

