#' saveData
#'
#' Function for saving new data in the folder drug_data.
#' If in the drug_data format, you can add new data from another drug.
#'
#' @param data A list of data that you want to save
#' @param new_add A list of new data you want to add to the existing xlsx files
#' @param new_drug A character vector which is the name of the new drug you're adding
#'
#'
#'
#' @export
#'
#'
saveData <- function(data, new_add = NULL, new_drug = NULL) {
  if(!file.exists("drug_data")){
    path <- getwd()
    dir.create(file.path(path, "drug_data"))

    # Write the file to the local system
    export_list(data,
                file.path(getwd(), "drug_data", paste0("%s", ".xlsx"))
    )
  }

  if(!is.null(new_add) & !is.null(new_drug)){
    if(length(new_add) == length(data)){
      for (i in names(data)){
        f <- file.path(getwd(), "drug_data", paste0(i, ".xlsx"))
        OUT <- openxlsx::loadWorkbook(f)

        openxlsx::addWorksheet(OUT, new_drug)
        openxlsx::writeData(OUT, sheet = new_drug, new_add[[i]])

        openxlsx::saveWorkbook(OUT, f, overwrite = TRUE)
      }
    }
  }

  message("Data saved with success !")
}



