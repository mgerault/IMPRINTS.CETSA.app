#' imprints_read_peptides
#'
#' Function to read your peptides files and concatenate them in one data frame.
#'
#' @param peptides_files The path to the peptides data. Can be one or more files.
#' @param treatment A character vector that contains the treatment applied for each channel like 'B1_Vehicle' for example.
#' @param temperatures A character vector that contains the temperatures used for each peptides file like '37C' for example.
#' @param proteins Either a data frame or a file that has the 2 columns 'id' and 'description' corresponding
#'                 to the proteins you want to keep (usually proteins after \code{imprints_rearrange} or \code{imprints_caldiff})
#' @param modification_tokeep A character vector that contains the peptides modifications you want to keep. Default is 'TMT'
#' @param dataset_name The name of your dataset
#'
#' @return The peptides data frame
#'
#' @export
#'
#'

imprints_read_peptides <- function(peptides_files, treatment, temperatures,
                                   proteins, modification_tokeep = c("TMT"),
                                   dataset_name = "Venetoclax6h"){
  if(!("Mix" %in% treatment)){
    message("Error: The channel 'Mix' should be in your treatment.")
    return()
  }
  n_peptides_files <- length(peptides_files)
  if(length(temperatures) != n_peptides_files){
    message("Error: The number of temperatures doesn't match the number of peptides files !")
    return()
  }
  if(is.character(temperatures)){
    if(!all(stringr::str_detect(temperatures, "^\\d{1,2}C$"))){
      message("Error: Your temperatures should be in the format '1 or 2 digits' followed by capital 'C', like '37C' for example.")
      return()
    }
  }
  else{
    message("Error: Your temperatures should be a character vector !")
    return()
  }
  peptides_files <- as.list(peptides_files)
  names(peptides_files) <- temperatures

  message("Reading peptides files")
  peptides_dataset <- list()
  for(p in 1:n_peptides_files){
    peptides_temperature <- names(peptides_files)[p]
    peptides <- peptides_files[[p]]

    peptides <- readr::read_tsv(peptides)
    abundance_columns <- stringr::str_subset(colnames(peptides), "^Abundances \\(Grouped\\):")
    if(!length(abundance_columns)){
      abundance_columns <- stringr::str_subset(colnames(peptides), "^Abundance:")
    }
    peptides <- peptides[,c("Master Protein Accessions", "Positions in Master Proteins",
                            "Annotated Sequence", "Modifications",
                            abundance_columns)]
    peptides <- peptides[order(peptides$`Master Protein Accessions`),]
    colnames(peptides)[5:ncol(peptides)] <- paste0(peptides_temperature, "_", treatment)
    peptides <- peptides[, -stringr::str_which(colnames(peptides), "_Mix$")]

    peptides_dataset[[peptides_temperature]] <- peptides
  }
  message("Concatening all peptides files")
  peptides_dataset <- plyr::join_all(peptides_dataset,
                                     by = c("Master Protein Accessions", "Positions in Master Proteins",
                                            "Annotated Sequence", "Modifications"),
                                     type = "full")


  # Take only the same proteins and adding protein description
  if(is.character(proteins)){
    if(file.exists(proteins)){
      proteins <- readr::read_tsv(proteins)
      if(all(c("id", "description") %in% colnames(proteins))){
        proteins <- proteins[,c("id", "description")]
      }
      else{
        message("Error: Your proteins dataset should contains the columns 'id' and 'description'")
        return()
      }
    }
    else{
      message("Error: Check your spelling, the protein file wasn't found")
      return()
    }
  }
  else{
    if(inherits(proteins, "data.frame")){
      if(all(c("id", "description") %in% colnames(proteins))){
        proteins <- proteins[,c("id", "description")]
      }
      else{
        message("Error: Your proteins dataset should contains the columns 'id' and 'description'")
        return()
      }
    }
    else if(is.null(proteins)){
      proteins <- data.frame(id = peptides_dataset$`Master Protein Accessions`,
                             description = NA)
    }
    else{
      message("Error: Your proteins data should either be NULL, a file or a data.frame that contains the columns 'id' and 'description'")
      return()
    }
  }
  colnames(proteins)[1] <- "Master Protein Accessions"
  peptides_dataset <- left_join(proteins, peptides_dataset, by = "Master Protein Accessions")

  if(all(stringr::str_length(modification_tokeep) > 0)){
    modification_tokeep <- paste(modification_tokeep, collapse = "|")

    ### remove some modifications, only keeping TMT modif
    modif <- unique(peptides_dataset$Modifications)
    modif  <- strsplit(modif, "; |(?>\\[.*?\\].*?\\K(; |$))", perl=T) # split according '; ' but not when inside []
    modif <- unique(unlist(modif))
    modif <- stringr::str_split(stringr::str_remove(modif, "\\d{1,}x"), " ")
    modif <- unique(unlist(lapply(modif, function(x) x[1])))
    modif <- modif[!stringr::str_detect(modif, modification_tokeep)]
    modif <- paste(modif, collapse = "|")

    peptides_dataset <- peptides_dataset[!stringr::str_detect(peptides_dataset$Modifications, modif),]
  }

  message("Saving files")
  # remove peptides with only NAs
  peptides_dataset <- peptides_dataset[apply(peptides_dataset, 1,
                                             function(x) sum(is.na(x)) < ncol(peptides_dataset) - 5),]


 readr::write_tsv(peptides_dataset,
                  file = paste0(format(Sys.time(), "%y%m%d_%H%M_"),
                                "peptides_", dataset_name, ".txt")
                  )
 return(peptides_dataset)
}
