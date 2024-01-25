#' imprints_read_peptides
#'
#' Function to read your peptides files and concatenate them in one data frame.
#'
#' @param peptides_files The path to the PeptidesGroups file from PD. Can be one or more files.
#' @param treatment A character vector that contains the treatment applied for each channel like 'B1_Vehicle' for example.
#' @param temperatures A character vector that contains the temperatures used for each peptides file like '37C' for example.
#' @param proteins Either a data frame or a file that has the 2 columns 'id' and 'description' corresponding
#'                 to the proteins you want to keep (usually proteins after \code{imprints_rearrange} or \code{imprints_caldiff}).
#'                 If not specified, no filtering is applied and you must have the column 'Master Protein Descriptions' in your PeptideGroup file(s).
#' @param modification_torm A character vector that contains the peptides modifications you want to remove. Default is NULL
#' @param dataset_name The name of your dataset
#'
#' @return The peptides data frame
#'
#' @export
#'
#'

imprints_read_peptides <- function(peptides_files, treatment, temperatures,
                                   proteins = NULL, modification_torm = NULL,
                                   dataset_name = "Venetoclax6h"){
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

    peptides <- readr::read_tsv(peptides, show_col_types = FALSE)
    abundance_columns <- grep("^Abundance:", colnames(peptides), value = TRUE)
    if(!length(abundance_columns)){
      abundance_columns <- grep("^Abundances \\(Grouped\\):", colnames(peptides), value = TRUE)
      if(!length(abundance_columns)){
        message(paste("Error: The file",  peptides_files[[p]],
                      "doesn't contain any abundance information ! It should contain columns starting with 'Abundance:' or 'Abundances (Grouped):'"))
        return()
      }
    }
    info_col <- c("Master Protein Accessions", "Master Protein Descriptions",
                  "Positions in Master Proteins",
                  "Annotated Sequence", "Modifications")
    miss_col <- info_col[!(info_col %in% colnames(peptides))]
    if(length(miss_col)){
      if("Master Protein Descriptions" %in% miss_col){
        message(paste("Warning: 'Master Protein Descriptions' wasn't found in", peptides_files[[p]]
                      )
                )
        if(is.null(proteins)){
          message(paste(peptides_files[[p]], "will not have the protein description information since you also didn't specify the 'protein' parameter."))
        }
        info_col <- info_col[-2]
        miss_col <- info_col[!(info_col %in% colnames(peptides))]
      }
      if(length(miss_col)){
        message(paste("Error:", paste(miss_col, collapse = ", "),
                      ifelse(length(miss_col) > 1, "are", "is"),
                      "missing in", peptides_files[[p]], "!")
                )
        return()
      }
    }

    peptides <- peptides[,c(info_col,
                            abundance_columns)]
    peptides <- peptides[order(peptides$`Master Protein Accessions`),]

    colnames(peptides)[(length(info_col)+1):ncol(peptides)] <- paste0(peptides_temperature, "_", treatment)

    if(length(grep("_Mix|_Empty", colnames(peptides)))){
      message(paste(paste(sub(".*_", "", grep("_Mix|_Empty", colnames(peptides), value = TRUE)),
                          collapse = ", "
                          ),
                    "channel(s) has(ve) been removed."
                    )
              )
      peptides <- peptides[, -grep("_Mix|_Empty", colnames(peptides))]
    }

    peptides_dataset[[peptides_temperature]] <- peptides
  }
  message("Concatening all peptides files")
  peptides_dataset <- plyr::join_all(peptides_dataset,
                                     by = info_col,
                                     type = "full")

  descr <- grep("Master Protein Descriptions", colnames(peptides_dataset))
  if(length(descr)){
    colnames(peptides_dataset)[descr] <- "description"
  }

  # Take only the same proteins and adding protein description
  if(!is.null(proteins)){
    peptides_dataset[["description"]] <- NULL
    if(is.character(proteins)){
      if(file.exists(proteins)){
        proteins <- readr::read_tsv(proteins)
        if(all(c("id", "description") %in% colnames(proteins))){
          proteins <- proteins[,c("id", "description")]
          colnames(proteins)[1] <- "Master Protein Accessions"
          peptides_dataset <- left_join(proteins, peptides_dataset, by = "Master Protein Accessions")
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
    else if(inherits(proteins, "data.frame")){
      if(all(c("id", "description") %in% colnames(proteins))){
        proteins <- proteins[,c("id", "description")]
        colnames(proteins)[1] <- "Master Protein Accessions"
        peptides_dataset <- left_join(proteins, peptides_dataset, by = "Master Protein Accessions")
      }
      else{
        message("Error: Your proteins dataset should contains the columns 'id' and 'description'")
        return()
      }
    }
    else{
      message("Error: Your proteins data should either be NULL, a file or a data.frame that contains the columns 'id' and 'description'")
      return()
    }
  }
  else if(!("description" %in% colnames(peptides_dataset))){
    cn <- colnames(peptides_dataset)
    peptides_dataset$description <- NA
    peptides_dataset <- peptides_dataset[,c(cn[1], "description", cn[-1])]
  }


  if(!is.null(modification_torm)){
    if(all(nchar(modification_torm))){
      modif_torm <- paste0("\\d{1}x", modification_torm, collapse = "|")
      modif_torm <- grep(modif_torm, peptides_dataset$Modifications)
      if(length(modif_torm)){
        if(length(modif_torm) != nrow(peptides_dataset)){
          message(paste("Removed", length(modif_torm), "modified peptides"))
          peptides_dataset <- peptides_dataset[-modif_torm,]
        }
        else{
          message("Warninig: Check the modifications you want to remove, it fitted with all the peptides from your dataset. Hence, no peptide were removed.")
        }
      }
      else{
        message("Warninig: Check the modifications you want to remove, it fitted with no peptides from your dataset. Hence, no peptide were removed.")
      }
    }
  }

  message("Saving files")
  # remove peptides with only NAs
  n <- nrow(peptides_dataset)
  peptides_dataset <- peptides_dataset[which(apply(peptides_dataset[,6:ncol(peptides_dataset)], 1,
                                             function(x) sum(is.na(x)) < ncol(peptides_dataset) - 5)),]
  n <- n - nrow(peptides_dataset)
  message(paste(n, "peptides without quantitative information were removed"))


 readr::write_tsv(peptides_dataset,
                  file = paste0(format(Sys.time(), "%y%m%d_%H%M_"),
                                "peptides_", dataset_name, ".txt")
                  )
 return(peptides_dataset)
}
