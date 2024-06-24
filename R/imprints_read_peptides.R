#' imprints_read_peptides
#'
#' Function to read your peptides files and concatenate them in one data frame.
#'
#' @param peptides_files The path to the PeptidesGroups file from PD. Can be one or more files.
#' @param treatment A character vector that contains the treatment applied for each channel like 'B1_Vehicle' for example.
#' @param temperatures A character vector that contains the temperatures used for each peptides file like '37C' for example.
#' @param prefixcontaminant Character corresponding to the prefix used to identify contaminants.
#' @param averagecount 	Whether to take the median of the abundance count numbers across the measured temperature
#'                      range and then use this value for filtering. Default set to TRUE.
#'                      Otherwise, filter the  peptides according to the associated count numbers at each temperature.
#' @param countthreshold The minimal threshold number of associated abundance count of peptides, default is 2.
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
                                   prefixcontaminant = "", averagecount = TRUE, countthreshold = 2,
                                   proteins = NULL, modification_torm = NULL,
                                   dataset_name = "imprints"){
  n_peptides_files <- length(peptides_files)
  if(length(temperatures) != n_peptides_files){
    message("Error: The number of temperatures doesn't match the number of peptides files !")
    return()
  }
  if(is.character(temperatures)){
    if(!all(grepl("^\\d{1,2}C$", temperatures))){
      message("Error: Your temperatures should be in the format '1 or 2 digits' followed by capital 'C', like '37C' for example.")
      return()
    }
  }
  else{
    message("Error: Your temperatures should be a character vector !")
    return()
  }

  if(countthreshold < 0){
    messsage("Parameter countthreshold can only be positive; value set to 0, hence no filtering will be applied.")
  }

  peptides_files <- as.list(peptides_files)
  names(peptides_files) <- temperatures

  message("Reading peptides files")
  peptides_dataset <- list()
  for(p in 1:n_peptides_files){
    peptides_temperature <- names(peptides_files)[p]
    peptides <- peptides_files[[p]]

    peptides <- readr::read_tsv(peptides, show_col_types = FALSE)
    abundance_columns <- grep("^Abundance: |^Abundances: ", colnames(peptides), value = TRUE)
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

    # extracting peptide abundance count
    countNum <- grep(pattern = paste0("^Abundances Count: [A-z0-9,. -]+: 126[A-z0-9,. -]+$"),
                     colnames(peptides), value = TRUE)
    if(length(countNum) == 1){
      peptides <- peptides[,c(info_col, abundance_columns, countNum)]
      colnames(peptides)[ncol(peptides)] <- paste0("countNum.", peptides_temperature)
    }
    else{
      peptides <- peptides[,c(info_col, abundance_columns)]
      message("Warning: There is no peptide abundance count information in the original data.\n
              Set a default number for the peptide abundance count to 2.")
      peptides[[paste0("countNum.", peptides_temperature)]] <- 2
    }

    # removing contaminants
    if (nchar(prefixcontaminant)) {
      ncont <- nrow(peptides)
      pattern <- grep(paste0("(^|(;|; ))", prefixcontaminant),
                      peptides$`Master Protein Accessions`)
      if (length(pattern) > 0) {
        peptides <- peptides[-pattern, ]
      }
      if (nrow(peptides) == 0) {
        stop("After removing Contaminant proteins, the dataset was empty.")
      }
      pattern <- grep("(B|b)os taurus", peptides$`Master Protein Descriptions`)
      if (length(pattern) > 0) {
        peptides <- peptides[-pattern, ]
      }
      if (nrow(peptides) == 0) {
        stop("After removing Bos taurus proteins, the dataset was empty.")
      }
      message(paste0(ncont - nrow(peptides), " contaminants were removed from the data at temperature ", peptides_temperature))
    }

    # removing mix and/or empty channels if any
    colnames(peptides)[(length(info_col)+1):(ncol(peptides) - 1)] <- paste0(peptides_temperature, "_", treatment)
    if(length(grep("_Mix|_Empty", colnames(peptides)))){
      message(paste(paste(sub(".*_", "",
                              grep("_Mix|_Empty", colnames(peptides), value = TRUE)
                              ),
                          collapse = ", "
                          ),
                    "channel(s) has(ve) been removed."
                    )
              )
      peptides <- peptides[, -grep("_Mix|_Empty", colnames(peptides))]
    }

    # filtering on abundance count if not on median abundance count but by temperature
    if(!averagecount){
      n_abc <- nrow(peptides)
      peptides <- peptides[which(peptides[[paste0("countNum.", peptides_temperature)]] >= countthreshold),]
      message(paste0(n_abc - nrow(peptides), " peptides didn't pass the cutoff of abundance count number\n",
                     countthreshold, " for the temperature ", peptides_temperature))
    }

    peptides <- peptides[order(peptides$`Master Protein Accessions`),]
    peptides_dataset[[peptides_temperature]] <- peptides
  }

  message("Concatening all peptides files")
  peptides_dataset <- plyr::join_all(peptides_dataset, by = info_col,
                                     type = "full")

  # computing median from abundance count
  peptides_dataset$countNum <- apply(peptides_dataset[,grep("^countNum\\.", colnames(peptides_dataset))],
                                     1, median, na.rm = TRUE)
  peptides_dataset <- peptides_dataset[,-grep("^countNum\\.", colnames(peptides_dataset))]

  # filtering on abundance count if filtering on median abundance count
  if(averagecount){
    n_abc <- nrow(peptides_dataset)
    peptides_dataset <- peptides_dataset[which(peptides_dataset$countNum >= countthreshold),]
    message(paste0(n_abc - nrow(peptides_dataset), " peptides didn't pass the cutoff of median abundance count number ", countthreshold))
  }

  # renaming column description if there
  descr <- grep("Master Protein Descriptions", colnames(peptides_dataset))
  if(length(descr)){
    colnames(peptides_dataset)[descr] <- "description"
  }

  # Take only the same proteins and adding protein description if added proteins argument
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


  # filtering modification
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
  peptides_dataset <- peptides_dataset[which(apply(peptides_dataset[,6:(ncol(peptides_dataset)-1)], 1,
                                             function(x) sum(is.na(x)) < ncol(peptides_dataset) - 6)),]
  n <- n - nrow(peptides_dataset)
  message(paste(n, "peptides without quantitative information were removed"))
  colnames(peptides_dataset)[1:5] <- gsub(" ", ".", colnames(peptides_dataset)[1:5])

 readr::write_tsv(peptides_dataset,
                  file = paste0(format(Sys.time(), "%y%m%d_%H%M_"),
                                "peptides_", dataset_name, ".txt")
                  )
 return(peptides_dataset)
}
