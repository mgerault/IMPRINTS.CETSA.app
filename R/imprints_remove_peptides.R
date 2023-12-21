#' imprints_remove_peptides
#'
#' Function to filter some peptides according there sequence. For example after using \code{imprints_sequence_peptides},
#' if you want to remove the peptides which corresponds to the cleaved sites.
#'
#' @param data A peptides dataset, typically after \code{imprints_sequence_peptides}.
#'             It needs to contains the columns 'Master Protein Accessions' and 'Positions in Master Proteins'.
#' @param proteins The proteins from which you want to remove the peptides.
#'                 If NULL and you put one sequence, it will remove all these peptides for all proteins from peptides data.
#'                 If NULL and you put several sequence, the number of sequence you put needs to match the number
#'                 of proteins you have in your peptides data.
#' @param sequence The peptide position you want to remove. If you put one sequence, it will remove it for all proteins
#'                 from your peptides data; otherwise it needs to match the number of proteins you put or have in your dataset.
#'                 The format needs to be a number followed by a dash and another number, like this : '208-221'.
#'
#' @return Your filtered dataset
#'
#' @export
#'
#'

imprints_remove_peptides <- function(data, proteins = NULL, sequence){
  if(length(sequence)){
    if(!all(stringr::str_detect(sequence, "^\\d{1,}-\\d{1,}"))){
      message("Error: 'sequence' needs to be a character vector and needs to be a number followed
              by a dash and another number like this '208-221'.")
      return()
    }
  }
  else{
    message("Error: You need to provide at least one sequence")
    return()
  }

  if(is.null(proteins)){
    proteins <- unique(stringr::str_remove_all(data$`Master Protein Accessions`, "\\s.*"))
  }
  else{
    if(!all(proteins %in% stringr::str_remove_all(data$`Master Protein Accessions`, "\\s.*"))){
      message("Error: Check the proteins you select, some are not in your data.")
      return()
    }
  }

  if(length(sequence) > 1){
    if(length(sequence) != length(proteins)){
      message("Error: The number of sequence doesn't match the number of proteins !")
      return()
    }
  }
  else{
    message(paste("The sequence", sequence, "will be removed from all proteins selected."))
  }

  for(p in 1:length(proteins)){
    protein_idx <- which(stringr::str_remove_all(data$`Master Protein Accessions`, "\\s.*") == proteins[p])
    protein_sequence <- unlist(stringr::str_extract_all(data$`Positions in Master Proteins`[protein_idx],
                                                        "(?<=\\[).+?(?=\\])")
                               )
    if(length(sequence) > 1){
      sequence_toremove <- sequence[p]
    }
    else{
      sequence_toremove <- sequence
    }
    sequence_toremove <- as.numeric(stringr::str_split(sequence_toremove, "-")[[1]])

    sequence_toremove <- lapply(stringr::str_split(protein_sequence, "-|~"),
                                function(y){y <- as.numeric(y)
                                            y <- y[1] >= sequence_toremove[1] & y[2] <= sequence_toremove[2];
                                            y
                                            }
                                )
    sequence_toremove <- which(unlist(sequence_toremove))

    if(length(sequence_toremove)){
      sequence_toremove <- protein_sequence[sequence_toremove]
      sequence_toremove <- paste0(proteins[p], " [", sequence_toremove, "]")

      message(paste("Removing", paste(sequence_toremove, collapse = ", ")))

      sequence_toremove <- which(!is.na(match(data$`Positions in Master Proteins`, sequence_toremove)))
      data <- data[-sequence_toremove,]
    }
    else{
      message(paste("No sequence has been removed for protein", proteins[p]))
    }
  }

  return(data)
}

