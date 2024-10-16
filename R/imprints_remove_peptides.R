#' imprints_remove_peptides
#'
#' Function to filter some peptides according there sequence. For example after using \code{imprints_sequence_peptides},
#' if you want to remove the peptides which corresponds to the cleaved sites.
#'
#' @param data A peptides dataset, typically after \code{imprints_sequence_peptides}.
#'             It needs to contains the columns 'Master Protein Accessions' and 'Positions in Master Proteins'.
#' @param proteins The proteins from which you want to remove/keep the peptides.
#'                 If NULL and you put one sequence, it will remove/keep all these peptides for all proteins from peptides data.
#'                 If NULL and you put several sequence, the number of sequence you put needs to match the number
#'                 of proteins you have in your peptides data.
#' @param sequence The peptide position you want to remove/keep. If you put one sequence, it will remove/keep it for all proteins
#'                 from your peptides data; otherwise it needs to match the number of proteins you put or have in your dataset.
#'                 The format needs to be a number followed by a dash and another number, like this : '208-221'.
#' @param margin Numeric to tell the margin number of amino acid added to the sequences selected. This avoid
#'   peptide selection ambigutiy when a peptide has only a small number of amino acid more than the sequence selected.
#'   Default is set to 2.
#' @param mode Character to specify if you want to remove or only the peptides selected; either 'remove' or 'keep'.
#'             Default is 'remove'.
#'
#' @return Your filtered dataset
#'
#' @export
#'
#'

imprints_remove_peptides <- function(data, proteins = NULL, sequence,
                                     margin = 2, mode = c("remove", "keep")){
  mode <- match.arg(mode)

  if(length(sequence)){
    if(!all(grepl("^\\d{1,}(-|~)\\d{1,}", sequence))){
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
    proteins <- unique(data$Master.Protein.Accessions)
  }
  else{
    if(!all(proteins %in% data$Master.Protein.Accessions)){
      message("Error: Check the proteins you selected, none are in your data.")
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

  sequence_tokeep <- c()
  for(p in 1:length(proteins)){
    protein_idx <- which(data$Master.Protein.Accessions == proteins[p])
    protein_sequence <- unlist(stringr::str_extract_all(sub("; .*", "",
                                                            data$Positions.in.Master.Proteins[protein_idx]), # if protein group just need first one
                                                        "(?<=\\[).+?(?=\\])")
                               )
    if(length(sequence) > 1){
      sequence_tofilter <- sequence[p]
    }
    else{
      sequence_tofilter <- sequence
    }
    sequence_tofilter <- sub("; .*", "", sequence_tofilter) # if protein group just need first one
    sequence_tofilter <- as.numeric(strsplit(sequence_tofilter, "-|~")[[1]])

    sequence_tofilter <- lapply(strsplit(protein_sequence, "-|~"),
                                function(y){y <- as.numeric(y)
                                            y <- y[1] >= sequence_tofilter[1] - margin & y[2] <= sequence_tofilter[2] + margin;
                                            y
                                            }
                                )
    sequence_tofilter <- which(unlist(sequence_tofilter))

    if(length(sequence_tofilter)){
      sequence_tofilter <- paste0("\\[", protein_sequence[sequence_tofilter], "\\]")
      pr_tofilter <- data$Positions.in.Master.Proteins[protein_idx]
      pr_tofilter_idx <- protein_idx

      sequence_tofilter_idx <- pr_tofilter_idx[grep(sequence_tofilter, sub("; .*", "", pr_tofilter))]
      sequence_tofilter <- pr_tofilter[grep(sequence_tofilter, sub("; .*", "", pr_tofilter))]


      message(paste(ifelse(mode == "remove", "Removing", "Keeping"),
                    paste(sequence_tofilter, collapse = ", ")))

      if(mode == "remove"){
        data <- data[-sequence_tofilter_idx,]
      }
      else if(mode == "keep"){
        sequence_tokeep <- c(sequence_tokeep, sequence_tofilter_idx)
      }
    }
    else{
      message(paste("No sequence has been", ifelse(mode == "remove", "removed", "kept"),
                    "for protein", proteins[p]))
    }
  }

  if(mode == "keep" & length(sequence_tokeep)){
    data <- data[sequence_tokeep,]
  }

  return(data)
}

