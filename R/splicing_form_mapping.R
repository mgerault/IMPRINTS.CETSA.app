#' imprints_isoform_peptides
#'
#' Function to check if the hits found by \code{\link{imprints_cleaved_peptides}} could be due to alternate
#' splicing forms being more expressed and not due to protein cleavage.
#'
#' @details
#' When a hit is returned by \code{\link{imprints_cleaved_peptides}}, it means that this protein has a peptide position
#' where the IMPRINTS profiles of the two obtained parts are significantly different. This difference can be caused
#' by protein modification and mainly proteolysis; but if a splicing form of protein is more expressed than its canonical
#' form, a significant difference in the profiles can also occur. The aim here is to refilter the hit list and give
#' the possible splicing forms which could be more expressed based on the output of \code{\link{imprints_cleaved_peptides}}
#' and the sequence alignments of the isoform sequence and its corresponding canonical form.
#'
#' @param data The normalized peptides data set, i.e. the outpout from \code{imprints_normalize_peptides}.
#' @param data_cleaved The cleavage hits data set, i.e. the outpout from \code{imprints_cleaved_peptides}.
#' @param control The control treatment from your dataset.
#' @param fasta The path to the FASTA file you used for the search
#' @param minimum_align Numeric to tell the minimum aligned sequence length. Default to 5.
#'   It is advised to put the same as the minimum peptide sequence length you selected in
#'   the proteomics search software you used.
#' @param save_xlsx Logical to tell if you want to save the categorized hits in an xlsx file.
#'   Default to TRUE.
#' @param xlsxname The name of your saved file.
#'
#'
#' @return A dataframe containing the potential splicing forms which could be more
#'         expressed instead of the protein being cleaved
#'
#' @export
#'
#' @seealso \code{\link{imprints_cleaved_peptides}}
#'

imprints_isoform_peptides <- function(data, data_cleaved, control, fasta,
                                      minimum_align = 5, save_xlsx = TRUE, xlsxname = "RESP_isoform_mapping"){
  if(!grepl("\\.fasta$", fasta)){
    message("the path to your FASTA file does not lead to a FASTA file !")
    return()
  }

  if(!(control %in% get_treat_level(data))){
    message(paste(control, "is not in your data !"))
    return()
  }

  if(control %in% unique(data_cleaved$treatment)){
    message(paste("Error:", control, "shouldn't be in your data_cleaved !"))
    return()
  }

  # fetch isoform from each hits from data_cleaved
  message("Fetching isoforms...")
  isoforms <- fetch_isoforms(unique(data_cleaved$id), fasta)

  if(nrow(isoforms) == 0){
    message("No isoform could have been found from your data")
    return()
  }

  # align isoform sequence with canonical
  message("Aligning isoform sequences...")
  isoforms$canonical_posalign <- apply(isoforms[,c("sequence", "canonical_sequence")], 1,
                                       function(x) compare_sequences(x[1], x[2])
                                       )

  # filter results
  isoforms$canonical_posalign <- sapply(isoforms$canonical_posalign,
                                        function(x){
                                          x <- strsplit(x, "; ")[[1]]
                                          x <- x[grep("-", x)] # removing similarity of only one AA
                                          if(length(x)){
                                            # removing sequence similarity of less than minimum_align AA
                                            xl <- sapply(strsplit(x, "-"),
                                                         function(z) diff(as.numeric(z))
                                            )
                                            x <- x[which(xl >= minimum_align)]
                                            if(length(x)){
                                              x <- paste(x, collapse = "; ")
                                            }
                                            else{
                                              x <- NA
                                            }
                                          }
                                          else{
                                            x <- NA
                                          };
                                          x
                                        }, USE.NAMES = FALSE)


  ### Map isoform alignment to RESP data
  message("Mapping alignments to RESP results and filtering candidates...")
  if(all(c("category", "details") %in% colnames(data_cleaved))){
    data_isoform_ambiguity <- data_cleaved[,c("id", "Gene", "description", "treatment", "cleaved_site",
                                              "RESP_score", "category", "details")]
  }
  else{
    data_isoform_ambiguity <- data_cleaved[,c("id", "Gene", "description", "treatment", "cleaved_site",
                                              "RESP_score")]
  }

  colnames(data_isoform_ambiguity)[1] <- "accession"
  data_isoform_ambiguity <- right_join(data_isoform_ambiguity, isoforms[,c("accession", "canonical_posalign",
                                                                           "isoforms",
                                                                           "length_isoform", "length_canonical")],
                                       by = "accession", relationship = "many-to-many")


  ## filtering
  # if no significant isoform alignment found, remove it (NAs)
  filtered_idx <- which(is.na(data_isoform_ambiguity$canonical_posalign))
  if(length(filtered_idx))
    data_isoform_ambiguity <- data_isoform_ambiguity[-filtered_idx,]

  if(nrow(data_isoform_ambiguity) == 0){
    message("No RESP hits could be due to alternative splicing forms being more expressed")
    return()
  }

  # if whole canonical sequence is in the isoform (simple addition at end of begining of sequence)
  # --> no way to detect it with RESP
  filtered_idx <- which(data_isoform_ambiguity$canonical_posalign ==
                          paste0("1-", data_isoform_ambiguity$length_canonical)
                        )
  if(length(filtered_idx))
    data_isoform_ambiguity <- data_isoform_ambiguity[-filtered_idx,]

  if(nrow(data_isoform_ambiguity) == 0){
    message("No RESP hits could be due to alternative splicing forms being more expressed")
    return()
  }

  # getting length of aligned sequence
  data_isoform_ambiguity$length_align <- sapply(strsplit(data_isoform_ambiguity$canonical_posalign, "; "),
                                                function(x){
                                                  x <- sapply(strsplit(x, "-"),
                                                              function(y){
                                                                y <- as.numeric(y)
                                                                y <- y[order(y, decreasing = TRUE)]
                                                                y[1] - y[2] + 1
                                                              })

                                                  sum(x)
                                                })

  # assign gap between common sequence
  # if length canonical 1-900 and alignment is 1-15; 80-900 --> gap of 65
  # if length canonical was 1-1000 --> gap of 65; 100
  data_isoform_ambiguity$gap_align <- mapply(function(align, canl){
    ord <- sapply(strsplit(align, "-"), function(y) sum(as.numeric(y)))
    align <- align[order(ord)]
    align <- unlist(sapply(strsplit(align, "-"), as.numeric,
                           simplify = FALSE))
    gap <- diff(align) - 1
    gap <- gap[2*c(1:(length(gap)-ceiling(length(gap)/2)))]
    if(align[1] != 1){
      gap <- c(align[1] - 1, gap[which(!is.na(gap))])
    }
    if(align[length(align)] != canl){
      gap <- c(gap[which(!is.na(gap))], canl - align[length(align)] - 1)
    };
    paste(gap, collapse = "; ")
  }, strsplit(data_isoform_ambiguity$canonical_posalign, "; "),
  data_isoform_ambiguity$length_canonical)

  # if the length of the aligned sequence is the same as the canonical and that the gap is only a 0 --> insertion
  # then, if the isoform is more present/stabilized than canonical, at worst one or two peptides around the insertion
  # will not be measured (and this no matter the size of the insetion) but the rest will shift the same way
  # --> impossible to see it with RESP analysis
  filtered_idx <- which(data_isoform_ambiguity$length_align == data_isoform_ambiguity$length_canonical &
                          data_isoform_ambiguity$gap_align == "0")
  if(length(filtered_idx))
    data_isoform_ambiguity <- data_isoform_ambiguity[-filtered_idx,]

  if(nrow(data_isoform_ambiguity) == 0){
    message("No RESP hits could be due to alternative splicing forms being more expressed")
    return()
  }

  # if the length of the aligned sequence is the same as the isoform --> deletions
  # if gaps too small in the middle, impossible to identify it with RESP
  # if gap around peptide length (above minimum_align) and at the end or the beginning, might come as hits in RESP (C-term or N-term)
  # not shifting or shifting negatively
  filtered_idx <- which(data_isoform_ambiguity$length_align == data_isoform_ambiguity$length_isoform &
                          sapply(strsplit(data_isoform_ambiguity$gap_align, "; "),
                                 function(x) max(sapply(x, as.numeric))) < minimum_align
                        )
  if(length(filtered_idx))
    data_isoform_ambiguity <- data_isoform_ambiguity[-filtered_idx,]

  if(nrow(data_isoform_ambiguity) == 0){
    message("No RESP hits could be due to alternative splicing forms being more expressed")
    return()
  }

  filtered_idx <- which(data_isoform_ambiguity$length_align == data_isoform_ambiguity$length_isoform &
                          grepl("; ", data_isoform_ambiguity$canonical_posalign) &
                          sapply(strsplit(data_isoform_ambiguity$gap_align, "; "),
                                 function(x) max(sapply(x, as.numeric))) < 15 &
                          !grepl("negatively", data_isoform_ambiguity$details)
                        )
  if(length(filtered_idx))
    data_isoform_ambiguity <- data_isoform_ambiguity[-filtered_idx,]

  if(nrow(data_isoform_ambiguity) == 0){
    message("No RESP hits could be due to alternative splicing forms being more expressed")
    return()
  }

  # if maximum gaps below minimum_align, impossible to impact on RESP analysis
  filtered_idx <- which(sapply(strsplit(data_isoform_ambiguity$gap_align, "; "),
                               function(x) max(sapply(x, as.numeric))) < minimum_align)
  if(length(filtered_idx))
    data_isoform_ambiguity <- data_isoform_ambiguity[-filtered_idx,]

  if(nrow(data_isoform_ambiguity) == 0){
    message("No RESP hits could be due to alternative splicing forms being more expressed")
    return()
  }

  # if sequence align contains start and end of canonical and more than 2 gaps, almost impossible to identify it with RESP
  filtered_idx <- which(grepl("(^|; )1-", data_isoform_ambiguity$canonical_posalign) &
                          mapply(function(p, canl) grepl(paste0("-", canl), p),
                                 data_isoform_ambiguity$canonical_posalign,
                                 data_isoform_ambiguity$length_canonical,
                                 USE.NAMES = FALSE) &
                          lengths(regmatches(data_isoform_ambiguity$canonical_posalign,
                                             gregexpr("; ", data_isoform_ambiguity$canonical_posalign))) >= 2 # more than 2 gaps
                        )
  if(length(filtered_idx))
    data_isoform_ambiguity <- data_isoform_ambiguity[-filtered_idx,]

  if(nrow(data_isoform_ambiguity) == 0){
    message("No RESP hits could be due to alternative splicing forms being more expressed")
    return()
  }

  # if align sequence has no gap and contains start or end  (like 1-792 or 280-950 for canonical of length 950)
  # and that the cleavage site is close from the align sequence, could be isoform
  # (--> will then need to shifting profiles to decide)
  filtered_idx <- which(lengths(regmatches(data_isoform_ambiguity$canonical_posalign,  # no gaps
                                           gregexpr("; ", data_isoform_ambiguity$canonical_posalign))) == 0 &
                          (grepl("(^|; )1-", data_isoform_ambiguity$canonical_posalign) |
                             mapply(function(p, canl) grepl(paste0("-", canl), p),
                                    data_isoform_ambiguity$canonical_posalign,
                                    data_isoform_ambiguity$length_canonical,
                                    USE.NAMES = FALSE)) & # contains start or end
                          !mapply(function(align, cl){
                            cl <- strsplit(cl, "; ")[[1]][1] # if protein group, only first protein is considered for isoform
                            cl <- as.numeric(strsplit(cl, "-|~")[[1]])
                            # add margin of 50 AA
                            cl[1] <- cl[1] - 50
                            cl[2] <- cl[2] + 50

                            align <- as.numeric(strsplit(align, "-")[[1]])
                            align <- ifelse(align[1] == 1, align[2], align[1]);
                            res <- cl[1] <= align & align <= cl[2]
                            ifelse(is.na(res), FALSE, res)
                          },
                          data_isoform_ambiguity$canonical_posalign,
                          data_isoform_ambiguity$cleaved_site,
                          USE.NAMES = FALSE)
                        )
  if(length(filtered_idx))
    data_isoform_ambiguity <- data_isoform_ambiguity[-filtered_idx,]

  if(nrow(data_isoform_ambiguity) == 0){
    message("No RESP hits could be due to alternative splicing forms being more expressed")
    return()
  }

  # biggest aligned part should shift the same way so cannot contain the cleavage site if isoform
  filtered_idx <- which(mapply(function(align, cl){
                                cl <- strsplit(cl, "; ")[[1]][1] # if protein group, only first protein is considered for isoform
                                cl <- as.numeric(strsplit(cl, "-|~")[[1]])
                                # add margin of 50 AA --> even if you take greater cleavge site still in largest aligned sequence ?
                                cl[1] <- cl[1] - 50
                                cl[2] <- cl[2] + 50

                                align <- strsplit(align, "; ")[[1]]
                                big_align <- sapply(align,
                                                    function(x){
                                                      diff(as.numeric(strsplit(x, "-")[[1]]))
                                                    }, USE.NAMES = FALSE)
                                big_align <- align[which.max(big_align)]
                                big_align <- as.numeric(strsplit(big_align, "-")[[1]])

                                big_align[1] <= cl[1] & cl[2] <= big_align[2]
                              },
                              data_isoform_ambiguity$canonical_posalign,
                              data_isoform_ambiguity$cleaved_site,
                              USE.NAMES = FALSE)
                        )
  if(length(filtered_idx))
    data_isoform_ambiguity <- data_isoform_ambiguity[-filtered_idx,]

  if(nrow(data_isoform_ambiguity) == 0){
    message("No RESP hits could be due to alternative splicing forms being more expressed")
    return()
  }


  # if no cleavage event and just an isoform is more abundant than canonical
  # --> either only the aligned part shift positively and/or only the non-aligned shift negatively
  # so if aligned belong to N-term (compare to clavage site) you expect in any case FCn > FCc --> RESP score < 0
  # conversely, if aligned belong to C-term (compare to clavage site) you expect in any case FCc > FCn --> RESP score > 0
  filtered_idx <- which(mapply(function(align, cl, resp_score){
                                cl <- strsplit(cl, "; ")[[1]][1] # if protein group, only first protein is considered for isoform
                                cl <- as.numeric(strsplit(cl, "-|~")[[1]])

                                align_before <- all(sapply(strsplit(align, "; ")[[1]],
                                                           function(x){
                                                             x <- as.numeric(strsplit(x, "-")[[1]])
                                                             x <= cl[1]
                                                           }))
                                if(all(align_before)){ # FCn > FCc
                                  return(resp_score < 0) # if FALSE, can't be isoform
                                }
                                else{
                                  align_after <- all(sapply(strsplit(align, "; ")[[1]],
                                                            function(x){
                                                              x <- as.numeric(strsplit(x, "-")[[1]])
                                                              x >= cl[2]
                                                            }))
                                  if(all(align_after)){ # FCc > FCn
                                    return(resp_score > 0) # if FALSE, can't be isoform
                                  }
                                  else{ # ambiguous --> need to keep
                                    return(TRUE)
                                  }
                                }
                              },
                              data_isoform_ambiguity$canonical_posalign,
                              data_isoform_ambiguity$cleaved_site,
                              data_isoform_ambiguity$RESP_score,
                              USE.NAMES = FALSE)
                        )
  if(length(filtered_idx))
    data_isoform_ambiguity <- data_isoform_ambiguity[filtered_idx,]

  if(nrow(data_isoform_ambiguity) == 0){
    message("No RESP hits could be due to alternative splicing forms being more expressed")
    return()
  }

  # if small insertion/deletion (~50 AA) in the canonical sequence, I don't see how it could be
  # a hit in RESP analysis
  filtered_idx <- which(sapply(strsplit(data_isoform_ambiguity$gap_align, "; "),
                               function(x) max(sapply(x, as.numeric))) < 50 &
                          lengths(regmatches(data_isoform_ambiguity$canonical_posalign,
                                             gregexpr("; ", data_isoform_ambiguity$canonical_posalign))) == 1
                        )
  if(length(filtered_idx))
    data_isoform_ambiguity <- data_isoform_ambiguity[-filtered_idx,]

  if(nrow(data_isoform_ambiguity) == 0){
    message("No RESP hits could be due to alternative splicing forms being more expressed")
    return()
  }

  # if majority aligned is toward N-term, you could expect RESP score < 0
  filtered_idx <- which(grepl("(^|; )1-", data_isoform_ambiguity$canonical_posalign) &
                          !mapply(function(p, canl) grepl(paste0("-", canl), p),
                                  data_isoform_ambiguity$canonical_posalign,
                                  data_isoform_ambiguity$length_canonical,
                                  USE.NAMES = FALSE) &
                          data_isoform_ambiguity$RESP_score > 0
                        )
  if(length(filtered_idx))
    data_isoform_ambiguity <- data_isoform_ambiguity[-filtered_idx,]

  if(nrow(data_isoform_ambiguity) == 0){
    message("No RESP hits could be due to alternative splicing forms being more expressed")
    return()
  }

  # conversely, if majority toward C-term, you could expect RESP score > 0
  filtered_idx <- which(!grepl("(^|; )1-", data_isoform_ambiguity$canonical_posalign) &
                          mapply(function(p, canl) grepl(paste0("-", canl), p),
                                 data_isoform_ambiguity$canonical_posalign,
                                 data_isoform_ambiguity$length_canonical,
                                 USE.NAMES = FALSE) &
                          data_isoform_ambiguity$RESP_score < 0
                        )
  if(length(filtered_idx))
    data_isoform_ambiguity <- data_isoform_ambiguity[-filtered_idx,]

  if(nrow(data_isoform_ambiguity) == 0){
    message("No RESP hits could be due to alternative splicing forms being more expressed")
    return()
  }

  # continue to follow this idea but check where's majority of aligned sequence compare to cleavage site
  filtered_idx <- which(mapply(function(align, cl, alil, resp_score){
                                cl <- strsplit(cl, "; ")[[1]][1] # if protein group, only first protein is considered for isoform
                                cl <- as.numeric(strsplit(cl, "-|~")[[1]])

                                nterm <- 1:cl[1]

                                nterm_align <- sum(sapply(strsplit(align, "; ")[[1]],
                                                          function(x){
                                                            x <- as.numeric(strsplit(x, "-")[[1]])
                                                            x <- x[1]:x[2]
                                                            sum(x %in% nterm)
                                                          }, USE.NAMES = FALSE))

                                nterm_align <- nterm_align/alil
                                if(nterm_align > 0.5){
                                  return(resp_score < 0)
                                }
                                else{
                                  return(resp_score > 0)
                                }
                              },
                              data_isoform_ambiguity$canonical_posalign,
                              data_isoform_ambiguity$cleaved_site,
                              data_isoform_ambiguity$length_align,
                              data_isoform_ambiguity$RESP_score,
                              USE.NAMES = FALSE))
  if(length(filtered_idx))
    data_isoform_ambiguity <- data_isoform_ambiguity[filtered_idx,]

  if(nrow(data_isoform_ambiguity) == 0){
    message("No RESP hits could be due to alternative splicing forms being more expressed")
    return()
  }

  ## data_isoform_ambiguity now contains only potential candidates
  ## only way to check now is to check peptides profiles
  # filtering protein of interest
  message("Comparing potential isoforms candiates with FC data...")
  data <- data[which(!is.na(match(data$Master.Protein.Accessions,
                                         data_isoform_ambiguity$accession
                                         ))),]

  # computing fold-changes
  directory_toremove <- list.files()
  data <- imprints_sequence_peptides(data, control = control) # will save txt file
  directory_toremove <- list.files()[!(list.files() %in% directory_toremove)]
  # remove this file
  directory_toremove <- grep("(?=^\\d{6}_\\d{4}_.*)(?=.*\\.txt$)", directory_toremove, value = TRUE, perl = TRUE)
  if(length(directory_toremove) == 1){
    unlink(directory_toremove, recursive = TRUE)
  }
                        
  data <- data[,-grep(paste0("_", control, "$"), colnames(data))]

  data <- data %>%
    select(-Annotated.Sequence, -countNum) %>%
    mutate(Gene = sapply(description, IMPRINTS.CETSA.app:::getGeneName, USE.NAMES = FALSE),
           description = sapply(description, IMPRINTS.CETSA.app:::getProteinName, USE.NAMES = FALSE)
           ) %>%
    tidyr::gather("key", "value", -Master.Protein.Accessions, -Gene, -description,
                  -Positions.in.Master.Proteins, -Modifications) %>%
    tidyr::separate(key, into = c("temperature", "biorep", "treatment"), sep = "_") %>%
    group_by(Master.Protein.Accessions, Gene, description,
             Positions.in.Master.Proteins, Modifications,
             temperature, treatment) %>%
    summarise(value = mean(value, na.rm = TRUE)) %>%
    ungroup() %>% group_by(Master.Protein.Accessions, Gene, description,
                           Positions.in.Master.Proteins, Modifications,
                           treatment) %>%
    summarise(value = value[which.max(abs(value))])

  colnames(data)[c(1,4)] <- c("accession", "pep_position")
  data$pep_position <- gsub(".* \\[|\\]", "", sub("; .*", "", data$pep_position))

  ## checking sign of aligned and non aligned sequence
  data_isoform_ambiguity <- left_join(data_isoform_ambiguity, data,
                                      by = c("accession", "Gene", "description", "treatment"),
                                      relationship = "many-to-many")


  message("Refiltering...")
  data_isoform_ambiguity <- data_isoform_ambiguity %>%
    group_by(accession, Gene, description, treatment, isoforms) %>%
    group_modify(~ {
      align <- .x$canonical_posalign[1]
      align <- strsplit(align, "; ")[[1]]
      align <- sapply(strsplit(align, "-"), as.numeric, simplify = FALSE)

      pep_pos <- sapply(strsplit(.x$pep_position, "-"), as.numeric, simplify = FALSE)
      pep_in_aligned <- sapply(pep_pos,
                               function(p){
                                 any(sapply(align,
                                            function(a){
                                              a[1] <= p[1] & p[2] <= a[2]
                                            }))
                               })
      pep_only_cano <- which(!pep_in_aligned)
      pep_in_aligned <- which(pep_in_aligned)

      # if only isoform, pep_in_aligned should shift positively and/or pep_only_cano should not shift or shift negatively
      # if only pep_aligned and non canonical measured, then if it was only an isoform more expressed, would come as RESP hit
      # same if it is only canonical
      if(length(pep_in_aligned) & length(pep_only_cano)){
        if(sum(.x$value[pep_in_aligned] > 0)/length(pep_in_aligned) > 0.5){
          if(sum(.x$value[pep_only_cano] <= 0)/length(pep_only_cano) > 0.5){
            isoform <- TRUE
          }
          else if(sum(.x$value[pep_in_aligned] >= 0.15)/length(pep_in_aligned) > 0.5){ # arbitrary cutoffs
            isoform <- TRUE
          }
          else{
            isoform <- FALSE
          }
        }
        else if(sum(.x$value[pep_only_cano] <= -0.15)/length(pep_only_cano) > 0.5){
          isoform <- TRUE
        }
        else{
          isoform <- FALSE
        }
      }
      else{
        isoform <- FALSE
      }

      if(isoform){
        # if potential isoform, meaning that the isoform is indeed more expressed
        # we should be confident that no peptides in the canonical one shift significantly more
        # than the others (the aligned ones); or more non negligently
        # --> comparing with median of aligned
        max_canonical <- max(.x$value[pep_only_cano], na.rm = TRUE)
        median_pos_align <- .x$value[pep_in_aligned][.x$value[pep_in_aligned] >= 0.15]
        if(!length(median_pos_align) & max_canonical >= 0.1){
          isoform <- FALSE
        }
        else{
          median_pos_align <- median(median_pos_align, na.rm = TRUE)
          if(max_canonical >= median_pos_align){
            isoform <- FALSE
          }
        }
      }

      .x$probable_isoform <- isoform

      return(.x)
    }) %>%
    filter(probable_isoform)

  if(nrow(data_isoform_ambiguity) == 0){
    message("No RESP hits could be due to alternative splicing forms being more expressed")
    return()
  }

  data_isoform_ambiguity$probable_isoform <- NULL
  data_isoform_ambiguity <- unique(data_isoform_ambiguity[,1:14])

  if(save_xlsx){
    openxlsx::write.xlsx(data_isoform_ambiguity, paste0(xlsxname, ".xlsx"))
  }
  message("Done !")

  return(data_isoform_ambiguity)
}


### Function to retrieve all isoform sequences from proteins
fetch_isoforms <- function(proteins, fasta){
  # get all canonical sequence from FASTA file
  fasta <- seqinr::read.fasta(fasta, seqtype = "AA")
  fasta <- lapply(fasta, as.character)
  names(fasta) <- sapply(strsplit(names(fasta), "\\|"), "[[", 2)
  fasta <- lapply(fasta, paste, collapse = "")

  if(!all(sub(";.*", "", proteins) %in% names(fasta))){
    missing_proteins <- sub(";.*", "", proteins)[!(sub(";.*", "", proteins) %in% names(fasta))]
    missing_proteins <- paste(missing_proteins, collapse = "; ")
    stop(paste(missing_proteins, "not in FASTA !"))
  }

  # fetching isoform sequence from UNIPROT database
  isoforms <- list()
  for(i in unique(proteins)){
    iso <- tryCatch(rbioapi::rba_uniprot_proteins(accession = i, isoforms = TRUE), # get isoform from uniprot
                    error=function(e) {
                      message(paste("No isoforms were found for", i))
                      return()
                    })

    if(!is.null(iso)){
      iso <- lapply(iso,
                    function(y){
                      y <- y[c("accession", "sequence")]
                      y$sequence <- y$sequence$sequence
                      as.data.frame(y)
                    })

      iso <- as.data.frame(Reduce(rbind, iso))
      iso$length_isoform <- nchar(iso$sequence)
      iso$canonical <- sapply(iso$sequence,
                              function(x) x == fasta[[sub(";.*", "", i)]], USE.NAMES = FALSE)

      iso$canonical_sequence <- iso$sequence[iso$canonical]
      iso <- iso[!iso$canonical,]
      if(nrow(iso)){
        iso$canonical <- NULL
        iso$length_canonical <- nchar(iso$canonical_sequence)

        iso$isoforms <- iso$accession
        iso$accession <- i
      }
      else{
        iso <- NULL
      }
    }

    isoforms[[i]] <- iso
  }

  isoforms <- as.data.frame(Reduce(rbind, isoforms))
  if(nrow(isoforms)){
    isoforms <- isoforms[,c("accession", "length_canonical", "canonical_sequence",
                            "isoforms", "length_isoform", "sequence")]
  }

  return(isoforms)
}

### Function to retrieve the sequence alignment from BLAST
# the function return position in the canonical form that aligns with the isoform
compare_sequences <- function(isoform, canonical) {
  # BLAST URL
  blast_url_post <- "https://blast.ncbi.nlm.nih.gov/BlastAlign.cgi?"
  blast_url_get <- "https://blast.ncbi.nlm.nih.gov/Blast.cgi?"

  # build query
  blast_response <- httr::POST(
    blast_url_post,
    body = list(
      CMD = "Put",
      PROGRAM = "blastp",
      BLAST_PROGRAMS = "blastp",
      BLAST_SPEC = "blast2seq",  # aligning two sequences
      QUERY = paste0(">Query\n", isoform),
      SUBJECTS = paste0(">Subject\n", canonical),
      FORMAT_TYPE = "XML"
    ),
    encode = "multipart"
  )

  # Check if job sent
  if(httr::status_code(blast_response) != 200){
    stop("Bad query")
  }

  # Get Request ID (RID) to retrieve results
  response_text <- httr::content(blast_response, "text")
  # cat("Response from BLAST API:\n", response_text, "\n")  --> for debugging
  rid <- sub(".*RID = (\\w+).*", "\\1", response_text)
  message(paste("RID:", rid))
  if(rid == response_text){
    stop("No RID, check query")
  }

  Sys.sleep(5) # wait for job to be done
  results <- NULL
  for(i in 1:10){ # try to fetch results
    result_response <- httr::GET(
      blast_url_get,
      query = list(
        CMD = "Get",
        RID = rid,
        FORMAT_TYPE = "XML"
      )
    )

    if(httr::status_code(result_response) == 200){
      results <- httr::content(result_response, "text")
      if(grepl("BlastOutput", results)) break
    }
    Sys.sleep(5) # break between request
  }

  if (is.null(results) || !grepl("BlastOutput", results)) {
    stop("Results not ready or bad query")
  }
  results_xml <- xml2::read_xml(results)

  # Get main alignment
  alignments <- xml2::xml_find_all(results_xml, ".//Hsp")
  all_positions <- c()

  for(i in seq_along(alignments)){
    hsp <- alignments[[i]]

    isoform_seq <- xml2::xml_text(xml2::xml_find_first(hsp, ".//Hsp_qseq"))
    canonical_seq <- xml2::xml_text(xml2::xml_find_first(hsp, ".//Hsp_hseq"))
    start_pos <- as.integer(xml2::xml_text(xml2::xml_find_first(hsp, ".//Hsp_hit-from")))

    # get positions
    positions_logical <- strsplit(isoform_seq, "")[[1]] == strsplit(canonical_seq, "")[[1]]
    positions_idx <- which(positions_logical)
    positions <- positions_idx + start_pos - 1
    for(i in 2:(length(positions) - 1)){
      if(positions_idx[i] - positions_idx[i-1] == 1 & positions_idx[i+1] - positions_idx[i] == 1){
        positions[i] <- 0
      }
    }
    positions <- gsub("  ", "-",
                      gsub("-", "",
                           strsplit(gsub(" 0 ", "-",
                                         paste(positions, collapse = "  ")
                           ),
                           "  ")[[1]])
    )

    # correct position if gaps to realign on canonical sequence
    if(start_pos == 1 | nchar(canonical_seq) + start_pos - 1 > nchar(canonical)){
      pos_diff <- which(!positions_logical)
      if(length(pos_diff)){
        gaps <- strsplit(canonical_seq, "")[[1]][pos_diff]
        gaps <- pos_diff[which(gaps == "-")]
        gaps <- gaps + start_pos - 1
        if(length(gaps)){
          positions <- sapply(strsplit(positions, "-"), as.numeric, simplify = FALSE)
          pos_inv <- list()
          for(i in 1:(length(positions) - 1)){
            pos_inv[[i]] <- c(positions[[i]][length(positions[[i]])], positions[[i+1]][1])
          }
          gaps <- sapply(pos_inv, function(x) sum(x[1] < gaps & gaps < x[2]))
          for(i in 1:length(gaps)){
            positions[(i+1):length(positions)] <- sapply(positions[(i+1):length(positions)],
                                                         function(x) x - gaps[i],
                                                         simplify = FALSE)
          }
          positions <- paste(sapply(positions, paste, collapse = "-"), collapse = "; ")
        }
      }
    }

    all_positions <- c(all_positions, positions)
  }

  all_positions <- unique(all_positions)
  all_positions <- paste(all_positions, collapse = "; ")
  return(all_positions)
}
