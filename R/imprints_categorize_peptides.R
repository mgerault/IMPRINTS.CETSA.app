#' imprints_categorize_peptides
#'
#' Function to categorize proteins found as hits by the function \code{imprints_cleaved_peptides}
#'
#' @details
#' This function will categorize the hits found by the function \code{imprints_cleaved_peptides} in 4
#' main categories: RESP for REgional Stabilization after Proteolysis, SP for Single Peptide,
#' MP for multiple peptides and FP for false positive. When a small 'm' is added to SP or MP it stands
#' for modified.
#'
#' @param data The normalized peptides data set, i.e. the outpout from \code{imprints_normalize_peptides}.
#' @param data_cleaved The cleavage hits data set, i.e. the outpout from \code{imprints_cleaved_peptides}.
#' @param control The control treatment from your dataset.
#' @param save_xlsx Logical to tell if you want to save the categorized hits in an xlsx file.
#'   Default to TRUE.
#' @param xlsxname The name of your saved file.
#'
#' @return A data.frame containing the categorized hits. It will add the two columns 'category' and 'details'
#' to the data.frame returned by the function \code{imprints_cleaved_peptides}.
#'
#' @export
#'


imprints_categorize_peptides <- function(data, data_cleaved, control,
                                         save_xlsx = TRUE, xlsxname = "RESP_summary_categorized"){
  if(control %in% unique(data_cleaved$treatment)){
    message(paste("Error: Your data_cleaved contain the treatment", control,
                  "which you selected as a control"))
    return(NULL)
  }
  if(!(control %in% get_treat_level(data))){
    message(paste("Error: Your data doesn't contain the treatment", control,
                  "which you selected as a control"))
    return(NULL)
  }
  if(!all(unique(data_cleaved$treatment) %in% get_treat_level(data))){
    tr_missing <- unique(data_cleaved$treatment)[!(unique(data_cleaved$treatment) %in% get_treat_level(data))]
    tr_missing <- paste(tr_missing, collapse = ", ")
    message(paste("Error: Your data doesn't contain the treatment(s)", tr_missing,
                  "which are in your data_cleaved"))
    return(NULL)
  }
  if("countNum" %in% colnames(data)){
    data$countNum <- NULL
  }

  # computing FC from cleavage hits
  message("Computing fold-changes...")
  directory_toremove <- list.files()
  data_diff <- imprints_sequence_peptides(data, control = control,
                                          proteins = unique(data_cleaved$id),
                                          dataset_name = "imprints_categorize_peptides_FC")
  directory_toremove <- list.files()[!(list.files() %in% directory_toremove)]
  directory_toremove <- grep("imprints_categorize_peptides_FC\\.txt$", directory_toremove, value = TRUE)
  if(length(directory_toremove) == 1){
    unlink(directory_toremove, recursive = TRUE)
  }

  data_diff <- data_diff[,-grep(paste0("_", control, "$"), colnames(data_diff))]

  # reshaping and filtering FC data
  message("Reshaping data...")
  data_diff <- data_diff %>%
    tidyr::gather("key", "value", -Master.Protein.Accessions, -description,
                  -Positions.in.Master.Proteins, -Annotated.Sequence, -Modifications) %>%
    tidyr::separate(key, into = c("temperature", "rep", "treatment"), sep = "_")
  data_diff <- data_diff[,c("Master.Protein.Accessions", "description", "Positions.in.Master.Proteins",
                            "Modifications", "temperature", "rep", "treatment", "value")]
  colnames(data_diff)[1] <- "id"
  data_diff <- dplyr::left_join(unique(data_cleaved[,c("id", "treatment")]),
                                data_diff,
                                by = c("id", "treatment"))

  # categorizing hits
  message("Categorizing cleaved proteins...")
  data_categorized <- data_diff %>%
    dplyr::group_by(id, description,
                    Positions.in.Master.Proteins, Modifications,
                    treatment, temperature) %>%
    dplyr::summarise(value = mean(value, na.rm = TRUE)) %>%
    tidyr::spread(temperature, value) %>%
    dplyr::group_by(id, description, treatment) %>%
    dplyr::group_modify(~{
      if(any(grepl("^36C", colnames(.x)))){ # removing QP
        .x <- .x[,-grep("^36C", colnames(.x))]
      }
      .x <- .x[which(apply(.x[,-c(1,2)], 1,  # removing peptides with only missing values
                           function(y){
                             any(!is.na(y))
                           })
                     ),]

      # order peptides from N-term to C-term
      .x <- .x[order(sapply(strsplit(gsub(".* \\[|\\]", "",
                                          sub(";.*", "", .x$Positions.in.Master.Proteins)
                                          ), "-"),
                            function(y) sum(as.numeric(y)))
                     ),]

      df <- .x
      df[,-c(1,2)] <- apply(df[,-c(1,2)], 2, tidyr::replace_na, replace = 0)
      df_d <- dist(apply(abs(df[,-c(1,2)]), 1, max)) # taking absolute value; only looking for single peptide for now so sign doesn't matter

      ### hierarchical clustering peptides the best way: taking mean from the method
      df_h <- data.frame(comp = cutree(hclust(df_d, method = "complete"), 2),
                         ward = cutree(hclust(df_d, method = "ward.D"), 2),
                         ward2 = cutree(hclust(df_d, method = "ward.D2"), 2),
                         ave = cutree(hclust(df_d, method = "average"), 2))
      df_h_name <- apply(df_h, 2,
                         function(y){
                           sapply(1:2, function(z) mean(apply(abs(df[,-c(1,2)]), 1, max)[y == z]))
                         })
      if(!all(apply(df_h_name, 2, which.max) == 1)){
        df_h[,which(apply(df_h_name, 2, which.max) == 2)] <-
          apply(df_h[,which(apply(df_h_name, 2, which.max) == 2),drop=FALSE], 2,
                function(y){
                  y[y == 1] <- 3
                  y[y == 2] <- 1
                  y[y == 3] <- 2;
                  y
                })
      }
      df_h <- apply(df_h, 1, mean)

      if(any(df_h == 1.5)){ # ambiguity, removing peptide
        p_torm <- which(df_h == 1.5)
        .x <- .x[-p_torm,]
        df <- df[-p_torm,]
        df_h <- df_h[-p_torm]
      }
      df_h <- round(df_h)
      ### final clustering of peptides done

      # saving cluster with less peptides
      df_mincl <- table(df_h)
      if(length(unique(df_mincl)) == 1){
        df_m <- apply(.x[,-c(1,2)], 1, function(y) y[which.max(abs(y))])
        df_m <- sapply(1:2, function(y) mean(abs(df_m[df_h == y])))
        df_mincl <- which.max(df_m)
      }
      else{
        df_mincl <- as.numeric(names(df_mincl)[which.min(df_mincl)])
      }
      pep <- .x[df_h == df_mincl,] # need to keep original FC values

      ### CATEGORIZATION ###
      if(nrow(pep) == 1){# if only one peptide was cluster in a group
        pep_s <- sign(apply(pep[,-c(1,2)], 1, mean, na.rm = TRUE))
        modif <- strsplit(pep$Modifications, "(?<=\\]); (?=\\d{1}x)", perl =  TRUE)[[1]]
        modif <- modif[-grep("\\d{1}xTMT", modif)]
        if(pep_s == 1){
          category <- paste("single peptide",
                            sub(".* ", "",
                                sub("; .*", "", pep$Positions.in.Master.Proteins)),
                            "positively shifting")
        }
        else{
          category <- paste("single peptide",
                            sub(".* ", "",
                                sub("; .*", "", pep$Positions.in.Master.Proteins)),
                            "negatively shifting")
        }
        if(length(modif)){
          category <- paste(category, paste(modif, collapse = "; ")
          )
        }
      }
      else{
        pep_seq <- gsub(".* \\[|\\]", "",
                        sub(";.*", "", pep$Positions.in.Master.Proteins)
                        )
        pep_seq <- lapply(strsplit(pep_seq, "-"), as.numeric)
        pep_mlength <- mean(sapply(pep_seq, diff))

        pep_seq_min <- min(sapply(pep_seq, "[[", 1))
        pep_seq_max <- max(sapply(pep_seq, "[[", 2))

        # check if peptides overlap
        pep_allin <- pep_seq_max - pep_seq_min <= pep_mlength + 4 # 4 is the AA margin
        if(pep_allin){ # if yes, most likely only one peptide shifting
          pep_s <- sign(apply(pep[,-c(1,2)], 1, mean, na.rm = TRUE))
          modif <- strsplit(pep$Modifications, "(?<=\\]); (?=\\d{1}x)", perl =  TRUE)
          modif <- lapply(modif, function(y) y[-grep("\\d{1}xTMT", y)])
          modif <- unique(unlist(modif))
          if(length(unique(pep_s)) != 1){
            category <- "single peptide with different profiles"
          }
          else{
            pep_s <- unique(pep_s)
            if(pep_s == 1){
              category <- paste("single peptide",
                                paste0("[", pep_seq_min, "-", pep_seq_max, "]"),
                                "positively shifting")
            }
            else{
              category <- paste("single peptide",
                                paste0("[", pep_seq_min, "-", pep_seq_max, "]"),
                                "negatively shifting")
            }
            if(length(modif)){
              category <- paste(category, paste(modif, collapse = ", ")
              )
            }
          }
        }
        else{
          # if non significant --> random --> most likely FP or isolated peptides shifting
          # if significant --> non random --> RESP or 'opposite' peptide shifting
          # significance can below 0.1
          pv_rnd <- randtests::runs.test(df_h, threshold = 1.5)$p.value
          if(pv_rnd <= 0.1){
            df_m <- apply(.x[,-c(1,2)], 1, function(y) y[which.max(abs(y))])
            pepm <- apply(.x[df_h == which.max(sapply(1:2, function(y) abs(mean(df_m[df_h == y]))
                                                      )
                                               ),
                             -c(1,2)],
                          1, function(y) y[which.max(abs(y))])
            pepm_p <- max(table(sign(pepm)))/length(pepm)
            if(length(pepm) <= 3){ # if only 2 or 3 peptides (can't be 1) and at least 2 of them are modified, most likely just a MPm
              new_pepm <- .x[df_m %in% pepm,]
              new_seq <- gsub(".* ", "", sub(";.*", "", new_pepm$Positions.in.Master.Proteins))
              modif <- strsplit(new_pepm$Modifications, "(?<=\\]); (?=\\d{1}x)", perl =  TRUE)
              modif <- lapply(modif, function(y){
                y <- y[-grep("\\d{1}xTMT", y)]
                if(length(y) == 0){
                  y <- ""
                };
                paste(y, collapse = " ")
              })
              if(sum(sapply(modif, nchar) > 0) >= 2){
                if(length(unique(sign(pepm))) == 1){
                  category <- paste(paste(new_seq, modif, collapse = ", "),
                                    ifelse(unique(sign(pepm)) == 1, "positively", "negatively"),
                                    "shifting")
                }
                else{
                  category <- paste(paste(paste(new_seq[which(sign(pepm) == 1)], modif[which(sign(pepm) == 1)], collapse = ", "),
                                          "positively shifting and"),
                                    paste(paste(new_seq[which(sign(pepm) == -1)], modif[which(sign(pepm) == -1)], collapse = ", "),
                                          "negatively shifting and")
                                    )
                }
              }
              else{
                category <- "RESP - low confidence"
              }
            }
            else if(pepm_p >= 0.8){ # checking that part with greater mean has overall same sign in shifting
              if(stats::shapiro.test(df_m)$p.value <= 0.1){ # checking if the max FC from each peptide are normally distributed
                category <- "RESP" # if not, most likely two parts --> RESP
              }
              else{
                pvm <- sapply(1:2, function(y) t.test(df_m[df_h == y])$p.value)
                pvm_m <- sapply(1:2, function(y) mean(df_m[df_h == y]))
                if(all(pvm <= 0.01)){ # both part are away from 0
                  if(sum(abs(pvm_m) <= 0.1) == 1){ # if only one part is below 0.1
                    category <- "RESP - low confidence"
                  }
                  else{
                    category <- "false positive"
                  }
                }
                else if(any(pvm <= 0.01) & any(pvm > 0.1)){ # only one is far from 0 and the other not
                  category <- "RESP"
                }
                else{ # more ambigous
                  category <- "to check"
                }
              }
            }
            else{ # checking each case when sign of shifting part are ~randomly distributed
              pepm_neg <- .x[df_m %in% pepm[sign(pepm) == -1],]
              pepm_pos <- .x[df_m %in% pepm[sign(pepm) == 1],]

              pepm_neg_seq <- lapply(strsplit(gsub(".* \\[|\\]", "",
                                                   sub(";.*", "", pepm_neg$Positions.in.Master.Proteins)
                                                   ), "-"),
                                     as.numeric)
              pepm_pos_seq <- lapply(strsplit(gsub(".* \\[|\\]", "",
                                                   sub(";.*", "", pepm_pos$Positions.in.Master.Proteins)
                                                   ), "-"),
                                     as.numeric)

              # checking if all peptides are same sequence, +/- 2 AA
              pepm_neg_seq_in <- sapply(pepm_neg_seq,
                                        function(y){
                                          all(sapply(pepm_neg_seq,
                                                     function(z){
                                                       z[1] - 2 <= y[1] & z[1] + 2 >= y[1] & z[2] - 2 <= y[2] & z[2] + 2 >= y[2]
                                                     }))
                                        })
              pepm_pos_seq_in <- sapply(pepm_pos_seq,
                                        function(y){
                                          all(sapply(pepm_pos_seq,
                                                     function(z){
                                                       z[1] - 2 <= y[1] & z[1] + 2 >= y[1] & z[2] - 2 <= y[2] & z[2] + 2 >= y[2]
                                                     }))
                                        })

              if(all(pepm_neg_seq_in)){ # if yes, merging them
                new_neg_seq <- paste0("[", min(sapply(pepm_neg_seq, "[[", 1)),
                                      "-", max(sapply(pepm_neg_seq, "[[", 2)), "]")
                modif_neg <- strsplit(pepm_neg$Modifications, "(?<=\\]); (?=\\d{1}x)", perl =  TRUE)
                modif_neg <- lapply(modif_neg, function(y) y[-grep("\\d{1}xTMT", y)])
                modif_neg <- unique(unlist(modif_neg))
                modif_neg <- ifelse(length(modif_neg), paste(modif_neg, collapse = ", "), "")

                category_neg <- paste(new_neg_seq, "negatively shifting", modif_neg)
              }
              else{
                new_neg_seq <- gsub(".* ", "", sub(";.*", "", pepm_neg$Positions.in.Master.Proteins))
                modif_neg <- strsplit(pepm_neg$Modifications, "(?<=\\]); (?=\\d{1}x)", perl =  TRUE)
                modif_neg <- lapply(modif_neg, function(y){
                  y <- y[-grep("\\d{1}xTMT", y)]
                  if(length(y) == 0){
                    y <- ""
                  };
                  paste(y, collapse = " ")
                })

                category_neg <- paste(paste(new_neg_seq, modif_neg, collapse = ", "), "negatively shifting")
              }
              if(all(pepm_pos_seq_in)){
                new_pos_seq <- paste0("[", min(sapply(pepm_pos_seq, "[[", 1)),
                                      "-", max(sapply(pepm_pos_seq, "[[", 2)), "]")
                modif_pos <- strsplit(pepm_pos$Modifications, "(?<=\\]); (?=\\d{1}x)", perl =  TRUE)
                modif_pos <- lapply(modif_pos, function(y) y[-grep("\\d{1}xTMT", y)])
                modif_pos <- unique(unlist(modif_pos))
                modif_pos <- ifelse(length(modif_pos), paste(modif_pos, collapse = ", "), "")

                category_pos <- paste(new_pos_seq, "positively shifting", modif_pos)
              }
              else{
                new_pos_seq <- gsub(".* ", "", sub(";.*", "", pepm_pos$Positions.in.Master.Proteins))
                modif_pos <- strsplit(pepm_pos$Modifications, "(?<=\\]); (?=\\d{1}x)", perl =  TRUE)
                modif_pos <- lapply(modif_pos, function(y){
                  y <- y[-grep("\\d{1}xTMT", y)]
                  if(length(y) == 0){
                    y <- ""
                  };
                  paste(y, collapse = " ")
                })

                category_pos <- paste(paste(new_pos_seq, modif_pos, collapse = ", "), "positively shifting")
              }

              category <- paste(category_neg, "and", category_pos)

              # if too many peptides to report, most likely a FALSE positive
              nb_reported_pep <- length(grep("\\[\\d{1,4}-\\d{1,4}\\]",
                                             strsplit(category, " ")[[1]])
                                        )
              if(nb_reported_pep > 3){
                category <- "false positive"
              }
            }
          }
          else{
            df_m <- apply(.x[,-c(1,2)], 1, function(y) y[which.max(abs(y))])
            df_pv <- sapply(1:2, function(y) t.test(df_m[df_h == y])$p.value)
            df_pvm <- sapply(1:2, function(y) mean(df_m[df_h == y]))
            df_pvma <- sapply(1:2, function(y) mean(abs(df_m[df_h == y])))

            if(sum(df_h == df_mincl) == 2){ # only 2 peptides in one cluster
              if(df_mincl == which.min(df_pvma)){ # if those are the minimum FC --> FP
                category <- "false positive"
              }
              else{# 2 peptides shifting at different places
                pepm <- df_m[df_h == df_mincl]

                pepm_neg <- .x[df_m %in% pepm[sign(pepm) == -1],]
                pepm_pos <- .x[df_m %in% pepm[sign(pepm) == 1],]

                category_pos <- ""
                category_neg <- ""

                if(nrow(pepm_neg)){
                  pepm_neg_seq <- lapply(strsplit(gsub(".* \\[|\\]", "",
                                                       sub(";.*", "", pepm_neg$Positions.in.Master.Proteins)
                                                       ), "-"),
                                         as.numeric)

                  # checking if all peptides are same sequence, +/- 2 AA
                  pepm_neg_seq_in <- sapply(pepm_neg_seq,
                                            function(y){
                                              all(sapply(pepm_neg_seq,
                                                         function(z){
                                                           z[1] - 2 <= y[1] & z[1] + 2 >= y[1] & z[2] - 2 <= y[2] & z[2] + 2 >= y[2]
                                                         }))
                                            })

                  if(all(pepm_neg_seq_in)){ # if yes, merging them
                    new_neg_seq <- paste0("[", min(sapply(pepm_neg_seq, "[[", 1)),
                                          "-", max(sapply(pepm_neg_seq, "[[", 2)), "]")
                    modif_neg <- strsplit(pepm_neg$Modifications, "(?<=\\]); (?=\\d{1}x)", perl =  TRUE)
                    modif_neg <- lapply(modif_neg, function(y) y[-grep("\\d{1}xTMT", y)])
                    modif_neg <- unique(unlist(modif_neg))
                    modif_neg <- ifelse(length(modif_neg), paste(modif_neg, collapse = ", "), "")

                    category_neg <- paste(new_neg_seq, "negatively shifting", modif_neg)
                  }
                  else{
                    new_neg_seq <- gsub(".* ", "", sub(";.*", "", pepm_neg$Positions.in.Master.Proteins))
                    modif_neg <- strsplit(pepm_neg$Modifications, "(?<=\\]); (?=\\d{1}x)", perl =  TRUE)
                    modif_neg <- lapply(modif_neg, function(y){
                      y <- y[-grep("\\d{1}xTMT", y)]
                      if(length(y) == 0){
                        y <- ""
                      };
                      paste(y, collapse = " ")
                    })

                    category_neg <- paste(paste(new_neg_seq, modif_neg, collapse = ", "), "negatively shifting")
                  }
                }

                if(nrow(pepm_pos)){
                  pepm_pos_seq <- lapply(strsplit(gsub(".* \\[|\\]", "",
                                                       sub(";.*", "", pepm_pos$Positions.in.Master.Proteins)
                                                       ), "-"),
                                         as.numeric)

                  # checking if all peptides are same sequence, +/- 2 AA
                  pepm_pos_seq_in <- sapply(pepm_pos_seq,
                                            function(y){
                                              all(sapply(pepm_pos_seq,
                                                         function(z){
                                                           z[1] - 2 <= y[1] & z[1] + 2 >= y[1] & z[2] - 2 <= y[2] & z[2] + 2 >= y[2]
                                                         }))
                                            })

                  if(all(pepm_pos_seq_in)){
                    new_pos_seq <- paste0("[", min(sapply(pepm_pos_seq, "[[", 1)),
                                          "-", max(sapply(pepm_pos_seq, "[[", 2)), "]")
                    modif_pos <- strsplit(pepm_pos$Modifications, "(?<=\\]); (?=\\d{1}x)", perl =  TRUE)
                    modif_pos <- lapply(modif_pos, function(y) y[-grep("\\d{1}xTMT", y)])
                    modif_pos <- unique(unlist(modif_pos))
                    modif_pos <- ifelse(length(modif_pos), paste(modif_pos, collapse = ", "), "")

                    category_pos <- paste(new_pos_seq, "positively shifting", modif_pos)
                  }
                  else{
                    new_pos_seq <- gsub(".* ", "", sub(";.*", "", pepm_pos$Positions.in.Master.Proteins))
                    modif_pos <- strsplit(pepm_pos$Modifications, "(?<=\\]); (?=\\d{1}x)", perl =  TRUE)
                    modif_pos <- lapply(modif_pos, function(y){
                      y <- y[-grep("\\d{1}xTMT", y)]
                      if(length(y) == 0){
                        y <- ""
                      };
                      paste(y, collapse = " ")
                    })

                    category_pos <- paste(paste(new_pos_seq, modif_pos, collapse = ", "), "positively shifting")
                  }
                }

                if(nchar(category_neg) & nchar(category_pos)){
                  category <- paste(category_neg, "and", category_pos)
                }
                else if(nchar(category_neg) == 0){
                  category <- category_pos
                }
                else if(nchar(category_pos) == 0){
                  category <- category_neg
                }
              }
            }
            else if(all(df_pv <= 0.01) & length(unique(sign(df_pvm))) == 1){ # two cluster significantly away from 0 and same sign
              if(max(df_pvma)/min(df_pvma) >= 3){ # if difference of maximum is quite high could be RESP of lower confidence
                category <- "RESP - low confidence"
              }
              else{
                category <- "false positive"
              }
            }
            else if(df_mincl == which.min(df_pvma)){
              df_lowpos <- strsplit(paste(df_h, collapse = ""), as.character(df_mincl))[[1]]
              if(mean(nchar(df_lowpos[nchar(df_lowpos) > 0])) >= 1){ # peptide with low FC spread out --> FP
                category <- "false positive"
              }
              else{ # very less likely
                category <- "to check"
              }
            }
            else{
              # this gives the number of amino acid between each group
              df_lowpos <- gsub(paste0("^", unique(df_h[df_h != df_mincl]), "*(?=",
                                       df_mincl, ")|(?<=", df_mincl, ")",
                                       unique(df_h[df_h != df_mincl]), "*$"),
                                "", paste(df_h, collapse = ""), perl =  TRUE)
              df_lowpos <- strsplit(df_lowpos, as.character(df_mincl))[[1]]
              df_lowpos <- nchar(df_lowpos[nchar(df_lowpos) > 0])
              # if only one separation (or two?) with only one or two amino acid
              # and that this group has peptides of same FC sign --> RESP (of lower confidence)
              # if lot of separation between the group, most likely a FP
              # if not lot of separation, but far appart, most likely single peptides shifting

              if(length(df_lowpos) <= 2){
                if(all(df_lowpos <= 2) & length(unique(sign(df_m[df_h == df_mincl]))) == 1){
                  if(sum(df_h == df_mincl) == 3){ # if only 3 peptides (can't be 2 or 1) and at least 2 of them are modified, most likely just a MPm
                    pepm <- .x[df_h == df_mincl,]
                    new_seq <- gsub(".* ", "", sub(";.*", "", pepm$Positions.in.Master.Proteins))
                    modif <- strsplit(pepm$Modifications, "(?<=\\]); (?=\\d{1}x)", perl =  TRUE)
                    modif <- lapply(modif, function(y){
                      y <- y[-grep("\\d{1}xTMT", y)]
                      if(length(y) == 0){
                        y <- ""
                      };
                      paste(y, collapse = " ")
                    })
                    if(sum(sapply(modif, nchar) > 0) >= 2){
                      category <- paste(paste(new_seq, modif, collapse = ", "),
                                        ifelse(unique(sign(df_m[df_h == df_mincl])) == 1, "positively", "negatively"),
                                        "shifting")
                    }
                    else{
                      category <- "RESP - low confidence"
                    }
                  }
                  else{
                    category <- "RESP - low confidence"
                  }
                }
                else{
                  ### single peptide routine
                  pepm <- df_m[df_h == df_mincl]

                  pepm_neg <- .x[df_m %in% pepm[sign(pepm) == -1],]
                  pepm_pos <- .x[df_m %in% pepm[sign(pepm) == 1],]

                  category_pos <- ""
                  category_neg <- ""

                  if(nrow(pepm_neg)){
                    pepm_neg_seq <- lapply(strsplit(gsub(".* \\[|\\]", "",
                                                         sub(";.*", "", pepm_neg$Positions.in.Master.Proteins)
                                                         ), "-"),
                                           as.numeric)

                    # checking if all peptides are same sequence, +/- 2 AA
                    pepm_neg_seq_in <- sapply(pepm_neg_seq,
                                              function(y){
                                                which(sapply(pepm_neg_seq,
                                                             function(z){
                                                               z[1] - 2 <= y[1] & z[1] + 2 >= y[1] & z[2] - 2 <= y[2] & z[2] + 2 >= y[2]
                                                             }))
                                              })

                    if(any(duplicated(pepm_neg_seq_in))){ # if yes, merging them
                      pepm_neg_seq_in <- pepm_neg_seq_in[-which(duplicated(pepm_neg_seq_in))]

                      category_neg <- sapply(pepm_neg_seq_in,
                                             function(y){
                                               new_neg_seq <- paste0("[", min(sapply(pepm_neg_seq[y], "[[", 1)),
                                                                     "-", max(sapply(pepm_neg_seq[y], "[[", 2)), "]"
                                               )
                                               modif_neg <- strsplit(pepm_neg$Modifications[y], "(?<=\\]); (?=\\d{1}x)", perl =  TRUE)
                                               modif_neg <- lapply(modif_neg, function(y) y[-grep("\\d{1}xTMT", y)])
                                               modif_neg <- unique(unlist(modif_neg))
                                               modif_neg <- ifelse(length(modif_neg), paste(modif_neg, collapse = "; "), "")

                                               paste(new_neg_seq, modif_neg)
                                             })

                      category_neg <- paste(paste(category_neg, collapse = ", "), "negatively shifting")
                    }
                    else{
                      new_neg_seq <- gsub(".* ", "", sub(";.*", "", pepm_neg$Positions.in.Master.Proteins))
                      modif_neg <- strsplit(pepm_neg$Modifications, "(?<=\\]); (?=\\d{1}x)", perl =  TRUE)
                      modif_neg <- lapply(modif_neg, function(y){
                        y <- y[-grep("\\d{1}xTMT", y)]
                        if(length(y) == 0){
                          y <- ""
                        };
                        paste(y, collapse = " ")
                      })

                      category_neg <- paste(paste(new_neg_seq, modif_neg, collapse = ", "), "negatively shifting")
                    }
                  }

                  if(nrow(pepm_pos)){
                    pepm_pos_seq <- lapply(strsplit(gsub(".* \\[|\\]", "",
                                                         sub(";.*", "", pepm_pos$Positions.in.Master.Proteins)
                                                         ), "-"),
                                           as.numeric)

                    # checking if all peptides are same sequence, +/- 2 AA
                    pepm_pos_seq_in <- sapply(pepm_pos_seq,
                                              function(y){
                                                which(sapply(pepm_pos_seq,
                                                             function(z){
                                                               z[1] - 2 <= y[1] & z[1] + 2 >= y[1] & z[2] - 2 <= y[2] & z[2] + 2 >= y[2]
                                                             }))
                                              })

                    if(any(duplicated(pepm_pos_seq_in))){
                      pepm_pos_seq_in <- pepm_pos_seq_in[-which(duplicated(pepm_pos_seq_in))]

                      category_pos <- sapply(pepm_pos_seq_in,
                                             function(y){
                                               new_pos_seq <- paste0("[", min(sapply(pepm_pos_seq[y], "[[", 1)),
                                                                     "-", max(sapply(pepm_pos_seq[y], "[[", 2)), "]")
                                               modif_pos <- strsplit(pepm_pos$Modifications[y], "(?<=\\]); (?=\\d{1}x)", perl =  TRUE)
                                               modif_pos <- lapply(modif_pos, function(y) y[-grep("\\d{1}xTMT", y)])
                                               modif_pos <- unique(unlist(modif_pos))
                                               modif_pos <- ifelse(length(modif_pos), paste(modif_pos, collapse = "; "), "")

                                               paste(new_pos_seq, modif_pos)
                                             })

                      category_pos <- paste(paste(category_pos, collapse = ", "), "positively shifting")
                    }
                    else{
                      new_pos_seq <- gsub(".* ", "", sub(";.*", "", pepm_pos$Positions.in.Master.Proteins))
                      modif_pos <- strsplit(pepm_pos$Modifications, "(?<=\\]); (?=\\d{1}x)", perl =  TRUE)
                      modif_pos <- lapply(modif_pos, function(y){
                        y <- y[-grep("\\d{1}xTMT", y)]
                        if(length(y) == 0){
                          y <- ""
                        };
                        paste(y, collapse = " ")
                      })

                      category_pos <- paste(paste(new_pos_seq, modif_pos, collapse = ", "), "positively shifting")
                    }
                  }

                  if(nchar(category_neg) & nchar(category_pos)){
                    category <- paste(category_neg, "and", category_pos)
                  }
                  else if(nchar(category_neg) == 0){
                    category <- category_pos
                  }
                  else if(nchar(category_pos) == 0){
                    category <- category_neg
                  }

                  # if too many peptides to report, most likely a FALSE positive
                  nb_reported_pep <- length(grep("\\[\\d{1,4}-\\d{1,4}\\]",
                                                 strsplit(category, " ")[[1]])
                                            )
                  if(nb_reported_pep > 3){
                    category <- "false positive"
                  }
                }
              }
              else{
                ### single peptide routine
                pepm <- df_m[df_h == df_mincl]

                pepm_neg <- .x[df_m %in% pepm[sign(pepm) == -1],]
                pepm_pos <- .x[df_m %in% pepm[sign(pepm) == 1],]

                category_pos <- ""
                category_neg <- ""

                if(nrow(pepm_neg)){
                  pepm_neg_seq <- lapply(strsplit(gsub(".* \\[|\\]", "",
                                                       sub(";.*", "", pepm_neg$Positions.in.Master.Proteins)
                                                       ), "-"),
                                         as.numeric)

                  # checking if all peptides are same sequence, +/- 2 AA
                  pepm_neg_seq_in <- sapply(pepm_neg_seq,
                                            function(y){
                                              which(sapply(pepm_neg_seq,
                                                           function(z){
                                                             z[1] - 2 <= y[1] & z[1] + 2 >= y[1] & z[2] - 2 <= y[2] & z[2] + 2 >= y[2]
                                                           }))
                                            })

                  if(any(duplicated(pepm_neg_seq_in))){ # if yes, merging them
                    pepm_neg_seq_in <- pepm_neg_seq_in[-which(duplicated(pepm_neg_seq_in))]

                    category_neg <- sapply(pepm_neg_seq_in,
                                           function(y){
                                             new_neg_seq <- paste0("[", min(sapply(pepm_neg_seq[y], "[[", 1)),
                                                                   "-", max(sapply(pepm_neg_seq[y], "[[", 2)), "]")
                                             modif_neg <- strsplit(pepm_neg$Modifications[y], "(?<=\\]); (?=\\d{1}x)", perl =  TRUE)
                                             modif_neg <- lapply(modif_neg, function(y) y[-grep("\\d{1}xTMT", y)])
                                             modif_neg <- unique(unlist(modif_neg))
                                             modif_neg <- ifelse(length(modif_neg), paste(modif_neg, collapse = "; "), "")

                                             paste(new_neg_seq, modif_neg)
                                           })

                    category_neg <- paste(paste(category_neg, collapse = ", "), "negatively shifting")
                  }
                  else{
                    new_neg_seq <- gsub(".* ", "", sub(";.*", "", pepm_neg$Positions.in.Master.Proteins))
                    modif_neg <- strsplit(pepm_neg$Modifications, "(?<=\\]); (?=\\d{1}x)", perl =  TRUE)
                    modif_neg <- lapply(modif_neg, function(y){
                      y <- y[-grep("\\d{1}xTMT", y)]
                      if(length(y) == 0){
                        y <- ""
                      };
                      paste(y, collapse = " ")
                    })

                    category_neg <- paste(paste(new_neg_seq, modif_neg, collapse = ", "), "negatively shifting")
                  }
                }

                if(nrow(pepm_pos)){
                  pepm_pos_seq <- lapply(strsplit(gsub(".* \\[|\\]", "",
                                                       sub(";.*", "", pepm_pos$Positions.in.Master.Proteins)
                                                       ), "-"),
                                         as.numeric)

                  # checking if all peptides are same sequence, +/- 2 AA
                  pepm_pos_seq_in <- sapply(pepm_pos_seq,
                                            function(y){
                                              which(sapply(pepm_pos_seq,
                                                           function(z){
                                                             z[1] - 2 <= y[1] & z[1] + 2 >= y[1] & z[2] - 2 <= y[2] & z[2] + 2 >= y[2]
                                                           }))
                                            })

                  if(any(duplicated(pepm_pos_seq_in))){
                    pepm_pos_seq_in <- pepm_pos_seq_in[-which(duplicated(pepm_pos_seq_in))]

                    category_pos <- sapply(pepm_pos_seq_in,
                                           function(y){
                                             new_pos_seq <- paste0("[", min(sapply(pepm_pos_seq[y], "[[", 1)),
                                                                   "-", max(sapply(pepm_pos_seq[y], "[[", 2)), "]")
                                             modif_pos <- strsplit(pepm_pos$Modifications[y], "(?<=\\]); (?=\\d{1}x)", perl =  TRUE)
                                             modif_pos <- lapply(modif_pos, function(y) y[-grep("\\d{1}xTMT", y)])
                                             modif_pos <- unique(unlist(modif_pos))
                                             modif_pos <- ifelse(length(modif_pos), paste(modif_pos, collapse = "; "), "")

                                             paste(new_pos_seq, modif_pos)
                                           })

                    category_pos <- paste(paste(category_pos, collapse = ", "), "positively shifting")
                  }
                  else{
                    new_pos_seq <- gsub(".* ", "", sub(";.*", "", pepm_pos$Positions.in.Master.Proteins))
                    modif_pos <- strsplit(pepm_pos$Modifications, "(?<=\\]); (?=\\d{1}x)", perl =  TRUE)
                    modif_pos <- lapply(modif_pos, function(y){
                      y <- y[-grep("\\d{1}xTMT", y)]
                      if(length(y) == 0){
                        y <- ""
                      };
                      paste(y, collapse = " ")
                    })

                    category_pos <- paste(paste(new_pos_seq, modif_pos, collapse = ", "), "positively shifting")
                  }
                }

                if(nchar(category_neg) & nchar(category_pos)){
                  category <- paste(category_neg, "and", category_pos)
                }
                else if(nchar(category_neg) == 0){
                  category <- category_pos
                }
                else if(nchar(category_pos) == 0){
                  category <- category_neg
                }

                # if too many peptides to report, most likely a FALSE positive
                nb_reported_pep <- length(grep("\\[\\d{1,4}-\\d{1,4}\\]",
                                               strsplit(category, " ")[[1]])
                                          )
                if(nb_reported_pep > 3){
                  category <- "false positive"
                }
              }
            }
          }
        }
      }

      df <- data.frame(details = category)
      if(grepl("RESP", category)){
        category <- "RESP"
      }
      else if(grepl("single", category)){
        category <- ifelse(grepl("\\d{1}x.* \\[", category),
                           "SPm", "SP")
      }
      else if(grepl("false positive", category)){
        category <- "FP"
      }
      else{
        category <- ifelse(grepl("\\d{1}x.* \\[", category),
                           "MPm", "MP")
      }

      df$category <- category
      return(df)
    }) %>%
    dplyr::mutate(Gene = IMPRINTS.CETSA.app:::getGeneName(description),
                  description = IMPRINTS.CETSA.app:::getProteinName(description))

  res <- dplyr::full_join(data_cleaved, data_categorized,
                          by = c("id", "Gene", "description", "treatment"))
  res <- res[,c("id", "Gene", "description", "treatment",
                "combined_pvalue", "RESP_score", "cleaved_site",
                "Nvalue_N-term", "Nvalue_C-term",
                "Npep_N-term", "Npep_C-term",
                "category", "details")]

  if(save_xlsx){
    openxlsx::write.xlsx(res, paste0(format(Sys.time(), "%y%m%d_%H%M_"), xlsxname, ".xlsx"))
  }

  message("Done !")
  return(res)
}
