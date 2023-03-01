#' imprints_sequence_peptides
#'
#' Function to select proteins and sequence from peptide data, sum up the peptide that are not selected
#' (before and after the sequence selected) and then plot the fold change on a pdf.
#'
#' @param data The input data set, the peptide data file filtered (non TMT modification removed) and non log2 transformed.
#' @param proteins The proteins you want to select. If NULL, select all.
#' @param sequence The peptide position you want to select. It needs to be a list that contains the same
#'                 number of element than the number of proteins or only one element. An element can be a string
#'                 like 191-199 or a vector that contains several strings like \code{c("52-71", "207-222")}.
#'                 Anyway, it needs to be in the format number-number.
#'                 If it's NULL, it will only select according the proteins so you can plot all the peptides from those.
#' @param control The condition corresponding to the control in your dataset.
#' @param dataset_name The name of your dataset to save your file.
#'
#' @return The caldiff output from the peptide position from the protein you selected; save it and also save the barplots.
#'
#' @export
#'
#'

imprints_sequence_peptides <- function(data, proteins = NULL, sequence = NULL, control = "Vehicle",
                                       dataset_name = ""){
  if(!is.null(sequence)){
    if(!(length(proteins) == length(sequence) | length(sequence) == 1) & !inherits(sequence, "list")){
      message("Error: The number of sequence needs to match the number of proteins or just provide one sequence")
      return()
    }
  }

  if(!is.null(proteins)){
    protein_in <- proteins %in% data$`Master Protein Accessions`

    if(all(!protein_in)){
      message("Error: None of the proteins you selected are in your data, please check your proteins.")
      return()
    }
    else if(any(!protein_in)){
      message(paste(paste(proteins[!protein_in], collapse = ", "), "were not found in your data and hence, has been removed."))
      proteins <- proteins[protein_in]

      if(!is.null(sequence) & length(sequence) > 1){
        sequence <- sequence[protein_in]
      }
    }

    data <- data[which(!is.na(match(data$`Master Protein Accessions`,
                                    proteins))
    ),]
  }
  else{
    proteins <- data$`Master Protein Accessions`
  }

  # reorder the data --> order according proteins and peptide sequence (~position)
  data$factor <- data$`Master Protein Accessions`
  data$`Master Protein Accessions` <- paste(data$`Positions in Master Proteins`,
                                            data$`Annotated Sequence`,
                                            data$Modifications, sep = " \n")
  data <- data[order(data$`Master Protein Accessions`),]
  ord <- lapply(stringr::str_split(data$`Master Protein Accessions`, "\\[|\\]"),
                function(x) {z = stringr::str_split(x[2], "-");
                z = as.numeric(z[[1]])
                z = sum(z);
                z
                }
  )
  ord <- unlist(ord)
  ord_f <- data.frame(factor = data$factor, x = ord)
  ord_f <- ord_f %>% dplyr::group_by(factor) %>%
    dplyr::mutate(y = order(x))
  b = (ord_f %>% dplyr::group_by(factor) %>% dplyr::count())$n
  if(length(b) > 1){
    b1 = b
    for(i in length(b):2){
      b1[i] <- sum(b1[(i-1):1])
    }
    b1[1] <- 0
  }
  else{
    b1 <- 0
  }
  b1 = rep(b1, b)
  ord_f$y <- ord_f$y + b1

  data <- data[ord_f$y,]
  data$`Master Protein Accessions` <- data$factor
  data$factor <- NULL

  cond_data <- unique(unlist(lapply(stringr::str_split(stringr::str_subset(colnames(data), "^\\d{1,}"), "_"), function(x) x[3])))
  cond_data <- cond_data[!(cond_data %in% control)]
  temp_data <- unique(unlist(lapply(stringr::str_split(stringr::str_subset(colnames(data), "^\\d{1,}"), "_"), function(x) x[1])))

  final_res <- list()
  if(!is.null(sequence)){
    message("Filtering and summing sequences...")
    for (i in 1:length(proteins)){
      message(paste0(proteins[i], " : ", i, "/", length(proteins)))

      tab_pr <- data[stringr::str_detect(data$`Master Protein Accessions`, proteins[i]),]
      tab_sequ = stringr::str_extract(tab_pr$`Positions in Master Proteins`, "(?<=\\[).+?(?=\\])")
      tab_sequ_n = stringr::str_split(tab_sequ, "-")

      if(length(sequence) == 1){
        seq = sequence[[1]]
      }
      else{
        seq = sequence[[i]]
      }

      if(length(seq) > 0){
        if(stringr::str_length(seq)[1] == 0){
          res <- tab_pr
        }
        else{
          seq <- seq[order(unlist(lapply(stringr::str_split(seq, "-"), function(x) as.numeric(x[2]))))]
          the_one <- lapply(stringr::str_split(seq, "-"), function(y){
            y <- as.numeric(y)
            x <- lapply(tab_sequ_n, function(x){
              x <- as.numeric(x)
              x[1] >= y[1] & x[2] <= y[2]  # if peptide sequence is 'in' the one selected, save it
            })
            x <- unlist(x);
            x
          })
          the_one <- lapply(the_one, which)
          the_one <- the_one[lapply(the_one, length) > 0] # only keep non empty results
          if(length(the_one) == 0){
            the_one <- lapply(stringr::str_split(seq, "-"), function(y){
              y <- as.numeric(y)
              x <- lapply(tab_sequ_n, function(x){
                x <- as.numeric(x)
                x[2] <= y[1]
              })
              x <- unlist(x);
              x
            })
            the_one <- lapply(the_one, which)
          }
          sequence_catego <- paste0(" ", paste(1:length(tab_sequ_n), collapse = " "), " ")
          for(l in 1:length(the_one)){
            to_insert <- paste0(" ", paste(the_one[[l]], collapse = " "), " ")
            to_add <- paste(sequence_catego[length(sequence_catego)], collapse = " ")
            if(l != 1){
              to_add <- paste0(" ", to_add)
            }
            to_add <- strsplit(to_add, to_insert)[[1]]
            sequence_catego <- sequence_catego[-length(sequence_catego)]
            sequence_catego <- append(sequence_catego, append(to_add, to_insert, after = 1))
          }

          sequence_catego <- stringr::str_remove_all(sequence_catego, "^ | $")
          sequence_catego <- lapply(stringr::str_split(sequence_catego, " "), as.numeric)
          sequence_catego <- sequence_catego[lapply(sequence_catego, length) > 0]
          sequence_catego <- sequence_catego[unlist(lapply(sequence_catego, function(x) all(!is.na(x))))]

          res <- list()
          for(k in 1:length(sequence_catego)){
            sumup_tab <- tab_pr[sequence_catego[[k]],]

            if(nrow(sumup_tab) > 1){
              annot <- sumup_tab %>% dplyr::select_if(~!is.numeric(.x))
              anno_seq <- as.numeric(unlist(tab_sequ_n[sequence_catego[[k]]]))
              annot$`Positions in Master Proteins`[1] <-  paste(annot$`Master Protein Accessions`[1],
                                                                paste0("[", min(anno_seq), "~", max(anno_seq), "]")
              )

              annot$`Annotated Sequence`[1] <- ""
              annot$Modifications[1] <- ""
              # we not gonna add modif informations, would be too much and only TMT modif
              annot <- annot[1,]

              num <- sumup_tab %>% dplyr::select(where(is.numeric)) %>%
                dplyr::summarise_all(sum, na.rm = TRUE)

              sumup_tab <- as.data.frame(cbind(annot, num))
            }

            res[[k]] <- sumup_tab
          }

          res <- as.data.frame(Reduce(rbind, res))
        }
      }
      else{
        annot <- tab_pr %>% dplyr::select_if(~!is.numeric(.x))
        anno_seq <- as.numeric(unlist(tab_sequ_n))
        annot$`Positions in Master Proteins`[1] <-  paste(annot$`Master Protein Accessions`[1],
                                                          paste0("[", min(anno_seq), "~", max(anno_seq), "]")
        )
        annot$`Annotated Sequence`[1] <- ""
        annot$Modifications[1] <- ""
        # we not gonna add modif informations, would be too much and only TMT modif
        annot <- annot[1,]

        num <- tab_pr %>% dplyr::select(where(is.numeric)) %>%
          dplyr::summarise_all(sum, na.rm = TRUE)

        res <- as.data.frame(cbind(annot, num))
      }

      final_res[[i]] <- res
    }

    final_res <- as.data.frame(Reduce(rbind, final_res))
  }
  else{
    final_res <- data
  }

  # getting back log2 fold change
  final_res[,6:ncol(final_res)] <-
    log2(final_res[,6:ncol(final_res)])

  # preparing data for cladiff function from IMPRINTS.CETSA and barplotting
  keep_n <- colnames(final_res)[1:5]
  colnames(final_res)[1:5] <- c("id", "description", "sumUniPeps", "sumPSMs", "countNum")
  final_res$id <- paste0(1:nrow(final_res), "_", final_res$id)

  message("Getting caldiff output") # will remove file saved by caldiff function
  diff_directory_toremove <- list.files()
  final_res_diff <- IMPRINTS.CETSA::imprints_caldiff(final_res, reftreatment = control)
  diff_directory_toremove <- list.files()[!(list.files() %in% diff_directory_toremove)]
  diff_directory_toremove <- stringr::str_subset(diff_directory_toremove, "final_res_\\d{6}_\\d{4}$")
  if(length(diff_directory_toremove) == 1){
    unlink(diff_directory_toremove, recursive = TRUE)
  }

  final_res_diff <- final_res_diff[order(as.numeric(unlist(
    lapply(stringr::str_split(final_res_diff$id, "_"),
           function(x) x[1]))
  )),]
  final_res_diff$id <- stringr::str_remove_all(final_res_diff$id, "^\\d{1,}_")

  final_res_diff$id <- paste(final_res_diff$sumUniPeps, final_res_diff$sumPSMs, final_res_diff$countNum, sep = " \n")

  # return data in same shape as input
  final_res_difffile <- final_res_diff[, c(1:2, (ncol(final_res_diff) - 2):ncol(final_res_diff),
                                           3:(ncol(final_res_diff) - 3))
                                       ]
  colnames(final_res_difffile)[1:5] <- keep_n
  final_res_difffile$`Master Protein Accessions` <- stringr::word(final_res_difffile$`Master Protein Accessions`, 1, 2)

  readr::write_tsv(final_res_difffile,
                   file = paste0(format(Sys.time(), "%y%m%d_%H%M_"),
                                 "CaldiffPeptides_", dataset_name, ".txt")
                   )


  message("Generates plot")
  imprints_barplotting_app(final_res_diff, ret_plot = FALSE, save_pdf = TRUE, layout = c(3,3),
                          pdfname = dataset_name)



  return(final_res_diff)
}


