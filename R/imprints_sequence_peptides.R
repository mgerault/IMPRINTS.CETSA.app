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
#' @param barplot Logical to tell if you want to plot the barplots from the obtained fold-changes. Default is FALSE
#' @param dataset_name The name of your dataset to save your file.
#'
#' @return The caldiff output from the peptide position from the protein you selected; save it and also save the barplots.
#'
#' @export
#'
#'

imprints_sequence_peptides <- function(data, proteins = NULL, sequence = NULL,
                                       control = "Vehicle", barplot = FALSE,
                                       dataset_name = "imprints"){
  if(!is.null(sequence)){
    if(!(length(proteins) == length(sequence) | length(sequence) == 1) & !inherits(sequence, "list")){
      message("Error: The number of sequence needs to match the number of proteins or just provide one sequence")
      return()
    }
  }

  if(!is.null(proteins)){
    protein_in <- proteins %in% data$Master.Protein.Accessions

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

    data <- data[which(!is.na(match(data$Master.Protein.Accessions,
                                    proteins))
                       ),]
  }
  else{
    proteins <- data$Master.Protein.Accessions
  }

  # reorder the data --> order according proteins and peptide sequence (~position)
  data$factor <- data$Master.Protein.Accessions
  data$Master.Protein.Accessions <- paste(data$Positions.in.Master.Proteins,
                                          data$Annotated.Sequence,
                                          data$Modifications, sep = " \n")
  data <- data[order(data$Master.Protein.Accessions),]
  ord <- gsub(".*\\[|\\]", "", gsub(";.*", "", data$Positions.in.Master.Proteins))
  ord <- sapply(strsplit(ord, "-|~"),
                function(x) {sum(as.numeric(x))})
  ord <- data.frame(factor = data$factor, idx = 1:nrow(data), sum.pos = ord) %>%
    dplyr::group_by(factor) %>%
    dplyr::mutate(ord.sum.pos = order(sum.pos),
                  nb.pep = length(ord.sum.pos),
                  final.order = idx[ord.sum.pos])

  data <- data[ord$final.order,]
  data$Master.Protein.Accessions <- data$factor
  data$factor <- NULL

  cond_data <- unique(unlist(lapply(strsplit(grep("^\\d{1,}", colnames(data), value = TRUE), "_"), "[[", 3)))
  cond_data <- cond_data[!(cond_data %in% control)]
  temp_data <- unique(unlist(lapply(strsplit(grep("^\\d{1,}", colnames(data), value = TRUE), "_"), "[[", 1)))

  final_res <- list()
  if(!is.null(sequence)){
    message("Filtering and summing sequences...")
    for (i in 1:length(proteins)){
      message(paste0(proteins[i], " : ", i, "/", length(proteins)))

      tab_pr <- data[which(data$Master.Protein.Accessions == proteins[i]),]
      tab_sequ = gsub(";.*", "", tab_pr$Positions.in.Master.Proteins) # if protein group, only take first one for simplicity
      tab_sequ = gsub(".* \\[|\\].*", "", tab_sequ)
      tab_sequ_n = strsplit(tab_sequ, "-")

      if(length(sequence) == 1){
        seq = sequence[[1]]
      }
      else{
        seq = sequence[[i]]
      }

      if(length(seq) > 0){
        if(nchar(seq)[1] == 0){
          res <- tab_pr
        }
        else{
          seq <- seq[order(unlist(lapply(strsplit(seq, "-"), function(x) as.numeric(x[2]))))]
          the_one <- lapply(strsplit(seq, "-"), function(y){
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
            the_one <- lapply(strsplit(seq, "-"), function(y){
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

          sequence_catego <- gsub("^ | $", "", sequence_catego)
          sequence_catego <- lapply(strsplit(sequence_catego, " "), as.numeric)
          sequence_catego <- sequence_catego[lapply(sequence_catego, length) > 0]
          sequence_catego <- sequence_catego[unlist(lapply(sequence_catego, function(x) all(!is.na(x))))]

          res <- list()
          for(k in 1:length(sequence_catego)){
            sumup_tab <- tab_pr[sequence_catego[[k]],]

            if(nrow(sumup_tab) > 1){
              if("countNum" %in% colnames(sumup_tab)){
                sumup_tab$countNum <- as.character(sumup_tab$countNum)
              }
              annot <- sumup_tab %>% dplyr::select_if(~!is.numeric(.x))
              anno_seq <- as.numeric(unlist(tab_sequ_n[sequence_catego[[k]]]))
              annot$Positions.in.Master.Proteins[1] <-  paste0(sub(";.*", "", annot$Master.Protein.Accessions[1]),
                                                               paste0(" [", min(anno_seq), "~", max(anno_seq), "]"),
                                                               ifelse(grepl(";", annot$Master.Protein.Accessions[1]),
                                                                      sub("^(.*?;)", ";", annot$Master.Protein.Accessions[1]),
                                                                      "")
                                                               )


              annot$Annotated.Sequence[1] <- ""
              annot$Modifications[1] <- ""
              if("countNum" %in% colnames(annot)){
                annot$countNum <- as.numeric(annot$countNum)
                annot$countNum[1] <- median(annot$countNum, na.rm = TRUE)
              }
              # we not gonna add modif informations, would be too much
              annot <- annot[1,]

              num <- sumup_tab %>% dplyr::select(where(is.numeric)) %>%
                dplyr::summarise_all(sum, na.rm = TRUE)

              sumup_tab <- as.data.frame(cbind(annot, num))
              if("countNum" %in% colnames(sumup_tab)){
                cn <- colnames(sumup_tab)[-which(colnames(sumup_tab) == "countNum")]
                sumup_tab <- sumup_tab[,c(cn, "countNum")]
              }
            }

            res[[k]] <- sumup_tab
          }

          res <- as.data.frame(Reduce(rbind, res))
        }
      }
      else{
        if("countNum" %in% colnames(tab_pr)){
          tab_pr$countNum <- as.character(tab_pr$countNum)
        }
        annot <- tab_pr %>% dplyr::select_if(~!is.numeric(.x))
        anno_seq <- as.numeric(unlist(tab_sequ_n))
        annot$Positions.in.Master.Proteins[1] <-  paste0(sub(";.*", "", annot$Master.Protein.Accessions[1]),
                                                         paste0(" [", min(anno_seq), "~", max(anno_seq), "]"),
                                                         ifelse(grepl(";", annot$Master.Protein.Accessions[1]),
                                                                sub("^(.*?;)", ";", annot$Master.Protein.Accessions[1]),
                                                                "")
                                                         )

        annot$Annotated.Sequence[1] <- ""
        annot$Modifications[1] <- ""
        # we not gonna add modif informations, would be too much and only TMT modif
        if("countNum" %in% colnames(annot)){
          annot$countNum <- as.numeric(annot$countNum)
          annot$countNum[1] <- median(annot$countNum, na.rm = TRUE)
        }
        annot <- annot[1,]

        num <- tab_pr %>% dplyr::select(where(is.numeric)) %>%
          dplyr::summarise_all(sum, na.rm = TRUE)

        res <- as.data.frame(cbind(annot, num))
        if("countNum" %in% colnames(res)){
          cn <- colnames(res)[-which(colnames(res) == "countNum")]
          res <- res[,c(cn, "countNum")]
        }
      }

      final_res[[i]] <- res
    }

    final_res <- as.data.frame(Reduce(rbind, final_res))
  }
  else{
    final_res <- data
  }

  countNum_in_peptides <- "countNum" %in% colnames(final_res)
  if(countNum_in_peptides){
    peptides_countNum <- final_res$countNum
    final_res$countNum <- NULL
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

  # checking if within rep can't be used
  nb_rep_ptreat <- unique(gsub("^\\d{2}C_", "", colnames(final_res)[grep("^\\d{2}C", colnames(final_res))]))
  nb_rep_ptreat <- table(gsub(".*_", "", nb_rep_ptreat))
  withinrep <- length(unique(nb_rep_ptreat)) == 1
  if(!withinrep){
    message("Warning: The treatments in your dataset doesn't have the same number of replicates ! The fold-changes weren't calculated within each replicate.")
  }

  final_res_diff <- IMPRINTS.CETSA::imprints_caldiff_f(final_res, reftreatment = control, withinrep = withinrep)
  diff_directory_toremove <- list.files()[!(list.files() %in% diff_directory_toremove)]
  diff_directory_toremove <- grep("final_res_\\d{6}_\\d{4}$", diff_directory_toremove, value = TRUE)
  if(length(diff_directory_toremove) == 1){
    unlink(diff_directory_toremove, recursive = TRUE)
  }

  final_res_diff <- final_res_diff[order(as.numeric(unlist(
                                          lapply(strsplit(final_res_diff$id, "_"),
                                                 function(x) x[1])
                                          )
                                        )),]
  final_res_diff$id <- gsub("^\\d{1,}_", "", final_res_diff$id)

  # return data in same shape as input
  final_res_difffile <- final_res_diff[, c(1:2, (ncol(final_res_diff) - 2):ncol(final_res_diff),
                                           3:(ncol(final_res_diff) - 3))
                                       ]
  colnames(final_res_difffile)[1:5] <- keep_n

  if(countNum_in_peptides){
    final_res_difffile$countNum <- peptides_countNum
  }

  readr::write_tsv(final_res_difffile,
                   file = paste0(format(Sys.time(), "%y%m%d_%H%M_"),
                                 "CaldiffPeptides_", dataset_name, ".txt")
                   )


  if(barplot){
    message("Generates plot")
    imprints_barplotting_peptides(final_res_difffile,
                                  ret_plot = FALSE, save_pdf = TRUE,
                                  layout = c(3,3), pdfname = dataset_name)
  }

  message("Done !")
  return(final_res_difffile)
}


