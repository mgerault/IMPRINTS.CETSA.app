#' join_drugdata
#'
#' Function to join a list of data.frames. It is the same as the join_all function
#' from the plyr package but allow to use the full_join function from the dplyr package
#' and so to add suffixes to duplicated names not selecting by the argument by.
#'
#' @param dfs A list of data frames
#' @param by A character vector of variables to join by.
#'           If NULL, the default, full_join() will perform a natural join, using all variables in common across the data frames.
#'           A message lists the variables so that you can check they're correct; suppress the message by supplying
#'           by explicitly.
#'
#'           To join by different variables, use a named vector. For example, by = c("a" = "b") will match x$a to y$b.
#'           To join by multiple variables, use a vector with length > 1. For example, by = c("a", "b") will match
#'           x$a to y$a and x$b to y$b. Use a named vector to match different variables in x and y.
#'           For example, by = c("a" = "b", "c" = "d") will match x$a to y$b and x$c to y$d.
#'           x being the data frame named 'joined' and y the data frame i from dfs.
#'
#'
#' @return The joined dataframe
#'
#' @seealso \code{\link{dplyr}}
#'
#' @export
#'


join_drugdata <- function (dfs, by = NULL){
  if(length(dfs) == 0){
    message("dfs length is 0, NULL has been returned.")
    return(NULL)
  }
  else if (length(dfs) == 1)
    return(dfs[[1]])

  if("description" %in% by){
    is_OX <- unlist(lapply(dfs,
                           function(z){
                             z <- z$description[1]
                             z <- grepl("OX=\\d{1,}", z);
                             z
                           })
                    )
    if(length(unique(is_OX)) > 1){   # if in description OX is precised and for other not, remove it
      dfs[is_OX] <- lapply(dfs[is_OX],
                           function(z){
                             z$description <- gsub("OX=\\d{1,} ", "", z$description);
                             z
                           })
    }

    id_diff_descr <- Reduce(rbind,
                            sapply(dfs, "[", 1:2, simplify = FALSE))
    id_diff_descr$dataset <- unlist(mapply(rep, 1:length(dfs), sapply(dfs, nrow),
                                           SIMPLIFY = FALSE))
    id_diff_descr <- id_diff_descr %>%
      group_by(id, description) %>%
      summarise(dataset = paste(dataset, collapse = " & ")) %>%
      ungroup() %>% group_by(id) %>%
      filter(length(description) > 1)

    if(nrow(id_diff_descr)){
      message("Warning: Some datasets were searched with different FASTA files ! Will update description")
      # finding best and most recent description for each id with different description
      id_diff_descr <- id_diff_descr %>% group_modify(~{
        description_best <- .x$description

        if(all(grepl(" PE=", .x$description))){
          pr_existence <- as.numeric(sub(" .*", "", sub(".* PE=", "", .x$description)))
          # 1. Experimental evidence at protein level
          # 2. Experimental evidence at transcript level
          # 3. Protein inferred from homology
          # 4. Protein predicted
          # 5. Protein uncertain
        }
        else{
          pr_existence <- 1
        }

        if(length(unique(pr_existence)) > 1){
          # the lower, the more sure of the protein existence
          description_best <- .x$description[which.min(pr_existence)]
        }
        else{
          if(all(grepl(" SV=", .x$description))){
            seq_version <- as.numeric(sub(" .*", "", sub(".* SV=", "", .x$description)))
          }
          else{
            seq_version <- 100
          }

          if(length(unique(seq_version)) > 1){
            # the greater, the most updated the sequence is
            description_best <- .x$description[which.max(seq_version)]
          }
          else{
            pr_descr <- sub(" OS=.*", "", .x$description)
            if(any(grepl("probable", pr_descr, ignore.case = TRUE)) &
               !all(grepl("probable", pr_descr, ignore.case = TRUE))){
              # keeping description not containing 'Probable'
              description_best <- .x$description[-grep("probable", pr_descr, ignore.case = TRUE)]
            }
            else if(any(grepl("uncharacterized", pr_descr, ignore.case = TRUE)) &
                    !all(grepl("uncharacterized", pr_descr, ignore.case = TRUE))){
              # keeping description not containing 'uncharacterized'
              description_best <- .x$description[-grep("uncharacterized", pr_descr, ignore.case = TRUE)]
            }
            else{
              gene <- sub(" .*", "", sub(".* GN=", "", .x$description))
              if(any(grepl("orf", gene)) & !all(grepl("orf", gene))){
                # keeping defined Gene, like C12orf55
                description_best <- .x$description[-grep("orf", gene)]
              }
            }
          }
        }


        if(length(description_best) > 1){
          description_best <- NA
        }
        .x$description_best <- description_best
        return(.x)
      })

      # order of most recent :
      fasta_most_updated <-  id_diff_descr$dataset[which(id_diff_descr$description == id_diff_descr$description_best)]
      fasta_most_updated <- c(1:length(dfs))[order(table(unlist(strsplit(fasta_most_updated,
                                                                               " & "))
      ), decreasing = TRUE)]

      id_diff_descr <- id_diff_descr %>%
        group_modify(~ {
          if(all(is.na(.x$description_best))){
            datasets_age <- sapply(strsplit(.x$dataset, " & "),
                                   function(y){
                                     y <- as.numeric(y)
                                     ages <- which(!is.na(match(fasta_most_updated, y)))
                                     min(ages)
                                   })

            .x$description_best <- .x$description[which.min(datasets_age)]
          }

          return(.x)
        }) %>%
        dplyr::filter(description != description_best) %>%
        dplyr::select(-description) %>%
        ungroup() %>% group_by(id, description_best) %>%
        dplyr::reframe(dataset = as.numeric(unlist(strsplit(dataset, " & "))))

      for(d in unique(id_diff_descr$dataset)){
        dfs[[d]] <- left_join(dfs[[d]], id_diff_descr[which(id_diff_descr$dataset == d),1:2], by = "id")
        dfs[[d]]$description[which(!is.na(dfs[[d]]$description_best))] <-
          dfs[[d]]$description_best[which(!is.na(dfs[[d]]$description_best))]

        dfs[[d]]$description_best <- NULL
      }
    }
  }

  joined <- dfs[[1]]
  for (i in 2:length(dfs)) {
    joined <- dplyr::full_join(joined, dfs[[i]], by = by)
  }
  return(joined)
}

