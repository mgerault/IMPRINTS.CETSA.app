#' join_cetsa
#'
#' Function to join a list of data frame from IMPRINTS-CETSA experiment and renamed the column.
#' (Allow to avoid the adding of ".x" or ".y" at the end of the column names)
#'
#' @param list_data A list of data frames. Usually the data frames are the output from the imprints_caldiff function.
#' @param new_names The new character element you want to add at the end of the columns name.
#'
#' @return The joined data frame (by 'id' and 'description')
#'
#' @export
#'

join_cetsa <- function(list_data, new_names = c("1h", "6h")){
  if(length(grep("_|/", new_names))){
    stop("Please enter a valid suffix. The character '_' and '/' are not allowed.")
  }
  if(inherits(list_data, "data.frame")){
    list_data <- list(list_data)
  }
  n <- length(list_data)

  if(length(new_names) != n){
    stop("You need to provide as many new_names as you have data.frame in your list !")
  }

  list_data <- mapply(function(x,y){ colnames(x)[!(colnames(x) %in% c("id", "description"))] <-
                                       paste0(colnames(x)[!(colnames(x) %in% c("id", "description"))], y); x},
                      list_data, new_names, SIMPLIFY = FALSE)

  is_OX <- unlist(lapply(list_data,
                         function(z){
                           z <- z$description
                           z <- grepl("OX=\\d{1,}", z);
                           any(z)
                         })
                  )
  if(any(is_OX)){ # if in description OX is precised and for other not, remove it
    list_data[is_OX] <- lapply(list_data[is_OX],
                                 function(z){
                                   z$description <- gsub("OX=\\d{1,} ", "", z$description);
                                   z
                              })
  }

  id_diff_descr <- Reduce(rbind,
                          sapply(list_data, "[", 1:2, simplify = FALSE))
  id_diff_descr$dataset <- unlist(mapply(rep, 1:length(list_data), sapply(list_data, nrow),
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
    fasta_most_updated <- c(1:length(list_data))[order(table(unlist(strsplit(fasta_most_updated,
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
      list_data[[d]] <- left_join(list_data[[d]], id_diff_descr[which(id_diff_descr$dataset == d),1:2], by = "id")
      list_data[[d]]$description[which(!is.na(list_data[[d]]$description_best))] <-
        list_data[[d]]$description_best[which(!is.na(list_data[[d]]$description_best))]

      list_data[[d]]$description_best <- NULL
    }
  }

  df <- plyr::join_all(list_data, by = c("id", "description"), type = "full")

  return(df)
}


