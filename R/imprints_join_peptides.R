#' imprints_join_peptides
#'
#' Function to join peptides dataset obtained after \code{imprints_sequence_peptides}.
#'
#' @param data A list of peptides data to join, typically after \code{imprints_sequence_peptides}
#'
#' @return Your filtered dataset
#'
#' @export
#'
#'

imprints_join_peptides <- function(data){
  data <- lapply(data, function(x){
    x$`Master Protein Accessions` <- unlist(lapply(stringr::str_split(x$`Master Protein Accessions`, " "), function(x) x[1]))
    x$`Annotated Sequence` <- ""
    x$Modifications <- ""
    x$description <- stringr::str_remove_all(x$description, " OX=\\d{1,}") # depending on the PD version, some dataset doesn't have this information --> remove it
    x <- x[order(x$`Positions in Master Proteins`),]
  })

  # return rank from a series of sequences
  get_sequence_order <- function(x){
    x <- stringr::str_split(x, "-|~")
    x <- lapply(x, function(x) as.numeric(x[2]))
    x <- unlist(x)
    x <- order(order(x)) # return rank from all sequence
    return(x)
  }

  # get position in master protein for each dataset and return as a dataframe
  data_protein_position <- lapply(data,
                                  function(x){
                                    x <- x$`Positions in Master Proteins`
                                    x_protein <- stringr::str_remove_all(x, "\\s.*")
                                    x_sequence <- unlist(stringr::str_extract_all(x, "(?<=\\[).+?(?=\\])"))
                                    x <- data.frame("protein" = x_protein,
                                                    "sequence" = x_sequence)
                                    x <- x %>% dplyr::group_by(protein) %>%
                                      dplyr::mutate(sequence_order = get_sequence_order(sequence)) %>%
                                      dplyr::ungroup();
                                    x
                                  })
  # rename 'sequence' column to differentiate dataset
  data_protein_position <- mapply(function(x,y){colnames(x)[2] <- paste0("sequence_", y);
                                                x},
                                  data_protein_position,
                                  1:length(data_protein_position),
                                  SIMPLIFY = FALSE
                                  )

  # join position in master protein from all data frames in one
  data_protein_position_number <- Reduce(function(x,y){dplyr::full_join(x,y, by = c("protein", "sequence_order"))},
                                         data_protein_position)

  # group sequence according there proximity
  group_sequence <- function(x, n_group){
    sequence <- lapply(stringr::str_split(x, "-|~"), as.numeric)
    sequence <- sequence[which(!is.na(sequence))]
    sequence <- Reduce(rbind, sequence)

    # get closest sequence
    sequence_dist <- dist(sequence)
    sequence_clust <- stats::hclust(sequence_dist)
    sequence_clust <- stats::cutree(sequence_clust, k = n_group)
    sequence_clust <- as.numeric(sequence_clust)

    if(length(sequence_clust) != length(x)){ # means that there are NAs --> need to assign a group to those
      x[which(!is.na(x))] <- sequence_clust
      sequence_clust <- as.numeric(x)
    }
    return(sequence_clust)
  }
  # change NA in missing group for sequence group
  group_sequence_changeNA <- function(group_sequenceNA, group_sequence){
    idx_NA <- which(is.na(group_sequenceNA))

    if(length(idx_NA)){
      group_toadd <- group_sequence[!(group_sequence %in% group_sequenceNA[-idx_NA])]
      group_sequenceNA[idx_NA] <- group_toadd
    }

    return(group_sequenceNA)
  }

  # group sequence for each dataset
  data_protein_position_number <- data_protein_position_number %>%
    tidyr::gather("sequence_dataset", "sequence", -protein, -sequence_order) %>%
    dplyr::group_by(protein) %>%
    dplyr::mutate(sequence_group = group_sequence(sequence, max(sequence_order))) %>%
    dplyr::group_by(protein, sequence_dataset) %>%
    dplyr::mutate(sequence_group = group_sequence_changeNA(sequence_group, sequence_order))
  data_protein_position_number$sequence_order <- NULL

  # function to get global sequence among datasets
  get_global_sequence <- function(x){
    sequence <- stringr::str_split(x, "-|~")
    sequence <- lapply(sequence, as.numeric)
    sequence <- unlist(sequence)

    global_sequence <- paste0(min(sequence, na.rm = TRUE), "~", max(sequence, na.rm = TRUE))
    return(global_sequence)
  }

  # extract global sequence from all dataset
  data_protein_position_number <- data_protein_position_number %>%
    dplyr::ungroup() %>%
    dplyr::group_by(protein, sequence_group) %>%
    dplyr::mutate(global_sequence = get_global_sequence(sequence))

  data_protein_position_number$global_sequence[is.na(data_protein_position_number$sequence)] <- NA # doesn't have the sequence --> replace by NA
  data_protein_position_number$sequence <- NULL
  data_protein_position_number$sequence_group <- NULL

  # get back the global sequence for each dataset --> return a list
  data_protein_position_number <- (data_protein_position_number %>% dplyr::filter(!is.na(global_sequence)) %>%
                                     dplyr::group_by(sequence_dataset) %>%
                                     dplyr::summarise(result = paste0(protein, " [", global_sequence, "]")) %>%
                                     dplyr::summarise(result = list(result[order(result)])) %>%
                                     dplyr::select(result))[[1]]

  # since we ordered postion, already matched
  # can now assign to corresponding dataset
  data <- mapply(function(x,y){x$`Positions in Master Proteins` <- y;
                               x
                              },
                 data,
                 data_protein_position_number,
                 SIMPLIFY = FALSE)

  # join datasets
  data <- Reduce(function(x,y) dplyr::full_join(x,y, by = c("Master Protein Accessions",
                                                            "description",
                                                            "Positions in Master Proteins",
                                                            "Annotated Sequence",
                                                            "Modifications")),
                 data)

  # ordering final dataset -> take into account sequence number
  data$factor <- data$`Master Protein Accessions`
  data <- data[order(data$`Master Protein Accessions`),]
  ord <- lapply(stringr::str_split(data$`Positions in Master Proteins`, "\\[|\\]"),
                function(x) {z = stringr::str_split(x[2], "~|-");
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
  data$factor <- NULL

  return(data)
}

