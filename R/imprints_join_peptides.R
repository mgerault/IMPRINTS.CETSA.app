#' imprints_join_peptides
#'
#' Function to join peptides dataset obtained after \code{imprints_sequence_peptides}.
#'
#' @param data A list of peptides data to join, typically after \code{imprints_sequence_peptides}
#' @param mode Either 'exact' or 'partial'. If 'exact', the datasets will be joined according their
#'   peptide sequence, i.e. to be joined they need to have the same annotated sequence, position and modification.
#'   If 'partial', peptides will only be joined on the protein position and will try to join close
#'   position together; i.e. it will not take into account the amino-acid sequence. It is mostly useful
#'   for joining dataset after summing different peptides parts together using the function \code{imprints_sequence_peptides}.
#'   Default is 'partial'.
#'
#' @return Your joined dataset
#'
#' @export
#'

imprints_join_peptides <- function(data, mode = c("partial", "exact")){
  mode <- match.arg(mode)

  data <- lapply(data, function(x){
    if(mode == "partial"){
      x$Annotated.Sequence <- ""
      x$Modifications <- ""
    }

    x$description <- gsub(" OX=\\d{1,}", "", x$description) # depending on the PD version, some dataset doesn't have this information --> remove it

    # ordering final dataset -> take into account sequence number
    x$factor <- x$Master.Protein.Accessions
    x <- x[order(x$Master.Protein.Accessions),]
    ord <- gsub(".*\\[|\\]", "",
                sub(";.*", "", x$Positions.in.Master.Proteins)
                )
    ord <- lapply(strsplit(ord, "-|~"),
                  function(x) {sum(as.numeric(x))}
    )
    ord <- unlist(ord)
    ord_f <- data.frame(factor = x$factor, x = ord)
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

    x <- x[ord_f$y,]
    x$factor <- NULL

    if("countNum" %in% colnames(x)){
      message("Removing countNum information from data for easier handle")
      x$countNum <- NULL
    };
    x
  })

  if(mode == "partial"){
    # return rank from a series of sequences
    get_sequence_order <- function(x){
      x <- strsplit(x, "-|~")
      x <- lapply(x, function(x) sum(as.numeric(x)))
      x <- unlist(x)
      x <- order(order(x)) # return rank from all sequence
      return(x)
    }

    # get position in master protein for each dataset and return as a dataframe
    data_protein_position <- lapply(data,
                                    function(x){
                                      x <- x$Positions.in.Master.Proteins
                                      x_protein <- sub(" .*", "", x)
                                      x_sequence <- unlist(stringr::str_extract(x, "(?<=\\[).+?(?=\\])"))
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
      sequence <- lapply(strsplit(x, "-|~"), as.numeric)
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
      sequence <- strsplit(x, "-|~")
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
    data_protein_position_number <- data_protein_position_number %>% dplyr::filter(!is.na(global_sequence)) %>%
      dplyr::group_by(sequence_dataset, protein) %>%
      mutate(global_sequence = global_sequence[get_sequence_order(global_sequence)])
    data_protein_position_number <- sapply(unique(data_protein_position_number$sequence_dataset),
                                           function(x){
                                             x <- data_protein_position_number[data_protein_position_number$sequence_dataset == x,]
                                             x$sequence_dataset <- NULL

                                             x$factor <- x$protein
                                             x <- x[order(x$protein),]
                                             ord <- x$global_sequence
                                             ord <- lapply(strsplit(ord, "-|~"),
                                                           function(x) {sum(as.numeric(x))}
                                             )
                                             ord <- unlist(ord)
                                             ord_f <- data.frame(factor = x$factor, x = ord)
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

                                             x <- x[ord_f$y,]
                                             x$factor <- NULL

                                             x <- paste0(x$protein, " [",  x$global_sequence, "]")
                                             return(x)
                                           }, simplify = FALSE)

    # since we ordered postion, already matched
    # can now assign to corresponding dataset
    data <- mapply(function(x,y){x$Positions.in.Master.Proteins <- y;
                                 x},
                   data,
                   data_protein_position_number,
                   SIMPLIFY = FALSE)
  }

  # join datasets
  data <- Reduce(function(x,y) dplyr::full_join(x,y, by = c("Master.Protein.Accessions",
                                                            "description",
                                                            "Positions.in.Master.Proteins",
                                                            "Annotated.Sequence",
                                                            "Modifications")),
                 data)

  # ordering final dataset -> take into account sequence number
  data$factor <- data$Master.Protein.Accessions
  data <- data[order(data$Master.Protein.Accessions),]
  ord <- gsub(".*\\[|\\]", "",
              sub(";.*", "", data$Positions.in.Master.Proteins)
              )
  ord <- lapply(strsplit(ord, "-|~"),
                function(x) {sum(as.numeric(x))}
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

