#' imprints_cleaved_peptides
#'
#' Function to find proteins that are potentially cleaved and their cleaved sites.
#' For more information, see Details section.
#'
#' @details
#' Let's say you found some hits in your dataset and you want to check if some of those
#' are cleaved. You first compute the fold changes for every peptides from each hit and then
#' search for potential cleaved sites.
#' The idea here is to compute the sum from all fold change for each peptide and then compute the cumulative sum of
#' these summed fold change of every peptides for each protein.
#' If the protein not cleaved, this cumulative sum should be constantly increasing/decreasing, i.e linear.
#' If the protein is cleaved, then this cumulative sum will either be convex or concave.
#' If convex, the last peptides are more abundant; if concave the first peptides are more abundant.
#'
#' The potential cleaved sites found, can only be peptides that are in the dataset.
#' Also, it's better to double check by plotting the peptides before and after the potential cleaved site found
#' by using \code{imprints_sequence_peptides}.
#'
#' @param data The input data set, i.e. the outpout from \code{imprints_sequence_peptides}.
#' @param R2 The R-squared cutoff. It characterize the linearity of the cumulative sum of the fold changes
#'           from every peptide from each protein. The higher it is, the less stringent you are on the decision
#'           of taking a protein as cleaved. Default is 0.9.
#' @param control The control treatment from your dataset. If the treatment is found in your data, then it is removed from it.
#' @param min_ValidValue The minimum proportion of valid values per peptides.
#'                       Default is 0.4; so if 7 temperatures need at least 3 valid values.
#'
#' @return The potential cleaved sites from the proteins considered as cleaved.
#'
#' @export
#'

imprints_cleaved_peptides <- function(data, R2 = 0.9, control = NULL, min_ValidValue = 0.4){
  if(!is.null(control)){
    if(any(grepl(paste0("_", control, "$"), colnames(data)))){
      data <- data[,-grep(paste0("_", control, "$"), colnames(data))]
    }
    else{
      message(paste("Warnings:", control, "hasn't been found in your data."))
    }
  }

  # data is output from 'imprints_sequence_peptides' --> remove unnecessary informations
  data$`Positions in Master Proteins` <- NULL
  data$`Annotated Sequence` <- NULL
  data$Modifications <- NULL
  colnames(data)[1] <- "id" # easier handle


  # get summed value from all temperature for every peptide
  data <- data %>%
    tidyr::gather("treatment", "value", -id, -description) %>%
    tidyr::separate(treatment, into = c("temp", "biorep", "treatment"), sep = "_") %>%
    dplyr::group_by(id, description, temp, treatment) %>%
    dplyr::reframe(mean_value = mean(value, na.rm = TRUE)) %>%
    dplyr::ungroup() %>% dplyr::group_by(id, description, treatment) %>%
    dplyr::filter(length(na.omit(mean_value))/length(mean_value) >= min_ValidValue) %>%  # keeping peptides with more than 40% of valid values
    dplyr::reframe(sum_profile = sum(mean_value, na.rm = TRUE))  # get sum value for each peptides --> sum all temperatures values


  # separate id in protein and sequence
  data$id <- stringr::word(data$id, 1, 2)
  data <- data %>%
    tidyr::separate(id, into = c("protein", "sequence"), sep = " ") %>%
    dplyr::mutate(sequence = gsub("\\[|\\]", "", sequence)) %>%
    dplyr::ungroup() %>% dplyr::group_by(protein, description, treatment) %>%
    dplyr::filter(length(sum_profile) > 3)  # only keeping proteins with more than 3 peptides


  # order data according proteins and sequence
  # ordering according sequence is capital
  data$factor <- data$protein
  data <- data[order(data$protein),]
  data$sequence <- gsub(";", "", data$sequence)
  ord <- lapply(strsplit(data$sequence, "-|~"),
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

  # compute cumulative sum for each protein
  data <- data %>% dplyr::group_by(protein, description, treatment) %>%
    dplyr::mutate(cumsum_profile = cumsum(tidyr::replace_na(sum_profile, 0))
           )

  # get cleaved site
  # if cumulative abundance concave --> first derivative is decreasing globally --> inflexion point is where second derivative is minimum
  # if cumulative abundance convex --> first derivative is increasing globally --> inflexion point is where second derivative is maximum
  inflex <- function(x, R2_cut = R2){
    R2_calc <- cor(1:length(x), x)  # if linear relationship --> ~ constantly increasing/decreasing --> no cleaved
    R2_calc <- R2_calc**2
    if(R2_calc <= R2_cut){
      x <- diff(x) # get first 'first derivative'
      # compute 'tendency' from 'first derivative', i.e. if increasing/decreasing
      df <- data.frame(x = 1:length(x),
                       y = x)
      res <- lm(y ~ x, data = df)
      res <- res$coeff[2]
      res <- sign(res)

      x <- diff(x) # get second 'first derivative'
      if(res == 1){
        inflex_point <- which.max(x) + 2 # double diff, so add 2 to index point
      }
      else{
        inflex_point <- which.min(x) + 2 # double diff, so add 2 to index point
      }

      if(inflex_point - 2 == length(x)){  # if inflex_point is last one, can't conclude that is cleaved
        inflex_point <- NaN
      }

      return(inflex_point)
    }
    else{
      return(NaN)
    }
  }

  data <- data %>% dplyr::ungroup() %>% dplyr::group_by(protein, description, treatment) %>%
    dplyr::reframe(cleaved_site = sequence[inflex(cumsum_profile)]) %>%  # get back potential cleaved site
    dplyr::filter(!is.na(cleaved_site))  # if NA --> not cleaved

  readr::write_tsv(data,
                   file = paste0(format(Sys.time(), "%y%m%d_%H%M_"),
                                 "Peptides_PotentialCleaved.txt")
                   )

  return(data)
}

