#' imprints_normalize_peptides
#'
#' Function to normalize your peptides IMPINTS-CETSA dataset.
#'
#' @details The function uses the function \code{imprints_normalization} from \code{IMPRINTS.CETSA} package.
#'          The function apply the VSN model to each temperature set and then call the function \code{removeBatchEffect}
#'          from \code{limma} package.
#'
#' @param peptides_data The data frame output from \code{imprints_read_peptides}. Can also be the path to that file.
#'
#' @return The normalized peptides data frame, non log2 transformed.
#'
#' @export
#'
#'

imprints_normalize_peptides <- function(peptides_data){
  if(is.character(peptides_data)){
    if(file.exists(peptides_data)){
      peptides_data <- readr::read_tsv(peptides_data)
      dataset_name <- "peptides_data"
    }
    else{
      message("Error: Your peptide file wasn't found, check your spelling.")
      return()
    }
  }
  else if(inherits(peptides_data, "data.frame")){
    necessary_columns <- c("Master.Protein.Accessions", "description",
                           "Positions.in.Master.Proteins",
                           "Annotated.Sequence", "Modifications")
    necessary_columns <- paste(colnames(peptides_data)[1:5], collapse = "|") == paste(necessary_columns, collapse = "|")
    if(!necessary_columns){
      message("Error: Your peptides data should start with the columns 'Master.Protein.Accessions',
                      'description', 'Positions.in.Master.Proteins', 'Annotated.Sequence' and 'Modifications'.
                      You should use the output from 'imprints_read_peptides'.")
      return()
    }
    dataset_name <- deparse(substitute(peptides_data))
  }
  else{
    message("Error: Your peptides data can only be the output file or the output data frame from 'imprints_read_peptides'.")
    return()
  }

  # prepare data shape in order to use imprints_normalization from IMPRINTS.CETSA package
  countNum_in_peptides <- "countNum" %in% colnames(peptides_data)
  if(countNum_in_peptides){
    peptides_countNum <- peptides_data$countNum
    peptides_data$countNum <- NULL
  }

  colnames(peptides_data)[1:5] <- c("sumUniPeps", "description", "id", "sumPSMs", "countNum")
  peptides_data$sumPSMs <- paste(peptides_data$id,  peptides_data$sumPSMs, sep = "_sep_")
  peptides_data$id <- as.character(1:nrow(peptides_data))
  peptides_data <- peptides_data[order(peptides_data$id),]
  peptides_data <- peptides_data[,c(3,2, 6:ncol(peptides_data), 1,4,5)]

  # normalize
  norm_directory_toremove <- list.files()
  peptides_data <- IMPRINTS.CETSA::imprints_normalization(peptides_data)
  norm_directory_toremove <- list.files()[!(list.files() %in% norm_directory_toremove)]
  norm_directory_toremove <- grep("peptides_data_\\d{6}_\\d{4}$", norm_directory_toremove, value = TRUE)
  if(length(norm_directory_toremove) == 1){
    unlink(norm_directory_toremove, recursive = TRUE)
  }

  # putting back the data in original shape
  peptides_data$id <- unlist(lapply(strsplit(peptides_data$sumPSMs, "_sep_"), function(x) x[1]))
  peptides_data$sumPSMs <- unlist(lapply(strsplit(peptides_data$sumPSMs, "_sep_"), function(x) x[2]))
  colnames(peptides_data)[c(1:2, (ncol(peptides_data) -  2):ncol(peptides_data))] <- c("Positions.in.Master.Proteins",
                                                                                       "description",
                                                                                       "Master.Protein.Accessions",
                                                                                       "Annotated.Sequence",
                                                                                       "Modifications" )
  peptides_data <- peptides_data[,c(ncol(peptides_data) -  2,2,1,
                                (ncol(peptides_data) -  1):ncol(peptides_data),
                                3:(ncol(peptides_data) -  3))]
  peptides_data[,6:ncol(peptides_data)] <- 2**peptides_data[,6:ncol(peptides_data)]

  if(countNum_in_peptides){
    peptides_data$countNum <- peptides_countNum
  }

  readr::write_tsv(peptides_data,
                   file = paste0(format(Sys.time(), "%y%m%d_%H%M_"),
                                 "NormPeptides_", dataset_name, ".txt")
                   )

  return(peptides_data)
}

