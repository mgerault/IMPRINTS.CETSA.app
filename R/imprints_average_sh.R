#' imprints_average_sh
#'
#' Function to calculate the averaged signals from the IMPRTINTS-CETSA result
#' It is totally based on the function ms_2D_average from the mineCETSA package.
#'
#' @param data dataset after calculating the relative protein abundance differences
#' @param savefile A logical to tell if you want save the results or not
#'
#' @return A dataframe
#'
#' @seealso \code{\link{ms_2D_caldiff}}
#'
#' @export
#'

imprints_average_sh <- function (data, savefile = TRUE)
{
  filename <- paste0(deparse(substitute(data)), "_average", ".txt")

  if (length(grep("description", names(data)))) {
    proteininfo <- unique(data[, c("id", "description")])
    data$description <- NULL
  }
  if (length(grep("countNum", names(data)))) {
    countinfo <- unique(data[,str_which(names(data), "^id$|^sumPSM|^countNum|^sumUniPeps")]) #allows to work with joined table
    data <- data[,-str_which(names(data), "^sumPSM|^countNum|^sumUniPeps")]
  }
  data <- tidyr::gather(data, condition, reading, -id)
  if (length(unlist(strsplit(data$condition[1], "_"))) ==
      4) {
    data <- tidyr::separate(data, condition, into = c("set",
                                                      "temperature", "replicate", "treatment"), sep = "_")

    data <- data %>% dplyr::group_by(id, set, temperature, treatment) %>%
      dplyr::summarise(reading.mean = mean(reading, na.rm = T))
    data <- tidyr::unite(data, condition, set, temperature,
                         treatment, sep = "_")
    data <- tidyr::spread(data, condition, reading.mean)
  }
  else if (length(unlist(strsplit(data$condition[1], "_"))) ==
           3) {
    data <- tidyr::separate(data, condition, into = c("temperature",
                                                      "replicate", "treatment"), sep = "_")
    data <- data %>% dplyr::group_by(id, temperature, treatment) %>%
      dplyr::summarise(reading.mean = mean(reading, na.rm = T))
    data <- tidyr::unite(data, condition, temperature, treatment,
                         sep = "_")
    data <- tidyr::spread(data, condition, reading.mean)
  }
  else {
    stop("make sure the namings of the columns of the dasaset are correct.")
  }
  data <- merge(data, countinfo)
  data <- merge(proteininfo, data)

  if(savefile){
    ms_filewrite(data, filename)
  }
  return(data)
}
