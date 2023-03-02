#' imprints_corr_to_ref_sh
#'
#' Function to generate 2D bar plot and pdf file with multipanel bar plots for 2D-CETSA data of
#' proteins which have similar profile from a selected protein.
#' It is totally based on the function ms_2D_barplotting from the mineCETSA package.
#'
#' @param data Dataset after calculating the relative protein abundance differences and deriving the average of the profile,
#'             ie imprints_average_sh, make sure the columns with readings are named in the format like "(Set_)37C_Treatment"
#' @param set A single character to indicate the sample name, if any
#' @param treatment A single character to indicate the sample name
#' @param reference A numeric vector with the profile readings
#' @param use_score A single character element that define the method score. Method available : 'euclidean' or 'pearson'
#' @param score_threshold A numeric value to indicate the threshold, default set to 0.9
#' @param include_neg Whether to include the negatively correlated proteins, default set to FALSE
#' @param max_na An integer indicating the maximun number of missing value for one protein, default is 0.
#'
#' @return A dataframe with the correlation from the profile and the proteins which passed the threshold
#'
#'
#' @export
#'

imprints_corr_to_ref_sh <- function (data = NULL, set = NULL, treatment = NULL, reference = NULL,
                                  use_score = c("euclidean", "pearson"),
                                  score_threshold = 0.9, include_neg = FALSE, max_na = 0)
{
  if (!length(set) %in% c(0, 1)) {
    stop("Pls provide one set name if any")
  }
  if (length(treatment) != 1) {
    stop("Pls provide one treatment name")
  }
  if (length(reference) == 0) {
    stop("Pls provide a numeric reading vector as the reference profile")
  }
  if (inherits(data, "data.frame")) {
    subset <- grep(paste0("_", treatment, "$"), names(data))
    if (length(set)) {
      subset_set <- grep(set, names(data))
      subset <- intersect(subset, subset_set)
    }
    if (length(subset) > 0 & length(subset) == length(reference)) {
      data1 = data[, c(1, 2, subset)]
    }
    else {
      stop("Pls provide the right set/treatment keyword character")
    }
  }
  if(!(use_score %in% c("euclidean", "pearson"))){
    stop("Please, choose a valid score method. For now, only 'euclidean' and 'pearson' are available.")
  }
  if (length(grep("countNum", names(data1)))) {
    countinfo1 <- unique(data1[,str_which(names(data1), "^id$|^sumPSM|^countNum|^sumUniPeps")])
    data1 <- data1[, -str_which(names(data1), "^sumPSM|^countNum|^sumUniPeps")]
  }
  if (length(grep("description", names(data1)))) {
    proteininfo <- unique(data1[, c("id", "description")])
    data1$description <- NULL
  }
  nb_na <- apply(data1, 1, function(x) sum(is.na(x)))
  if(max_na < 0){
    stop("Please enter a positive value for the maximum NA per row")
  }
  na_thresh <- which(nb_na <= max_na)
  data1 <- data1[na_thresh,]

  datal <- tidyr::gather(data1, treatment, reading, -id)
  if (length(unlist(strsplit(datal$treatment[1], "_"))) == 3) {
    datal <- tidyr::separate(datal, treatment, into = c("set",
                                                 "temperature", "condition"))
    datal <- tidyr::unite(datal, "condition", c("set", "condition"))
  }
  else if (length(unlist(strsplit(datal$treatment[1], "_"))) == 2) {
    datal <- tidyr::separate(datal, treatment, into = c("temperature",
                                                 "condition"))
  }

  treatment = unique(datal$condition)
  datal$condition <- NULL

  dataw <- tidyr::spread(datal, id, reading)
  dataw <- t(dataw)
  colnames(dataw) <- dataw[1,]
  dataw <- dataw[-1,]
  rname <- rownames(dataw)
  dataw <- apply(dataw, 2, as.numeric)
  rownames(dataw) <- rname
  dataw <- as.data.frame(dataw)


  if(use_score == "euclidean"){
    dataw$distance <- apply(dataw, 1, function(x) dist(rbind(x, reference))[1])
    dataw$score <- 1/(1 + dataw$distance)
    dataw$id <- rownames(dataw)
    rownames(dataw) <- 1:nrow(dataw)

    cortable <- dataw[, c("id", "score")]

    hist(cortable$score, breaks = 100, xlab = "", main = "Euclidean distance score")
  }
  else if(use_score == "pearson"){
    cortable <- cor(t(dataw[,1:ncol(dataw)]), reference, use = "pairwise.complete.obs")
    cortable <- as.data.frame(cortable)
    colnames(cortable) <- "score"
    cortable <- tibble::rownames_to_column(cortable, var = "id")

    hist(cortable$score, breaks = 100, xlab = "", main = "Distribution of Pearson correlations")
  }


  if (score_threshold > 0) {
    if (include_neg) {
      cortable <- subset(cortable, abs(score) >=
                           score_threshold)
    }
    else {
      cortable <- subset(cortable, score >= score_threshold)
    }
  }
  if (nrow(cortable) == 0) {
    g <- ggplot(data.frame(x = c(0,1), y = c(0,1)), aes(x,y, label = "s")) +
      geom_text(x=0.5, y=0.5, label = paste0("No proteins pass the score threshold of ",
                                             score_threshold, "\n try to lower the threshold",
                                             "\nor change the maximum number",
                                             "\nof missing values"), size = 6) +
      theme_cowplot() +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank())

    return(g)
  }
  else {
    message(paste0(nrow(cortable), " proteins pass the score threshold of ",
                   score_threshold))
  }

  cortable <- merge(proteininfo, cortable)

  return(cortable[order(cortable$score, decreasing = TRUE),])
}

