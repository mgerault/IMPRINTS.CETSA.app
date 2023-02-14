#' imprints_heatmap
#'
#' Function to get the heatmap from your data.
#'
#' @param data Output data from imprints_average_sh
#' @param hit_summary The summary file from the hitlist output
#' @param NN_data The NN file from the hitlist output
#' @param PRcomplex_data Output data from imprints_complex_mapping_sh. If not NULL and hit_summary NULL, will
#'                 print different heatmaps according to the protein complex.
#' @param treatment A character telling the condition from which you want to see the heatmap
#' @param max_na An integer indicating the maximum number of missing values per row (per protein)
#' @param response A character to tell if you want to the destabilized proteins, stabilized or both.
#'                 Accepted value are 'D' for destabilize, 'S' for stabilized or 'both'.
#' @param select_cat A character vector indicating the categories from which you want to see the heatmap
#' @param saveHeat Logical to tell if you want to save the heatmap
#' @param file_type The format file, either 'png' or 'pdf'.
#' @param file_name The file name
#' @param titleH The title for your heatmap
#' @param gradient_color The color for the gardient of the heatmap. Can only be of length three.
#' @param cat_color A list which contains the colors for each categories you selected.
#' @param back_color The color from the background of the heatmap (can be NULL)
#' @param border_color The color from the border of the plot (can be NULL)
#'
#' @return A grob object, the heatmap.
#'
#' @seealso \code{\link{imprints_average_sh}} , \code{\link{imprints_complex_mapping_sh}}
#'
#' @export
#'

imprints_heatmap <- function(data, hit_summary = NULL, NN_data = NULL,
                          PRcomplex_data = NULL,
                          treatment, max_na = 0,
                          response = c("both", "S", "D"),
                          select_cat = c("CC", "CN", "NC"),
                          saveHeat = FALSE, file_type = c("png", "pdf"), file_name = "Heatmap",
                          titleH = "Elutriation heatmap", gradient_color = c("#005EFF", "#FFFFFF", "#FF0000"),
                          cat_color = list("CC" = "#FB4F0B", "CN" = "#0FAEB9", "NC" = "#E7B700"),
                          back_color = "#FFFFFF", border_color = NULL){
  response <- match.arg(response)
  file_type <- match.arg(file_type)
  if (length(treatment) != 1) {
    stop("Please provide one treatment name")
  }
  if (inherits(data, "data.frame")) {
    subset <- stringr::str_which(names(data), paste0("_", treatment, "$", collapse = "|"))
    if (length(subset) > 0) {
      data1 = data[, c(1, 2, subset)]
    }
    else {
      stop("Please provide the right treatment keyword character")
    }
  }
  if (length(grep("countNum", names(data1)))) {
    countinfo1 <- unique(data1[,str_which(names(data1), "^id$|^sumPSM|^countNum|^sumUniPeps")])
    data1 <- data1[, -str_which(names(data1), "^sumPSM|^countNum|^sumUniPeps")]
  }
  if (length(grep("description", names(data1)))) {
    proteininfo <- unique(data1[, c("id", "description")])
    data1$description <- NULL
  }
  if(max_na < 0){
    stop("Please enter a positive value for the maximum NA per row")
  }
  nb_na <- apply(data1, 1, function(x) sum(is.na(x)))
  na_thresh <- which(nb_na <= max_na)
  data1 <- data1[na_thresh,]

  if(length(gradient_color) < 3){
    stop("Please select at least 3 colors")
  }

  message("Clustering data")
  d <- dist(data1[,-1])
  if(length(which(is.na(d)))){
    d[which(is.na(d))] <- 0
  }
  prot_dend <- hclust(d)

  message("Prepare data for plotting")
  datal <- gather(data1, treatment, reading, -id)
  if (length(unlist(strsplit(datal$treatment[1], "_"))) == 3) {
    datal <- separate(datal, treatment, into = c("set",
                                                 "temperature", "condition"))
    datal <- unite(datal, "condition", c("set", "condition"))
  }
  else if (length(unlist(strsplit(datal$treatment[1], "_"))) == 2) {
    datal <- separate(datal, treatment, into = c("temperature",
                                                 "condition"))
  }


  message("Filtering your data")
  treatment = unique(datal$condition)
  datal$condition <- NULL


  pr_axis <- list(element_blank(), element_blank())
  face_sw <- "y"
  if(!is.null(hit_summary)){
    df_hits <- hit_summary %>% dplyr::filter(Condition == treatment)
    if(!is.null(NN_data)){
      df_NN <- NN_data %>% dplyr::filter(Condition == treatment) %>%
        dplyr::group_by(id,Condition,category) %>%  dplyr::summarise()

      df_hits <- rbind(df_hits, df_NN)
    }
    df_hits$Condition <- NULL

    datah <- dplyr::inner_join(datal, df_hits, by = "id")

    if(sum((select_cat %in% unique(datah$category))) != length(select_cat)){
      stop("Please choose valid categories (the ones present in your data).")
    }
    else if(length(select_cat) > 0 & length(select_cat) != length(unique(datah$category))){
      datah <- datah %>% dplyr::filter(!is.na(match(category, select_cat)))
    }
  }
  else if(!is.null(PRcomplex_data)){
    PRcomplex_data <- PRcomplex_data[,c("id", "ComplexName")]
    colnames(PRcomplex_data)[2] <- "category"

    datah <- dplyr::inner_join(datal, PRcomplex_data, by = "id")
    datah$category <- unlist(lapply(datah$category, function(x){
                                                     if(str_length(x) > 25){
                                                       s <- str_split(x, " ")[[1]]
                                                       s <- str_remove_all(s, "^ ")
                                                       s <- s[str_length(s) != 0]
                                                       s <- paste(s, collapse = " \n")
                                                     }
                                                     else{
                                                       s <- x
                                                     };
                                                    s
                                                    })
                             )

    pr_axis <- list(element_text(), element_line())
    face_sw <- NULL
  }
  else{
    datah <- datal
    datah$category <- rep("", nrow(datah))
  }

  if(response == "S"){
    datah <- datah %>% group_by(id) %>% mutate(M = mean(reading, na.rm = TRUE))
    datah <- datah %>% filter(M >= 0)
  }
  else if(response == "D"){
    datah <- datah %>% group_by(id) %>% mutate(M = mean(reading, na.rm = TRUE))
    datah <- datah %>% filter(M <= 0)
  }
  else if(response != "both"){
    stop("Please enter a valide reponse : 'S' for stabilization, 'D' for destabilization or 'both'.")
  }

  lcol <- max(c(abs(round(min(datah$reading, na.rm = TRUE))),
                ceiling(max(datah$reading, na.rm = TRUE))))
  br <- c(-lcol, 0, lcol)
  cl <- gradient_color

  datah$id <- factor(datah$id, levels = data1$id[prot_dend$order])
  datah$temperature <- factor(datah$temperature, levels = unique(datah$temperature), ordered = TRUE)
  if(!is.null(hit_summary)){
    datah$category <- factor(datah$category, levels = c("CN", "NC", "CC", "ND", "NN"), ordered = TRUE)
  }

  message("Getting plot")
  H <- ggplot(datah, aes(temperature, id, fill = reading)) + geom_tile() +
    scale_fill_gradientn(breaks = br,
                         colors = cl,
                         limits = c(br[1], br[length(br)])) +
    facet_grid(rows = vars(datah$category),
               scales = "free", space="free_y", switch = face_sw) +
    labs(x="", y="", title = titleH, subtitle = paste("Condition :", treatment),
         fill = "Protein \nabundance \ndifference") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, family = "serif", size = 20),
          plot.subtitle = element_text(size = 10, family = "serif", face = "italic"),
          legend.title = element_text(size = 9, family = "serif"),
          axis.text.y = pr_axis[[1]],
          axis.ticks.y = pr_axis[[2]],
          legend.background = element_rect(
            fill = back_color,
            size = 1
          ),
          plot.background = element_rect(
            fill = back_color,
            colour = border_color,
            size = 1
          ))


  H <- ggplot_gtable(ggplot_build(H))

  if(!is.null(hit_summary)){
    stripr <- which(grepl('strip-l', H$layout$name)) #strip-l correspond to the left facet label
    fills <- c("CN" = "#0FAEB9", "NC" = "#E7B700", "CC" = "#FB4F0B", "ND" = "#8F3A8461", "NN" = "#ABABAB")
    if(!is.null(cat_color)){
      for(i in names(cat_color)){
        fills[[i]] <- cat_color[[i]]
      }
    }
    fills_ord <- c("CN" = 1, "NC" = 2, "CC" = 3, "ND" = 4, "NN" = 5)
    fills <- fills[select_cat][order(fills_ord[select_cat])]
  }
  else if(!is.null(PRcomplex_data)){
    stripr <- which(grepl('strip-r', H$layout$name)) #strip-r correspond to the right facet label
    fills <- PaletteWithoutGrey(unique(datah$category))
  }
  else{
    stripr <- which(grepl('strip-l', H$layout$name)) #strip-l correspond to the left facet label
    fills <- "#5691FC"
  }
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', H$grobs[[i]]$grobs[[1]]$childrenOrder))
    H$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }

  if(saveHeat){
    if(file_type == "png"){
      png(paste0(file_name, ".png"))
      plot(H)
      dev.off()
      message("Heatmap saved !")
    }
    else if(file_type == "pdf"){
      pdf(paste0(file_name, ".pdf"))
      plot(H)
      dev.off()
      message("Heatmap saved !")
    }
    else{
      message("The plot couldn't have been saved. Please choose a valide format file : 'png' of 'pdf'.")
    }
  }

  return(H)
}


### PaletteWithoutGrey function ###
#generates a color list depending on the number of element of a character vector
PaletteWithoutGrey <- function(treatment){

  n = length(unique(treatment))
  x <- grDevices::colors(distinct = TRUE)                           #all the color from R
  mycol <- x[which(is.na(stringr::str_extract(x, "gr(e|a)y")))]   #keep only colors that are not grey

  listcolor <- c()
  for (i in 0:(n-1))
    listcolor <- append(listcolor, mycol[i*20 + 9])      #save a color from the list (the number 20 and 9 were chosen in order to have distincts colors, this is empirical, can be changed)

  return(listcolor)
}



