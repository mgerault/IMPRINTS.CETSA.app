#' imprints_barplotting_peptides
#'
#' Function to generate IMPRINTS bar plot and pdf file with multipanel bar plots for IMPRINTS-CETSA data.
#' This function is based on the function imprints_barplotting from the IMPRINTS.CETSA package.
#'
#' @param data dataset after \code{imprints_sequence_peptides} to plot.
#' @param treatmentlevel a vector of treatment labels, such as c("DMSO","TNFa","AT26533")
#'                       the order determines the arrangement, so in this case DMSO
#'                       group would be the first group
#' @param RESP Logical to tell if you want to plot IMPRINTS profile in the RESP format (\code{\link{imprints_cleaved_peptides}}).
#'             In order to be able to do it, you'll need exactly 2 peptide sequence per protein.
#' @param printBothName A logical to tell if you want to print the both protein names on the plot
#' @param printGeneName A logical to tell if you want to print the gene names on the plot
#' @param witherrorbar A logical to print or not the error bar on the plot
#' @param withpoint A logical to print or not the data point of each replicate on the plot on top of the bars
#' @param colorpanel a vector of color scheme provided by default with the function PaletteWithoutGrey
#' @param usegradient whether the barplot should be draw in color gradient format
#' @param colorgradient the color scheme of gradient applied, default value c("#4575B4","ivory", "#D73027")
#' @param linegraph whether to plot the graph in a line graph format, default set to FALSE
#' @param log2scale whether the yscales should be in log2 scale, default set to TRUE
#' @param ratio aspect ratio of the plot, default set to 0.6
#' @param ret_plot Logical to tell if you want to return the last plot
#' @param save_pdf A logical to tell if you want to save plots in a pdf file
#' @param layout a vector indicating the panel layout for multi-panel plots per page,
#'               default value is c(2,3) for set containing data, otherwise c(4,3), use when save_pdf = TRUE
#' @param toplabel textual label at the top part of the page
#' @param leftlabel textual label at the left side of the page
#' @param bottomlabel textual label at the bottom part of the page
#' @param pdfname textual label of the pdf file
#' @param pdfheight a number indicate the height of pdf file, default value 12
#' @param pdfwidth a number indicate the width of pdf file, default value 12
#'
#'
#' @return The imprints barplot
#'
#' @seealso \code{\link{imprints_barplotting_app}}
#'
#' @export
#'

imprints_barplotting_peptides <- function(data, treatmentlevel = get_treat_level(data), RESP = FALSE,
                                          printBothName = TRUE, printGeneName = FALSE,
                                          witherrorbar = TRUE, withpoint = FALSE, layout = NULL,
                                          colorpanel = PaletteWithoutGrey(treatmentlevel),
                                          usegradient = FALSE, colorgradient = c("#4575B4", "ivory", "#D73027"),
                                          linegraph = FALSE, log2scale = TRUE, ratio = 1.2,
                                          ret_plot = TRUE, save_pdf = FALSE,
                                          toplabel = "IMPRINTS-CETSA bar plotting", leftlabel = "", bottomlabel = "",
                                          pdfname = "barplot",  pdfheight = 12, pdfwidth = 12){
  needed_columns <- c("Master.Protein.Accessions", "description",
                      "Positions.in.Master.Proteins", "Annotated.Sequence",
                      "Modifications")
  missing_columns <- needed_columns[!(needed_columns %in% colnames(data))]
  if(length(missing_columns)){
    message(paste("Error:",
                  paste(missing_columns, collapse = ", "),
                  ifelse(length(missing_columns) > 1, "are", "is"),
                  "missing in your data ! Are you sure you selected a peptides dataset ?"))
    return()
  }

  if("countNum" %in% colnames(data)){
    data$countNum <- NULL
  }

  if(RESP){
    check_resp <- data %>%
      group_by(Master.Protein.Accessions) %>%
      count()
    check_resp <- all(check_resp$n == 2)
    if(!check_resp){
      message("Error: to be able to plot your IMPRINTS peptides dataset in the RESP format, you need excatly two peptides sequence per protein !")
      return()
    }
  }

  if(save_pdf){
    dataname <- deparse(substitute(data))
  }

  ### function to plot IMPRINTS profiles
  barplotting <- function(d1, RESP = TRUE){
    if (!log2scale) {
      minreading = 0.5
      maxreading = 2
      if(withpoint){
        pts <- as.numeric(unlist(strsplit(d1$pts, "; ")))
        legendscale = c(min(max(min(pts, na.rm = T) - 0.5, 0), minreading),
                        max(max(pts, na.rm = T) + 0.5, maxreading))
      }
      else{
        legendscale = c(min(max(min(d1$mean, na.rm = T) -  0.5, 0), minreading),
                        max(max(d1$mean, na.rm = T) + 0.5, maxreading))
      }
    }
    else {
      minreading = -0.5
      maxreading = 0.5
      if(withpoint){
        pts <- as.numeric(unlist(strsplit(d1$pts, "; ")))
        legendscale = c(min(min(pts, na.rm = T) - 0.1, minreading),
                        max(max(pts, na.rm = T) + 0.1, maxreading))
      }
      else{
        legendscale = c(min(min(d1$mean, na.rm = T) - 0.1, minreading),
                        max(max(d1$mean, na.rm = T) + 0.1, maxreading))
      }
    }

    if(RESP){
      ord_position <- unique(d1$global.position)
      ord_position <- gsub(".*\\[|\\]", "", ord_position)
      ord_position <- lapply(strsplit(ord_position, "-|~"),
                             function(x){
                               sum(as.numeric(x))
                             })
      ord_position <- unique(d1$global.position)[order(unlist(ord_position))]
      d1$global.position <- factor(d1$global.position, levels = ord_position)
    }

    d1$QP <- FALSE
    if("36C" %in% d1$temperature){
      d1$QP[which(d1$temperature == "36C")] <- TRUE
      lvl_tokeep <- levels(d1$condition)
      lvl_tokeep <- gsub("36C", "QP", lvl_tokeep)
      d1$condition <- as.character(d1$condition)
      d1$condition[which(d1$temperature == "36C")] <- gsub("36C", "QP", d1$condition[which(d1$temperature == "36C")])
      d1$condition <- factor(d1$condition, levels = lvl_tokeep)
    }
    if (linegraph) {
      colorpanel <- PaletteWithoutGrey(temperature)
      q <- ggplot(d1, aes(x = treatment, y = mean,
                          group = temperature, color = temperature)) +
        geom_line() + geom_point() + coord_cartesian(ylim = legendscale) +
        scale_color_manual(drop = FALSE, values = colorpanel)
    }
    else if (!usegradient) {
      q <- ggplot(d1, aes(x = condition, y = mean,
                          fill = treatment)) +
        geom_bar(stat = "identity", aes(color = QP), size = rel(0.85)) +
        coord_cartesian(ylim = legendscale) +
        scale_fill_manual(drop = FALSE, values = colorpanel) +
        scale_color_manual(values = c("TRUE" = "#656565", "FALSE" = "#FFFFFF00")) +
        guides(color = "none") +
        scale_x_discrete(labels = gsub("_.{1,}", "", levels(d1$condition)))
    }
    else {
      q <- ggplot(d1, aes(x = condition, y = mean,
                          fill = mean)) +
        geom_bar(stat = "identity", aes(color = QP), size = rel(0.85)) +
        coord_cartesian(ylim = legendscale) +
        scale_fill_gradient2(limits = legendscale,
                             low = colorgradient[1], mid = colorgradient[2],
                             high = colorgradient[3], midpoint = 0, na.value = "gray90",
                             guide = guide_colorbar("")) +
        scale_color_manual(values = c("TRUE" = "#656565", "FALSE" = "#FFFFFF00")) +
        guides(color = "none") +
        scale_x_discrete(labels = gsub("_.{1,}", "", levels(d1$condition)))
    }
    if (witherrorbar) {
      if (linegraph) {
        q <- q + geom_errorbar(aes(ymin = mean - se,
                                   ymax = mean + se), width = 0.1)
      }
      else {
        q <- q + geom_errorbar(aes(ymin = mean - se,
                                   ymax = mean + se), width = 0.2, position = position_dodge(0.9))
      }
    }
    if (log2scale) {
      q <- q + ylab("fold change(log2)") + ggtitle(as.character(unique(d1$id)))
    }
    else {
      q <- q + ylab("fold change") + ggtitle(as.character(unique(d1$id)))
    }
    if(withpoint){
      if(RESP){
        d1_pts <- d1 %>%
          group_by(id, temperature, treatment, condition, global.position) %>%
          group_modify(~ {
            pts <- as.numeric(unlist(strsplit(.x$pts, "; ")))
            rep <- unlist(strsplit(.x$biorep, "; "))
            df <- .x
            df$pts <- NULL
            df$biorep <- NULL
            df <- Reduce(rbind, lapply(1:length(pts), function(x) df))
            df$pts <- pts
            df$replicate <- rep

            return(df)
          })
      }
      else{
        d1_pts <- d1 %>%
          group_by(id, temperature, treatment, condition) %>%
          group_modify(~ {
            pts <- as.numeric(unlist(strsplit(.x$pts, "; ")))
            rep <- unlist(strsplit(.x$biorep, "; "))
            df <- .x
            df$pts <- NULL
            df$biorep <- NULL
            df <- Reduce(rbind, lapply(1:length(pts), function(x) df))
            df$pts <- pts
            df$replicate <- rep

            return(df)
          })
      }

      q <- q +
        geom_point(data = d1_pts, aes(x = condition, y = pts, shape = replicate),
                   size = rel(1.5), fill = NA) +
        scale_shape_manual(values = c(1,2,4,5,6,7,8))
    }

    if(RESP){
      q <- q + facet_wrap(~global.position) +
        cowplot::theme_cowplot() +
        theme(text = element_text(size = 10),
              strip.text.x = element_text(size = 8),
              strip.background.x = element_rect(fill = "white", color = "black"),
              plot.title = element_text(hjust = 0.5,size = rel(0.8)),
              legend.background = element_rect(fill = NULL),
              legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.15,"cm"),
              legend.title = element_text(face = "bold"),
              legend.text = element_text(size = rel(0.7)),
              legend.justification = "center", panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), strip.background = element_blank(),
              axis.line.x = element_line(), axis.line.y = element_line(),
              axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.7)),
              aspect.ratio = ratio)
    }
    else{
      q <- q  +
        cowplot::theme_cowplot() +
        theme(text = element_text(size = 10),
              strip.text.x = element_text(size = 8),
              plot.title = element_text(hjust = 0.5,size = rel(0.8)),
              legend.background = element_rect(fill = NULL),
              legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.15,"cm"),
              legend.title = element_text(face = "bold"),
              legend.text = element_text(size = rel(0.7)),
              legend.justification = "center", panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), strip.background = element_blank(),
              axis.line.x = element_line(), axis.line.y = element_line(),
              axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.7)),
              aspect.ratio = ratio)
    }

    return(q)
  }

  message("Preparing data for plotting...")
  nrowdata <- nrow(data)
  if (printBothName) {
    data <- data %>% dplyr::rowwise() %>%
      dplyr::mutate(description1 = getProteinName(description),
                    description2 = getGeneName(description),
                    Master.Protein.Accessions = paste(Master.Protein.Accessions,
                                                      description1, description2,
                                                      sep = "\n"))
      data$description1 <- NULL
      data$description2 <- NULL
  }
  else if (printGeneName) {
    data <- data %>% dplyr::rowwise() %>%
      dplyr::mutate(description = getGeneName(description),
                    Master.Protein.Accessions = paste(Master.Protein.Accessions,
                                                      description, sep = "\n"))
  }
  else {
    data <- data %>% dplyr::rowwise() %>%
      dplyr::mutate(description = getProteinName(description),
                    Master.Protein.Accessions = paste(Master.Protein.Accessions,
                                                      description, sep = "\n"))
  }
  data$description <- NULL

  if(any(is.na(data$Modifications))){
    data$Modifications[which(is.na(data$Modifications))] <- ""
  }
  data$id <- paste(data$Master.Protein.Accessions,
                   paste(data$Positions.in.Master.Proteins, data$Modifications),
                   data$Annotated.Sequence, sep = "\n")
  data <- data[,-grep(paste(needed_columns, collapse = "|"),
                      colnames(data))
               ]

  data1 <- data %>%
    tidyr::gather("condition", "reading", -id)

  if (!log2scale) {
    data1 <- data1 %>%
      dplyr::mutate(reading = 2^reading)
  }
  if(length(treatmentlevel) != length(get_treat_level(data))){
    data1 <- data1[grep(paste0("_", treatmentlevel, "($|_)", collapse = "|"), data1$condition),]
  }
  a <- data1$condition[1]
  if(length(unlist(strsplit(a, "_"))) == 3){
    data1 <- tidyr::separate(data1, condition, into = c("temperature", "replicate", "treatment"), sep = "_")
    temperature <- sort(unique(data1$temperature))
    temp_idx <- grep("^[0-9]", temperature)

    if(length(temp_idx) != length(temperature)){
      temperature <- c(sort(temperature[-temp_idx]), sort(temperature[temp_idx]))
    }
    data1$id <- factor(data1$id, levels = unique(data1$id), ordered = TRUE) #preserve order

    if(withpoint){
      cdata <- plyr::ddply(data1, c("id", "temperature", "treatment"),
                           summarise,
                           N = length(na.omit(reading)),
                           mean = mean(reading,na.rm = T),
                           sd = sd(reading, na.rm = T),
                           se = sd/sqrt(N),
                           pts = paste0(reading, collapse = "; "),
                           biorep = paste0(replicate, collapse = "; ")
                           )
    }
    else{
      cdata <- plyr::ddply(data1, c("id", "temperature", "treatment"),
                           summarise, N = length(na.omit(reading)),
                           mean = mean(reading,na.rm = T),
                           sd = sd(reading, na.rm = T),
                           se = sd/sqrt(N))
    }

    cdata$id <- as.character(cdata$id)

    if(length(layout) == 0) {
      layout <- c(4, 3)
    }

  }
  else {
    stop("make sure the namings of the columns of the dasaset are correct.")
  }

  cdata <- cdata %>%
    dplyr::rowwise() %>%
    dplyr::mutate(condition = paste(temperature, treatment, sep = "_"))
  cdata$id <- factor(cdata$id, levels = unique(cdata$id), ordered = TRUE)

  cdata$treatment <- factor(as.character(cdata$treatment),
                            levels = treatmentlevel)
  cdata$condition <- factor(as.character(cdata$condition),
                            levels = apply(expand.grid(temperature, treatmentlevel),
                                           1, paste, collapse = "_"))

  # if data with different temperatures, prevent from creating non sense factors
  cdata$condition <- factor(as.character(cdata$condition),
                            levels = levels(cdata$condition)[levels(cdata$condition)  %in% as.character(cdata$condition)]
                            )

  message("Generating fitted plot, pls wait.")
  if(RESP){
    cdata$id <- as.character(cdata$id)

    cdata$global.position <- unlist(lapply(strsplit(cdata$id, "\n"),
                                           "[[", ifelse(printBothName, 4, 3))
                                    )
    cdata$global.position <- sub(" {1,}$", "", cdata$global.position)

    cdata$id <- unlist(lapply(lapply(strsplit(cdata$id, "\n"),
                                     "[", 1:ifelse(printBothName, 3, 2)),
                              paste, collapse = "\n")
                       )
    cdata$id <- factor(cdata$id, ordered = TRUE)

    plots <- plyr::dlply(cdata, plyr::.(id), .fun = barplotting, RESP = RESP)
  }
  else{
    plots <- plyr::dlply(cdata, plyr::.(id), .fun = barplotting, RESP = RESP)
  }

  if(save_pdf){
    message("Start saving plot")

    params <- list(nrow = layout[1], ncol = layout[2])
    n <- with(params, nrow * ncol)
    pages <- length(plots)%/%n + as.logical(length(plots)%%n)
    groups <- split(seq_along(plots), gl(pages, n, length(plots)))
    n_p <- length(names(groups))
    pl <- lapply(names(groups), function(i) {
      message(paste("Saving page", i, "/", n_p))
      do.call(gridExtra::arrangeGrob,
              c(plots[groups[[i]]], params,
                top = toplabel, left = leftlabel,
                bottom = bottomlabel)
              )
      })

    class(pl) <- c("arrangelist", "ggplot", class(pl))
    pdfname <- paste0(pdfname, ".pdf")
    ggsave(file = paste0(format(Sys.time(), "%y%m%d_%H%M_"),
                         dataname, "_", pdfname),
           pl, height = pdfheight, width = pdfwidth)

    message("IMPRINTS-CETSA bar plot file generated successfully.")
  }

  if(ret_plot){
    message("IMPRINTS-CETSA bar plot generated successfully.")
    return(plots)
  }
  else{
    g <- ggplot(data.frame(x = c(0,1), y = c(0,1)), aes(x,y, label = "s")) +
      geom_text(x=0.5, y=0.5, label = "All the barplots has been saved succesfully !
                                       \nGo check your files", size = 6) +
      cowplot::theme_cowplot() +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank())

    return(g)
  }
}


### PaletteWithoutGrey function ###
#generates a color list depending on the number of element of a character vector
PaletteWithoutGrey <- function(treatment){
  n = length(unique(treatment))
  x <- grDevices::colors(distinct = TRUE)                           #all the color from R
  mycol <- x[-grep("gr(e|a)y", x)]   #keep only colors that are not grey

  listcolor <- c()
  for (i in 0:(n-1)){
    listcolor <- append(listcolor, mycol[((i*20 + 9) %% length(mycol)) + 1])      #save a color from the list (the number 20 and 9 were chosen in order to have distincts colors, this is empirical, can be changed)
  }

  return(listcolor)
}

getGeneName <- function(x){
  gene = strsplit(strsplit(x, "GN=")[[1]][2], " ")[[1]][1]
  if (length(gene) == 0) {
    return(" ")
  }
  else {
    return(gene)
  }
}

getProteinName <- function(x){
  protein = strsplit(x, " OS=")[[1]][1]
  if (length(protein) == 0) {
    return(" ")
  }
  else {
    return(protein)
  }
}

