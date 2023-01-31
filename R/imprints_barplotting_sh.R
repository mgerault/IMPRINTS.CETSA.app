#' imprints_barplotting_sh
#'
#' Function to generate 2D bar plot and pdf file with multipanel bar plots for 2D-CETSA data.
#' It is totally based on the function ms_2D_barplotting from the mineCETSA package.
#' It was retaken for easier use in the shiny app.
#'
#' @param data dataset after ms_2D_caldiff to plot. Can also be a list of this dataset.
#' @param treatmentlevel a vector of treatment labels, such as c("DMSO","TNFa","AT26533")
#'                       the order determines the arrangement, so in this case DMSO
#'                       group would be the first group
#' @param corrtable A correlation table
#' @param setlevel a vector of set information if any, such as c("M13","M16")
#' @param plotseq a vector of plots arragement sequence (in composite ID)
#' @param printBothName A logical to tell if you want to print the both protein names on the plot
#' @param printGeneName A logical to tell if you want to print the gene names on the plot
#' @param pfdatabase A logical for using pdf database or not
#' @param witherrorbar A logical to print or not the error bar on the plot
#' @param colorpanel a vector of customizable color scheme provided by default with the function PaletteWithoutGrey
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
#' @return The ms 2D barplot
#'
#' @examples
#' library(mineCETSA)
#' library(mineCETSAapp)
#'
#' elutriation_wVeh <- elutriation[,-grep("G1",names(elutriation))]
#' O43776_elu <- ms_subsetting(elutriation_wVeh, isfile = FALSE, hitidlist = c("O43776"))
#' imprints_barplotting_sh(O43776_elu)
#'
#' @seealso \code{\link{ms_2D_barplotting}}
#'
#' @export
#'

imprints_barplotting_sh <- function (data, treatmentlevel = get_treat_level(data), setlevel = NULL, corrtable = NULL,
                                  plotseq = NULL, printBothName = TRUE, printGeneName = FALSE,
                                  pfdatabase = FALSE, witherrorbar = TRUE, layout = NULL,
                                  colorpanel = PaletteWithoutGrey(treatmentlevel),
                                  usegradient = FALSE, colorgradient = c("#4575B4", "ivory", "#D73027"),
                                  linegraph = FALSE, log2scale = TRUE, ratio = 0.6,
                                  ret_plot = TRUE,
                                  save_pdf = FALSE, toplabel = "IMPRINTS-CETSA bar plotting",
                                  leftlabel = "", bottomlabel = "", pdfname = "barplot",
                                  pdfheight = 12, pdfwidth = 12)
{

  if(save_pdf){
    dataname <- deparse(substitute(data))
  }

  barplotting <- function(d1, withset = FALSE) {
    if (withset) {
      d1 <- droplevels(d1)
      d1_list <- split(d1, d1$set)
      q_list <- list()
      n_loop <- 1
      for (j in names(d1_list)) {
        if (nrow(d1_list[[j]]) > 0) {
          d2 <- d1_list[[j]]
          if (!log2scale) {
            minreading = 0.5
            maxreading = 2
            legendscale = c(min(max(min(d2$mean, na.rm = T) -
                                      0.5, 0), minreading), max(max(d2$mean,
                                                                    na.rm = T) + 0.5, maxreading))
          }
          else {
            minreading = -0.5
            maxreading = 0.5
            legendscale = c(min(min(d2$mean, na.rm = T) -
                                  0.1, minreading), max(max(d2$mean, na.rm = T) +
                                                          0.1, maxreading))
          }
          q <- ggplot(d2, aes(x = condition, y = mean,
                              fill = treatment)) + geom_bar(stat = "identity") +
            coord_cartesian(ylim = legendscale) + scale_fill_manual(drop = FALSE,
                                                                    values = colorpanel)
          if (witherrorbar) {
            q <- q + geom_errorbar(aes(ymin = mean -
                                         se, ymax = mean + se), width = 0.2, position = position_dodge(0.9))
          }
          if (log2scale) {
            q <- q + ylab("fold change(log2)") + ggtitle(paste(j,
                                                               as.character(unique(d2$id)), sep = "\n"))
          }
          else {
            q <- q + ylab("fold change") + ggtitle(paste(j,
                                                         as.character(unique(d2$id)), sep = "\n"))
          }
          q <- q + labs(subtitle = subt$category[n_loop]) +
            theme_cowplot() + theme(text = element_text(size = 10),
                                    strip.text.x = element_text(size = 5),
                                    plot.title = element_text(hjust = 0.5,
                                                              size = rel(0.8)),
                                    legend.background = element_rect(fill = NULL),
                                    legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.15,"cm"),
                                    legend.title = element_text(face = "bold"),
                                    legend.text = element_text(size = rel(0.7)),
                                    legend.justification = "center", panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), strip.background = element_blank(),
                                    axis.line.x = element_line(), axis.line.y = element_line(),
                                    axis.text.x = element_text(angle = 45, hjust = 1,
                                                               size = rel(0.7)), aspect.ratio = 0.6)
          q_list[[j]] <- q
          n_loop <- n_loop + 1
        }
        else {
          q <- ggplot()
          q_list[[j]] <- q
        }
      }
      q_list <- gridExtra::grid.arrange(grobs = q_list,
                                        ncol = 1)
      return(q_list)
    }

    else {
      if (!log2scale) {
        minreading = 0.5
        maxreading = 2
        legendscale = c(min(max(min(d1$mean, na.rm = T) -
                                  0.5, 0), minreading), max(max(d1$mean, na.rm = T) +
                                                              0.5, maxreading))
      }
      else {
        minreading = -0.5
        maxreading = 0.5
        legendscale = c(min(min(d1$mean, na.rm = T) -
                              0.1, minreading), max(max(d1$mean, na.rm = T) +
                                                      0.1, maxreading))
      }
      d1$QP <- FALSE
      if("36C" %in% d1$temperature){
        d1$QP[which(d1$temperature == "36C")] <- TRUE
        lvl_tokeep <- levels(d1$condition)
        lvl_tokeep <- stringr::str_replace_all(lvl_tokeep, "36C", "QP")
        d1$condition <- as.character(d1$condition)
        d1$condition[which(d1$temperature == "36C")] <- stringr::str_replace_all(d1$condition[which(d1$temperature == "36C")],
                                                                                 "36C", "QP")
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
                            fill = treatment)) + geom_bar(stat = "identity", aes(color = QP),
                                                          size = rel(0.85)) +
          coord_cartesian(ylim = legendscale) + scale_fill_manual(drop = FALSE,
                                                                  values = colorpanel) +
          scale_color_manual(values = c("TRUE" = "#656565", "FALSE" = "#FFFFFF00")) +
          guides(color = "none") +
          scale_x_discrete(labels = stringr::str_remove_all(levels(d1$condition), "_.{1,}"))
      }
      else {
        q <- ggplot(d1, aes(x = condition, y = mean,
                            fill = mean)) + geom_bar(stat = "identity", aes(color = QP),
                                                     size = rel(0.85)) +
          coord_cartesian(ylim = legendscale) +
          scale_fill_gradient2(limits = legendscale,
                               low = colorgradient[1], mid = colorgradient[2],
                               high = colorgradient[3], midpoint = 0, na.value = "gray90",
                               guide = guide_colorbar("")) +
          scale_color_manual(values = c("TRUE" = "#656565", "FALSE" = "#FFFFFF00")) +
          guides(color = "none") +
          scale_x_discrete(labels = stringr::str_remove_all(levels(d1$condition), "_.{1,}"))
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
      q <- q + labs(subtitle = subt[as.character(unique(d1$id)), "category"]) +
        theme_cowplot() + theme(text = element_text(size = 10),
                                strip.text.x = element_text(size = 5),
                                plot.title = element_text(hjust = 0.5,size = rel(0.8)),
                                legend.background = element_rect(fill = NULL),
                                legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.15,"cm"),
                                legend.title = element_text(face = "bold"),
                                legend.text = element_text(size = rel(0.7)),
                                legend.justification = "center", panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(), strip.background = element_blank(),
                                axis.line.x = element_line(), axis.line.y = element_line(),
                                axis.text.x = element_text(angle = 45, hjust = 1,
                                                           size = rel(0.7)),
                                aspect.ratio = ratio)
      return(q)
    }
  }
  if(("list" %in% class(data))){
    if(!save_pdf){
      stop("Your input data is a list. The aim is to save in the same pdf with diffrent data sets.
           Retry with setting 'save_pdf' to TRUE.")
    }
    else{
      pl <- list()
      sav_data <- data
      for(k in names(data)){
        data <- sav_data[[k]]
        treatmentlevel <- get_treat_level(data)

        nrowdata <- nrow(data)
        if (nrowdata == 0) {
          message("Make sure there are more than one experimental condition in dataset.")
          stop("Otherwise specify remsinglecondprot==FALSE !")
        }
        if (printBothName & !pfdatabase) {
          data <- data %>% dplyr::rowwise() %>% dplyr::mutate(description1 = getProteinName(description,
                                                                                          pfdatabase)) %>%
            dplyr::mutate(description2 = getGeneName(description)) %>%
            dplyr::mutate(id = paste(id, description1, description2,
                              sep = "\n"))
          data$description1 <- NULL
          data$description2 <- NULL
        }
        else if (printGeneName & !pfdatabase) {
          data <- data %>% dplyr::rowwise() %>%
            dplyr::mutate(description = getGeneName(description)) %>%
            dplyr::mutate(id = paste(id, description, sep = "\n"))
        }
        else {
          data <- data %>% dplyr::rowwise() %>%
            dplyr::mutate(description = getProteinName(description, pfdatabase)) %>%
            dplyr::mutate(id = paste(id, description, sep = "\n"))
        }

        if(length(grep("^category", names(data)))){
          if(length(grep("^score", names(data)))){
            subt <- data[, c(1, grep("^category", names(data)), grep("^score", names(data)))]
            subt <- as.data.frame(subt)
            colnames(subt) <- c("id", "category", "score")
            subt$category <- paste("Category :", subt$category, ", Score :", round(subt$score,4))
            subt$score <- NULL
            rownames(subt) <- subt$id
            data <- data[,-grep("^score", names(data))]
            data[grep("^category", names(data))] <- subt$category
          }
          else{
            subt <- data[, c(1, grep("^category", names(data)))]
            subt <- as.data.frame(subt)
            colnames(subt) <- c("id", "category")
            subt$category <- paste("Category :", subt$category)
            rownames(subt) <- subt$id
            data[grep("^category", names(data))] <- subt$category
          }

          ord_data <- data[NULL,]
          for(i in c("CN", "NC", "CC", "ND", "NN")){
            cat_idx <- str_which(data$category, paste0("^Category : ", i))
            if(length(cat_idx) > 0){
              cat_idx <- data[cat_idx,]
              w <- cat_idx[order(cat_idx$category, decreasing = TRUE),]
              ord_data <- rbind(ord_data, w)
            }
          }

          if(nrow(ord_data) != 0)
            data <- ord_data
        }
        else if(length(grep("^score", names(data)))){
          subt <- data[, c(1, grep("^score", names(data)))]
          subt <- as.data.frame(subt)
          colnames(subt) <- c("id", "score")
          subt$category <- paste("Score :", round(subt$score,4))  # keep same name for simplicity
          subt$score <- NULL
          rownames(subt) <- subt$id
          data <- data[order(data[,grep("^score", names(data))], decreasing = TRUE),]
          data <- data[,-grep("^score", names(data))]
        }
        else{
          subt <- NULL
        }

        data$description <- NULL
        data1 <- tidyr::gather(data[, -str_which(names(data), "^sumPSM|^countNum|^sumUniPeps|^drug$|^category")],
                               condition, reading, -id)

        if (!log2scale) {
          data1 <- dplyr::mutate(data1, reading = 2^reading)
        }
        a <- data1$condition[1]
        if (length(unlist(strsplit(a, "_"))) == 4) {
          withset <- TRUE
          data1 <- tidyr::separate(data1, condition, into = c("set",
                                                              "temperature", "replicate", "treatment"), sep = "_")
          temperature <- sort(unique(data1$temperature))
          temp_idx <- stringr::str_which(temperature, "^[0-9]")
          if(length(temp_idx) != length(temperature)){
            temperature <- c(sort(temperature[-temp_idx]), sort(temperature[temp_idx]))
          }
          data1$id <- factor(data1$id, levels = unique(data1$id), ordered = TRUE) #preserve order
          cdata <- plyr::ddply(data1, c("id", "set", "temperature",
                                        "treatment"),
                               summarise, N = length(na.omit(reading)),
                               mean = mean(reading, na.rm = T), sd = sd(reading, na.rm = T), se = sd/sqrt(N))
          cdata$id <- as.character(cdata$id)
          if (length(layout) == 0) {
            layout <- c(2, 3)
          }

        }
        else if (length(unlist(strsplit(a, "_"))) == 3) {
          withset <- FALSE
          data1 <- tidyr::separate(data1, condition, into = c("temperature",
                                                              "replicate", "treatment"), sep = "_")
          temperature <- sort(unique(data1$temperature))
          temp_idx <- stringr::str_which(temperature, "^[0-9]")
          if(length(temp_idx) != length(temperature)){
            temperature <- c(sort(temperature[-temp_idx]), sort(temperature[temp_idx]))
          }
          data1$id <- factor(data1$id, levels = unique(data1$id), ordered = TRUE) #preserve order
          cdata <- plyr::ddply(data1, c("id", "temperature", "treatment"),
                               summarise, N = length(na.omit(reading)), mean = mean(reading,na.rm = T),
                               sd = sd(reading, na.rm = T), se = sd/sqrt(N))
          cdata$id <- as.character(cdata$id)
          if (length(layout) == 0) {
            layout <- c(4, 3)
          }

        }
        else {
          stop("make sure the namings of the columns of the dasaset are correct.")
        }
        cdata <- cdata %>% dplyr::rowwise() %>% dplyr::mutate(condition = paste(temperature,
                                                                  treatment, sep = "_"))
        if (withset) {
          cdata$set <- factor(as.character(cdata$set), levels = setlevel)
        }
        if (class(corrtable) != "NULL") {
          corrtable <- corrtable[order(corrtable$correlation,
                                       decreasing = T), ]
          if (printBothName & !pfdatabase) {
            corrtable <- ms_composite_ID_Gene_Protein(corrtable,
                                                      pfdatabase)
          }
          else if (printGeneName & !pfdatabase) {
            corrtable <- ms_composite_ID_Gene(corrtable)
          }
          else {
            corrtable <- ms_composite_ID_Protein(corrtable,
                                                 pfdatabase)
          }
          cdata$id <- factor(cdata$id, levels = corrtable$id)
        }
        if (length(plotseq)) {
          cdata$id <- factor(cdata$id, levels = plotseq)
        }
        else {
          cdata$id <- factor(cdata$id, levels = unique(cdata$id), ordered = TRUE)
        }
        cdata$treatment <- factor(as.character(cdata$treatment),
                                  levels = treatmentlevel)
        cdata$condition <- factor(as.character(cdata$condition),
                                  levels = apply(expand.grid(temperature, treatmentlevel),
                                                 1, paste, collapse = "_"))
        # if data with different temperatures, prevent from creating non sense factors
        cdata$condition <- factor(as.character(cdata$condition),
                                  levels = levels(cdata$condition)[levels(cdata$condition)
                                                                   %in% as.character(cdata$condition)
                                                                   ]
                                  )
        message("Generating fitted plot, pls wait.")



        plots <- plyr::dlply(cdata, plyr::.(id), .fun = barplotting,
                             withset = withset)

        params <- list(nrow = layout[1], ncol = layout[2])
        n <- with(params, nrow * ncol)
        pages <- length(plots)%/%n + as.logical(length(plots)%%n)
        groups <- split(seq_along(plots), gl(pages, n, length(plots)))
        n_p <- length(names(groups))
        pl[[k]] <- lapply(names(groups), function(i) {
          message(paste("Saving page", i, "/", n_p))
          do.call(gridExtra::arrangeGrob,
                  c(plots[groups[[i]]], params,
                    top = paste(toplabel, k), left = leftlabel,
                    bottom = bottomlabel)
                  )
          })

        class(pl[[k]]) <- c("arrangelist", "ggplot", class(pl[[k]]))
      }

      message("Start saving plot")
      class(pl) <- c("arrangelist", "ggplot", class(pl))
      pdfname <- paste0(pdfname, ".pdf")
      ggpubr::ggexport(filename = paste0(format(Sys.time(), "%y%m%d_%H%M_"), dataname, "_", pdfname),
                       plotlist = pl,
                       height = pdfheight, width = pdfwidth)

      message("IMPRINTS-CETSA bar plot file generated successfully.")

      if(ret_plot){
        return(plots)
      }
      else{
        g <- ggplot(data.frame(x = c(0,1), y = c(0,1)), aes(x,y, label = "s")) +
          geom_text(x=0.5, y=0.5, label = "All the barplots has been saved succesfully !
                                         \nGo check your files", size = 6) +
          theme_cowplot() +
          theme(axis.text.x = element_blank(),
                axis.title.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text.y = element_blank(),
                axis.title.y = element_blank(),
                axis.ticks.y = element_blank())

        return(g)
      }
    }

  }
  else{
    nrowdata <- nrow(data)
    if (nrowdata == 0) {
      message("Make sure there are more than one experimental condition in dataset.")
      stop("Otherwise specify remsinglecondprot==FALSE !")
    }
    if (printBothName & !pfdatabase) {
      data <- data %>% dplyr::rowwise() %>% dplyr::mutate(description1 = getProteinName(description,
                                                                                      pfdatabase)) %>%
        dplyr::mutate(description2 = getGeneName(description)) %>%
        dplyr::mutate(id = paste(id, description1, description2,
                          sep = "\n"))
      data$description1 <- NULL
      data$description2 <- NULL
    }
    else if (printGeneName & !pfdatabase) {
      data <- data %>% dplyr::rowwise() %>%
        dplyr::mutate(description = getGeneName(description)) %>%
        dplyr::mutate(id = paste(id, description, sep = "\n"))
    }
    else {
      data <- data %>% dplyr::rowwise() %>%
        dplyr::mutate(description = getProteinName(description, pfdatabase)) %>%
        dplyr::mutate(id = paste(id, description, sep = "\n"))
    }

    if(length(grep("^category", names(data)))){
      if(length(grep("^score", names(data)))){
        subt <- data[, c(1, grep("^category", names(data)), grep("^score", names(data)))]
        subt <- as.data.frame(subt)
        colnames(subt) <- c("id", "category", "score")
        subt$category <- paste("Category :", subt$category, ", Score :", round(subt$score,4))
        subt$score <- NULL
        rownames(subt) <- subt$id
        data <- data[,-grep("^score", names(data))]
        data[grep("^category", names(data))] <- subt$category
      }
      else{
        subt <- data[, c(1, grep("^category", names(data)))]
        subt <- as.data.frame(subt)
        colnames(subt) <- c("id", "category")
        subt$category <- paste("Category :", subt$category)
        rownames(subt) <- subt$id
        data[grep("^category", names(data))] <- subt$category
      }

      ord_data <- data[NULL,]
      for(i in c("CN", "NC", "CC", "ND", "NN")){
        cat_idx <- str_which(data$category, paste0("^Category : ", i))
        if(length(cat_idx) > 0){
          cat_idx <- data[cat_idx,]
          w <- cat_idx[order(cat_idx$category, decreasing = TRUE),]
          ord_data <- rbind(ord_data, w)
        }
      }

      if(nrow(ord_data) != 0)
        data <- ord_data
    }
    else if(length(grep("^score", names(data)))){
      subt <- data[, c(1, grep("^score", names(data)))]
      subt <- as.data.frame(subt)
      colnames(subt) <- c("id", "score")
      subt$category <- paste("Score :", round(subt$score,4))  # keep same name for simplicity
      subt$score <- NULL
      rownames(subt) <- subt$id
      data <- data[order(data[,grep("^score", names(data))], decreasing = TRUE),]
      data <- data[,-grep("^score", names(data))]
    }
    else{
      subt <- NULL
    }

    data$description <- NULL
    data1 <- tidyr::gather(data[, -str_which(names(data), "^sumPSM|^countNum|^sumUniPeps|^drug$|^category")],
                           condition, reading, -id)
    if (!log2scale) {
      data1 <- dplyr::mutate(data1, reading = 2^reading)
    }
    if(length(treatmentlevel) != length(get_treat_level(data))){
      data1 <- data1 %>% filter(stringr::str_detect(condition, paste(treatmentlevel, collapse = "|")))
    }
    a <- data1$condition[1]
    if (length(unlist(strsplit(a, "_"))) == 4) {
      withset <- TRUE
      data1 <- tidyr::separate(data1, condition, into = c("set",
                                                          "temperature", "replicate", "treatment"), sep = "_")
      temperature <- sort(unique(data1$temperature))
      temp_idx <- stringr::str_which(temperature, "^[0-9]")
      if(length(temp_idx) != length(temperature)){
        temperature <- c(sort(temperature[-temp_idx]), sort(temperature[temp_idx]))
      }
      data1$id <- factor(data1$id, levels = unique(data1$id), ordered = TRUE) #preserve order
      cdata <- plyr::ddply(data1, c("id", "set", "temperature",
                                    "treatment"),
                           summarise, N = length(na.omit(reading)),
                           mean = mean(reading, na.rm = T), sd = sd(reading, na.rm = T), se = sd/sqrt(N))
      cdata$id <- as.character(cdata$id)
      if (length(layout) == 0) {
        layout <- c(2, 3)
      }

    }
    else if (length(unlist(strsplit(a, "_"))) == 3) {
      withset <- FALSE
      data1 <- tidyr::separate(data1, condition, into = c("temperature",
                                                          "replicate", "treatment"), sep = "_")
      temperature <- sort(unique(data1$temperature))
      temp_idx <- stringr::str_which(temperature, "^[0-9]")
      if(length(temp_idx) != length(temperature)){
        temperature <- c(sort(temperature[-temp_idx]), sort(temperature[temp_idx]))
      }
      data1$id <- factor(data1$id, levels = unique(data1$id), ordered = TRUE) #preserve order
      cdata <- plyr::ddply(data1, c("id", "temperature", "treatment"),
                           summarise, N = length(na.omit(reading)), mean = mean(reading,na.rm = T),
                           sd = sd(reading, na.rm = T), se = sd/sqrt(N))
      cdata$id <- as.character(cdata$id)

      if (length(layout) == 0) {
        layout <- c(4, 3)
      }

    }
    else {
      stop("make sure the namings of the columns of the dasaset are correct.")
    }
    cdata <- cdata %>% dplyr::rowwise() %>% dplyr::mutate(condition = paste(temperature,
                                                              treatment, sep = "_"))
    if (withset) {
      cdata$set <- factor(as.character(cdata$set), levels = setlevel)
    }
    if (class(corrtable) != "NULL") {
      corrtable <- corrtable[order(corrtable$correlation,
                                   decreasing = T), ]
      if (printBothName & !pfdatabase) {
        corrtable <- ms_composite_ID_Gene_Protein(corrtable,
                                                  pfdatabase)
      }
      else if (printGeneName & !pfdatabase) {
        corrtable <- ms_composite_ID_Gene(corrtable)
      }
      else {
        corrtable <- ms_composite_ID_Protein(corrtable,
                                             pfdatabase)
      }
      cdata$id <- factor(cdata$id, levels = corrtable$id)
    }
    if (length(plotseq)) {
      cdata$id <- factor(cdata$id, levels = plotseq)
    }
    else {
      cdata$id <- factor(cdata$id, levels = unique(cdata$id), ordered = TRUE)
    }
    cdata$treatment <- factor(as.character(cdata$treatment),
                              levels = treatmentlevel)
    cdata$condition <- factor(as.character(cdata$condition),
                              levels = apply(expand.grid(temperature, treatmentlevel),
                                             1, paste, collapse = "_"))
    # if data with different temperatures, prevent from creating non sense factors
    cdata$condition <- factor(as.character(cdata$condition),
                              levels = levels(cdata$condition)[levels(cdata$condition)
                                                               %in% as.character(cdata$condition)
                                                               ]
                              )

    message("Generating fitted plot, pls wait.")



    plots <- plyr::dlply(cdata, plyr::.(id), .fun = barplotting,
                         withset = withset)


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
                           dataname, "_", pdfname), pl, height = pdfheight,
             width = pdfwidth)

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
        theme_cowplot() +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank())

      return(g)
    }
  }

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

getGeneName <- function (x){
  gene = stringr::str_split(stringr::str_split(x, "GN=")[[1]][2], " ")[[1]][1]
  if (length(gene) == 0) {
    return(" ")
  }
  else {
    return(gene)
  }
}

getProteinName <- function (x, pfdatabase = FALSE)
{
  if (pfdatabase) {
    protein = gsub("gene_product=", "", stringr::str_split(x, "\\|")[[1]][4])
  }
  else {
    protein = stringr::str_split(x, " OS=")[[1]][1]
  }
  if (length(protein) == 0) {
    return(" ")
  }
  else {
    return(protein)
  }
}




