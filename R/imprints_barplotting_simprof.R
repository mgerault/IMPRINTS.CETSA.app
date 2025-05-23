#' imprints_barplotting_simprof
#'
#' Function to generate IMPRINTS bar plot and pdf file with multipanel bar plots for IMPRINTS-CETSA data of
#' proteins which have similar profile from a selected protein.
#'
#' @param data Dataset after imprints_caldiff to plot.
#' @param data_average Dataset after imprints_average. If null, will get it from data.
#' @param treatmentlevel A single character element which corresponds to one of the condition from the data
#' @param protein_profile A single character element which corresponds to one
#'                        of the protein from the data that you want the similar profile
#' @param setlevel a vector of set information if any, such as c("M13","M16")
#' @param printBothName A logical to tell if you want to print the both protein names on the plot
#' @param score_threshold A numeric value to indicate the threshold, default set to 0.9
#' @param use_score A single character element that define the method score. Method available : 'euclidean' or 'pearson'
#' @param max_na_prow An integer indicating the maximum number of missing value for one protein
#' @param printGeneName A logical to tell if you want to print the gene names on the plot
#' @param witherrorbar A logical to print or not the error bar on the plot
#' @param colorpanel a vector of color scheme provided by default with the function PaletteWithoutGrey
#' @param usegradient whether the barplot should be draw in color gradient format
#' @param colorgradient the color scheme of gradient applied, default value c("#4575B4","ivory", "#D73027")
#' @param linegraph whether to plot the graph in a line graph format, default set to FALSE
#' @param log2scale whether the yscales should be in log2 scale, default set to TRUE
#' @param ratio aspect ratio of the plot, default set to 0.6
#' @param ret_plot Logical to tell if you want to return the bar plot from the protein selected
#' @param withpopup A logical, only use in shiny context. If TRUE, will call a modal dialog in a shiny app;
#'                  else, it will ask a question directly in the R console (readline function)
#' @param modvar A character used when withpopup is TRUE, so in shiny context. Y or YES and the function goes on;
#'               N or NO and the function stop (doesn't save and return the similar protein profiles)
#' @param continue A logical to tell if you wanna continue and get the bar plots after getting the similar profile.
#'                 This was thought when withpopup is TRUE. If continue is FALSE will only stop after getting similar profiles.
#' @param got_it A logical to tell if you already have your data filtered. If TRUE, will directly starting to get the bar plot.
#' @param save_prlist A logical to tell if you want to save or not the protein list with a similar profile
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
#' @return The similar barplots
#'
#' @seealso \code{\link{imprints_barplotting_app}} , \code{\link{imprints_corr_to_ref}}
#'
#' @export
#'

imprints_barplotting_simprof <- function(data, data_average = NULL,
                                         treatmentlevel = "Buparlisib6h", protein_profile = "P85037",
                                         setlevel = NULL, printBothName = TRUE, printGeneName = FALSE,
                                         score_threshold = 0.9, max_na_prow = 0,
                                         use_score = "euclidean",
                                         witherrorbar = TRUE, layout = NULL,
                                         colorpanel = "#18FF00",
                                         usegradient = FALSE, colorgradient = c("#4575B4", "ivory", "#D73027"),
                                         linegraph = FALSE, log2scale = TRUE, ratio = 0.6,
                                         ret_plot = FALSE,
                                         withpopup = FALSE, continue = TRUE, modvar = "", got_it = FALSE,
                                         save_prlist = TRUE,
                                         save_pdf = TRUE, toplabel = "IMPRINTS-CETSA bar plotting",
                                         leftlabel = "", bottomlabel = "", pdfname = "barplot",
                                         pdfheight = 12, pdfwidth = 12){
  if(is.null(shiny::getDefaultReactiveDomain()) & withpopup){
    stop("withpopup is for a shiny context only. Please set it to FALSE.")
  }

  if(length(treatmentlevel) != 1){
    stop("You must select only one condition !")
  }
  else if(!(treatmentlevel %in% get_treat_level(data))){
    stop("You must select a condition present in your data !")
  }

  if(length(protein_profile) != 1){
    stop("You must select only one protein !")
  }

  if(save_pdf){
    dataname <- deparse(substitute(data))
  }

  ### function to plot IMPRINTS profiles
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
          q <- q + labs(subtitle = subt$score[n_loop]) +
            cowplot::theme_cowplot() + theme(text = element_text(size = 10),
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
                            fill = treatment)) + geom_bar(stat = "identity", aes(color = QP),
                                                          size = rel(0.85)) +
          coord_cartesian(ylim = legendscale) + scale_fill_manual(drop = FALSE,
                                                                  values = colorpanel) +
          scale_color_manual(values = c("TRUE" = "#656565", "FALSE" = "#FFFFFF00")) +
          guides(color = "none") +
          scale_x_discrete(labels = gsub("_.{1,}", "", levels(d1$condition)))
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
      q <- q + labs(subtitle = subt[as.character(unique(d1$id)), "score"]) +
        cowplot::theme_cowplot() + theme(text = element_text(size = 10),
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



  if(!got_it){
    if(is.null(data_average)){
      message("Start average calculation")
      data_ave <- imprints_average(data, savefile = TRUE)
      message("Average calculation done !")
    }
    else{
      data_ave <- data_average
    }


    target_profile <- data_ave[which(data_ave$id == protein_profile),]
    if(nrow(target_profile) == 0){
      stop("You must select a  protein present in your data")
    }
    idx_cond <- get_treat_level(target_profile)[(get_treat_level(target_profile) %in% treatmentlevel)]
    target_profile <- as.numeric(target_profile[,grep(paste0("_", idx_cond, "$"), names(target_profile))])

    if(sum(is.na(target_profile)) == length(target_profile)){
      g <- ggplot(data.frame(x = c(0,1), y = c(0,1)), aes(x,y, label = "s")) +
        geom_text(x=0.5, y=0.5, label = "The profile you selected
                                       \ncontains only missing values !", size = 6) +
        cowplot::theme_cowplot() +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank())

      return(g)
    }

    message("Getting similar profile")
    data_simi <- imprints_corr_to_ref(data_ave, treatment = treatmentlevel,
                                      reference = target_profile,
                                      use_score = use_score,
                                      score_threshold = score_threshold,
                                      max_na = max_na_prow)

    if(!inherits(data_simi, "data.frame")){
      return(data_simi)
    }

    data <- data[which(!is.na(match(data$id, data_simi$id))),]
    tr_data <- get_treat_level(data)[!(get_treat_level(data) %in% treatmentlevel)]
    tr_data <- paste0("_", tr_data, "$", collapse = "|")
    data <- data[,-grep(tr_data, names(data))]
    data_simi <- data_simi[,c("id", "score")]
    data <- dplyr::left_join(data, data_simi, by = "id")
    message("Filtering done !")
  }

  nrowdata <- nrow(data)

  if(!withpopup){
    go <- ''
    while(!(go %in% c('YES','NO','Y','N')) ){
      go <- toupper(readline(prompt =
                               paste(nrowdata, "proteins with similar profiles have been found. Do you want to continue ? (Yes/No): ")))
      if (go %in% c('YES','Y')){
        message("Let's get this profiles then !")
      }
      else if (go %in% c('NO','N')) {
        g <- ggplot(data.frame(x = c(0,1), y = c(0,1)), aes(x,y, label = "s")) +
          geom_text(x=0.5, y=0.5, label = "Try to change the threshold or
                                         \nthe score method then", size = 6) +
          cowplot::theme_cowplot() +
          theme(axis.text.x = element_blank(),
                axis.title.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text.y = element_blank(),
                axis.title.y = element_blank(),
                axis.ticks.y = element_blank())

        return(g)
      }
      else {
        message("Invalid choice")
      }
    }
  }
  else if(!continue){
    popupModal <- function() {
      modalDialog(
        HTML(paste("<h3>", nrowdata - 1, "proteins with similar profiles with a score of", score_threshold,
                        "have been found. <br>
                   You can continue by clicking on 'OK' or cancel and change the paramater.</h3>")),

        footer = tagList(
          actionButton("cancel", "Cancel"),
          actionButton("ok", "OK")
        )
      )
    }
    if (modvar %in% c('YES','Y')){
      message("Let's get this profiles then !")

      res <- imprints_barplotting_simprof(data, data_average = data_ave,
                                       treatmentlevel = treatmentlevel, protein_profile = protein_profile,
                                       setlevel = setlevel, printBothName = printBothName, printGeneName = printGeneName,
                                       score_threshold = score_threshold, max_na_prow = max_na_prow,
                                       use_score = use_score,
                                       witherrorbar = witherrorbar, layout = layout,
                                       colorpanel = colorpanel,
                                       usegradient = usegradient, colorgradient = colorgradient,
                                       linegraph = linegraph, log2scale = log2scale, ratio = ratio,
                                       ret_plot = ret_plot,
                                       withpopup = TRUE, continue = TRUE, modvar = "", got_it = TRUE,
                                       save_prlist = save_prlist,
                                       save_pdf = save_pdf, toplabel = toplabel,
                                       leftlabel = leftlabel, bottomlabel = bottomlabel, pdfname = pdfname,
                                       pdfheight = pdfheight, pdfwidth = pdfwidth)

      return(res)
    }
    else if (modvar %in% c('NO','N')) {
      g <- ggplot(data.frame(x = c(0,1), y = c(0,1)), aes(x,y, label = "s")) +
        geom_text(x=0.5, y=0.5, label = "Try to change the threshold or
                                         \nthe score method then", size = 6) +
        cowplot::theme_cowplot() +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank())

      return(g)
    }
    else{
      showModal(popupModal())
      return(data)
    }
  }


  if(continue){
    if(save_prlist){
      tab_sim <- data[,c("id", "description", "score")]
      tab_sim$Gene.name <- as.character(lapply(tab_sim$description, getGeneName))
      tab_sim$Protein.name <- as.character(lapply(tab_sim$description, function(x) getProteinName(x)))
      tab_sim$description <- NULL
      colnames(tab_sim)[1] <- "UniprotID"
      tab_sim <- tab_sim[order(tab_sim$score, decreasing = TRUE),]
      tab_sim <- tab_sim[,c("UniprotID", "Protein.name", "Gene.name", "score")]

      openxlsx::write.xlsx(tab_sim, paste0(format(Sys.time(), "%y%m%d_%H%M_"), dataname, "_ProteinList.xlsx"))
    }
    if (nrowdata == 0) {
      message("Make sure there are more than one experimental condition in dataset.")
      stop("Otherwise specify remsinglecondprot==FALSE !")
    }
    if (printBothName) {
      data <- data %>% dplyr::rowwise() %>% dplyr::mutate(description1 = getProteinName(description)) %>%
        dplyr::mutate(description2 = getGeneName(description)) %>%
        dplyr::mutate(id = paste(id, description1, description2,
                                 sep = "\n"))
      data$description1 <- NULL
      data$description2 <- NULL
    }
    else if (printGeneName) {
      data <- data %>% dplyr::rowwise() %>%
        dplyr::mutate(description = getGeneName(description)) %>%
        dplyr::mutate(id = paste(id, description, sep = "\n"))
    }
    else {
      data <- data %>% dplyr::rowwise() %>%
        dplyr::mutate(description = getProteinName(description)) %>%
        dplyr::mutate(id = paste(id, description, sep = "\n"))
    }

    #get subtitle
    subt <- data[, c(1, grep("^score", names(data)))]
    subt <- as.data.frame(subt)
    colnames(subt) <- c("id", "score")
    subt$score <- paste(use_score, "score with", protein_profile, ":", round(subt$score, 4))
    rownames(subt) <- subt$id


    data$description <- NULL
    data <- data[order(data$score, decreasing = TRUE),]
    data1 <- tidyr::gather(data[, -grep("^sumPSM|^countNum|^sumUniPeps|^drug$|^score", names(data))],
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
    cdata$id <- factor(cdata$id, levels = unique(cdata$id), ordered = TRUE)

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

    cdata_ <- cdata[-grep(paste0("^", protein_profile, "\\n"), cdata$id),]
    plots <- plyr::dlply(cdata_, plyr::.(id), .fun = barplotting,
                         withset = withset)

    subt <- NULL
    main_prof <- cdata[grep(protein_profile, cdata$id),]
    plotsmain <- plyr::dlply(main_prof, plyr::.(id), .fun = barplotting,
                             withset = withset)

    if(save_pdf){
      message("Start saving plot")

      pl <- list()

      groups <- split(seq_along(plotsmain), gl(1, 1, 1))
      pl[["main"]] <- lapply(names(groups), function(i) {
        do.call(gridExtra::arrangeGrob,
                c(plotsmain[groups[[i]]], list(nrow = 1, ncol = 1),
                  top = toplabel, left = leftlabel,
                  bottom = bottomlabel)
                )
      })


      params <- list(nrow = layout[1], ncol = layout[2])
      n <- with(params, nrow * ncol)
      pages <- length(plots)%/%n + as.logical(length(plots)%%n)
      groups <- split(seq_along(plots), gl(pages, n, length(plots)))
      n_p <- length(names(groups))
      pl[["all"]] <- lapply(names(groups), function(i) {
        message(paste("Saving page", i, "/", n_p))
        do.call(gridExtra::arrangeGrob,
                c(plots[groups[[i]]], params,
                  top = toplabel, left = leftlabel,
                  bottom = bottomlabel)
                )
      })

      class(pl[["main"]]) <- c("arrangelist", "ggplot", class(pl[["main"]]))
      class(pl[["all"]]) <- c("arrangelist", "ggplot", class(pl[["all"]]))
      class(pl) <- c("arrangelist", "ggplot", class(pl))
      pdfname <- paste0(pdfname, ".pdf")

      message("Saving final pdf file")
      ggpubr::ggexport(filename = paste0(format(Sys.time(), "%y%m%d_%H%M_"), dataname, "_", pdfname),
                       plotlist = pl,
                       height = pdfheight, width = pdfwidth)


      message("IMPRINTS-CETSA bar plot file generated successfully.")
    }


    if(ret_plot){
      return(plotsmain)
    }
    else{
      g <- ggplot(data.frame(x = c(0,1), y = c(0,1)), aes(x,y, label = "s")) +
        geom_text(x=0.5, y=0.5, label = "All the barplots have been saved succesfully !
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

}
