#' ms_2D_barplotting_sh
#'
#' Function to generate 2D bar plot and pdf file with multipanel bar plots for 2D-CETSA data.
#' It is totally based on the function ms_2D_barplotting from the mineCETSA package.
#' It was retaken for easier use in the shiny app.
#'
#' @param data dataset after ms_2D_caldiff to plot
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
#' PI3K1h6h_wVeh <- PI3K1h6h_file[,-grep("Vehicle",names(PI3K1h6h_file))]
#' P85037_1h6h <- ms_subsetting(PI3K1h6h_wVeh, isfile = FALSE, hitidlist = c("P85037"))
#' ms_2D_barplotting_sh(P85037_1h6h)
#'
#' @seealso \code{\link{ms_2D_barplotting}}
#'
#' @export
#'

ms_2D_barplotting_sh <- function (data, treatmentlevel = get_treat_level(data), setlevel = NULL, corrtable = NULL,
                                  plotseq = NULL, printBothName = TRUE, printGeneName = FALSE,
                                  pfdatabase = FALSE, witherrorbar = TRUE, layout = NULL,
                                  colorpanel = PaletteWithoutGrey(treatmentlevel),
                                  usegradient = FALSE, colorgradient = c("#4575B4", "ivory", "#D73027"),
                                  linegraph = FALSE, log2scale = TRUE, ratio = 0.6,
                                  save_pdf = FALSE, toplabel = "IMPRINTS-CETSA bar plotting",
                                  leftlabel = "", bottomlabel = "", pdfname = "bar_ggplotting",
                                  pdfheight = 12, pdfwidth = 12)
{

  if(save_pdf){
    dataname <- deparse(substitute(data))
    outdir <- mineCETSA:::ms_directory(data, dataname)$outdir
  }


  nrowdata <- nrow(data)
  if (nrowdata == 0) {
    message("Make sure there are more than one experimental condition in dataset.")
    stop("Otherwise specify remsinglecondprot==FALSE !")
  }
  if (printBothName & !pfdatabase) {
    data <- data %>% rowwise() %>% mutate(description1 = mineCETSA:::getProteinName(description,
                                                                        pfdatabase)) %>%
      mutate(description2 = mineCETSA:::getGeneName(description)) %>%
      mutate(id = paste(id, description1, description2,
                        sep = "\n"))
    data$description1 <- NULL
    data$description2 <- NULL
  }
  else if (printGeneName & !pfdatabase) {
    data <- data %>% rowwise() %>%
      mutate(description = getGeneName(description)) %>%
      mutate(id = paste(id, description, sep = "\n"))
  }
  else {
    data <- data %>% rowwise() %>%
      mutate(description = getProteinName(description, pfdatabase)) %>%
      mutate(id = paste(id, description, sep = "\n"))
  }
  data$description <- NULL
  data1 <- tidyr::gather(data[, -str_which(names(data), "^sumPSM|^countNum|^sumUniPeps|^drug$")],
                         condition, reading, -id)
  if (!log2scale) {
    data1 <- mutate(data1, reading = 2^reading)
  }
  a <- data1$condition[1]
  if (length(unlist(strsplit(a, "_"))) == 4) {
    withset <- TRUE
    data1 <- tidyr::separate(data1, condition, into = c("set",
                                                        "temperature", "replicate", "treatment"), sep = "_")
    temperature <- sort(unique(data1$temperature))
    cdata <- plyr::ddply(data1, c("id", "set", "temperature",
                                  "treatment"),
                         summarise, N = length(na.omit(reading)),
                         mean = mean(reading, na.rm = T), sd = sd(reading, na.rm = T), se = sd/sqrt(N))
    if (length(layout) == 0) {
      layout <- c(2, 3)
    }

  }
  else if (length(unlist(strsplit(a, "_"))) == 3) {
    withset <- FALSE
    data1 <- tidyr::separate(data1, condition, into = c("temperature",
                                                        "replicate", "treatment"), sep = "_")
    temperature <- sort(unique(data1$temperature))
    cdata <- plyr::ddply(data1, c("id", "temperature", "treatment"),
                         summarise, N = length(na.omit(reading)), mean = mean(reading,na.rm = T),
                         sd = sd(reading, na.rm = T), se = sd/sqrt(N))
    if (length(layout) == 0) {
      layout <- c(4, 3)
    }

  }
  else {
    stop("make sure the namings of the columns of the dasaset are correct.")
  }
  cdata <- cdata %>% rowwise() %>% mutate(condition = paste(temperature,
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
    cdata$id <- factor(cdata$id)
  }
  cdata$treatment <- factor(as.character(cdata$treatment),
                            levels = treatmentlevel)
  cdata$condition <- factor(as.character(cdata$condition),
                            levels = apply(expand.grid(temperature, treatmentlevel),
                                           1, paste, collapse = "_"))
  message("Generating fitted plot, pls wait.")
  mineCETSA:::external_graphs(F)
  barplotting <- function(d1, withset = FALSE) {
    if (withset) {
      d1 <- droplevels(d1)
      d1_list <- split(d1, d1$set)
      q_list <- list()
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
          q <- q + theme_cowplot() + theme(text = element_text(size = 10),
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
      if (linegraph) {
        colorpanel <- PaletteWithoutGrey(temperature)
        q <- ggplot(d1, aes(x = treatment, y = mean,
                            group = temperature, color = temperature)) +
          geom_line() + geom_point() + coord_cartesian(ylim = legendscale) +
          scale_color_manual(drop = FALSE, values = colorpanel)
      }
      else if (!usegradient) {
        q <- ggplot(d1, aes(x = condition, y = mean,
                            fill = treatment)) + geom_bar(stat = "identity") +
          coord_cartesian(ylim = legendscale) + scale_fill_manual(drop = FALSE,
                                                                  values = colorpanel)
      }
      else {
        q <- ggplot(d1, aes(x = condition, y = mean,
                            fill = mean)) + geom_bar(stat = "identity") +
          coord_cartesian(ylim = legendscale) +
          scale_fill_gradient2(limits = legendscale,
                               low = colorgradient[1], mid = colorgradient[2],
                               high = colorgradient[3], midpoint = 0, na.value = "gray90",
                               guide = guide_colorbar(""))
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
      q <- q + theme_cowplot() + theme(text = element_text(size = 10),
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
                                                                  size = rel(0.7)), aspect.ratio = ratio)
      return(q)
    }
  }
  plots <- plyr::dlply(cdata, plyr::.(id), .fun = barplotting,
                       withset = withset)


  if(save_pdf){
    message("Start saving plot")

    params <- list(nrow = layout[1], ncol = layout[2])
    n <- with(params, nrow * ncol)
    pages <- length(plots)%/%n + as.logical(length(plots)%%n)
    groups <- split(seq_along(plots), gl(pages, n, length(plots)))
    pl <- lapply(names(groups), function(i) {
      gridExtra::grid.arrange(do.call(gridExtra::arrangeGrob,
                                      c(plots[groups[[i]]], params, top = toplabel, left = leftlabel,
                                        bottom = bottomlabel)))
    })
    class(pl) <- c("arrangelist", "ggplot", class(pl))
    pdfname <- paste0("/", pdfname, ".pdf")
    if (length(outdir)) {
      ggsave(file = paste0(outdir, "/", format(Sys.time(),
                                               "%y%m%d_%H%M_"), dataname, "_", pdfname), pl, height = pdfheight,
             width = pdfwidth)
    }
    else {
      ggsave(file = paste0(format(Sys.time(), "%y%m%d_%H%M_"),
                           dataname, "_", pdfname), pl, height = pdfheight,
             width = pdfwidth)
    }

    message("IMPRINTS-CETSA bar plot file generated successfully.")
  }

 return(plots)

}



### PaletteWithoutGrey function ###
#generates a color list depending on the number of different treatment
PaletteWithoutGrey <- function(treatment){

  n = length(unique(treatment))
  x <- grDevices::colors(distinct = TRUE)                           #all the color from R
  mycol <- x[which(is.na(stringr::str_extract(x, "gr(e|a)y")))]   #keep only colors that are not grey

  listcolor <- c()
  for (i in 0:(n-1))
    listcolor <- append(listcolor, mycol[i*20 + 9])      #save a color from the list (the number 20 and 9 were chosen in order to have distincts colors, this is empirical, can be changed)

  return(listcolor)
}


















