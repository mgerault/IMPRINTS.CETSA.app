#' imprints_barplotting_categorize_peptides
#'
#' Function to plot the categorized proteins found as hits by the function \code{imprints_cleaved_peptides}
#'
#'
#' @param data The normalized peptides data set, i.e. the outpout from \code{imprints_normalize_peptides}.
#' @param data_cleaved The categorized cleavage hits data set, i.e. the outpout from \code{imprints_categorize_peptides}.
#' @param treatment The treatment from which you want to see the plots. Can only be one.
#' @param control The control treatment from your dataset.
#' @param format Format of the plot; either \code{RESP_peptide} which will plot the RESP plot and then all peptides
#'   separately from the protein or \code{peptide_one} which will plot all peptides in one plot for each protein.
#'   Default is \code{RESP_peptide}.
#' @param color The color of the bar plots.
#' @param pdfname textual label of the pdf file
#'
#' @return NULL. Will save a pdf file with ordered and catgorized plots
#'
#' @seealso \code{\link{imprints_barplotting_peptides}}
#' @seealso \code{\link{imprints_categorize_peptides}}
#'
#' @export
#'
imprints_barplotting_categorize_peptides <- function(data, data_cleaved, treatment, control,
                                                     format = c("RESP_peptide", "peptide_one"), color = "red",
                                                     pdfname = "categorized_RESP_barplots"){
  format <- match.arg(format)

  # data verifications
  if(!(treatment %in% get_treat_level(data))){
    message(paste("Error: Your data doesn't contain the treatment", treatment))
    return(NULL)
  }
  if(!(control %in% get_treat_level(data))){
    message(paste("Error: Your data doesn't contain the control", control))
    return(NULL)
  }
  if(!(treatment %in% unique(data_cleaved$treatment))){
    message(paste("Error: Your data_cleaved doesn't contain the treatment", treatment))
    return(NULL)
  }
  if(control %in% unique(data_cleaved$treatment)){
    message(paste("Error: Your data_cleaved contain the treatment", control,
                  "which you selected as a control"))
    return(NULL)
  }

  # checking data
  cn_missing <- c("Master.Protein.Accessions", "description", "Positions.in.Master.Proteins",
                  "Annotated.Sequence", "Modifications", "countNum")
  cn_missing <- cn_missing[!(cn_missing %in% colnames(data))]
  if(length(cn_missing)){
    message(paste("Error: Your data miss the column(s)",
                  paste(cn_missing, collapse = ", "))
            )
    return(NULL)
  }

  # checking data_cleaved
  cn_missing <- c("id", "Gene", "description",
                  "cleaved_site", "category", "details")
  cn_missing <- cn_missing[!(cn_missing %in% colnames(data_cleaved))]
  if(length(cn_missing)){
    if(any(c("category", "details") %in% cn_missing) & length(cn_missing) <= 2){
      message("Categorization missing, performing categorization...")
      data_cleaved <- imprints_categorize_peptides(data, data_cleaved, control,
                                                   save_xlsx = FALSE)
    }
    else{
      message(paste("Error: Your data_cleaved miss the column(s)",
                    paste(cn_missing, collapse = ", "))
              )
      return(NULL)
    }
  }

  message("Computing necessary FC...")
  # filtering data
  data <- data[,c("Master.Protein.Accessions", "description", "Positions.in.Master.Proteins",
                  "Annotated.Sequence", "Modifications",
                  grep(paste0("_", treatment, "$|_", control, "$"),
                       colnames(data), value = TRUE),
                  "countNum")]
  data_cleaved <- unique(data_cleaved[data_cleaved$treatment == treatment,
                                      c("id", "Gene", "description",
                                        "cleaved_site", "category", "details")])
  # all peptides
  directory_toremove <- list.files()
  data_diff_peptides <- imprints_sequence_peptides(data, control = control,
                                                   proteins = data_cleaved$id,
                                                   dataset_name = "imprints_categorize_peptides_FC")
  directory_toremove <- list.files()[!(list.files() %in% directory_toremove)]
  directory_toremove <- grep("imprints_categorize_peptides_FC\\.txt$", directory_toremove, value = TRUE)
  if(length(directory_toremove) == 1){
    unlink(directory_toremove, recursive = TRUE)
  }
  data_diff_peptides <- data_diff_peptides[,-grep(paste0("_", control, "$"), colnames(data_diff_peptides))]

  if(format == "RESP_peptide"){
    # summarized peptides for RESP plot
    directory_toremove <- list.files()
    data_diff_summary <- imprints_sequence_peptides(data, control = control,
                                                    proteins = data_cleaved$id,
                                                    sequence = sub("~", "-", data_cleaved$cleaved_site),
                                                    dataset_name = "imprints_categorize_peptides_FC")
    directory_toremove <- list.files()[!(list.files() %in% directory_toremove)]
    directory_toremove <- grep("imprints_categorize_peptides_FC\\.txt$", directory_toremove, value = TRUE)
    if(length(directory_toremove) == 1){
      unlink(directory_toremove, recursive = TRUE)
    }
    data_diff_summary <- imprints_remove_peptides(data_diff_summary, proteins = data_cleaved$id,
                                                  sequence = sub("~", "-", data_cleaved$cleaved_site))
    data_diff_summary <- data_diff_summary[,-grep(paste0("_", control, "$"), colnames(data_diff_summary))]

  }

  message("Plotting...\n")
  # ordering data_cleaved according category and Gene
  data_cleaved <- data_cleaved[order(data_cleaved$Gene),]
  data_cleaved <- data_cleaved[order(as.numeric(factor(data_cleaved$category,
                                                       levels = c("RESP", "SP", "SPm", "MP", "MPm", "FP"))
                                                )
                                     ),]

  # actual plotting
  if(format == "RESP_peptide"){
    l <- list()
    for(p in data_cleaved$id){
      # plotting summary plot
      X <- data_diff_summary[which(data_diff_summary$Master.Protein.Accessions == p),]
      l[[paste0(p, "_summary")]] <- imprints_barplotting_peptides(X, format = RESP_peptide, colorpanel = color)
      # plotting individual peptides
      X <- data_diff_peptides[which(data_diff_peptides$Master.Protein.Accessions == p),]
      l[[paste0(p, "_peptides")]] <- imprints_barplotting_peptides(X, colorpanel = color)
    }

    for(n in names(l)){
      message(n)
      # adding category to plots
      p <- unique(sub("_.*", "", n))
      p_category <- data_cleaved$category[data_cleaved$id == p]
      p_category <- paste(paste0("<span style='font-size:16pt; color:",
                                 ifelse(p_category == "RESP", "blue",
                                        ifelse(p_category == "FP",
                                               "red",
                                               "orange")), "'>**",
                                 p_category,
                                 "**</span>"),
                          "<span style='color:black'>-", data_cleaved$details[data_cleaved$id == p], "</span>")
      if(grepl("_summary$", n)){
        l[[n]] <- gridExtra::marrangeGrob(l[[n]],
                                          nrow = 1, ncol = 1,
                                          bottom = "Summary plot",
                                          top = gridtext::richtext_grob(p_category,
                                                                        gp = grid::gpar(col = c("black", "blue", "orange", "red"),
                                                                                        fontfamily = c("Times"))
                                          )
        )
      }
      else{
        np <- 9
        pages <- length(l[[n]])%/%np + as.logical(length(l[[n]])%%np)
        message(paste("Writing", pages, "pages..."))

        l[[n]] <- gridExtra::marrangeGrob(l[[n]], layout_matrix = matrix(seq_len(9), nrow = 3, ncol = 3,
                                                                         byrow = TRUE),
                                          nrow = 3, ncol = 3,
                                          bottom = "IMPRINTS peptides plots",
                                          top = gridtext::richtext_grob(p_category,
                                                                        gp = grid::gpar(col = c("black", "blue", "orange", "red"),
                                                                                        fontfamily = c("Times"))
                                          )
        )
      }
      class(l[[n]]) <- c("arrangelist", "ggplot", class(l[[n]]))
    }

    message("\nSaving final pdf file...")
    ggpubr::ggexport(filename =  paste0(format(Sys.time(), "%y%m%d_%H%M_"), pdfname, ".pdf"),
                     plotlist = l, height = 12, width = 14)
    message("Done !")
  }
  else if(format == "peptide_one"){
    l <- list()
    for(p in data_cleaved$id){
      X <- data_diff_peptides[which(data_diff_peptides$Master.Protein.Accessions == p),]
      l[[p]] <- plyr::dlply(X, plyr::.(Master.Protein.Accessions), .fun = barplotting_peptides, color = color)
    }

    np = length(l)
    npi = 1
    for(n in names(l)){
      message(paste("Writing page", npi, "/", np))

      # adding category to plots
      p_category <- data_cleaved$category[data_cleaved$id == n]
      p_category <- paste(paste0("<span style='font-size:16pt; color:",
                                 ifelse(p_category == "RESP", "blue",
                                        ifelse(p_category == "FP",
                                               "red",
                                               "orange")), "'>**",
                                 p_category,
                                 "**</span>"),
                          "<span style='color:black'>-", data_cleaved$details[data_cleaved$id == n], "</span>")


        l[[n]] <- gridExtra::marrangeGrob(l[[n]],
                                          nrow = 1, ncol = 1,
                                          bottom = "IMPRINTS peptides plots",
                                          top = gridtext::richtext_grob(p_category,
                                                                        gp = grid::gpar(col = c("black", "blue", "orange", "red"),
                                                                                        fontfamily = c("Times"))
                                                                        )
                                          )

      class(l[[n]]) <- c("arrangelist", "ggplot", class(l[[n]]))

      npi = npi + 1
    }

    message("\nSaving final pdf file...")
    ggpubr::ggexport(filename =  paste0(format(Sys.time(), "%y%m%d_%H%M_"), pdfname, ".pdf"),
                     plotlist = l, height = 8, width = 24)
    message("Done !")
  }
}


barplotting_peptides <- function(x, color){
  x$countNum <- NULL

  # preparing data - mean and se
  x <- x %>%
    tidyr::gather("key", "value", -Master.Protein.Accessions, -description,
                  -Positions.in.Master.Proteins,
                  -Annotated.Sequence, -Modifications) %>%
    tidyr::separate(key, into = c("temperature", "rep", "treatment"), sep = "_") %>%
    group_by(Master.Protein.Accessions, description,
             Positions.in.Master.Proteins,
             Annotated.Sequence, Modifications,
             treatment, temperature) %>%
    summarise(se = sd(value, na.rm = TRUE)/sqrt(length(na.omit(value))),
              value = mean(value, na.rm = TRUE))

  # preparing title
  x$Gene <- sapply(x$description, getGeneName, USE.NAMES = FALSE)
  x$description <- sapply(x$description, getProteinName, USE.NAMES = FALSE)
  x$temperature <- factor(x$temperature)
  x$Positions.in.Master.Proteins <- sub(".* ", "", x$Positions.in.Master.Proteins)
  x$Positions.in.Master.Proteins <- gsub("\\[|\\]", "", x$Positions.in.Master.Proteins)

  # formating modification
  x$Modifications <- unlist(lapply(strsplit(x$Modifications, "\\];"),
                                    function(y){
                                      y <- y[-grep("TMT", y)]
                                      if(length(y)){
                                        y <- sub("\\d{1}x", "", y)
                                        y <- sub("^ ", "", y)
                                        y <- gsub("\\[|\\]", "", y)
                                        y <- paste(y, collapse = "\n")
                                      }
                                      else{
                                        y <- ""
                                      };
                                      y
                                    })
                            )

  # formating labels/title
  x$Positions.in.Master.Proteins <- paste0(x$Positions.in.Master.Proteins, "\n",
                                            x$Modifications)
  x$Positions.in.Master.Proteins <- factor(x$Positions.in.Master.Proteins,
                                            unique(x$Positions.in.Master.Proteins)[order(sapply(strsplit(gsub("\\[|\\]|\n.*", "",
                                                                                                               unique(x$Positions.in.Master.Proteins)
                                            ),
                                            "-"),
                                            function(y) sum(as.numeric(y)))
                                            )])

  # handling if QP in temperature
  x$temperature <- factor(x$temperature)
  x$QP <- FALSE
  lvl_temperature <- levels(x$temperature)
  if("36C" %in% x$temperature){
    x$QP[which(x$temperature == "36C")] <- TRUE
    x$temperature <- factor(sub("36C", "QP", x$temperature),
                            levels = c("QP", lvl_temperature[-grep("36C", lvl_temperature)]))
  }

  # plot
  g <- ggplot(x, aes(temperature, value, fill = treatment)) +
    geom_bar(stat = "identity", aes(color = QP)) +
    geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 0.1) +
    facet_wrap(~Positions.in.Master.Proteins, strip.position = "bottom",
               nrow = 1) +
    scale_fill_manual(values = color) +
    scale_color_manual(values = c("TRUE" = "#656565", "FALSE" = "#FFFFFF00")) +
    guides(color = "none") +
    labs(y = "fold-change (log2)", x = "peptides",
         title = paste(x$Master.Protein.Accessions[1], "\n",
                       x$description[1], "\n",
                       x$Gene[1])) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5,size = rel(1.5), face = "bold"),
          panel.grid.major.x = element_blank(),
          axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.7)),
          axis.text.y = element_text(size = rel(1.4)),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(face = "bold", angle = 90, size = rel(1.1)))

  return(g)
}

getGeneName <- function (x){
  gene = strsplit(strsplit(x, "GN=")[[1]][2], " ")[[1]][1]
  if (length(gene) == 0) {
    return(" ")
  }
  else {
    return(gene)
  }
}

getProteinName <- function (x){
  protein = strsplit(x, " OS=")[[1]][1]
  if (length(protein) == 0) {
    return(" ")
  }
  else {
    return(protein)
  }
}

