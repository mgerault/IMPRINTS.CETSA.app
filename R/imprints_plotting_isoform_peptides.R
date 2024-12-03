#' imprints_plotting_isoform_peptides
#'
#' Function to visualize the isoform sequence alignment with its canonical form obtained from
#' the function \code{\link{imprints_isoform_peptides}} to compare with the cleavage site obtained
#' by \code{\link{imprints_cleaved_peptides}}.
#'
#' @details
#' The plot will show you how each peptides aligns with the canonical sequence. Each peptide will be clustered
#' in two groups, the one with greatest log2 fold-change and the ones with the lowest.
#' This has been thought for a final control of the results from \code{\link{imprints_isoform_peptides}} and
#' \code{\link{imprints_cleaved_peptides}} in order to see if a protein having a significant RESP effect could
#' be due to the expression of an alternative splicing form.
#'
#' @param data The normalized peptides data set, i.e. the outpout from \code{imprints_normalize_peptides}.
#' @param data_isoform The cleavage hits data aligned with the potential isoforms, i.e. the outpout from \code{imprints_isoform_peptides}.
#' @param control The control treatment from your dataset.
#' @param treatment The treatment from which you want to see the peptides fold-changes. Can only be one.
#' @param min_ValidValue The minimum proportion of non-missing values per peptides.
#'                       Default is 0.4; so if 6 temperatures need at least 3 non missing values.
#' @param ret_plot Logical to tell if you want to return the last plot
#' @param save_pdf A logical to tell if you want to save plots in a pdf file
#' @param layout A vector indicating the panel layout for multi-panel plots per page,
#'               default value is c(2,2), use when save_pdf = TRUE
#' @param toplabel textual label at the top part of the page
#' @param leftlabel textual label at the left side of the page
#' @param bottomlabel textual label at the bottom part of the page
#' @param pdfname Character to label the pdf file
#' @param pdfheight Numeric indicating the height of pdf file, default value 12
#' @param pdfwidth Numeric indicating the width of pdf file, default value 12
#'
#'
#' @return The plots of the isoform sequence alignment
#'
#' @seealso \code{\link{imprints_isoform_peptides}}
#'
#' @export
#'

imprints_plotting_isoform_peptides <- function(data, data_isoform, control, treatment,
                                               min_ValidValue = 0.4,
                                               ret_plot = TRUE, save_pdf = FALSE, layout = c(1,2),
                                               toplabel = "", leftlabel = "", bottomlabel = "",
                                               pdfname = "isoform_alignement_plots",
                                               pdfheight = 9, pdfwidth = 18){
  ### checking inputs
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

  needed_columns <- c("accession", "description", "Gene", "treatment", "cleaved_site",
                      "isoforms", "canonical_posalign", "length_canonical")
  missing_columns <- needed_columns[!(needed_columns %in% colnames(data_isoform))]
  if(length(missing_columns)){
    message(paste("Error:",
                  paste(missing_columns, collapse = ", "),
                  ifelse(length(missing_columns) > 1, "are", "is"),
                  "missing in your data_isoform !\nAre you sure you selected the output from imprints_isoform_peptides ?"))
    return()
  }

  if(!(control %in% get_treat_level(data))){
    message(paste("Error:", control, "is not in your data !"))
    return()
  }
  if(!(treatment %in% get_treat_level(data))){
    message(paste("Error:", treatment, "is not in your data !"))
    return()
  }

  if(control %in% unique(data_isoform$treatment)){
    message(paste("Error:", control, "shouldn't be in your data_isoform !"))
    return()
  }
  if(!(treatment %in% unique(data_isoform$treatment))){
    message(paste("Error:", treatment, "is not in your data_isoform !"))
    return()
  }

  missing_proteins <- unique(data_isoform$accession)[!(unique(data_isoform$accession) %in% unique(data$Master.Protein.Accessions))]
  if(length(missing_proteins)){
    message(paste("Error:",
                  paste(missing_proteins, collapse = ", "),
                  ifelse(length(missing_proteins) > 1, "are", "is"),
                  "missing in your data!\nAre you sure you selected the right dataset ?"))
    return()
  }

  ### filtering treatments and proteins
  data_isoform <- data_isoform[which(data_isoform$treatment == treatment),]
  treat_torm <- get_treat_level(data)[!(get_treat_level(data) %in% c(control, treatment))]
  data <- data[which(!is.na(match(data$Master.Protein.Accessions, unique(data_isoform$accession)))),
               -grep(paste0("_", treat_torm, "$", collapse = "|"), colnames(data))]

  ### computing FC
  directory_toremove <- list.files()
  data <- imprints_sequence_peptides(data, control = control)
  directory_toremove <- list.files()[!(list.files() %in% directory_toremove)]
  # remove this file
  directory_toremove <- grep("(?=^\\d{6}_\\d{4}_.*)(?=.*\\.txt$)", directory_toremove, value = TRUE, perl = TRUE)
  if(length(directory_toremove) == 1){
    unlink(directory_toremove, recursive = TRUE)
  }
  data <- data[,-grep(paste0("_", control, "$"), colnames(data))]


  ### filter data and compute max FC
  message("Preparing data")
  data$Positions.in.Master.Proteins <- gsub(".* \\[|\\]", "",
                                                    sub(";.*", "", data$Positions.in.Master.Proteins)
                                            )
  data$Annotated.Sequence <- NULL
  data$countNum <- NULL
  data$Gene <- sapply(data$description, getGeneName, USE.NAMES = FALSE)
  data$description <- sapply(data$description, getProteinName, USE.NAMES = FALSE)
  colnames(data)[c(1,3)] <- c("id", "sequence")

  data <- data %>%
    tidyr::gather("treatment", "value", -id, -sequence, -description, -Gene, -Modifications) %>%
    tidyr::separate(treatment, into = c("temp", "biorep", "treatment"), sep = "_") %>%
    dplyr::group_by(id, description, Gene, sequence, Modifications, temp, treatment, biorep) %>% # if peptide modified, mean values from same peptide sequence
    dplyr::summarise(value = mean(value, na.rm = TRUE)) %>%
    dplyr::ungroup() %>% dplyr::group_by(id, description, Gene, sequence, Modifications, temp, treatment) %>%
    dplyr::summarise(mean_value = mean(value, na.rm = TRUE)) %>%
    dplyr::ungroup() %>% dplyr::group_by(id, description, Gene, sequence, Modifications, treatment) %>%
    # filter out peptides with less than min_ValidValue
    dplyr::filter(length(na.omit(mean_value))/length(mean_value) >= min_ValidValue) %>%
    # compute greatest FC from all temperatures
    dplyr::reframe(maxFC = mean_value[which.max(abs(mean_value))]) %>%
    # getting back start and end of peptide sequence
    mutate(sequence_start = as.numeric(sub("-.*", "", sequence)),
           sequence_end = as.numeric(sub(".*-", "", sequence))) %>%
    select(-sequence) %>%
    # cluster peptides as in imprints_categorize_peptides but not on absolute value
    dplyr::ungroup() %>% dplyr::group_by(id, description, Gene, treatment) %>%
    dplyr::group_modify(~ {
      h <- data.frame(comp = cutree(hclust(dist(.x$maxFC), method = "complete"), 2),
                      ward = cutree(hclust(dist(.x$maxFC), method = "ward.D"), 2),
                      ward2 = cutree(hclust(dist(.x$maxFC), method = "ward.D2"), 2),
                      ave = cutree(hclust(dist(.x$maxFC), method = "average"), 2))
      h_name <- apply(h, 2,
                      function(y){
                        sapply(1:2, function(z) mean(.x$maxFC[y == z]))
                      })
      if(!all(apply(h_name, 2, which.max) == 1)){
        h[,which(apply(h_name, 2, which.max) == 2)] <-
          apply(h[,which(apply(h_name, 2, which.max) == 2),drop=FALSE], 2,
                function(y){
                  y[y == 1] <- 3
                  y[y == 2] <- 1
                  y[y == 3] <- 2;
                  y
                })
      }
      h <- apply(h, 1, mean)

      if(any(h == 1.5)){ # ambiguous peptides will be in cluster named 3
        p_ambi <- which(h == 1.5)
        h[p_ambi] <- 3
      }
      h <- round(h)

      .x$cluster <- h

      return(.x)
    })

  ### join isoform alignment information
  colnames(data_isoform)[grep("^accession$", colnames(data_isoform))] <- "id"
  data <- inner_join(data, data_isoform[,c("id", "description", "Gene", "treatment", "cleaved_site",
                                        "isoforms", "canonical_posalign", "length_canonical")],
                     by = c("id", "description", "Gene", "treatment"),
                     relationship = "many-to-many")

  ### set variable for plot
  data$cluster <- factor(data$cluster, levels = c(1,2,3))
  data$y_start <- 1.1
  data$y_end <- 1.2

  # plotting
  plots <- plyr::dlply(data, plyr::.(isoforms), .fun = isoform_seqalign_plotting)

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
    message("Saving final pdf file")
    ggsave(file = paste0(format(Sys.time(), "%y%m%d_%H%M_"), pdfname),
           pl, height = pdfheight, width = pdfwidth)

    message("Sequence alignment plots generated successfully.")
  }

  # handling output of function
  if(ret_plot){
    message("Sequence alignment plots generated successfully.")
    return(plots)
  }
  else{
    g <- ggplot(data.frame(x = c(0,1), y = c(0,1)), aes(x,y, label = "s")) +
      geom_text(x=0.5, y=0.5,
                label = "All the sequence alignment plots\nhave been saved succesfully !\nGo check your files",
                size = 6) +
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


### actual plotting function
isoform_seqalign_plotting <- function(data){
  message(paste("Preparing plot for", data$isoforms[1]))

  ### Define Data for the whole canonical sequence and its alignment with the selected isoform
  positions <- unique(c(1,
                        sort(unlist(sapply(strsplit(strsplit(data$canonical_posalign, "; ")[[1]],
                                                    "-"), as.numeric,
                                           USE.NAMES = FALSE, simplify = FALSE))
                        ),
                        data$length_canonical))

  positions_mapping <- ifelse(sapply(paste0(positions[-length(positions)], "-", positions[-1]),
                                     function(x){
                                       align <- strsplit(strsplit(data$canonical_posalign, "; ")[[1]],
                                                         "-")
                                       align <- lapply(align, as.numeric)

                                       x <- as.numeric(strsplit(x, "-")[[1]])
                                       any(sapply(align, function(y) y[1] == x[1] & y[2] == x[2], USE.NAMES = FALSE))
                                     }, USE.NAMES = FALSE),
                              "aligned", "canonical")

  whole_sequence <- data.frame(
    x = positions[-length(positions)],  # x-coordinates (start and end points)
    xend = positions[-1],
    y = 0.9,                        # y-coordinates (fixed radius for the circle)
    yend = 1,
    type = factor(positions_mapping, c("canonical", "aligned"))
  )

  ## saving cleavage site
  cleavage_site <- as.numeric(strsplit(sub(";.*", "", data$cleaved_site[1]), # if protein group, only taking the first portein
                                       "-|~")[[1]])


  g <- ggplot() +
    # draw canonical + aligned sequence
    geom_rect(data = whole_sequence,
              aes(xmin = x, ymin = y, xmax = xend, ymax = yend, fill = type),
              color = "black", linewidth = 0.4) +
    scale_fill_manual(values = c("canonical" = "grey25", "aligned" = "#365A8FCC"),
                      labels = c("Sequence unique\nto canonical",
                                 "Isoform sequence\naligned with canonical"),
                      name = "") +
    guides(fill = guide_legend(order = 1)) +
    ggnewscale::new_scale_fill() + # define new fill scale to have color gradient according fold-change
    # plot aligned peptides
    geom_rect(aes(xmin = sequence_start, xmax = sequence_end,
                  ymin = y_start, ymax = y_end, fill = maxFC, color = cluster),
              data = data, linewidth = 1)


  ### draw links between peptides --> first with clustered in group 1
  if(sum(data$cluster == 1) > 1){
    grp1 <- which(data$cluster == 1)
    grp1 <- expand.grid(grp1, grp1)
    grp1 <- grp1[-which(apply(grp1, 1, function(y) any(duplicated(y)))),]
    grp1 <- grp1[-which(duplicated(t(apply(grp1, 1, sort))))]

    for(p in 1:nrow(grp1)){
      # Define Links between peptides of 'same cluster'
      # (need parabol link since geom_curve doesn't work in coord_polar)
      p1_med_pos <- median(c(data$sequence_start[grp1[p, 1]],
                             data$sequence_end[grp1[p, 1]])
                           )
      p2_med_pos <- median(c(data$sequence_start[grp1[p, 2]],
                             data$sequence_end[grp1[p, 2]])
                           )

      p_med_pos <- sort(c(p1_med_pos, p2_med_pos))

      peptides_cluster_link <- data.frame(
        x = seq(p_med_pos[1], p_med_pos[2], 0.01),
        # parabol equation linking two points
        y= ((1.1-0.3)/(diff(p_med_pos)/2)**2)*(seq(-diff(p_med_pos)/2, diff(p_med_pos)/2, 0.01))**2 + 0.3,
        cluster = factor(1, levels = levels(data$cluster))
      )


      g <- g +
        # draw links between peptides of same cluster
        geom_line(data = peptides_cluster_link,
                  aes(x = x, y = y, color = cluster), linewidth = 0.6,
                  show.legend = FALSE, alpha = 0.8)
    }
  }

  ### then in group 2
  if(sum(data$cluster == 2) > 1){
    grp2 <- which(data$cluster == 2)
    grp2 <- expand.grid(grp2, grp2)
    grp2 <- grp2[-which(apply(grp2, 1, function(y) any(duplicated(y)))),]
    grp2 <- grp2[-which(duplicated(t(apply(grp2, 1, sort)))),]

    for(p in 1:nrow(grp2)){
      # Define Links between peptides of 'same cluster'
      # (need parabol link since geom_curve doesn't work in coord_polar)
      p1_med_pos <- median(c(data$sequence_start[grp2[p, 1]],
                             data$sequence_end[grp2[p, 1]])
                           )
      p2_med_pos <- median(c(data$sequence_start[grp2[p, 2]],
                             data$sequence_end[grp2[p, 2]])
                           )

      p_med_pos <- sort(c(p1_med_pos, p2_med_pos))

      peptides_cluster_link <- data.frame(
        x = seq(p_med_pos[1], p_med_pos[2], 0.01),
        # parabol equation linking two points
        y= ((1.1-0.3)/(diff(p_med_pos)/2)**2)*(seq(-diff(p_med_pos)/2, diff(p_med_pos)/2, 0.01))**2 + 0.3,
        cluster = factor(2, levels = levels(data$cluster))
      )


      g <- g +
        # draw links between peptides of same cluster
        geom_line(data = peptides_cluster_link,
                  aes(x = x, y = y, color = cluster), linewidth = 0.6,
                  show.legend = FALSE, alpha = 0.8)
    }
  }

  g <- g +
    # adding labels
    geom_text(data = data.frame(x = c(-2, positions[length(positions)] + 9,  # start and end
                                      median(cleavage_site),  # middle of cleavage site
                                      cleavage_site[1], cleavage_site[2], # bounds of cleavage site
                                      positions[-c(1,length(positions))]
                                      ), # where sequence aligns
                                y = c(0.95,0.95,
                                      0.7,
                                      1.1, 1.1,
                                      rep(0.8, length(positions) - 2)),
                                label = c("1", as.character(positions[length(positions)]),
                                          "cleavage\nsite",
                                          as.character(cleavage_site),
                                          as.character(positions[-c(1,length(positions))])
                                          )
                                ),
              aes(x = x, y = y, label = label),
              color = "black", size = 5) +

    # adding cleavage site
    geom_line(data = data.frame(x = cleavage_site[1], y = seq(0.85,1.05, 0.01)),
              aes(x, y),
              linetype = "dashed") +
    geom_line(data = data.frame(x = cleavage_site[2], y = seq(0.85,1.05, 0.01)),
              aes(x, y),
              linetype = "dashed") +
    geom_line(data = data.frame(x = seq(cleavage_site[1], cleavage_site[2], 0.01), y = 0.85),
              aes(x, y),
              linetype = "dashed") +
    geom_line(data = data.frame(x = seq(cleavage_site[1], cleavage_site[2], 0.01), y = 1.05),
              aes(x, y),
              linetype = "dashed") +

    # set colors and aestethic
    scale_fill_gradientn(breaks = c(-2,-1.5,seq(-1,1,0.25), 1.5, 2),
                         colors = c("#4A86FF", "#17B5FD", "#00E49D", "#40FF5B", "#10B527", "#00570C",
                                    "black",
                                    "#880000", "#F01414", "#FF4646", "#E21257", "#FF2EA1", "#C456FF"),
                         values = scales::rescale(c(-2,-1.5,seq(-1,1,0.25), 1.5, 2)),
                         labels = c("-2", "", "-1", "", "-0.5", "-0.25", "0", "0.25", "0.5", "", "1", "", "2"),
                         limits = c(-2,2), guide = "colorbar") +
    scale_color_manual(values = c("1" = "#5C42FF", "2" = "#1284C0", "3" = "grey40"),
                       labels = c("Cluster with \ngreatest FC", "Cluster with \nlowest FC", "Ambiguous peptides"),
                       name = "Cluster") +
    ylim(c(0,1.2)) + xlim(c(-5, positions[length(positions)]*(1+1/9))) +
    labs(title = paste0("Peptide sequence alignments of ", data$Gene[1],
                        "\n",data$id[1], " ", data$description[1],
                        "\nwith its isoform ", data$isoforms[1]),
         fill = "Greatest log2\nfold-change") +
    guides(color = guide_legend(override.aes = list(fill = "transparent"), order = 3),
           fill = guide_colorbar(barheight = unit(0.75, "npc"))) +
    coord_polar() + # Set up polar coordinates to draw circles
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "right",
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 12, face = "bold"))


  return(g)
}


### Functions to extract gene and protein name from description
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
