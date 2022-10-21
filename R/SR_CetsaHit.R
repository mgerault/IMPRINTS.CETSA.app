#' SR_CetsaHit
#'
#' Function to categorize protein according to their Expression/Stability change
#' based on their stability rate and their combined p-value on their two biggest fold changes.
#'
#' @details This function is still on test. Your data should contain the columns named
#'          Protein.Group, Genes and peptides_counts_all.
#'
#' @param data The input data set from which categorization is performed on and hitlist is produced from.
#' @param data_diff The output from ms_2D_caldiff; can be NULL and if so will compute it; can also be a path to the data.
#' @param ctrl The name of the control.
#' @param valid_val The percentage of valid values you want per condition. If less, score will be set to NA,
#'                  i.e. the protein will not be a hit
#' @param SR_cutoff The stability rate cutoff
#' @param pval_cutoff The combined p-value cutoff
#' @param folder_name The name of the folder in which you want to save the results.
#' @param peptide_count_col The name of the column that contain the unique peptide count.
#'                          If it is sumUniPeps, don't bother with this parameter.
#' @param species The species on which you did the experiment (not necessary if already present in your datas).
#'                Default is 'Homo Sapiens'.
#'
#' @return A dataframe which contains the hits.
#'
#' @export

SR_CetsaHit <- function(data, data_diff = NULL, ctrl, valid_val = NULL,
                        SR_cutoff = 1, pval_cutoff = 0.05,
                        folder_name = "Hits_analysis",
                        peptide_count_col = "peptides_counts_all",
                        species = "Homo Sapiens"){
  wd <- getwd()
  outdir <- paste0(wd, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", folder_name)
  if (dir.exists(outdir) == FALSE) {
    dir.create(outdir)
  }

  n <- nrow(data)
  name_cond <- stringr::str_subset(colnames(data), "^\\d{2}")
  if(ncol(data) - length(name_cond) != 5){
    stop("Your data must have 5 information columns,
          like 'id', 'description', 'genes', 'peptide_count' and 'protein_names' for example.")
  }
  cond <- stringr::str_split(name_cond, "_")
  cond <- unique(unlist(lapply(cond, function(x) x[3])))
  if(!(ctrl %in% cond)){
    stop("'ctrl' is not in the conditions. Please, check spelling or column names.")
  }
  cond <- cond[-which(cond == ctrl)]

  temp <- stringr::str_split(name_cond, "_")
  temp <- unique(unlist(lapply(temp, function(x) x[1])))

  min_validval <- 0
  if(!is.null(valid_val)){
    min_validval <- round(valid_val*length(temp))
    message(paste("Will remove protein with less than", min_validval, "valid values."))
  }

  if(is.null(data_diff)){
    data_diff <- data
    if(!("description" %in% colnames(data_diff))){
      to_describe <- ""
      while(length(to_describe) != 2 & all(!(to_describe %in% colnames(data)))){
        to_describe <- readline(prompt = "Type the column name of the description you want to keep in first and the gene column name in second; separated by 2 spaces"
        )
        to_describe <- stringr::str_split(to_describe, "  ")[[1]]

        if(length(to_describe) != 2){
          print("Separate the columns names by 2 spaces")
        }
        if(all(!(to_describe %in% colnames(data)))){
          print("Check columns names; it doesn't match")
        }
      }

      data_diff$description <- paste(data[[to_describe[1]]], "OS=Homo sapiens", paste0("GN=", data[[to_describe[2]]]))
      data_diff[[to_describe[1]]] <- NULL
    }
    if(!("sumUniPeps" %in% colnames(data_diff))){
      colnames(data_diff)[which(colnames(data_diff) == peptide_count_col)] <- "sumUniPeps"
    }
    message("Getting fold change...")
    colnames(data_diff)[-stringr::str_which(colnames(data_diff), "^\\d{1,}|description|sumUniPeps")] <- c("id", "sumPSMs", "countNum")
    data_diff <- data_diff[, c("id", "description", name_cond, "sumUniPeps", "sumPSMs", "countNum")]
    data_diff <- mineCETSA::ms_2D_caldiff(data_diff, treatmentlevel = c(ctrl, cond))
  }
  else if(class(data_diff) == "character"){
    if(stringr::str_detect(data_diff, "\\.tsv$|\\.csv$|\\.txt$"))
      data_diff <- readr::read_tsv(data_diff)
    else
      stop("Format isn't recognize")
  }
  else if(!("data.frame" %in% class(data_diff))){
    stop("data_diff is neither a file, neither a data frame !")
  }


  message("Computing mean values...")
  # get average value among bioreplicate for each protein
  diff_SR <- data_diff[,-stringr::str_which(colnames(data_diff), paste0("_", ctrl, "$"))]
  diff_SR <- tidyr::gather(diff_SR, condition, reading, -id, -description, -sumUniPeps, -sumPSMs, -countNum)
  diff_SR <- tidyr::separate(diff_SR, condition, into = c("temperature",
                                                          "replicate", "treatment"), sep = "_")
  diff_SR <- diff_SR %>% dplyr::group_by(id, description, temperature, treatment) %>%
    dplyr::summarise(reading.mean = mean(reading, na.rm = T))
  diff_SR <- tidyr::unite(diff_SR, condition, temperature, treatment,
                          sep = "_")
  diff_SR <- tidyr::spread(diff_SR, condition, reading.mean)


  message("Computing p-values...")
  pval <- list()
  for(k in cond){
    message(k)
    res <- list()
    M <- data[,stringr::str_which(colnames(data), paste0("_", ctrl, "$|_", k, "$"))]
    for(i in temp){  # 6 temperatures = 6 fractions
      message(i)
      X <- M[,stringr::str_which(colnames(M), paste0("^", i, "_"))]
      grp <- unname(sapply(colnames(X), function(x) stringr::str_split(x, "_")[[1]][3]))

      res[[i]] <- MKmisc::mod.t.test(as.matrix(X),
                                     group = factor(grp),
                                     adjust.method = "BH")$p.value
    }
    pval[[k]] <- res
  }


  message("Computing stability rate and Fisher p-value on top 2 fractions")
  diff_SR <- diff_SR[,-c(1:2)]
  for(k in cond){# get top 2 biggest diff_SR
    message(k)

    message("Get top2 fraction")
    diff_SR[[paste0("Top2_", k)]] <- apply(diff_SR[,stringr::str_which(colnames(diff_SR), paste0("_", k, "$"))],
                                           1,
                                           function(x) paste(order(as.numeric(abs(x)),
                                                                   decreasing = TRUE)[1:2],
                                                             collapse = ";")
    )


    message("Get stability rate")
    diff_SR[[paste0("SR_", k)]] <-
      apply(diff_SR[,stringr::str_which(colnames(diff_SR), paste0("_", k, "$"))], 1,
            function(d){
              top = d[length(d)]
              top <- stringr::str_split(top, ";")[[1]]
              top <- as.numeric(top)
              d = d[-length(d)]
              d = as.numeric(d)
              d = abs(d)

              crit = FALSE
              if(!is.null(valid_val)){
                if(sum(!is.na(d)) < min_validval){
                  crit = TRUE
                }
              }
              else if(sum(is.na(d)) == length(d)){
                crit = TRUE
              }

              if(crit){
                res = NA
              }
              else{
                ord = order(d, decreasing = TRUE)
                d = d[ord]
                dat = data.frame(x = 1:length(temp),
                                 y = d)
                wght = rep(0.1, length(d))
                wght[1:2] <- 0.5
                res = lm(y ~ x,
                         data = dat,
                         weights = wght
                )$coefficients[1]
              };
              res = abs(res)
            }
      )

    message("Get Fisher p-value")
    diff_SR[[paste0("pvalF1_", k)]] <- rep(1,n)
    diff_SR[[paste0("pvalF2_", k)]] <- rep(1,n)
    diff_SR[[paste0("Fisher_", k)]] <- rep(1,n)
    for(i in 1:n){
      topF <- diff_SR[[paste0("Top2_", k)]][i]
      topF <- stringr::str_split(topF, ";")[[1]]
      topF <- as.numeric(topF)

      pF1 <- pval[[k]][[temp[topF[1]]]][i]
      pF2 <- pval[[k]][[temp[topF[2]]]][i]
      if(is.na(pF1)){
        pF1 <- 1
      }
      if(is.na(pF2)){
        pF2 <- 1
      }

      diff_SR[[paste0("pvalF1_", k)]][i] <- pF1
      diff_SR[[paste0("pvalF2_", k)]][i] <- pF2

      diff_SR[[paste0("Fisher_", k)]][i] <- metap::sumlog(c(pF1, pF2))$p
    }

    diff_SR[[paste0("GlobalScore_", k)]] <- tidyr::replace_na(diff_SR[[paste0("SR_", k)]]*(-log10(diff_SR[[paste0("Fisher_", k)]])), 0)

    diff_SR[[paste0("Stabilisation_", k)]] <- apply(diff_SR[,paste0(temp, "_", k)], 1, function(x){m = mean(x, na.rm = TRUE)
    m = -1*(m < 0) + 1*(m >=0);
    m}
    )
  }
  info <- data_diff[,c("id", "description", "sumUniPeps")]
  info$description <- stringr::str_extract(paste0(info$description, " "), "(?<=GN=).+?(?= )") # extract genes
  colnames(info) <- c("id", "Genes", "peptide_count")
  diff_SR <- as.data.frame(cbind(info, diff_SR))

  diff_SR$GlobalScore <- as.numeric(apply(diff_SR[,stringr::str_which(colnames(diff_SR), "^GlobalScore_")], 1, prod))
  higher_GS <- order(diff_SR$GlobalScore, decreasing = TRUE)[1:30]


  diff_SR_plot <- diff_SR[,c("id", "Genes", stringr::str_subset(colnames(diff_SR), "^SR_|^Fisher_"))]
  diff_SR_plot <- tidyr::gather(diff_SR_plot, condition, reading, -id, -Genes)
  diff_SR_plot <- tidyr::separate(diff_SR_plot, condition, into = c("Value", "condition"))
  diff_SR_plot <- tidyr::spread(diff_SR_plot, Value, reading)
  diff_SR_plot$condition <- factor(diff_SR_plot$condition)
  diff_SR_plot$criteria <- -log10(diff_SR_plot$Fisher) >= -log10(pval_cutoff) & diff_SR_plot$SR >= SR_cutoff
  diff_SR_plot$criteria <- tidyr::replace_na(diff_SR_plot$criteria, FALSE)

  message("Creating and saving plot")
  g_h <- ggplot(diff_SR_plot, aes(SR, -log10(Fisher), color = criteria)) +
    geom_point(show.legend = FALSE) +
    geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed") +
    geom_vline(xintercept = SR_cutoff, linetype = "dashed") +
    labs(title = "Stability rate plot",
         y = "-log10(Fisher p-value)",
         x = "Stability rate") +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70")) +
    gghighlight::gghighlight(criteria,
                label_key = Genes,
                label_params = list(color = "blue"),
                use_direct_label = TRUE,
                max_highlight = 100,
                calculate_per_facet = TRUE) +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~condition)

  ggsave(paste0(format(Sys.time(), "%y%m%d_%H%M"), "_", "hits_plot.png"),
         plot = g_h,
         device = "png",
         path = outdir,
         width = 14,
         height = 8)


  g_I <- ggplot(diff_SR_plot, aes(SR, -log10(Fisher), color = criteria,
                                  group = Genes, group2 = id)) +
    geom_point(show.legend = FALSE) +
    geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed") +
    geom_vline(xintercept = SR_cutoff, linetype = "dashed") +
    labs(title = "Stability rate plot",
         y = "-log10(Fisher p-value)",
         x = "Stability rate") +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70"))  +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~condition)


  g_I <- plotly::ggplotly(g_I, width = 1080, height = 560)
  htmltools::save_html(g_I, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "hits_plotInt.html"))

  if(length(cond) == 3){
    a <- diff_SR$peptide_count
    a_c <- rep("", length(a))
    a_c[which(a <= 3)] <- 1
    a_c[which(a > 3 & a <= 10)] <- 2
    a_c[which(a > 10 & a <= 20)] <- 3
    a_c[which(a > 20 & a <= 50)] <- 4
    a_c[which(a > 50)] <- 5
    diff_SR$n_pep <- factor(a_c, ordered = TRUE)
    ## Global Score plot
    for(i in 1:3){
      GS <- ggplot(diff_SR, aes(.data[[paste0("GlobalScore_", cond[((3 + i) %% 3) + 1])]],
                                .data[[paste0("GlobalScore_", cond[((2 + i) %% 3) + 1])]],
                                color = .data[[paste0("GlobalScore_", cond[((1 + i) %% 3) + 1])]],
                                size = n_pep,
                                alpha = 0.85,
                                label = Genes,
                                text = paste("Genes :", Genes, "\n",
                                             "Protein Groups :", id, "\n",
                                             paste("Global Score", cond[3], ":"), round(.data[[paste0("GlobalScore_",
                                                                                                      cond[3])]],
                                                                                        2), "\n",
                                             paste("Global Score", cond[2], ":"), round(.data[[paste0("GlobalScore_",
                                                                                                      cond[2])]],
                                                                                        2), "\n",
                                             paste("Global Score", cond[1], ":"), round(.data[[paste0("GlobalScore_",
                                                                                                      cond[1])]],
                                                                                        2), "\n",
                                             "Number of peptides :", peptide_count)
      )
      ) +
        geom_point() +
        scale_size_manual(values = c(1, 2, 4, 6, 8),
                          labels = c("<= 3", "> 3 and <= 10", "> 10 and <= 20", "> 20 and <= 50", "> 50")) +
        ggrepel::geom_label_repel(data = diff_SR[higher_GS,],
                                  aes(color = diff_SR[higher_GS,][[paste0("GlobalScore_", cond[((1 + i) %% 3) + 1])]]),
                                  size = 5,
                                  point.padding = 0,
                                  min.segment.length = 0,
                                  box.padding = 0.3,
                                  alpha = 1,
                                  max.overlaps = 50) +
        scale_color_gradient2(low = "#0027FF", mid = "#FFD800", high = "#FF0000",
                              midpoint = 10, oob = scales::squish, limits = c(0,20)
        ) +
        labs(x = paste("Global score", cond[((3 + i) %% 3) + 1]),
             y = paste("Global score", cond[((2 + i) %% 3) + 1]),
             color = paste("Global score", cond[((1 + i) %% 3) + 1]),
             title = "Global Score plot",
             size = "NUmber of peptides") +
        theme(plot.title = element_text(hjust = 0.5),
              panel.background = element_rect(fill = "#080808",
                                              colour = "#080808",
                                              size = 0.5,
                                              linetype = "solid"),
              plot.background = element_rect(fill = "#D4F0FF"),
              legend.background = element_rect(fill = "#D4F0FF")
        ) +
        guides(alpha = "none")

      ggsave(paste0(format(Sys.time(), "%y%m%d_%H%M"), "_", "GS_plot", i, ".png"),
             plot = GS,
             device = "png",
             path = outdir,
             width = 14,
             height = 8)

      GS_I <- plotly::ggplotly(GS, tooltip = c("text"), width = 1080, height = 560)
      htmltools::save_html(GS_I, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "GS_plotInt", i, ".html"))
    }
  }

  message("Saving datas...")
  diff_SR_plot <- diff_SR_plot[diff_SR_plot$criteria,]
  diff_SR_plot$criteria = NULL
  diff_SR_plot$GlobalScore <- -log10(diff_SR_plot$Fisher)*diff_SR_plot$SR
  openxlsx::write.xlsx(diff_SR, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "hits_analysis_tab.xlsx"))
  openxlsx::write.xlsx(diff_SR_plot, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "hits_summary.xlsx"))

  if(nrow(diff_SR_plot) > 1){
    vennlist <- (diff_SR_plot %>% dplyr::group_by(condition) %>%
                   dplyr::summarize(vennlist = list(id)) %>%
                   dplyr::select(vennlist))[[1]]
    names(vennlist) <- levels(diff_SR_plot$condition)
    VennDiagram::venn.diagram(vennlist, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "VennDiagram.png"),
                              output = TRUE,
                              fill = RColorBrewer::brewer.pal(length(names(vennlist)), "Pastel2")[1:length(vennlist)],
                              fontface = "bold",
                              fontfamiliy = "sans",
                              cat.cex = 1.6,
                              cat.fontface = "bold",
                              cat.fontfamily = "sans")
    vennlist <- mineCETSAapp::com_protein_loop(vennlist)
    for(i in 1:length(vennlist)){
      prot <- data.frame("id" = vennlist[[i]],
                         "Genes" = info$Genes[which(!is.na(match(info$id,
                                                                 vennlist[[i]])))]
      )
      score_info <- diff_SR[,c(1, stringr::str_which(colnames(diff_SR), "^SR_|^Fisher_|^GlobalScore"))]

      vennlist[[i]] <- dplyr::left_join(prot, score_info, by = "id")
    }
    openxlsx::write.xlsx(vennlist, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "Venn tab.xlsx"))
  }

  return(diff_SR_plot)
}
