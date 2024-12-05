#' imprints_IS
#'
#' Function to categorize protein according to their Expression/Stability change
#' based on their I-score and their combined p-value on their two biggest fold changes.
#'
#' @param data The input data set from which categorization is performed on and hitlist is produced from.
#' @param data_diff The output from imprints_caldiff; can be NULL and if so will compute it; can also be a path to the data.
#' @param ctrl The name of the control.
#' @param valid_val The percentage of non-missing values you want per treatment. If less, score will be set to NA,
#'                  i.e. the protein will not be a hit
#' @param IS_cutoff The I-score cutoff. Default is 1.5.
#' @param fixed_score_cutoff Logical to tell if you want to use a fixed cutoff for the I-score.
#'   If TRUE, the value IS_cutoff will directly be used as the cutoff and for all treatments. If FALSE,
#'   the I-score cutoff will be calculated as the value selected for IS_cutoff plus the median of the
#'   I-scores of the proteins which have a p-value lower than the median of all p-values for a given treatment.
#'   Default is FALSE.
#' @param curvature The curvature used for the curve on the volcano plot
#' @param FDR The FDR used for the BH corrected combined p-value
#' @param FDR_category The FDR used for the BH corrected  p-value at 37°C used in order to categorize the hits
#' @param folder_name The name of the folder in which you want to save the results.
#' @param peptide_count_col The name of the column that contain the unique peptide count.
#'                          If it is sumUniPeps, don't bother with this parameter.
#' @param species The species on which you did the experiment (not necessary if already present in your datas).
#'                Default is 'Homo Sapiens'.
#'
#' @return A dataframe which contains the hits.
#'
#' @export

imprints_IS <- function(data, data_diff = NULL, ctrl, valid_val = NULL,
                        IS_cutoff = 1.5, fixed_score_cutoff = FALSE,
                        FDR = 0.01, curvature = 0.05,
                        FDR_category = 0.1,
                        folder_name = "Hits_analysis",
                        peptide_count_col = "sumUniPeps",
                        species = "Homo Sapiens"){
  wd <- getwd()
  outdir <- paste0(wd, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", folder_name)
  if (dir.exists(outdir) == FALSE) {
    dir.create(outdir)
  }

  QP_column <- grep("^36C_", colnames(data))
  if(length(QP_column)){
    data <- data[,-QP_column]
  }

  data <- data[order(data$id),] # to prevent bad matching between data and data_diff
  n <- nrow(data)
  name_cond <- grep("^\\d{2}", colnames(data), value = TRUE)
  if(ncol(data) - length(name_cond) != 5){
    stop("Your data must have 5 information columns,
          like 'id', 'description', 'genes', 'peptide_count' and 'protein_names' for example.")
  }
  cond <- strsplit(name_cond, "_")
  cond <- unique(unlist(lapply(cond, function(x) x[3])))
  if(!(ctrl %in% cond)){
    stop("'ctrl' is not in the treatments. Please, check spelling or column names.")
  }
  cond <- cond[-which(cond == ctrl)]

  temp <- strsplit(name_cond, "_")
  temp <- unique(unlist(lapply(temp, function(x) x[1])))

  min_validval <- 0
  if(!is.null(valid_val)){
    if(valid_val == 0 | valid_val > 1){
      valid_val <- NULL
    }
    else{
      min_validval <- round(valid_val*length(temp))
      message(paste("Will remove protein with less than", min_validval, "valid values."))
    }
  }

  if(is.null(data_diff)){
    data_diff <- data
    if(!("description" %in% colnames(data_diff))){
      to_describe <- ""
      while(length(to_describe) != 2 & all(!(to_describe %in% colnames(data)))){
        to_describe <- readline(prompt = "Type the column name of the description you want to keep in first and the gene column name in second; separated by 2 spaces"
        )
        to_describe <- strsplit(to_describe, "  ")[[1]]

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
    colnames(data_diff)[-grep("^\\d{1,}|description|sumUniPeps", colnames(data_diff))] <- c("id", "sumPSMs", "countNum")
    data_diff <- data_diff[, c("id", "description", name_cond, "sumUniPeps", "sumPSMs", "countNum")]
    data_diff <- IMPRINTS.CETSA::imprints_caldiff_f(data_diff, reftreatment = ctrl)
  }
  else if(inherits(data_diff, "character")){
    if(grepl("\\.tsv$|\\.csv$|\\.txt$", data_diff))
      data_diff <- readr::read_tsv(data_diff)
    else
      stop("Format isn't recognize")
  }
  else if(!inherits(data_diff, "data.frame")){
    stop("data_diff is neither a file, neither a data frame !")
  }
  data_diff <- data_diff[order(data_diff$id),]


  QP_column <- grep("^36C_", colnames(data_diff))
  if(length(QP_column)){
    data_diff <- data_diff[,-QP_column]
  }

  message("Computing mean values...")
  # get average value among bioreplicate for each protein
  diff_IS <- data_diff[,-grep(paste0("_", ctrl, "$"), colnames(data_diff))]
  diff_IS <- tidyr::gather(diff_IS, treatment, reading, -id, -description, -sumUniPeps, -sumPSMs, -countNum)
  diff_IS <- tidyr::separate(diff_IS, treatment, into = c("temperature",
                                                          "replicate", "treatment"), sep = "_")
  diff_IS <- diff_IS %>% dplyr::group_by(id, description, temperature, treatment) %>%
    dplyr::reframe(reading.mean = mean(reading, na.rm = T))
  diff_IS <- tidyr::unite(diff_IS, treatment, temperature, treatment,
                          sep = "_")
  diff_IS <- tidyr::spread(diff_IS, treatment, reading.mean)


  message("Computing p-values...")
  pval <- list()
  for(k in cond){
    message(k)
    res <- list()
    M <- data[,grep(paste0("_", ctrl, "$|_", k, "$"), colnames(data))]
    for(i in temp){  # 6 temperatures = 6 fractions
      message(i)
      X <- M[,grep(paste0("^", i, "_"), colnames(M))]
      grp <- unname(sapply(colnames(X), function(x) strsplit(x, "_")[[1]][3]))

      res[[i]] <- MKmisc::mod.t.test(as.matrix(X),
                                     group = factor(grp),
                                     adjust.method = "BH")$p.value
    }
    pval[[k]] <- res
  }


  message("Computing I-score and Fisher p-value on top 2 fractions")
  diff_IS <- diff_IS[,-c(1:2)]
  for(k in cond){# get top 2 biggest diff_IS
    message(k)

    message("Getting top2 fraction")
    diff_IS[[paste0("Top2_", k)]] <- apply(diff_IS[,grep(paste0("_", k, "$"), colnames(diff_IS))],
                                           1,
                                           function(x) paste(order(as.numeric(abs(x)),
                                                                   decreasing = TRUE)[1:2],
                                                             collapse = ";")
    )


    message("Getting I-score")
    diff_IS[[paste0("IS_", k)]] <-
      apply(diff_IS[,grep(paste0("_", k, "$"), colnames(diff_IS))], 1,
            function(d){
              top = d[length(d)]
              top <- strsplit(top, ";")[[1]]
              top <- as.numeric(top)
              d = d[-length(d)]
              d = as.numeric(d)
              stabilization = sign(mean(d, na.rm = TRUE))
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
              res = abs(res)*stabilization
            }
      )
    diff_IS[[paste0("IS_", k)]] <- (diff_IS[[paste0("IS_", k)]] -
                                      mean(diff_IS[[paste0("IS_", k)]], na.rm=TRUE))/sd(diff_IS[[paste0("IS_", k)]], na.rm=TRUE)

    message("Getting Fisher p-value")
    diff_IS[[paste0("pvalF1_", k)]] <- rep(1,n)
    diff_IS[[paste0("pvalF2_", k)]] <- rep(1,n)
    diff_IS[[paste0("Fisher_", k)]] <- rep(1,n)
    for(i in 1:n){
      topF <- diff_IS[[paste0("Top2_", k)]][i]
      topF <- strsplit(topF, ";")[[1]]
      topF <- as.numeric(topF)

      pF1 <- pval[[k]][[temp[topF[1]]]][i]
      pF2 <- pval[[k]][[temp[topF[2]]]][i]
      if(is.na(pF1)){
        pF1 <- 1
      }
      if(is.na(pF2)){
        pF2 <- 1
      }

      diff_IS[[paste0("pvalF1_", k)]][i] <- pF1
      diff_IS[[paste0("pvalF2_", k)]][i] <- pF2

      diff_IS[[paste0("Fisher_", k)]][i] <- metap::sumlog(c(pF1, pF2))$p
    }
    diff_IS[[paste0("pval", temp[1],  "_", k)]] <- pval[[k]][[temp[1]]]

    diff_IS[[paste0("GlobalScore_", k)]] <- tidyr::replace_na(diff_IS[[paste0("IS_", k)]]*(-log10(diff_IS[[paste0("Fisher_", k)]])), 0)
  }
  info <- data_diff[,c("id", "description", "sumUniPeps")]
  info$description <- stringr::str_extract(paste0(info$description, " "), "(?<=GN=).+?(?= )") # extract genes
  colnames(info) <- c("id", "Gene", "sumUniPeps")
  diff_IS <- as.data.frame(cbind(info, diff_IS))

  diff_IS_plot <- diff_IS[,c("id", "Gene", grep("^IS_|^Fisher_", colnames(diff_IS), value = TRUE))]
  diff_IS_plot <- tidyr::gather(diff_IS_plot, treatment, reading, -id, -Gene)
  diff_IS_plot <- tidyr::separate(diff_IS_plot, treatment, into = c("Value", "treatment"), sep = "_")
  diff_IS_plot <- tidyr::spread(diff_IS_plot, Value, reading)
  diff_IS_plot$treatment <- factor(diff_IS_plot$treatment)

  if(fixed_score_cutoff){
    cutoff <- diff_IS_plot %>% dplyr::group_by(treatment) %>%
      dplyr::mutate(BH = (order(order(Fisher))/length(Fisher))*FDR) %>%
      dplyr::reframe(pval = find_cutoff(Fisher, BH),
                     IS_pos = IS_cutoff,
                     IS_neg = -IS_cutoff)
  }
  else{
    cutoff <- diff_IS_plot %>% dplyr::group_by(treatment) %>%
      dplyr::mutate(BH = (order(order(Fisher))/length(Fisher))*FDR) %>%
      dplyr::reframe(pval = find_cutoff(Fisher, BH),
                     IS_pos = IS_cutoff + median(abs(IS)[which(Fisher < quantile(Fisher, 0.5, na.rm = TRUE))], na.rm = TRUE),
                     IS_neg = -IS_cutoff - median(abs(IS)[which(Fisher < quantile(Fisher, 0.5, na.rm = TRUE))], na.rm = TRUE))
  }

  cutoff_file <- paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "cutoff.txt")
  readr::write_tsv(cutoff, cutoff_file)
  extra_info <- paste0("\nParameters: \nIS cutoff=", IS_cutoff, ", FDR=", FDR*100, "%, curvature=", curvature)
  write(extra_info, cutoff_file, sep = "\n", append = TRUE)


  diff_IS_plot <- diff_IS_plot %>% dplyr::group_by(id, Gene, treatment) %>%
    dplyr::mutate(criteria = Fisher <= cutoff$pval[which(cutoff$treatment == treatment)] &
                    (IS >= cutoff$IS_pos[which(cutoff$treatment == treatment)] | IS <= cutoff$IS_neg[which(cutoff$treatment == treatment)]),
                  curve = curve(IS, cutoff$IS_neg[which(cutoff$treatment == treatment)],
                                cutoff$IS_pos[which(cutoff$treatment == treatment)],
                                cutoff$pval[which(cutoff$treatment == treatment)],
                                curvature = curvature
                  ),
                  criteria_curve = -log10(Fisher) >= curve
    )
  diff_IS_plot$criteria_curve <- tidyr::replace_na(diff_IS_plot$criteria_curve, FALSE)

  cond <- unique(diff_IS_plot$treatment)
  n_cond <- length(cond)
  df_curve <- data.frame(IS = rep(seq(-max(abs(diff_IS_plot$IS), na.rm = TRUE),
                                      max(abs(diff_IS_plot$IS), na.rm = TRUE), 0.005),
                                  n_cond))
  df_curve$treatment <- rep(cond, each = nrow(df_curve)/n_cond)
  df_curve <- df_curve %>% dplyr::group_by(treatment, rownames(df_curve)) %>%
    dplyr::mutate(curve = curve(IS, cutoff$IS_neg[which(cutoff$treatment == treatment)],
                                cutoff$IS_pos[which(cutoff$treatment == treatment)], cutoff$pval[which(cutoff$treatment == treatment)],
                                curvature = curvature)
    )

  message("Creating and saving plot")
  g_h <- ggplot(diff_IS_plot, aes(IS, -log10(Fisher), color = criteria_curve)) +
    geom_point() +
    geom_line(data = df_curve, aes(x = IS, y = curve), linetype = "dashed", color = "black") +
    facet_wrap(~treatment) +
    ggrepel::geom_label_repel(data = diff_IS_plot[diff_IS_plot$criteria_curve,],
                              aes(IS, -log10(Fisher),
                                  label = Gene), show.legend = FALSE) +
    ylim(c(0, max(-log10(diff_IS_plot$Fisher)))) +
    xlim(c(-max(abs(diff_IS_plot$IS), na.rm = TRUE),
           max(abs(diff_IS_plot$IS), na.rm = TRUE))
         ) +
    labs(title = "I-score plot",
         y = "-log10(Fisher p-value)",
         x = "I-score") +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(face = "bold"))

  ggsave(paste0(format(Sys.time(), "%y%m%d_%H%M"), "_", "hits_plot.png"),
         plot = g_h, device = "png",  path = outdir,
         width = 16, height = 9)

  g_I <- ggplot(diff_IS_plot, aes(IS, -log10(Fisher), color = criteria_curve)) +
    geom_point(aes(text = paste0("Fisher p-value: ", diff_IS_plot$Fisher, "\n",
                                 "I-score: ", diff_IS_plot$IS, "\n",
                                 "Protein: ", diff_IS_plot$id, "\n",
                                 "Gene: ", diff_IS_plot$Gene))
               ) +
    geom_line(data = df_curve, aes(x = IS, y = curve), linetype = "dashed", color = "black") +
    facet_wrap(~treatment) +
    ylim(c(0, max(-log10(diff_IS_plot$Fisher), na.rm = TRUE))) +
    xlim(c(-max(abs(diff_IS_plot$IS), na.rm = TRUE),
           max(abs(diff_IS_plot$IS), na.rm = TRUE))
         ) +
    labs(title = "I-score plot",
         y = "-log10(Fisher p-value)",
         x = "I-score") +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(face = "bold"))


  g_I <- plotly::ggplotly(g_I, tooltip = "text", width = 1080, height = 560)
  htmltools::save_html(g_I, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "hits_plotInt.html"))

  message("Categorization...")
  diff_IS_plot$criteria <- NULL
  diff_IS_plot <- diff_IS_plot[diff_IS_plot$criteria_curve,]
  diff_IS_plot$criteria_curve <- NULL
  diff_IS_plot$curve <- NULL

  for_categorize <- diff_IS[,grep("^id|Gene|^\\d{1,}|^pval37", colnames(diff_IS))]
  for_categorize <- for_categorize %>% tidyr::gather("key", "value", -id, -Gene) %>%
    tidyr::separate(key, into = c("key", "treatment"), sep = "_") %>%
    tidyr::spread(key, value)

  cutoff <- for_categorize %>% dplyr::group_by(treatment) %>%
    dplyr::mutate(BH = (order(order(pval37C))/length(pval37C))*FDR_category) %>% # FDR of 10% for significant 37
    dplyr::reframe(pval = find_cutoff(pval37C, BH))

  diff_IS_plot <- left_join(diff_IS_plot, for_categorize, by = c("id", "Gene", "treatment"))
  diff_IS_plot$category <- NA
  diff_IS_plot$category[which(is.na(diff_IS_plot$pval37C))] <- "ND"

  diff_IS_plot <- diff_IS_plot %>% group_by(id, Gene, treatment) %>%
    mutate(is_C = pval37C <= cutoff$pval[which(cutoff$treatment == treatment)])

  diff_IS_plot$category[which(!diff_IS_plot$is_C)] <- "NC"
  diff_IS_plot$category[which(is.na(diff_IS_plot$category) & is.na(diff_IS_plot[[temp[length(temp)]]]))] <- "CC" #denatured in the end but significant 37

  diff_IS_plot$coeff <- apply(diff_IS_plot[,grep("^\\d{1,}", colnames(diff_IS_plot))], 1,
                              function(x){
                                x <- abs(x)
                                n <- 1/length(na.omit(x))
                                x <- na.omit(x)
                                x <- x/sum(x)
                                x <- sum((x - n)**2);
                                x
                              })
  diff_IS_plot$category[which(is.na(diff_IS_plot$category) & diff_IS_plot$coeff <= 1e-2)] <- "CN"
  diff_IS_plot$category[which(is.na(diff_IS_plot$category))] <- "CC"
  diff_IS_plot <- diff_IS_plot[,c("id", "Gene", "treatment", "Fisher", "IS", "category")]

  for_categorize <- diff_IS_plot[,c("id", "Gene", "treatment", "category")] %>%
    tidyr::spread(treatment, category)
  colnames(for_categorize)[-c(1:2)] <- paste0("category_", colnames(for_categorize)[-c(1:2)])
  diff_IS <- diff_IS %>% dplyr::left_join(for_categorize, by = c("id", "Gene"))

  diff_IS[,grep("category", colnames(diff_IS))] <- apply(as.data.frame(diff_IS[,grep("category", colnames(diff_IS))]),
                                                         2, function(x) tidyr::replace_na(x, "NN")
                                                         )

  message("Saving datas...")
  openxlsx::write.xlsx(diff_IS, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "hits_analysis_tab.xlsx"))
  openxlsx::write.xlsx(diff_IS_plot, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "hits_summary.xlsx"))

  if(nrow(diff_IS_plot) > 1 & length(cond) > 1){
    diff_IS_plot$treatment <- factor(diff_IS_plot$treatment)
    vennlist <- diff_IS_plot
    vennlist$category <- NULL
    vennlist <- (vennlist %>% dplyr::ungroup() %>% dplyr::group_by(treatment) %>%
                   dplyr::summarize(vennlist = list(id)) %>%
                   dplyr::select(vennlist))[[1]]
    names(vennlist) <- levels(diff_IS_plot$treatment)
    VennDiagram::venn.diagram(vennlist, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "VennDiagram.png"),
                              output = TRUE, disable.logging = TRUE, imagetype = "png",
                              fill = RColorBrewer::brewer.pal(length(names(vennlist)), "Pastel2")[1:length(vennlist)],
                              fontface = "bold",
                              fontfamiliy = "sans",
                              cat.cex = 1.6,
                              cat.fontface = "bold",
                              cat.fontfamily = "sans")
    vennlist <- com_protein_loop(vennlist)
    too_long <- which(sapply(names(vennlist), nchar) > 31)
    if(length(too_long)){
      for(n in too_long){
        name_toolong <- names(vennlist)[n]
        in_common <- Reduce(intersect, strsplit(strsplit(name_toolong, " & ")[[1]], ""))
        if(length(in_common)){
          name_toolong <- gsub(paste(in_common, collapse = "|"), "", name_toolong)
        }
        if(nchar(name_toolong) > 31){
          name_toolong <- gsub(" ", "", name_toolong)
        }
        if(nchar(name_toolong) > 31){
          name_toolong <- paste0("&", name_toolong, "&")
          name_toolong <- stringr::str_remove_all(name_toolong, "(?<=&[a-zA-Z]).+?(?=&)")
          name_toolong <-  gsub("^&|&$", "", name_toolong)
        }
        names(vennlist)[n] <- name_toolong
      }
    }
    for(i in 1:length(vennlist)){
      prot <- data.frame("id" = vennlist[[i]],
                         "Gene" = info$Gene[which(!is.na(match(info$id,
                                                                 vennlist[[i]])))]
      )
      score_info <- diff_IS[,c(1, grep("^IS_|^Fisher_", colnames(diff_IS)))]

      vennlist[[i]] <- dplyr::left_join(prot, score_info, by = "id")
    }
    openxlsx::write.xlsx(vennlist, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "Venn tab.xlsx"))
  }

  message("Calculation done !")
  return(diff_IS_plot)
}


### function to compute cutoff curve
curve <- function(x, cut_neg, cut_pos, cut_p, curvature = curvature){
  y <- rep(NA, length(x))
  neg <- which(x <= cut_neg)
  pos <- which(x >= cut_pos)
  y[neg] <- curvature/abs(x[neg] - cut_neg) + -log10(cut_p)
  y[pos] <- curvature/abs(x[pos] - cut_pos) + -log10(cut_p)

  return(y)
}

### function to find p-value cutoff
find_cutoff <- function(x,y){
  id <- order(x)
  x <- x[id]
  y <- y[id]

  x <- x[which(x < y)]
  return(x[length(x)])
}
