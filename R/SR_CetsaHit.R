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

SR_CetsaHit <- function(data, data_diff = NULL, ctrl, valid_val = NULL,
                        SR_cutoff = 1.5,
                        FDR = 0.01, curvature = 0.5,
                        FDR_category = 0.1,
                        folder_name = "Hits_analysis",
                        peptide_count_col = "peptides_counts_all",
                        species = "Homo Sapiens"){
  wd <- getwd()
  outdir <- paste0(wd, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", folder_name)
  if (dir.exists(outdir) == FALSE) {
    dir.create(outdir)
  }

  QP_column <- stringr::str_which(colnames(data), "^36C_")
  if(length(QP_column)){
    data <- data[,-QP_column]
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
  else if("character" %in% class(data_diff)){
    if(stringr::str_detect(data_diff, "\\.tsv$|\\.csv$|\\.txt$"))
      data_diff <- readr::read_tsv(data_diff)
    else
      stop("Format isn't recognize")
  }
  else if(!("data.frame" %in% class(data_diff))){
    stop("data_diff is neither a file, neither a data frame !")
  }


  QP_column <- stringr::str_which(colnames(data_diff), "^36C_")
  if(length(QP_column)){
    data_diff <- data_diff[,-QP_column]
  }

  curve <- function(x, cut_neg, cut_pos, cut_p, curvature = curvature){
    y <- rep(NA, length(x))
    neg <- which(x <= cut_neg)
    pos <- which(x >= cut_pos)
    y[neg] <- curvature/abs(x[neg] - cut_neg) + -log10(cut_p)
    y[pos] <- curvature/abs(x[pos] - cut_pos) + -log10(cut_p)

    return(y)
  }
  find_cutoff <- function(x,y){
    id <- order(x)
    x <- x[id]
    y <- y[id]

    x <- x[which(x > y)]
    return(x[1])
  }


  message("Computing mean values...")
  # get average value among bioreplicate for each protein
  diff_SR <- data_diff[,-stringr::str_which(colnames(data_diff), paste0("_", ctrl, "$"))]
  diff_SR <- tidyr::gather(diff_SR, Condition, reading, -id, -description, -sumUniPeps, -sumPSMs, -countNum)
  diff_SR <- tidyr::separate(diff_SR, Condition, into = c("temperature",
                                                          "replicate", "treatment"), sep = "_")
  diff_SR <- diff_SR %>% dplyr::group_by(id, description, temperature, treatment) %>%
    dplyr::summarise(reading.mean = mean(reading, na.rm = T))
  diff_SR <- tidyr::unite(diff_SR, Condition, temperature, treatment,
                          sep = "_")
  diff_SR <- tidyr::spread(diff_SR, Condition, reading.mean)


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

    message("Getting top2 fraction")
    diff_SR[[paste0("Top2_", k)]] <- apply(diff_SR[,stringr::str_which(colnames(diff_SR), paste0("_", k, "$"))],
                                           1,
                                           function(x) paste(order(as.numeric(abs(x)),
                                                                   decreasing = TRUE)[1:2],
                                                             collapse = ";")
    )


    message("Getting stability rate")
    diff_SR[[paste0("SR_", k)]] <-
      apply(diff_SR[,stringr::str_which(colnames(diff_SR), paste0("_", k, "$"))], 1,
            function(d){
              top = d[length(d)]
              top <- stringr::str_split(top, ";")[[1]]
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
    diff_SR[[paste0("SR_", k)]] <- (diff_SR[[paste0("SR_", k)]] -
                                      mean(diff_SR[[paste0("SR_", k)]], na.rm=TRUE))/sd(diff_SR[[paste0("SR_", k)]], na.rm=TRUE)

    message("Getting Fisher p-value")
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
    diff_SR[[paste0("pval", temp[1],  "_", k)]] <- pval[[k]][[temp[1]]]

    diff_SR[[paste0("GlobalScore_", k)]] <- tidyr::replace_na(diff_SR[[paste0("SR_", k)]]*(-log10(diff_SR[[paste0("Fisher_", k)]])), 0)
  }
  info <- data_diff[,c("id", "description", "sumUniPeps")]
  info$description <- stringr::str_extract(paste0(info$description, " "), "(?<=GN=).+?(?= )") # extract genes
  colnames(info) <- c("id", "Genes", "peptide_count")
  diff_SR <- as.data.frame(cbind(info, diff_SR))

  diff_SR_plot <- diff_SR[,c("id", "Genes", stringr::str_subset(colnames(diff_SR), "^SR_|^Fisher_"))]
  diff_SR_plot <- tidyr::gather(diff_SR_plot, Condition, reading, -id, -Genes)
  diff_SR_plot <- tidyr::separate(diff_SR_plot, Condition, into = c("Value", "Condition"))
  diff_SR_plot <- tidyr::spread(diff_SR_plot, Value, reading)
  diff_SR_plot$Condition <- factor(diff_SR_plot$Condition)

  cutoff <- diff_SR_plot %>% dplyr::group_by(Condition) %>%
    dplyr::mutate(BH = (order(order(Fisher))/length(Fisher))*FDR) %>%
    dplyr::summarise(pval = find_cutoff(Fisher, BH),
                     SR_pos = SR_cutoff + median(SR[which(Fisher < quantile(Fisher, 0.5))], na.rm = TRUE),
                     SR_neg = -SR_cutoff - median(SR[which(Fisher < quantile(Fisher, 0.5))], na.rm = TRUE))

  print(cutoff)

  diff_SR_plot <- diff_SR_plot %>% dplyr::group_by(id, Genes, Condition) %>%
    dplyr::mutate(criteria = Fisher <= cutoff$pval[which(cutoff$Condition == Condition)] &
                    (SR >= cutoff$SR_pos[which(cutoff$Condition == Condition)] | SR <= cutoff$SR_neg[which(cutoff$Condition == Condition)]),
                  curve = curve(SR, cutoff$SR_neg[which(cutoff$Condition == Condition)],
                                cutoff$SR_pos[which(cutoff$Condition == Condition)],
                                cutoff$pval[which(cutoff$Condition == Condition)],
                                curvature = curvature
                  ),
                  criteria_curve = -log10(Fisher) >= curve
    )
  diff_SR_plot$criteria_curve <- tidyr::replace_na(diff_SR_plot$criteria_curve, FALSE)

  cond <- unique(diff_SR_plot$Condition)
  n_cond <- length(cond)
  df_curve <- data.frame(SR = rep(seq(min(diff_SR_plot$SR, na.rm = TRUE), max(diff_SR_plot$SR, na.rm = TRUE), 0.01), n_cond))
  df_curve$Condition <- rep(cond, each = nrow(df_curve)/n_cond)
  df_curve <- df_curve %>% dplyr::group_by(Condition, rownames(df_curve)) %>%
    dplyr::mutate(curve = curve(SR, cutoff$SR_neg[which(cutoff$Condition == Condition)],
                                cutoff$SR_pos[which(cutoff$Condition == Condition)], cutoff$pval[which(cutoff$Condition == Condition)],
                                curvature = curvature)
    )

  message("Creating and saving plot")
  g_h <- ggplot(diff_SR_plot, aes(SR, -log10(Fisher), color = criteria_curve)) +
    geom_point() +
    geom_line(data = df_curve, aes(x = SR, y = curve), linetype = "dashed", color = "black") +
    ylim(c(0, max(-log10(diff_SR_plot$Fisher)))) +
    labs(title = "Stability rate plot",
         y = "-log10(Fisher p-value)",
         x = "Stability rate") +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70"))  +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~Condition) +
    ggrepel::geom_label_repel(data = diff_SR_plot[diff_SR_plot$criteria_curve,],
                              aes(SR, -log10(Fisher),
                                  label = Genes), show.legend = FALSE)

  ggsave(paste0(format(Sys.time(), "%y%m%d_%H%M"), "_", "hits_plot.png"),
         plot = g_h,
         device = "png",
         path = outdir,
         width = 14,
         height = 8)


  g_I <-ggplot(diff_SR_plot, aes(SR, -log10(Fisher), color = criteria_curve)) +
    geom_point(aes(text = paste0("Fisher p-value: ", diff_SR_plot$Fisher, "\n",
                                 "SR: ", diff_SR_plot$SR, "\n",
                                 "Protein: ", diff_SR_plot$id, "\n",
                                 "Genes: ", diff_SR_plot$Genes))
               ) +
    geom_line(data = df_curve, aes(x = SR, y = curve), linetype = "dashed", color = "black") +
    ylim(c(0, max(-log10(diff_SR_plot$Fisher)))) +
    labs(title = "Stability rate plot",
         y = "-log10(Fisher p-value)",
         x = "Stability rate") +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70"))  +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~Condition)

  g_I <- plotly::ggplotly(g_I, tooltip = "text", width = 1080, height = 560)
  htmltools::save_html(g_I, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "hits_plotInt.html"))

  message("Categorization...")
  diff_SR_plot$criteria <- NULL
  diff_SR_plot <- diff_SR_plot[diff_SR_plot$criteria_curve,]
  diff_SR_plot$criteria_curve <- NULL
  diff_SR_plot$curve <- NULL

  for_categorize <- diff_SR[,str_which(colnames(diff_SR), "^id|Genes|^\\d{1,}|^pval37")]
  for_categorize <- for_categorize %>% tidyr::gather("key", "value", -id, -Genes) %>%
    tidyr::separate(key, into = c("key", "Condition"), sep = "_") %>%
    tidyr::spread(key, value)

  cutoff <- for_categorize %>% dplyr::group_by(Condition) %>%
    dplyr::mutate(BH = (order(order(pval37C))/length(pval37C))*FDR_category) %>% # FDR of 10% for significant 37
    dplyr::summarise(pval = find_cutoff(pval37C, BH))

  diff_SR_plot <- left_join(diff_SR_plot, for_categorize, by = c("id", "Genes", "Condition"))
  diff_SR_plot$category <- NA
  diff_SR_plot$category[which(is.na(diff_SR_plot$pval37C))] <- "ND"

  diff_SR_plot <- diff_SR_plot %>% group_by(id, Genes, Condition) %>%
    mutate(is_C = pval37C <= cutoff$pval[which(cutoff$Condition == Condition)])

  diff_SR_plot$category[which(!diff_SR_plot$is_C)] <- "NC"
  diff_SR_plot$category[which(is.na(diff_SR_plot$category) & is.na(diff_SR_plot[[temp[length(temp)]]]))] <- "CC" #denatured in the end but significant 37

  diff_SR_plot$coeff <- apply(diff_SR_plot[,str_which(colnames(diff_SR_plot), "^\\d{1,}")], 1,
                              function(x){
                                x <- abs(x)
                                n <- 1/length(na.omit(x))
                                x <- na.omit(x)
                                x <- x/sum(x)
                                x <- sum((x - n)**2);
                                x
                              })
  diff_SR_plot$category[which(is.na(diff_SR_plot$category) & diff_SR_plot$coeff <= 1e-2)] <- "CN"
  diff_SR_plot$category[which(is.na(diff_SR_plot$category))] <- "CC"
  diff_SR_plot <- diff_SR_plot[,c("id", "Genes", "Condition", "Fisher", "SR", "category")]

  for_categorize <- diff_SR_plot[,c("id", "Genes", "Condition", "category")] %>%
    tidyr::spread(Condition, category)
  colnames(for_categorize)[-c(1:2)] <- paste0("category_", colnames(for_categorize)[-c(1:2)])
  diff_SR <- diff_SR %>% dplyr::left_join(for_categorize, by = c("id", "Genes"))
  diff_SR[,str_which(colnames(diff_SR), "category")] <- apply(as.data.frame(diff_SR[,str_which(colnames(diff_SR), "category")]),
                                                              2, function(x) tidyr::replace_na(x, "NN")
  )

  message("Saving datas...")
  openxlsx::write.xlsx(diff_SR, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "hits_analysis_tab.xlsx"))
  openxlsx::write.xlsx(diff_SR_plot, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "hits_summary.xlsx"))

  if(nrow(diff_SR_plot) > 1){
    diff_SR_plot$Condition <- factor(diff_SR_plot$Condition)
    vennlist <- diff_SR_plot
    vennlist$category <- NULL
    vennlist <- (vennlist %>% dplyr::ungroup() %>% dplyr::group_by(Condition) %>%
                   dplyr::summarize(vennlist = list(id)) %>%
                   dplyr::select(vennlist))[[1]]
    names(vennlist) <- levels(diff_SR_plot$Condition)
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
      score_info <- diff_SR[,c(1, stringr::str_which(colnames(diff_SR), "^SR_|^Fisher_"))]

      vennlist[[i]] <- dplyr::left_join(prot, score_info, by = "id")
    }
    openxlsx::write.xlsx(vennlist, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "Venn tab.xlsx"))
  }

  return(diff_SR_plot)
}

