#' imprints_cleaved_peptides
#'
#' Function to find proteins that are potentially cleaved and their cleaved sites.
#' For more information, see Details section.
#'
#' @details
#' The idea here is to compute the sum from all fold change for each peptide and then compute the cumulative sum of
#' these summed fold change of every peptides for each protein.
#' If the protein is not cleaved, this cumulative sum should be constantly increasing/decreasing, i.e linear.
#' If the protein is cleaved, then this cumulative sum will either be convex or concave.
#' If convex, the last peptides are more abundant; if concave the first peptides are more abundant.
#' In the end a volcano plot will be plotted where the p-value corresponds to the combined p-value from
#' the six temperatures comparison between peptides located towards the N-terminal and the ones located
#' towards the C-terminal position.
#'
#' You can then use the function \code{imprints_sequence_peptides} to plot the peptides before and after the potential cleaved site found.
#'
#' @param data The normalized peptides data set, i.e. the outpout from \code{imprints_normalize_peptides}.
#' @param data_diff The log2 fold-changes peptides data set, i.e. the outpout from \code{imprints_sequence_peptides}.
#'                  If NULL, it will be computed.
#' @param control The control treatment from your dataset.
#' @param min_ValidValue The minimum proportion of valid values per peptides.
#'                       Default is 0.4; so if 6 temperatures need at least 3 non missing values.
#' @param FDR The FDR used to obtained the final p-value cutoff and maximum fold-change cutoff. Default is 1%
#' @param diff_cutoff The maximum fold-change cutoff. Default is 0.25
#' @param curvature The curvature used for the curve on the volcano plot
#' @param folder_name The name of the folder in which you want to save the results.
#'
#' @return The potential cleaved sites from the proteins considered as cleaved.
#'
#' @export
#'

imprints_cleaved_peptides <- function(data, data_diff = NULL,
                                      control = "Vehicle", min_ValidValue = 0.4,
                                      FDR = 0.01, diff_cutoff = 0.25,
                                      curvature = 0.1, folder_name = ""){
  if(!any(grepl(paste0("_", control, "$"), colnames(data)))){
    message(paste("Error:", control, "is wasn't found in your data !"))
    return()
  }

  wd <- getwd()
  outdir <- paste0(wd, "/", format(Sys.time(), "%y%m%d_%H%M"), "_RESP_analysis",
                   ifelse(nchar(folder_name) > 0, "_", ""), folder_name)
  if (dir.exists(outdir) == FALSE) {
    dir.create(outdir)
  }

  require(limma)

  # if any, removing QP columns
  if(any(grepl("^36C_", colnames(data)))){
    data <- data[,-grep("^36C_", colnames(data))]
  }

  if(is.null(data_diff)){
    message("No fold-change dataset specified, computing...")
    data_diff <- imprints_sequence_peptides(data, control = control,
                                            dataset_name = deparse(substitute(data))
                                            )
  }

  if(any(grepl("^36C_", colnames(data_diff)))){
    data_diff <- data_diff[,-grep("^36C_", colnames(data_diff))]
  }

  if("countNum" %in% colnames(data)){
    data$countNum <- NULL
  }
  if("countNum" %in% colnames(data_diff)){
    data_diff$countNum <- NULL
  }

  data_cleaved <- data_diff
  if(any(grepl(paste0("_", control, "$"), colnames(data_cleaved)))){
      data_cleaved <- data_cleaved[,-grep(paste0("_", control, "$"), colnames(data_cleaved))]
  }

  # data_cleaved is output from 'imprints_sequence_peptides' --> remove unnecessary informations
  data_cleaved$Master.Protein.Accessions <- data_cleaved$Positions.in.Master.Proteins
  data_cleaved$Positions.in.Master.Proteins <- NULL
  data_cleaved$Annotated.Sequence <- NULL
  data_cleaved$Modifications <- NULL
  colnames(data_cleaved)[1] <- "id" # easier handle


  message("Averaging value, summing profiles...")
  # get summed value from all temperature for every peptide + # keeping track of number of non-missing values
  data_cleaved <- data_cleaved %>%
    tidyr::gather("treatment", "value", -id, -description) %>%
    tidyr::separate(treatment, into = c("temp", "biorep", "treatment"), sep = "_") %>%
    dplyr::group_by(id, description, temp, treatment, biorep) %>% # if peptide modified, mean values from same peptide sequence
    dplyr::summarise(value = mean(value, na.rm = TRUE)) %>%
    dplyr::mutate(N = paste(biorep, as.numeric(!is.na(value)), sep = "_")) %>%  # keeping track of number of non-missing values
    dplyr::ungroup() %>% dplyr::group_by(id, description, temp, treatment) %>%
    dplyr::summarise(mean_value = mean(value, na.rm = TRUE),
                     Nvalue = paste(N, collapse = ";")) %>%  # keeping track of number of non-missing values
    dplyr::mutate(Nvalue = paste(temp, Nvalue, sep = "-")) %>%
    dplyr::ungroup() %>% dplyr::group_by(id, description, treatment) %>%
    dplyr::filter(length(na.omit(mean_value))/length(mean_value) >= min_ValidValue) %>%  # keeping peptides with more than 40% of valid values
    dplyr::reframe(sum_profile = sum(mean_value, na.rm = TRUE),  # get sum value for each peptides --> sum all temperatures values
                   Nvalue = paste(Nvalue, collapse = "|"))   # keeping track of number of non-missing values


  message("Filtering and ordering peptides...")
  # separate id in protein and sequence (need to take into account protein groups)
  data_cleaved$sequence <- gsub("(?<= |^).{5,6}(?= \\[)", "", data_cleaved$id, perl = TRUE)
  data_cleaved$sequence <- gsub(" ", "", data_cleaved$sequence)
  data_cleaved$sequence <- gsub(";", "; ", data_cleaved$sequence)

  data_cleaved$id <- gsub("(?<=\\[)\\d{1,4}-\\d{1,4}(?=\\])", "", data_cleaved$id, perl = TRUE)
  data_cleaved$id <- gsub(" \\[(\\]|\\];)", "", data_cleaved$id)
  data_cleaved$id <- gsub(" ", "; ", data_cleaved$id)

  data_cleaved <- data_cleaved %>%
    dplyr::mutate(sequence = gsub("\\[|\\]", "", sequence)) %>%
    dplyr::ungroup() %>% dplyr::group_by(id, description, treatment) %>%
    dplyr::filter(length(sum_profile) > 3)  # only keeping proteins with more than 3 peptides


  # order data according proteins and sequence
  # ordering according sequence is capital
  data_cleaved$factor <- data_cleaved$id
  data_cleaved <- data_cleaved[order(data_cleaved$id),]
  ord <- sapply(strsplit(gsub(";.*", "", data_cleaved$sequence), "-|~"), # onmy taking first protein if proteinn group
                function(x) {sum(as.numeric(x))}
                )
  ord <- data.frame(factor = data_cleaved$factor, idx = 1:nrow(data_cleaved), sum.pos = ord) %>%
    dplyr::group_by(factor) %>%
    dplyr::mutate(ord.sum.pos = order(sum.pos),
                  nb.pep = length(ord.sum.pos),
                  final.order = idx[ord.sum.pos])

  data_cleaved <- data_cleaved[ord$final.order,]
  data_cleaved$Master.Protein.Accessions <- data_cleaved$factor
  data_cleaved$factor <- NULL

  message("Computing cumulative sum to find potential cleaved sites...")
  # compute cumulative sum for each protein
  data_cleaved <- data_cleaved %>%
    dplyr::group_by(id, description, treatment) %>%
    dplyr::mutate(cumsum_profile = cumsum(tidyr::replace_na(sum_profile, 0)))

  # get cleaved site
  # if cumulative abundance concave --> first derivative is decreasing globally --> inflexion point is where second derivative is minimum
  # if cumulative abundance convex --> first derivative is increasing globally --> inflexion point is where second derivative is maximum
  data_cleaved <- data_cleaved %>% dplyr::ungroup() %>%
    dplyr::group_by(id, description, treatment) %>%
    dplyr::mutate(cleaved_site = sequence[inflex(cumsum_profile)]) %>%  # get back potential cleaved site
    dplyr::filter(!is.na(cleaved_site))  # if NA --> not cleaved

  message("Checking cleaved sites, formating...")
  # set peptide position --> before cleaved site = N-terminal side / after cleaved site = C-terminal side
  # if protein group, only taking first protein
  data_cleaved$position_global <- unlist(lapply(strsplit(gsub(";.*", "", data_cleaved$sequence),
                                                 "~|-"),
                                        function(s) sum(as.numeric(s)))
                                 )
  data_cleaved$cleaved.global <- unlist(lapply(strsplit(gsub(";.*", "", data_cleaved$cleaved_site),
                                                "~|-"),
                                       function(s) sum(as.numeric(s)))
                                )

  data_cleaved <- data_cleaved %>%
    dplyr::mutate(position_global = ifelse(position_global < cleaved.global,
                                           "N", ifelse(sequence == cleaved_site,
                                                       "cleaved", "C"))
           ) %>%
    dplyr::select(-cleaved.global)

  # checking cleaved site, if not assign right position
  data_cleaved <- data_cleaved %>% dplyr::group_by(id, description, treatment) %>%
    dplyr::group_modify(~ {
      val_N <- .x$sum_profile[which(.x$position_global == "N")]
      val_C <- .x$sum_profile[which(.x$position_global == "C")]
      cle <- .x$sum_profile[which(.x$position_global == "cleaved")]
      non_cleaved <- .x$cleaved_site[1]

      bef_cleaved <- .x$sequence[which(.x$position_global == "cleaved") -1]
      aft_cleaved <- .x$sequence[which(.x$position_global == "cleaved") +1]

      bef_cleaved <- sapply(strsplit(bef_cleaved, "; ")[[1]],
                            function(x){
                              x <- strsplit(x, "-|~")[[1]]
                              as.numeric(x[2]) + 1
                            }, USE.NAMES = FALSE)
      aft_cleaved <- sapply(strsplit(aft_cleaved, "; ")[[1]],
                            function(x){
                              x <- strsplit(x, "-|~")[[1]]
                              as.numeric(x[1]) - 1
                            }, USE.NAMES = FALSE)

      # proportion augmentation/diminution of sd
      pN <- (sd(c(val_N, cle), na.rm = TRUE)/sd(val_N, na.rm = TRUE))
      pC <- (sd(c(val_C, cle), na.rm = TRUE)/sd(val_C, na.rm = TRUE))

      if(any(is.na(c(pN, pC)))){ # if only one peptide in N or C --> sd is NA; only one of them can be NA (at least 4 peptides in total)
        # if only one, cleaved site is assign to this group
        .x$position_global[which(.x$position_global == "cleaved")] <- ifelse(is.na(pN), "N", "C")
        non_cleaved <- sapply(strsplit(non_cleaved, "; ")[[1]],
                              function(x){
                                x <- strsplit(x, "-|~")[[1]]
                                as.numeric(x)
                              }, simplify = FALSE)
        if(is.na(pN)){
          non_cleaved <- sapply(non_cleaved, "[[", 2) + 1
          non_cleaved <- paste(non_cleaved, aft_cleaved, sep = "~")
        }
        else{
          non_cleaved <- sapply(non_cleaved, "[[", 1) - 1
          non_cleaved <- paste(bef_cleaved, non_cleaved, sep = "~")
        }
        non_cleaved <- paste(non_cleaved, collapse = "; ")
      }
      else if(cle >= 0){ # to be cleaved, peptide should be less abundant, i.e. negative fold changes
        .x$position_global[which(.x$position_global == "cleaved")] <- ifelse(pN <= pC, "N", "C")
        non_cleaved <- sapply(strsplit(non_cleaved, "; ")[[1]],
                              function(x){
                                x <- strsplit(x, "-|~")[[1]]
                                as.numeric(x)
                              }, simplify = FALSE)
        if(pN <= pC){
          non_cleaved <- sapply(non_cleaved, "[[", 2) + 1
          non_cleaved <- paste(non_cleaved, aft_cleaved, sep = "~")
        }
        else{
          non_cleaved <- sapply(non_cleaved, "[[", 1) - 1
          non_cleaved <- paste(bef_cleaved, non_cleaved, sep = "~")
        }
        non_cleaved <- paste(non_cleaved, collapse = "; ")
      }
      else{
        if(any(c(pN, pC) < 1)){ # arbitrary, sd augmentation should be more than 100%
          .x$position_global[which(.x$position_global == "cleaved")] <- ifelse(pN <= pC, "N", "C")
          non_cleaved <- sapply(strsplit(non_cleaved, "; ")[[1]],
                                function(x){
                                  x <- strsplit(x, "-|~")[[1]]
                                  as.numeric(x)
                                }, simplify = FALSE)

          if(pN <= pC){
            non_cleaved <- sapply(non_cleaved, "[[", 2) + 1
            non_cleaved <- paste(non_cleaved, aft_cleaved, sep = "~")
          }
          else{
            non_cleaved <- sapply(non_cleaved, "[[", 1) - 1
            non_cleaved <- paste(bef_cleaved, non_cleaved, sep = "~")
          }
          non_cleaved <- paste(non_cleaved, collapse = "; ")
        }
      }
      .x$cleaved_site <- non_cleaved

      return(.x)
    }) %>%
    dplyr::filter(position_global != "cleaved")


  ### sum number of values per peptide N and peptide C
  data_cleaved <- data_cleaved %>%
    dplyr::group_by(id, description, treatment, cleaved_site, position_global) %>%
    dplyr::summarise(Nvalue = nb_pepval(Nvalue),
                     Npep = length(sequence)
                     )


  message("Summing peptides profiles from N-terminal side and C-terminal side from each treatment")
  treat <- unique(data_cleaved$treatment)
  new_data_diff <- list()
  res <- list()
  for(t in treat){
    treat_data <- unique(data_cleaved[data_cleaved$treatment == t, c("id", "cleaved_site")])
    treat_data_diff <- data
    treat_torm <- treat[!(treat %in% t)]
    if(length(treat_torm)){
      treat_torm <- paste0("_", treat_torm, collapse = "|")
      treat_data_diff <- treat_data_diff[,-grep(treat_torm, colnames(treat_data_diff))]
    }

    directory_toremove <- list.files()
    treat_data_diff <- imprints_sequence_peptides(treat_data_diff,
                                                  proteins = treat_data$id,
                                                  sequence = sub("~", "-", treat_data$cleaved_site),
                                                  control = control, dataset_name = t
                                                  )
    directory_toremove <- list.files()[!(list.files() %in% directory_toremove)]
    directory_toremove <- grep(paste0(t, "\\.txt$"), directory_toremove, value = TRUE)
    if(length(directory_toremove) == 1){
      unlink(directory_toremove, recursive = TRUE)
    }


    res[[t]] <- treat_data
    ### removing peptides that could be cleaved site
    treat_data_diff <- imprints_remove_peptides(treat_data_diff,
                                                treat_data$id,
                                                treat_data$cleaved_site
                                                )
    treat_data_diff <- treat_data_diff[,-grep(paste0("_", control), colnames(treat_data_diff))]
    treat_data_diff$Master.Protein.Accessions <- treat_data_diff$Positions.in.Master.Proteins
    treat_data_diff$Modifications <- NULL
    treat_data_diff$Annotated.Sequence <- NULL
    treat_data_diff$Positions.in.Master.Proteins <- NULL
    colnames(treat_data_diff)[1] <- "id"

    treat_data_diff <- treat_data_diff %>%
      tidyr::gather("key", "value", -id, -description) %>%
      tidyr::separate(key, into = c("temperature", "biorep", "treatment"), sep = "_") %>%
      mutate(protein = gsub(" ", "; ", gsub(" \\[(\\]|\\];)", "", gsub("(?<=\\[)\\d{1,4}(-|~)\\d{1,4}(?=\\])", "", id, perl = TRUE))),
             position = gsub(";", "; ", gsub(" ", "", gsub("(?<=; |^).{5,6}(?= \\[)", "", id, perl = TRUE)))
             ) %>%
      select(-id)
    treat_data_diff <- treat_data_diff[,c("protein", "position", "description", "temperature", "biorep", "treatment", "value")]
    new_data_diff[[t]] <- treat_data_diff
  }
  new_data_diff <- as.data.frame(Reduce(rbind, new_data_diff))

  ##### computing p-values for each protein --> H0 = protein is not cleaved
  message("Reshaping data and computing experimental design")
  new_data_diff$position_global <- unlist(lapply(strsplit(gsub("\\[|\\]", "", gsub(";.*", "", new_data_diff$position)),
                                                          "~|-"),
                                                 function(s) sum(as.numeric(s)))
                                          )
  new_data_diff <- new_data_diff %>%
    dplyr::group_by(protein, description, treatment) %>%
    dplyr::mutate(position_global = ifelse(position_global == min(position_global),
                                           "N", "C")) %>%
    dplyr::select(-position)

  new_data_diff$temperature <- factor(new_data_diff$temperature)
  new_data_diff$biorep <- factor(new_data_diff$biorep)
  new_data_diff$position_global <- factor(new_data_diff$position_global,
                                          levels = c("N", "C"))

  ntemp <- nlevels(new_data_diff$temperature)
  nrep <- nlevels(new_data_diff$biorep)
  # npos = 2 ; N or C

  ### limma design
  # design
  position <- rep(levels(new_data_diff$position_global), ntemp*nrep)
  temperature <- rep(levels(new_data_diff$temperature), each=nrep*2)
  replicate <- rep(rep(levels(new_data_diff$biorep), each = 2), ntemp)
  subj <- paste0(position, replicate)
  # g features points of interest to build contrasts
  g <- factor(paste0(position, temperature),
              levels = c(paste0(levels(new_data_diff$position_global)[1],
                                levels(new_data_diff$temperature)
                                ),
                         paste0(levels(new_data_diff$position_global)[2],
                                levels(new_data_diff$temperature)
                                )
                         )
              )

  design <- model.matrix(~ 0 + g + replicate)
  colnames(design)[1:nlevels(g)] <- levels(g)

  # reshape data and order columns
  new_data_diff_matrix <- new_data_diff %>%
    tidyr::unite("key", position_global, biorep, temperature, sep = "_") %>%
    tidyr::spread(key, value) %>%
    dplyr::mutate(description = unname(sapply(description, IMPRINTS.CETSA.app:::getGeneName))) %>%
    tidyr::unite("id", protein, description, sep = "_")
  new_data_diff_matrix <- new_data_diff_matrix[,c("id", "treatment", paste(position, replicate, temperature, sep = "_"))]

  ##### Get weight matrix for lmFit
  message("Computing weight matrix from number of measurements")
  weights <- data_cleaved
  weights$description <- unname(sapply(weights$description, IMPRINTS.CETSA.app:::getGeneName))
  weights$cleaved_site <- NULL
  weights$Npep <- NULL

  ## extract number of measurements
  weights <- weights %>%
    dplyr::group_by(id, description, treatment, position_global) %>%
    dplyr::group_modify(~ {
      x <- lapply(strsplit(.x$Nvalue, "\\|")[[1]], function(n) strsplit(n, "-")[[1]])

      temp <- unlist(lapply(x, "[[", 1))
      biorep <- strsplit(unlist(lapply(x, "[[", 2)), ";")

      x <- mapply(function(b, t){
        b <- strsplit(b, "_")
        b <- t(as.data.frame(b))
        rownames(b) <- 1:nrow(b)
        colnames(b) <- c("biorep", "Nv")
        b <- as.data.frame(b)
        b$Nv <- as.numeric(b$Nv)
        b$temperature <- t;
        b
      },
      biorep, temp, SIMPLIFY = FALSE)

      x <- as.data.frame(Reduce(rbind, x))
      return(x)
    })

  weights <- weights %>%
    tidyr::unite("id", id, description, sep = "_") %>%
    tidyr::unite("key", position_global, biorep, temperature, sep = "_") %>%
    tidyr::spread(key, Nv)
  weights <- weights[,colnames(new_data_diff_matrix)]

  message("Computing p-values for each temperatures for each treatment...")
  res <- list()
  for(tr in treat){
    X <- new_data_diff_matrix[new_data_diff_matrix$treatment == tr,] %>%
      dplyr::ungroup() %>% dplyr::select(-treatment) %>%
      tibble::column_to_rownames("id") %>%
      as.matrix()
    Xw <- weights[weights$treatment == tr,] %>%
      dplyr::ungroup() %>% dplyr::select(-treatment) %>%
      tibble::column_to_rownames("id") %>%
      as.matrix()

    ### Add blocking on subject in duplicateCorrelation().
    corfit <- duplicateCorrelation(X, design, block=subj)
    ### Fit our model with weigths (not normalized)
    fitW <- lmFit(X, design, weights = Xw,
                  block=subj, correlation=corfit$consensus)

    temperature <- levels(new_data_diff$temperature)
    res[[tr]] <- list()
    for(tmp in temperature){
      # How do the two profiles differ in their response over temperature ?
      con.1 <- makeContrasts(contrasts = paste0("C", tmp, " - N", tmp),
                             levels = design)
      fitcon.1 <- contrasts.fit(fitW, con.1)
      fitcon.1 <- eBayes(fitcon.1)
      res[[tr]][[tmp]] <- topTable(fitcon.1, number=nrow(X), adjust.method="BH")
    }

    res[[tr]] <- as.data.frame(Reduce(rbind,
                                     mapply(function(x, n){
                                       x$temperature <- n
                                       x <- x %>%
                                         tibble::rownames_to_column("id") %>%
                                         tidyr::separate(id, into = c("id", "Gene"), sep = "_");
                                       x
                                     },
                                     res[[tr]], names(res[[tr]]),
                                     SIMPLIFY = FALSE)
                                     )
                               )

    res[[tr]]$treatment <- tr
  }
  res <- as.data.frame(Reduce(rbind, res))

  message("Finding most likely proteins with RESP effect")
  res <- res %>%
    dplyr::mutate(adj.P.Val = tidyr::replace_na(adj.P.Val, 1)) %>%
    dplyr::group_by(id, Gene, treatment) %>%
    dplyr::summarise(combined_pvalue = ifelse(length(na.omit(logFC)),
                                              metap::logitp(adj.P.Val)$p, # Goerge's method
                                              NA),
                     maxFC = ifelse(length(na.omit(logFC)),
                                    logFC[which(abs(logFC) == max(abs(logFC), na.rm = T))],
                                    NA)
                     )
  # FDR 1%
  cutoff <- res %>% dplyr::ungroup() %>% dplyr::group_by(treatment) %>%
    dplyr::mutate(BH = (order(order(combined_pvalue))/length(combined_pvalue))*FDR) %>%
    dplyr::summarise(pval = find_cutoff(combined_pvalue, BH),
                     maxFC_pos = diff_cutoff + median(maxFC[which(combined_pvalue < quantile(combined_pvalue, 0.5, na.rm = TRUE))], na.rm = TRUE),
                     maxFC_neg = -diff_cutoff - median(maxFC[which(combined_pvalue < quantile(combined_pvalue, 0.5, na.rm = TRUE))], na.rm = TRUE))

  res <- res %>% dplyr::group_by(id, Gene, treatment) %>%
    dplyr::mutate(criteria = combined_pvalue <= cutoff$pval[which(cutoff$treatment == treatment)] &
                    (maxFC >= cutoff$maxFC_pos[which(cutoff$treatment == treatment)] | maxFC <= cutoff$maxFC_neg[which(cutoff$treatment == treatment)]),
                  curve = curve(maxFC,
                                cutoff$maxFC_neg[which(cutoff$treatment == treatment)],
                                cutoff$maxFC_pos[which(cutoff$treatment == treatment)],
                                cutoff$pval[which(cutoff$treatment == treatment)],
                                curvature = curvature),
                  criteria_curve = -log10(combined_pvalue) >= curve
                  )
  res$criteria_curve <- tidyr::replace_na(res$criteria_curve, FALSE)

  n_treat <- length(treat)
  df_curve <- data.frame(maxFC = rep(seq(min(res$maxFC, na.rm = TRUE),
                                         max(res$maxFC, na.rm = TRUE), 0.01),
                                     n_treat)
                         )
  df_curve$treatment <- rep(treat, each = nrow(df_curve)/n_treat)
  df_curve <- df_curve %>% dplyr::group_by(treatment, rownames(df_curve)) %>%
    dplyr::mutate(curve = curve(maxFC,
                                cutoff$maxFC_neg[which(cutoff$treatment == treatment)],
                                cutoff$maxFC_pos[which(cutoff$treatment == treatment)],
                                cutoff$pval[which(cutoff$treatment == treatment)]
                                )
                  )

  message("Creating and saving plot")
  g_h <- ggplot(res, aes(maxFC, -log10(combined_pvalue), color = criteria_curve)) +
    geom_point() +
    geom_line(data = df_curve, aes(x = maxFC, y = curve), linetype = "dashed", color = "black") +
    ylim(c(0, max(-log10(res$combined_pvalue), na.rm = TRUE))) +
    xlim(c(-max(abs(res$maxFC), na.rm = TRUE),
           max(abs(res$maxFC), na.rm = TRUE))
         ) +
    labs(title = "RESP volcano plot",
         y = "-log10(combined p-value)",
         x = "maximum log2(FC)") +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70"))  +
    facet_wrap(~treatment) +
    ggrepel::geom_label_repel(data = res[res$criteria_curve,],
                              aes(maxFC, -log10(combined_pvalue), label = Gene),
                              show.legend = FALSE, min.segment.length = 0) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")

  ggsave(paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "RESP_hits_plot.png"),
         plot = g_h, device = "png",
         width = 16, height = 9)

  message("Preparing and saving results")
  res <- res[res$criteria_curve,c("id", "Gene", "treatment", "combined_pvalue", "maxFC")]
  data_cleaved$description <- unname(sapply(data_cleaved$description, IMPRINTS.CETSA.app:::getProteinName))
  data_cleaved$position_global <- paste0(data_cleaved$position_global, "-term")
  data_cleaved <- data_cleaved %>%
    tidyr::pivot_wider(id_cols = c("id", "description", "treatment", "cleaved_site"),
                       names_from = position_global,
                       values_from = c("Nvalue", "Npep"))
  res <- dplyr::left_join(res, data_cleaved, by = c("id", "treatment"),
                          relationship = "one-to-many")
  res <- res[order(res$combined_pvalue),]
  openxlsx::write.xlsx(res, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "RESP_hits_summary.xlsx"))

  message("Done !")
  return(res)
}


### function to compute discrete inflexion point
inflex <- function(x){
  x <- diff(x) # get first 'first derivative'
  # compute 'tendency' from 'first derivative', i.e. if increasing/decreasing
  df <- data.frame(x = 1:length(x),
                   y = x)
  res <- lm(y ~ x, data = df)
  res <- res$coeff[2]
  res <- sign(res)

  x <- diff(x) # get second 'first derivative'
  if(res == 1){
    inflex_point <- which.max(x) + 2 # double diff, so add 2 to index point
  }
  else{
    inflex_point <- which.min(x) + 2 # double diff, so add 2 to index point
  }

  if(inflex_point - 2 == length(x)){  # if inflex_point is last one, can't conclude that is cleaved
    inflex_point <- NA
  }

  return(inflex_point)
}

### function to sum number of values per peptide N and peptide C
nb_pepval <- function(x){
  x <- lapply(strsplit(x, "\\|"), strsplit, split = "-")

  x <- lapply(x,
              function(t){
                temp <- unlist(lapply(t, "[[", 1))
                biorep <- strsplit(unlist(lapply(t, "[[", 2)), ";")
                names(biorep) <- temp
                biorep <- lapply(biorep,
                                 function(b){
                                   b <- strsplit(b, "_")
                                   b_names <- unlist(lapply(b, "[[", 1))
                                   nv <- lapply(b, "[[", 2)
                                   nv <- lapply(nv, as.numeric)
                                   names(nv) <- b_names;
                                   nv
                                 });
                biorep
              })

  x <- as.data.frame(x, check.names = FALSE)
  x <- sapply(unique(names(x)), function(nv) sum(x[,names(x) == nv]))

  temp <- unique(sub("\\..*", "", names(x)))
  x <- sapply(temp,
              function(nv){
                b <- x[grep(nv, names(x))]
                names(b) <- sub(".*\\.", "", names(b))
                b <- paste(paste(names(b), b, sep = "_"), collapse = ";")
                b <- paste(nv, b, sep = "-");
                b
              })

  x <- paste(x, collapse = "|")
  return(x)
}

### function to find p-value cutoff
find_cutoff <- function(x,y){
  id <- order(x)
  x <- x[id]
  y <- y[id]

  x <- x[which(x > y)]
  return(x[1])
}

### function to compute cutoff curve
curve <- function(x, cut_neg, cut_pos, cut_p, curvature = 0.1){
  y <- rep(NA, length(x))
  neg <- which(x <= cut_neg)
  pos <- which(x >= cut_pos)
  y[neg] <- curvature/abs(x[neg] - cut_neg) + -log10(cut_p)
  y[pos] <- curvature/abs(x[pos] - cut_pos) + -log10(cut_p)

  return(y)
}


