#' imprints_barplotting_categorize_peptides
#'
#' Function to plot the categorized proteins found as hits by the function \code{imprints_cleaved_peptides}
#'
#'
#' @param data The normalized peptides data set, i.e. the outpout from \code{imprints_normalize_peptides}.
#' @param data_cleaved The categorized cleavage hits data set, i.e. the outpout from \code{imprints_categorize_peptides}.
#' @param treatment The treatment from which you want to see the plots. Can only be one.
#' @param control The control treatment from your dataset.
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
                                                     color = "red", pdfname = "categorized_RESP_barplots"){
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

  # summarized peptides
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


  message("Plotting...\n")
  # ordering data_cleaved according category and Gene
  data_cleaved <- data_cleaved[order(data_cleaved$Gene),]
  data_cleaved <- data_cleaved[order(as.numeric(factor(data_cleaved$category,
                                                       levels = c("RESP", "SP", "SPm", "MP", "MPm", "FP"))
                                                )
                                     ),]

  # actual plotting
  l <- list()
  for(p in data_cleaved$id){
    # plotting summary plot
    X <- data_diff_summary[which(data_diff_summary$Master.Protein.Accessions == p),]
    l[[paste0(p, "_summary")]] <- imprints_barplotting_peptides(X, RESP = TRUE, colorpanel = color)
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

  message("\nSaving pdf file...")
  ggpubr::ggexport(filename =  paste0(format(Sys.time(), "%y%m%d_%H%M_"), pdfname, ".pdf"),
                   plotlist = l, height = 12, width = 14)
  message("Done !")
}

