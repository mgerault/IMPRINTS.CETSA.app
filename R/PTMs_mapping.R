#' imprints_ptms_peptides
#'
#' Function to check if the hits found by \code{\link{imprints_cleaved_peptides}} could be due to PTMs and
#'  not  to protein cleavage.
#'
#' @details
#' When a hit is returned by \code{\link{imprints_cleaved_peptides}}, it means that this protein has a peptide position
#' where the IMPRINTS profiles of the two obtained parts are significantly different. This difference can be caused
#' by protein modification and mainly proteolysis but other more common PTMs can cause such a difference like
#' phosphorylation or ubiquitination. The aim here is to refilter the hit list and give the possible proteins which
#' could be modified based on the output of \code{\link{imprints_cleaved_peptides}}.
#'
#' @param data The normalized peptides data set, i.e. the outpout from \code{imprints_normalize_peptides}.
#' @param data_cleaved The cleavage hits data set, i.e. the outpout from \code{imprints_cleaved_peptides}.
#' @param control The control treatment from your dataset.
#' @param minref The minimum number of references reporting the PTM in PhosphoSitePlus; Default is 2
#' @param PTM_FC_cutoff The minimum value for the maximum log2 fold-change of a peptide to be considered
#'   significantly modified
#' @param PTM_direction Character to set if you look for PTMs in higher or lower abundance or both you can
#'   select respectively hyper, hypo or both. Default is hyper.
#'   Since the function assume you didn't add the PTM as a dynamic modification during your protein identification search,
#'   if a peptide is "hypermodified' then peptide with negative value will be selected and conversely if a peptide is
#'   'hypomodified' peptide with positive value will be kept.
#' @param save_xlsx Logical to tell if you want to save the categorized hits in an xlsx file.
#'   Default to TRUE.
#' @param xlsxname The name of your saved file.
#'
#'
#' @return A dataframe containing the potential proteins having PTMs
#'  instead of the protein being cleaved
#'
#' @export
#'
#' @seealso \code{\link{imprints_cleaved_peptides}}
#'

imprints_ptms_peptides <- function(data, data_cleaved, control, minref = 2,
                                   PTM_FC_cutoff = 0.25, PTM_direction = c("hyper", "hypo", "both"),
                                   save_xlsx = TRUE, xlsxname = "RESP_PTMs_mapping"){
  PTM_direction <- match.arg(PTM_direction)

  if(!(control %in% get_treat_level(data))){
    message(paste(control, "is not in your data !"))
    return()
  }

  if(control %in% unique(data_cleaved$treatment)){
    message(paste("Error:", control, "shouldn't be in your data_cleaved !"))
    return()
  }

  # fetch isoform from each hits from data_cleaved
  message("Fetching PTMs...")
  ptms <- fetch_ptms(unique(data_cleaved$id))
  ptms <- ptms[which(ptms$nbref >= minref),]

  if(nrow(ptms) == 0){
    message("No PTMs could have been found for any of your proteins")
    return()
  }

  message(paste(length(unique(ptms$id)), "proteins were found to have known PTMs"))

  ### Map fetched PTMs to RESP data
  if(all(c("category", "details") %in% colnames(data_cleaved))){
    data_ptms <- data_cleaved[,c("id", "Gene", "description", "treatment", "cleaved_site",
                                 "RESP_score", "category", "details")]
  }
  else{
    data_ptms <- data_cleaved[,c("id", "Gene", "description", "treatment", "cleaved_site", "RESP_score")]
  }
  data_ptms <- full_join(data_ptms, ptms, by = "id", relationship = "many-to-many")

  ### computing fold-changes of potential PTMs
  data <- data[which(!is.na(match(data$Master.Protein.Accessions, unique(data_ptms$id)))),]

  directory_toremove <- list.files()
  data <- imprints_sequence_peptides(data, control = control) # will save txt file
  directory_toremove <- list.files()[!(list.files() %in% directory_toremove)]
  # remove this file
  directory_toremove <- grep("(?=^\\d{6}_\\d{4}_.*)(?=.*\\.txt$)", directory_toremove, value = TRUE, perl = TRUE)
  if(length(directory_toremove) == 1){
    unlink(directory_toremove, recursive = TRUE)
  }
  data <- data[,-grep(paste0("_", control, "$"), colnames(data))]

  message("Comparing PTMs candiates with FC data...")
  # summarizing peptide with their maximum FC value
  data <- data %>%
    dplyr::select(-Annotated.Sequence, -countNum) %>%
    dplyr::mutate(Gene = sapply(description, IMPRINTS.CETSA.app:::getGeneName, USE.NAMES = FALSE),
                  description = sapply(description, IMPRINTS.CETSA.app:::getProteinName, USE.NAMES = FALSE)
                  ) %>%
    tidyr::gather("key", "value", -Master.Protein.Accessions, -Gene, -description,
                  -Positions.in.Master.Proteins, -Modifications) %>%
    tidyr::separate(key, into = c("temperature", "biorep", "treatment"), sep = "_") %>%
    dplyr::group_by(Master.Protein.Accessions, Gene, description,
                    Positions.in.Master.Proteins, Modifications,
                    temperature, treatment) %>%
    dplyr::summarise(value = mean(value, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Master.Protein.Accessions, Gene, description,
                    Positions.in.Master.Proteins, Modifications,
                    treatment) %>%
    dplyr::summarise(value = value[which.max(abs(value))]) %>%
    dplyr::ungroup() %>% dplyr::group_by(Master.Protein.Accessions, Gene, description, treatment) %>%
    dplyr::mutate(npep = length(unique(Positions.in.Master.Proteins)))
  colnames(data)[grep("Master.Protein.Accessions", colnames(data))] <- "id"

  # summarizing quantified PTMs
  data$Modifications <- mapply(function(m, p){
    m <- strsplit(m, "\\]; ")[[1]]
    m <- m[-grep("TMT(pro|6plex)", m)]
    if(length(m)){
      p <- sub(";.*", "", p)
      p <- gsub(".* \\[|-.*", "", p)
      p <- as.numeric(p) - 1

      m_change <- gsub(".*\\[|\\]", "", m)
      m_change <- unlist(strsplit(m_change, "; "))
      m_change <- sapply(m_change,
                         function(y){
                           yp <- sub("^.", "", y)
                           y <- sub(yp, "", y)
                           if(nchar(yp)){
                             yp <- as.numeric(yp) + p
                           }
                           else{
                             yp <- p
                           };
                           paste0(y, yp)
                         })

      m <- paste(m, collapse = "]; ")
      if(!grepl("\\]$", m)){
        m <- paste0(m, "]")
      }

      for(i in names(m_change)){
        m <- sub(paste0("(?<=(\\[| ))", i), m_change[[i]], m, perl = TRUE)
      }
    }
    else{
      m <- NA
    };
    m
  }, data$Modifications, data$Positions.in.Master.Proteins, USE.NAMES = FALSE)

  data$Modifications.Position <- sapply(data$Modifications,
                                                 function(x){
                                                   if(!is.na(x)){
                                                     x <- strsplit(x, "\\]; ")[[1]]
                                                     x <- gsub(".* \\[|\\]", "", x)
                                                     x <- paste(x, collapse = "; ")
                                                     x <- strsplit(x, "; ")[[1]]
                                                     x <- x[order(as.numeric(sub("^.", "", x)))]
                                                     x <- paste(x, collapse = "; ")
                                                   };
                                                   x
                                                 }, USE.NAMES = FALSE)
  data$Modifications <- sapply(data$Modifications,
                                        function(x){
                                          if(!is.na(x)){
                                            x <- strsplit(x, "\\]; ")[[1]]
                                            x <- sub(" .*", "", x)
                                            x <- sub("^\\d{1}x", "", x)
                                            x <- x[order(x)]
                                            x <- paste(x, collapse = "; ")
                                          };
                                          x
                                        }, USE.NAMES = FALSE)

  # mapping and filtering with PTMs data
  data_ptms <- dplyr::inner_join(data, data_ptms, by = c("id", "Gene", "description", "treatment"), relationship = "many-to-many")
  data_ptms <- data_ptms[which(abs(data_ptms$value) >= PTM_FC_cutoff),]

  if(nrow(data_ptms)){
    # filtering peptides containing PTMs
    ### if PTMs caused by treatment and that the PTM is not added in the protein identification search,
    # then the unmodified peptide should be found to be significantly less abundant than the modified one
    # i.e. the FC from the corresponding treatment should be significantly negative and inversely
    data_ptms <- data_ptms %>%
      dplyr::ungroup() %>% dplyr::group_by(id, Gene, description, Positions.in.Master.Proteins,
                                           Modifications, Modifications.Position, treatment, PTM_position) %>%
      mutate(pep_pos = gsub(".* \\[|\\]", "", sub(";.*", "",  Positions.in.Master.Proteins)),
             pep_pos1 = as.numeric(strsplit(pep_pos, "-")[[1]][1]),
             pep_pos2 = as.numeric(strsplit(pep_pos, "-")[[1]][2]),
             pep_in = pep_pos1 - 2 <= PTM_position & PTM_position <= pep_pos2 + 2) # add margin of 2 --> useful when modification on K

    if(PTM_direction == "hyper"){
      data_ptms <- data_ptms %>%
        dplyr::filter((pep_in & value < 0) | !is.na(Modifications))
    }
    else if(PTM_direction == "hypo"){
      data_ptms <- data_ptms %>%
        dplyr::filter((pep_in & value > 0) | !is.na(Modifications))
    }
    else if(PTM_direction == "both"){
      data_ptms <- data_ptms %>%
        dplyr::filter(pep_in | !is.na(Modifications))
    }

    if(nrow(data_ptms)){
      data_ptms <- data_ptms %>%
        dplyr::select(-pep_pos, -pep_pos1, -pep_pos2) %>%
        dplyr::ungroup() %>% dplyr::group_by(id, Gene, description, Positions.in.Master.Proteins,
                                             Modifications, Modifications.Position, treatment) %>%
        dplyr::group_modify(~{
          if(all(!.x$pep_in, na.rm = TRUE)){ # NA means there were no mapped PTMs but quantified modified peptide
            df <- .x[1,]
            df$pep_in <- NULL
            df[,c("PTM", "PTM_position", "PTM_aa", "nbref")] <- NA
          }
          else{
            df <- .x
            if(any(!df$pep_in | is.na(df$pep_in))){
              df[which(!df$pep_in | is.na(df$pep_in)),
                 c("PTM", "PTM_position", "PTM_aa", "nbref")] <- NA
            }

            if(PTM_direction == "hyper"){
              if(any(df$pep_in & df$value > 0)){
                df[which(df$pep_in & df$value > 0),c("PTM", "PTM_position", "PTM_aa", "nbref")] <- NA
              }
            }
            else if(PTM_direction == "hypo"){
              if(any(df$pep_in & df$value < 0)){
                df[which(df$pep_in & df$value < 0),c("PTM", "PTM_position", "PTM_aa", "nbref")] <- NA
              }
            }

            if(!(all(is.na(df$PTM)))){ # avoid repetition of information
              df <- df[-which(is.na(df$PTM)),]
            }

            df$pep_in <- NULL
            df <- unique(df)
          }
          return(df)
        })  %>%
      dplyr::ungroup() %>% dplyr::group_by(id, Gene, description, treatment) %>%
        mutate(npep_modif = length(unique(Positions.in.Master.Proteins)),
               ppep_modif = npep_modif/npep)

      if(save_xlsx){
        openxlsx::write.xlsx(data_ptms, paste0(format(Sys.time(), "%y%m%d_%H%M"), "_", xlsxname, ".xlsx"))
      }
      message("Done !")
      return(data_ptms)
    }
    else{
      message("No RESP hits could be due to PTMs")
      return()
    }

  }
  else{
    message("No RESP hits could be due to PTMs")
    return()
  }
}

### Function to retrieve all PTMS from proteins from the PhosphoSitePlus database
fetch_ptms <- function(proteins){
  nprot <- length(proteins)
  ptms <- sapply(proteins,
                 function(p){
                   # get PTMs reported by PhosphoSitePlus from the protein p
                   ptm_p <- NULL

                   message(paste0("Requesting PTMs from ", p, " ", which(proteins == p), "/", nprot))
                   for(i in 1:15){ # try to fetch results
                     h <- httr::handle("https://www.phosphosite.org/")
                     x <- httr::GET(handle = h, url = NULL,
                                    path = paste0("uniprotAccAction?id=", sub(";.*","",p)
                                                  )) # if protein group, only taking first one
                     status_request <- httr::status_code(x)

                     if(status_request == 200){
                       ptm_p <- httr::content(x, "text")
                       if(grepl("PTMsites", ptm_p)){
                         rm(h, x) # allow to close connection
                         break
                       }
                     }
                     else if(status_request == 429){
                       message(paste("Too many requests, waiting 5 minute(s) before making another request for", p))
                       rm(h, x) # allow to close connection
                       Sys.sleep(300)
                       message(paste0("Start requesting PTMs again for ", p, " ", which(proteins == p), "/", nprot))
                     }
                     else if(status_request == 500){
                       message(paste("Internal server error, waiting 5 minute(s) before making another request for", p))
                       rm(h, x) # allow to close connection
                       Sys.sleep(300)
                       message(paste0("Start requesting PTMs again for ", p, " ", which(proteins == p), "/", nprot))
                     }
                     else{
                       # break between requests
                       rm(h, x) # allow to close connection
                       Sys.sleep(sample(1:5, 1))
                     }
                   }

                   if(is.null(ptm_p) || !grepl("PTMsites", ptm_p)){
                     ptm_p <- NULL
                     m <- paste(p, "didn't retrieve any PTMs:")
                     m <- paste(m, ifelse(status_request != 200,
                                          paste("Request error", status_request),
                                          "no PTMs were found for this protein"))
                     message(m)
                   }
                   else{
                     ptm_p <- sub(".*PTMsites", "", ptm_p)
                     ptm_p <- sub("\\].*", "", ptm_p)
                     ptm_p <- sub(".*value=.\\[", "", ptm_p)
                     ptm_p <- gsub("^\\{|\\}$", "", ptm_p)
                     ptm_p <- strsplit(ptm_p, "\\},\\{")[[1]]
                     ptm_p <- gsub(":&quot;", "", ptm_p)

                     if(all(!grepl(";ID&quot;", ptm_p))){
                       ptm_p <- NULL
                       message(paste("No PTMs are reported in PhosphoSitePlus for", p))
                     }
                     else{
                       ptm_p <- sapply(ptm_p,
                                       function(y){
                                         pos <- sub("&quot;.*", "", sub(".*;ID&quot;", "", y));
                                         data.frame(PTM = sub("&quot;.*", "", sub(".*;MODIFICATION&quot;", "", y)),
                                                    PTM_position = as.numeric(sub("^.", "", pos)),
                                                    PTM_aa = sub("\\d{1,}", "", pos),
                                                    nbref = as.numeric(sub(",&quot;.*", "", sub(".*;REF&quot;:", "", y)))
                                                    )
                                         }, simplify = FALSE, USE.NAMES = FALSE)
                       ptm_p <- as.data.frame(Reduce(rbind, ptm_p))

                       ptm_p$id <- p
                       ptm_p <- ptm_p[,c(5,1:4)]
                     }
                   };
                   ptm_p
                   }, simplify = FALSE, USE.NAMES = FALSE)

  ptms <- as.data.frame(Reduce(rbind, ptms))

  return(ptms)
}
