#' find_in_pubmed
#'
#' Function to check if there is already some publications in Pubmed on your subject.
#' If so, export word files of the publications
#'
#' @param data The word you want to search a publication on. It can be a file path, a character vector, a list or
#'             the dataframe hitlist output from the Hitlist function.
#' @param feat A second word that will be match to every word from data, can be NULL.
#' @param imp_by_hitlist Logical to tell if the data are in the hitlist format.
#' @param condition A character vector for selecting specific condition. Use when data are in hitlist format.
#' @param language A character to tell if you want to select publications in a specific language like 'english'.
#' @param year_rg A character to tell a range of year of publication. The format is '2000:2020'
#' @param your_API A character taht is your NCBI API key. If you don't have an account, use NULL.
#' @param newfolder_name The name of the folder that will be created for saving the word documents
#'
#'
#' @return Dataframe which contains the proteins and if they give a result or not
#'
#' @import pubmedR
#' @import bibliometrix
#' @import officer
#'
#' @export
#'
#' @examples
#' library(mineCETSAapp)
#' find_in_pubmed("Lamin B1", imp_by_hitlist = FALSE)
#'
#' @seealso \code{\link{hitlist}} and \code{\link{pmApi2df}} for more details
#'

find_in_pubmed <- function(data, feat = "PI3K", imp_by_hitlist = FALSE, condition = "",
                           language = NULL, year_rg = NULL, your_API = NULL,
                           newfolder_name = "PI3K_pubmed_search"){

  path <- getwd()
  if(!file.exists(file.path(path, newfolder_name))){
    dir.create(file.path(path, newfolder_name))
  }

  no_res <- c()

  if (typeof(data) == "character"){
    if(length(data) == 1 & all(str_detect(data, "\\.txt$|\\.csv$|\\.xlsx$")))
      data <- import_list(data)[[1]]
  }

  if (imp_by_hitlist){
    condi <- unique(data$Condition)
    if (str_length(condition) != 0){
      condi <- condition
    }
    prot_description_list <- unique(as.character(sapply(data$description[which(!is.na(match(data$Condition, condi)))],
                                                        mineCETSA:::getProteinName)))

    gene_list <- unique(as.character(sapply(data$description[which(!is.na(match(data$Condition, condi)))],
                                            mineCETSA:::getGeneName)))
    prot_list <- unique(data$id[which(!is.na(match(data$Condition, condi)))])
  }
  else{
    prot_description_list <- data
  }
  prot_description_list <- str_replace_all(prot_description_list, "-", " ")
  prot_description_list <- str_replace_all(prot_description_list, "/", " ")

  LA <- NULL
  DP <- NULL
  if(!is.null(language)){
    LA <- paste0(language, "[LA]")
  }
  if(!is.null(year_rg)){
    DP <- paste0(year_rg, "[DP]")
  }
  T_AB2 <- paste0(feat, "*[Title/Abstract]")

  have_paper <- c()
  for(i in prot_description_list){
    T_AB1 <- paste0(i, "*[Title/Abstract]")

    query <- paste(c(T_AB1, T_AB2, LA, DP), collapse = " AND ")

    res <- pmQueryTotalCount(query = query, api_key = your_API)

    D <- pmApiRequest_m(query = query, limit = res$total_count, api_key = your_API)

    doc_1 <- read_docx()
    if(D$total_count != 0){
      M <- pmApi2df_m(D)
      M <- convert2df_m(D, dbsource = "pubmed", format = "api")

      M_interest <- M[, c("TI", "PY", "AU", "AB")]

      for (k in 1:nrow(M_interest)){
        AB <- tolower(M_interest$AB[k])
        AB <- str_split(AB, "\\.")[[1]]
        AB <- str_replace_all(AB, "^ ", "")
        AB <- str_replace(AB, "^.", toupper(str_sub(AB, 1,1)))
        AB <- paste(AB, collapse = ". ")

        doc_1 %>%
          body_add_par(paste(M_interest$TI[k], M_interest$PY[k], sep = ", "), style = "heading 1") %>%
          body_add_par("", style = "Normal") %>%
          body_add_par(M_interest$AU[k], style = "Normal") %>%
          body_add_par("", style = "Normal") %>%
          body_add_par("Abstract", style = "heading 2") %>%
          body_add_par("", style = "Normal") %>%
          body_add_par(AB, style = "Normal")

      }

      print(doc_1, target = paste0(newfolder_name, "/", i, ".docx"))

      have_paper <- append(have_paper, TRUE)
    }
    else{
      no_res <- append(no_res, i)
      doc_1 %>%
        body_add_par("No results found", style = "heading 1")

      have_paper <- append(have_paper, FALSE)
    }

  }

  if (!is.null(no_res)){
    message(paste(length(no_res), "/",
                  length(prot_description_list),
                  "proteins didn't give any results"))
  }
  message("All results have been saved !")

  if (imp_by_hitlist){
    return(data.frame("UniprotID" = prot_list,
                      "Gene.name" = gene_list,
                      "Protein_name" =  prot_description_list,
                      "Have_publication" = have_paper))
  }
  else{
    return(data.frame("Protein_name" =  prot_description_list,
                      "Have_publication" = have_paper))
  }
}




### exact same function as pmApiRequest from pubmedR package,
### just had to change 's <- 0' in order to handle if only one result is found
pmApiRequest_m <-  function (query, limit, api_key = NULL){
  res <- pmQueryTotalCount(query = query, api_key = NULL)
  n <- min(res$total_count, res$total_count)
  step <- 200
  step <- min(res$total_count, step)
  metadata <- list()
  stop <- FALSE
  s <- 0
  while (!isTRUE(stop)) {
    cat("Documents ", s + step, " of ", n, "\n")
    multi_summs <- rentrez::entrez_fetch(db = "pubmed", web_history = res$web_history,
                                         retstart = s, retmax = step, rettype = "xml", parsed = T,
                                         api_key = NULL)
    multi_summs <- XML::xmlToList(multi_summs, simplify = F)
    metadata <- c(metadata, multi_summs)
    if (n <= (s + step)) {
      stop <- TRUE
    }
    else {
      s <- s + step
      if ((s + step) > limit) {
        step <- (n - s + 1)
      }
    }
  }
  P <- list(data = metadata, query = query, query_translation = res$query_translation,
            records_downloaded = n, total_count = res$total_count)
  return(P)
}

### exact same function as pmApi2df from pubmedR package,
### just had to start from 0 in the txtProgressBar for the same reason as above
pmApi2df_m <- function (P, format = "bibliometrix")
{
  P <- P$data
  n <- length(P)
  df <- data.frame(AU = rep(NA, n), AF = "NA", TI = "NA",
                   SO = "NA", SO_CO = NA, LA = NA, DT = NA, DE = NA, ID = NA,
                   MESH = NA, AB = "NA", C1 = NA, CR = "NA", TC = NA, SN = NA,
                   J9 = NA, JI = NA, PY = NA, PY_IS = NA, VL = NA, DI = NA,
                   PG = NA, GRANT_ID = NA, GRANT_ORG = NA, UT = NA, PMID = NA,
                   DB = "PUBMED", AU_UN = NA, stringsAsFactors = FALSE)
  pb <- utils::txtProgressBar(min = 0, max = n, initial = 0,
                              char = "=")
  for (i in 1:n) {
    utils::setTxtProgressBar(pb, i)
    a <- pubmedR:::list2char(P[[i]])
    items <- names(a)
    df$LA[i] <- a["MedlineCitation.Article.Language"]
    df$DT[i] <- a["MedlineCitation.Article.PublicationTypeList.PublicationType.text"]
    df$TI[i] <- a["MedlineCitation.Article.ArticleTitle"]
    ind <- which(items == "PubmedData.History.PubMedPubDate.Year")
    if (length(ind) > 0) {
      df$PY[i] <- min(as.numeric(a[ind]), na.rm = TRUE)
    }
    else {
      df$PY[i] = NA
    }
    df$PY_IS[i] <- a["MedlineCitation.Article.Journal.JournalIssue.PubDate.Year"]
    AU_last_ind <- which(items == "MedlineCitation.Article.AuthorList.Author.LastName")
    AU_first_ind <- which(items == "MedlineCitation.Article.AuthorList.Author.ForeName")
    AU_init_ind <- which(items == "MedlineCitation.Article.AuthorList.Author.Initials")
    nameAF <- paste(a[AU_last_ind], a[AU_first_ind], sep = ", ")
    nameAU <- paste(a[AU_last_ind], a[AU_init_ind], sep = " ")
    df$AF[i] <- paste(nameAF, collapse = ";")
    df$AU[i] <- paste(nameAU, collapse = ";")
    Aff_name_ind <- which(items == "MedlineCitation.Article.AuthorList.Author.AffiliationInfo.Affiliation")
    Affiliations <- a[Aff_name_ind]
    Affiliations <- lapply(Affiliations, function(l) {
      l <- unlist(strsplit(l, ", "))
      l <- paste(l[!(regexpr("\\@", l) > -1)], collapse = ", ")
    })
    df$C1[i] <- df$AU_UN[i] <- paste(Affiliations, collapse = ";")
    DE_ind <- which(items == "MedlineCitation.KeywordList.Keyword.text")
    df$DE[i] <- paste(a[DE_ind], collapse = ";")
    ID_ind <- which(items == "MedlineCitation.MeshHeadingList.MeshHeading.DescriptorName.text")
    df$ID[i] <- df$MESH[i] <- paste(a[ID_ind], collapse = ";")
    ind <- which(items %in% "MedlineCitation.Article.Abstract.AbstractText.text")
    if (length(ind) > 0) {
      df$AB[i] <- paste(a[ind], collapse = " ")
    }
    else {
      ind <- which(items %in% "MedlineCitation.Article.Abstract.AbstractText")
      if (length(ind) > 0) {
        df$AB[i] <- a[ind]
      }
    }
    df$SO[i] <- a["MedlineCitation.Article.Journal.Title"]
    df$JI[i] <- df$J9[i] <- a["MedlineCitation.Article.Journal.ISOAbbreviation"]
    df$SO_CO[i] <- a["MedlineCitation.MedlineJournalInfo.Country"]
    doi_ind <- which(items == "PubmedData.ArticleIdList.ArticleId..attrs.IdType")
    ind <- which(a[doi_ind] == "doi")
    if (length(ind) > 0) {
      doi_ind <- doi_ind[ind] - 1
      df$DI[i] <- a[doi_ind]
    }
    df$SN[i] <- a["MedlineCitation.Article.Journal.ISSN.text"]
    df$PG[i] <- a["MedlineCitation.Article.Pagination.MedlinePgn"]
    df$VL[i] <- a["MedlineCitation.Article.Journal.JournalIssue.Volume"]
    df$UT[i] <- df$PMID[i] <- a["MedlineCitation.PMID.text"]
    GR_ID <- which(items %in% "MedlineCitation.Article.GrantList.Grant.GrantID")
    df$GRANT_ID[i] <- paste(a[GR_ID], collapse = ";")
    GR_ORG <- which(items %in% "MedlineCitation.Article.GrantList.Grant.Agency")
    df$GRANT_ORG[i] <- paste(a[GR_ORG], collapse = ";")
  }
  if (format == "bibliometrix") {
    DI <- df$DI
    df <- data.frame(lapply(df, toupper), stringsAsFactors = FALSE)
    df$DI <- DI
    df$AU_CO = "NA"
    df$AU1_CO = "NA"
  }
  df$PY <- as.numeric(df$PY)
  df$TC <- as.numeric(df$TC)
  df$TC[is.na(df$TC)] <- 0
  df = df[!is.na(df$DT), ]
  close(pb)
  return(df)
}

### exact same function as convert2df from bibliometrix package,
### just had to start from 0 in the settxtProgressbar for the same reason as above
convert2df_m <- function (file, dbsource = "wos", format = "plaintext")
{
  allowed_formats <- c("api", "bibtex", "csv", "endnote",
                       "excel", "plaintext", "pubmed")
  allowed_db <- c("cochrane", "dimensions", "generic", "isi",
                  "pubmed", "scopus", "wos", "lens")
  cat("\nConverting your", dbsource, "collection into a bibliographic dataframe\n\n")
  if (length(setdiff(dbsource, allowed_db)) > 0) {
    cat("\n 'dbsource' argument is not properly specified")
    cat("\n 'dbsource' argument has to be a character string matching one among:",
        allowed_db, "\n")
  }
  if (length(setdiff(format, allowed_formats)) > 0) {
    cat("\n 'format' argument is not properly specified")
    cat("\n 'format' argument has to be a character string matching one among:",
        allowed_formats, "\n")
  }
  if (length(setdiff(format, c("api", "plaintext", "bibtex",
                               "csv", "excel", "endnote"))) > 0) {
    D <- importFiles(file)
    D <- iconv(D, "latin1", "ASCII", sub = "")
  }
  if (dbsource == "wos")
    dbsource <- "isi"
  if (format == "endnote")
    format <- "plaintext"
  if (format == "lens")
    format <- "csv"
  switch(dbsource, isi = {
    switch(format, bibtex = {
      D <- importFiles(file)
      M <- bib2df(D, dbsource = "isi")
    }, plaintext = {
      D <- importFiles(file)
      M <- isi2df(D)
    })
  }, scopus = {
    switch(format, bibtex = {
      D <- importFiles(file)
      M <- bib2df(D, dbsource = "scopus")
    }, csv = {
      M <- csvScopus2df(file)
    })
  }, generic = {
    D <- importFiles(file)
    M <- bib2df(D, dbsource = "generic")
  }, lens = {
    M <- csvLens2df(file)
  }, pubmed = {
    switch(format, api = {
      M <- pmApi2df_m(file)
      M$DB <- "PUBMED"
    }, {
      D <- importFiles(file)
      M <- pubmed2df(D)
    })
  }, cochrane = {
    D <- importFiles(file)
    M <- cochrane2df(D)
  }, dimensions = {
    switch(format, api = {
      M <- dsApi2df(file)
      M$DB <- "DIMENSIONS"
    }, {
      M <- dimensions2df(file, format = format)
    })
  })
  if ("PY" %in% names(M)) {
    M$PY = as.numeric(M$PY)
  }
  else {
    M$PY <- NA
  }
  if ("TC" %in% names(M)) {
    M$TC = as.numeric(M$TC)
    M$TC[is.na(M$TC)] <- 0
  }
  else {
    M$TC <- 0
  }
  if (!("CR" %in% names(M))) {
    M$CR = "none"
  }
  else {
    M$CR <- trim.leading(trimES(gsub("\\[,||\\[||\\]|| \\.\\. || \\. ",
                                     "", M$CR)))
  }
  if (dbsource != "cochrane") {
    M$AU = gsub(intToUtf8(8217), intToUtf8(39), M$AU)
  }
  cat("Done!\n\n")
  if (!(dbsource %in% c("dimensions", "pubmed", "lens"))) {
    if ("C1" %in% names(M)) {
      cat("\nGenerating affiliation field tag AU_UN from C1:  ")
      M <- metaTagExtraction(M, Field = "AU_UN")
      cat("Done!\n\n")
    }
    else {
      M$C1 = NA
      M$AU_UN = NA
    }
    M$AU = unlist(lapply(strsplit(M$AU, ";"), function(x) {
      x = trimws(trimES(gsub("[^[:alnum:][-]']", " ",
                             x)))
      x = paste(x, collapse = ";")
    }))
  }
  if ((dbsource == "pubmed") & (format == "pubmed")) {
    if ("C1" %in% names(M)) {
      cat("\nGenerating affiliation field tag AU_UN from C1:  ")
      M <- metaTagExtraction(M, Field = "AU_UN")
      cat("Done!\n\n")
    }
    else {
      M$C1 = NA
      M$AU_UN = NA
    }
  }
  suppressWarnings(M <- metaTagExtraction(M, Field = "SR"))
  d <- duplicated(M$SR)
  if (sum(d) > 0)
    cat("\nRemoved ", sum(d), "duplicated documents\n")
  M <- M[!d, ]
  row.names(M) <- M$SR
  class(M) <- c("bibliometrixDB", "data.frame")
  return(M)
}




