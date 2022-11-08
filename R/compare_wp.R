#' compare_wp
#'
#' Function to run enrichment analysis on your hits from the wiki-pathway database and return
#' a plot that compare the pathways found between your conditions.
#'
#' @param hits The hitlist; a data.frame containing the genes id and preferably a condition column but not necessary.
#' @param gene_column The name of the coulumn that contains the genes. Default is 'Genes'.
#' @param condition_column The name of the column that contains the conditions. Default is NULL.
#' @param n_pathway Number of pathway to show on plot. Default is 5.
#'                  For more info you, see \code{\link{compareCluster}}.
#'
#' @return A list that contains the results and the plot.
#'
#' @details For now, only Homo Sapiens is supported.
#'
#' @export
#'
#' @seealso \code{\link{clusterProfiler}}

compare_wp <- function(hits, gene_column = "Genes", condition_column = NULL,
                       n_pathway = 5){
  ### load database
  wp <- get_wikipath() # wiki pathway
  ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl") # genes_id / gene symbols

  # get genes id from gene symbol
  hits_gene_id <- biomaRt::getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
                                 filters = "hgnc_symbol",
                                 curl = curl::handle_setopt(new_handle(), timeout = 30000),
                                 values = unique(sort(hits[[gene_column]])),
                                 bmHeader = TRUE,
                                 mart = ensembl)
  colnames(hits_gene_id) <- c(gene_column, "Genes_id")

  hits <- dplyr::left_join(hits, hits_gene_id, by = gene_column)

  if(is.null(condition_column)){
    hits$Condition <- "condition"
  }
  else{
    colnames(hits)[stringr::str_which(colnames(hits), condition_column)] <- "Condition"
  }

  # cluster compare enrichment analysis
  hits_enrich <- clusterProfiler::compareCluster(Genes_id~Condition,
                                                 data = hits, fun = clusterProfiler::enricher,
                                                 TERM2GENE=wp[,c("wpid", "gene")],
                                                 TERM2NAME=wp[,c("wpid", "name")])


  res <- fortify(hits_enrich, showCategory = 6)
  res$Condition <- factor(res$Condition, levels = unique(res$Condition))
  graph <- ggplot(res, aes(Condition, Description,
                color = p.adjust,
                size = GeneRatio)) +
    geom_point() +
    scale_y_discrete(labels = enrichplot:::default_labeller(30)) +
    scale_color_continuous(low = "#B2ECBF", high = "#3821A3", name = "p.adjust",
                           guide = guide_colorbar(reverse = TRUE)) +
    scale_size(range = c(3,8)) +
    DOSE::theme_dose(12) +
    theme(axis.title.x = element_blank())

  res$geneSymbol <- unlist(lapply(strsplit(res$geneID, "/"), function(x){x <- as.numeric(x)
  g <- hits_gene_id$Genes[which(!is.na(match(hits_gene_id$Genes_id, x)))]
  g <- paste(g, collapse = "/");
  g
  }
  ))

  return(list("res" = res, "graph" = graph))
}


# function to always get most recent wiki pathway database (update every 10 of month)
get_wikipath <- function(wp = TRUE){
  date <- stringr::str_split(Sys.Date(), "-")[[1]]
  date <- as.numeric(date)
  if(date[3] < 10){
    date[2] <- date[2] - 1
  }
  date[3] <- 10
  if(date[2] == 0){
    date[2] <- 12
    date[1] <- date[1] - 1
  }
  else if(date[2] < 10){
    date[2] <- paste0("0", date[2])
  }
  date <- paste0(date, collapse = "")

  url_wiki <- paste0("https://wikipathways-data.wmcloud.org/current/gmt/wikipathways-",
                     date, "-gmt-Homo_sapiens.gmt")
  url_wiki <- url(url_wiki)

  if(wp){
    wikipath <- clusterProfiler::read.gmt.wp(url_wiki)
  }
  else{
    wikipath <- clusterProfiler::read.gmt(url_wiki)
  }
  close(url_wiki)
  return(wikipath)
}
