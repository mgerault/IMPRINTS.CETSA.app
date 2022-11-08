#' compare_wp
#'
#' Function to run enrichment analysis on your hits from the wiki-pathway database and return
#' a plot that compare the pathways found between your conditions.
#'
#' @param hits A data.frame containing the genes id, a score value and preferably a condition column but not necessary.
#' @param gene_column The name of the coulumn that contains the genes. Default is 'Genes'.
#' @param score_column The name of the coulumn that contains the score values. Default is 'SR'.
#' @param condition_column The name of the column that contains the conditions. Default is NULL.
#' @param condition The name of the condition you ant to keep. Default is NULL.
#' @param pval_cutoff The p-value cutoff for the gene concept network.
#'
#' @return The gene concept network plot.
#'
#' @details For now, only Homo Sapiens is supported.
#'
#' @export
#'
#' @seealso \code{\link{clusterProfiler}}

gene_conceptNet <- function(hits, gene_column = "Genes", score_column = "SR",
                            condition_column = NULL, condition = NULL,
                            pval_cutoff = 0.01){
  ### load database
  wp <- get_wikipath(FALSE) # wiki pathway
  ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl") # genes_id / gene symbols

  if(!is.null(condition_column) & !is.null(condition)){
    hits <- hits[which(hits[[condition_column]] == condition),]
  }
  if(any(is.na(hits[[score_column]]))){
    hits <- hits[which(!is.na(hits[[score_column]])),]
  }

  # get genes id from gene symbol
  hits_gene_id <- biomaRt::getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
                                 filters = "hgnc_symbol",
                                 values = unique(sort(hits[[gene_column]])),
                                 bmHeader = TRUE,
                                 mart = ensembl)
  colnames(hits_gene_id) <- c(gene_column, "Genes_id")

  hits <- dplyr::left_join(hits, hits_gene_id, by = gene_column)

  # start enrichment analysis
  hits_enrich <- clusterProfiler::enricher(hits$Genes_id,
                          TERM2GENE = wp,
                          pvalueCutoff = 0.05)
  hits_enrich@result$geneSymbol <- unlist(lapply(strsplit(hits_enrich@result$geneID, "/"),
                                                 function(x){x <- as.numeric(x)
                                                 g <- hits[[gene_column]][which(!is.na(match(hits$Genes_id, x)))]
                                                 g <- paste(g, collapse = "/");
                                                 g
                                                 }
  ))
  hits_enrich@result$Description <- stringr::str_remove_all(hits_enrich@result$Description, "%.{1,}")
  colnames(hits_enrich@result)[c(8,10)] <- c("geneNumber", "geneID")
  n_toshow <- length(which(hits_enrich@result$p.adjust <= pval_cutoff))

  fold <- hits[,c(gene_column, score_column)]
  fold <- as.data.frame(fold)
  if(any(duplicated(fold[[gene_column]]))){
    fold <- fold[-which(duplicated(fold[[gene_column]])),]
  }
  rownames(fold) <- fold[[gene_column]]
  fold[[gene_column]] <- NULL
  fold <- unlist(as.list(as.data.frame(t(fold))))


  graph <- enrichplot::cnetplot(hits_enrich, showCategory = n_toshow, foldChange = fold,
                    color_category = "#0059FE") +
    scale_color_gradient(low = "#87FE00", high = "#FE0000") +
    labs(title = paste("Gene-concept network", condition, "hits"),
         subtitle = paste("p-value cutoff:", pval_cutoff),
         color = score_column,
         size = "Count")

  return(graph)
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
