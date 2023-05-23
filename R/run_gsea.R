#' run_gsea
#'
#' Function to run enrichment analysis and plot a GSEA plot.
#'
#' @param hits A data.frame containing the genes id, a score value and preferably a treatment column but not necessary.
#' @param gene_column The name of the coulumn that contains the genes. Default is 'Gene'.
#' @param score_column The name of the coulumn that contains the score values. Default is 'SR'.
#' @param treatment_column The name of the column that contains the treatments. Default is NULL.
#' @param treatment The name of the treatment you want to keep. Default is NULL.
#' @param species Specify the species. Currently, only 'human' and 'mouse' are available.
#' @param pos_enrichment Logical to tell if you want to only look at positive enrichment score.
#'                       If FALSE, show only negative.
#' @param pval_cutoff The p-value cutoff for the enrichment analysis.
#' @param database Specify the database. Currently, WikiPathway, KEGG and GO ar available.
#'
#' @return A list that contains the results and the plot.
#'
#' @export
#'
#' @seealso \code{\link{clusterProfiler}}

run_gsea <- function(hits, gene_column = "Gene", score_column = "SR",
                     treatment_column = NULL, treatment = NULL,
                     species = c("human", "mouse"), pos_enrichment = TRUE,
                     pval_cutoff = 0.01,
                     database = c("WikiPathway", "KEGG", "GO")){
  require(clusterProfiler)
  if(!("KEGGREST" %in% installed.packages())){
    message("Installing KEGGREST package")
    BiocManager::install("KEGGREST")
  }
  if(!("biomaRt" %in% installed.packages())){
    message("Installing biomaRt package")
    BiocManager::install("biomaRt")
  }

  species <- tolower(species)
  species <- match.arg(species)

  database <- match.arg(database)

  biomart_data <- ifelse(species == "human", "hsapiens_gene_ensembl", "mmusculus_gene_ensembl")

  ### load database
  ensembl <- NULL
  while(is.null(ensembl)){
    ensembl <- tryCatch(biomaRt::useMart("ensembl", dataset = biomart_data, port = ""),
                        error = function(e) message("Timeout reached for getting gene ensemble,
                                                    fetching it again.")) # genes_id / gene symbols
  }

  if(!is.null(treatment_column) & !is.null(treatment)){
    hits <- hits[which(hits[[treatment_column]] == treatment),]
  }
  if(any(is.na(hits[[score_column]]))){
    hits <- hits[which(!is.na(hits[[score_column]])),]
  }
  if(any(is.na(hits[[gene_column]]))){
    hits <- hits[which(!is.na(hits[[gene_column]])),]
  }

  # get genes id from gene symbol
  hits_gene_id <- biomaRt::getBM(attributes = c("uniprot_gn_symbol", "entrezgene_id"),
                                 filters = "uniprot_gn_symbol",
                                 values = unique(sort(hits[[gene_column]])),
                                 bmHeader = TRUE,
                                 mart = ensembl)
  colnames(hits_gene_id) <- c(gene_column, "Gene_id")

  hits <- dplyr::left_join(hits, hits_gene_id, by = gene_column, multiple = "all")
  if(any(is.na(hits$Gene_id))){
    hits <- hits[which(!is.na(hits$Gene_id)),]
  }
  hits$score <- hits[[score_column]]
  hits <- hits[, c("Gene_id", "score")]
  if(any(duplicated(hits$Gene_id))){
    hits <- hits[-which(duplicated(hits$Gene_id)),]
  }
  hits <- hits[order(hits$score, decreasing = TRUE),]
  hits <- as.data.frame(hits)
  rownames(hits) <- hits$Gene_id
  hits$Gene_id <- NULL

  hits <- unlist(as.list(as.data.frame(t(hits))))

  if(database == "WikiPathway"){
    wp <- get_wikipath(species = species) # wiki pathway
    gsea_res <- clusterProfiler::GSEA(hits,
                                      TERM2GENE=wp[,c("wpid", "gene")],
                                      TERM2NAME=wp[,c("wpid", "name")],
                                      scoreType = "pos", pvalueCutoff = pval_cutoff)
  }
  else if(database == "KEGG"){
    if(species == "human"){
      gsea_res <- clusterProfiler::gseKEGG(hits, organism = "hsa", pvalueCutoff = pval_cutoff)
    }
    else if(species == "mouse"){
      gsea_res <- clusterProfiler::gseKEGG(hits, organism = "mmu", pvalueCutoff = pval_cutoff)
    }
    rm(.KEGG_clusterProfiler_Env, envir=sys.frame()) # hidden object from clusterprofiler prevent dbplyr to load when in the environment
  }
  else if(database == "GO"){
    if(species == "human"){
      if(!("org.Hs.eg.db" %in% installed.packages())){
        message("Installing org.Hs.eg.db package")
        BiocManager::install("org.Hs.eg.db")
      }
      gsea_res <- clusterProfiler::gseGO(hits, OrgDb = "org.Hs.eg.db", pvalueCutoff = pval_cutoff)
    }
    else if(species == "mouse"){
      if(!("org.Mm.eg.db" %in% installed.packages())){
        message("Installing org.Mm.eg.db package")
        BiocManager::install("org.Mm.eg.db")
      }
      gsea_res <- clusterProfiler::gseGO(hits, OrgDb = "org.Mm.eg.db", pvalueCutoff = pval_cutoff)
    }
    rm(.GO_clusterProfiler_Env, .GOTERM_Env, envir=sys.frame()) # hidden object from clusterprofiler prevent dbplyr to load when in the environment
  }

  if(nrow(gsea_res@result) == 0){
    #no term enriched under specific pvalueCutoff...
    graph <- ggplot(data.frame(x = c(0,1), y = c(0,1)), aes(x,y, label = "s")) +
      geom_text(x=0.5, y=0.5, label = paste("No term enriched \nunder p-value of", pval_cutoff), size = 10) +
      cowplot::theme_cowplot() +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank())

    return(list("res" = NULL, "graph" = graph))
  }
  else{
    gsea_res@result$geneSymbol <- unlist(lapply(strsplit(gsea_res@result$core_enrichment, "/"),
                                                function(x){x <- as.numeric(x)
                                                g <- hits_gene_id$Gene[which(!is.na(match(hits_gene_id$Gene_id, x)))]
                                                g <- paste(g, collapse = "/");
                                                g
                                                }
                                                )
                                         )

    if(pos_enrichment){
      if(length(which(gsea_res@result$enrichmentScore > 0))){
        graph <- enrichplot::gseaplot2(gsea_res,
                                       geneSetID = which(gsea_res@result$enrichmentScore > 0))
      }
      else{
        message("All enrichment score are negative, showing those instead.")
        graph <- enrichplot::gseaplot2(gsea_res,
                                       geneSetID = which(gsea_res@result$enrichmentScore < 0))
      }
    }
    else{
      if(length(which(gsea_res@result$enrichmentScore < 0))){
        graph <- enrichplot::gseaplot2(gsea_res,
                                       geneSetID = which(gsea_res@result$enrichmentScore < 0))
      }
      else{
        message("All enrichment score are positive, showing those instead.")
        graph <- enrichplot::gseaplot2(gsea_res,
                                       geneSetID = which(gsea_res@result$enrichmentScore > 0))
      }
    }
    graph[[1]]$labels$subtitle <- paste(database, " pvalueCutoff:", pval_cutoff)

    return(list("res" = gsea_res@result, "graph" = graph))
  }
}


# function to always get most recent wiki pathway database (update every 10 of month)
get_wikipath <- function(wp = TRUE, species = "human"){
  if(species == "human"){
    species <- "Homo_sapiens"
  }
  else if(species == "mouse"){
    species <- "Mus_musculus"
  }

  date <- stringr::str_split(Sys.Date(), "-")[[1]]
  date <- as.numeric(date)
  if(date[3] < 11){
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

  url_wiki <- paste0("https://wikipathways-data.wmcloud.org/", date, "/gmt/wikipathways-",
                     date, "-gmt-", species, ".gmt")
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
