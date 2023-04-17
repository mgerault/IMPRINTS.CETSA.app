#' compare_enrich
#'
#' Function to run enrichment analysis on your hits and return
#' a plot that compare the pathways found between your treatments.
#'
#' @param hits The hitlist; a data.frame containing the genes id and preferably a treatment column but not necessary.
#' @param gene_column The name of the coulumn that contains the genes. Default is 'Genes'.
#' @param treatment_column The name of the column that contains the treatments. Default is NULL.
#' @param species Specify the species. Currently, only 'human' and 'mouse' are available.
#' @param n_pathway Number of pathway to show on plot. Default is 5.
#'                  For more info you, see \code{\link{compareCluster}}.
#' @param pval_cutoff The p-value cutoff for the enrichment analysis.
#' @param database Specify the database. Currently, WikiPathway, KEGG and GO ar available.
#'
#' @return A list that contains the results and the plot.
#'
#' @export
#'
#' @seealso \code{\link{clusterProfiler}}

compare_enrich <- function(hits, gene_column = "Genes", treatment_column = NULL,
                       species = c("human", "mouse"), n_pathway = 5,
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
                        error = function(e) message("Timeout reached for getting gene ensemble, fetching it again.")) # genes_id / gene symbols
  }

  if(any(is.na(hits[[gene_column]]))){
    hits <- hits[which(!is.na(hits[[gene_column]])),]
  }

  # get genes id from gene symbol
  hits_gene_id <- biomaRt::getBM(attributes = c("uniprot_gn_symbol", "entrezgene_id"),
                                 filters = "uniprot_gn_symbol",
                                 curl = curl::handle_setopt(curl::new_handle(), timeout = 30000),
                                 values = unique(sort(hits[[gene_column]])),
                                 bmHeader = TRUE,
                                 mart = ensembl)
  colnames(hits_gene_id) <- c(gene_column, "Genes_id")

  hits <- dplyr::left_join(hits, hits_gene_id, by = gene_column, multiple = "all")

  if(is.null(treatment_column)){
    hits$treatment <- "treatment"
  }
  else{
    colnames(hits)[stringr::str_which(colnames(hits), treatment_column)] <- "treatment"
  }

  # cluster compare enrichment analysis
  if(database == "WikiPathway"){
    wp <- get_wikipath(species = species) # wiki pathway
    hits_enrich <- clusterProfiler::compareCluster(Genes_id~treatment,
                                                   data = hits, fun = "enricher",
                                                   TERM2GENE=wp[,c("wpid", "gene")],
                                                   TERM2NAME=wp[,c("wpid", "name")],
                                                   pvalueCutoff = pval_cutoff)
  }
  else if(database == "KEGG"){
    if(species == "human"){
      hits_enrich <- clusterProfiler::compareCluster(Genes_id~treatment,
                                                     data = hits, fun = "enrichKEGG",
                                                     organism = "hsa", pvalueCutoff = pval_cutoff)
    }
    else if(species == "mouse"){
      hits_enrich <- clusterProfiler::compareCluster(Genes_id~treatment,
                                                     data = hits, fun = "enrichKEGG",
                                                     organism = "mmu", pvalueCutoff = pval_cutoff)
    }
    rm(.KEGG_clusterProfiler_Env, envir=sys.frame()) # hidden object from clusterprofiler prevent dbplyr to load when in the environment
  }
  else if(database == "GO"){
    if(species == "human"){
      if(!("org.Hs.eg.db" %in% installed.packages())){
        message("Installing org.Hs.eg.db package")
        BiocManager::install("org.Hs.eg.db")
      }
      hits_enrich <- clusterProfiler::compareCluster(Genes_id~treatment,
                                                     data = hits, fun = "enrichGO",
                                                     OrgDb = "org.Hs.eg.db", pvalueCutoff = pval_cutoff)
    }
    else if(species == "mouse"){
      if(!("org.Mm.eg.db" %in% installed.packages())){
        message("Installing org.Mm.eg.db package")
        BiocManager::install("org.Mm.eg.db")
      }
      hits_enrich <- clusterProfiler::compareCluster(Genes_id~treatment,
                                                     data = hits, fun = "enrichGO",
                                                     OrgDb = "org.Mm.eg.db", pvalueCutoff = pval_cutoff)
    }
    rm(.GO_clusterProfiler_Env, .GOTERM_Env, envir=sys.frame()) # hidden object from clusterprofiler prevent dbplyr to load when in the environment
  }

  if(is.null(hits_enrich)){
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
    res <- fortify(hits_enrich, showCategory = 6)
    res$treatment <- factor(res$treatment, levels = unique(res$treatment))
    graph <- ggplot(res, aes(treatment, Description,
                             color = p.adjust,
                             size = GeneRatio)) +
      geom_point() +
      scale_y_discrete(labels = enrichplot:::default_labeller(30)) +
      scale_color_continuous(low = "#B2ECBF", high = "#3821A3", name = "p.adjust",
                             guide = guide_colorbar(reverse = TRUE)) +
      scale_size(range = c(3,8)) +
      DOSE::theme_dose(12) +
      labs(subtitle = paste(database, " pvalueCutoff:", pval_cutoff)) +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.85)),
            axis.text.y = element_text(angle = 30, size = rel(0.85)))

    res$geneSymbol <- unlist(lapply(strsplit(res$geneID, "/"), function(x){x <- as.numeric(x)
    g <- hits_gene_id$Genes[which(!is.na(match(hits_gene_id$Genes_id, x)))]
    g <- paste(g, collapse = "/");
    g
    }
    ))

    return(list("res" = res, "graph" = graph))
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




