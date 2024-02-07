#' gene_concept_net
#'
#' Function to run over representation analysis on your hits and return
#' a network plot.
#'
#' @param hits A data.frame containing the genes id, a score value and preferably a treatment column but not necessary.
#' @param gene_column The name of the column that contains the genes. Default is 'Gene'.
#' @param score_column The name of the column that contains the score values. Default is 'IS'.
#' @param treatment_column The name of the column that contains the treatments. Default is NULL.
#' @param treatment The name of the treatment you ant to keep. Default is NULL.
#' @param species Specify the species. Currently, only 'human' and 'mouse' are available.
#' @param pval_cutoff The p-value cutoff for the gene concept network.
#' @param minGSSize minimal size of each gene set for analyzing. default here is 3
#' @param database Specify the database. Currently, WikiPathway, KEGG, GO and CETSA are available.
#'
#' @return The gene concept network plot.
#'
#' @export
#'
#' @seealso \code{\link{clusterProfiler}}

gene_concept_net <- function(hits, gene_column = "Gene", score_column = "IS",
                            treatment_column = NULL, treatment = NULL,
                            species = c("human", "mouse"),
                            pval_cutoff = 0.01, minGSSize = 3,
                            database = c("WikiPathway", "KEGG", "GO", "CETSA")){
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

  if(!is.null(treatment_column) & !is.null(treatment)){
    hits <- hits[which(hits[[treatment_column]] == treatment),]
  }
  if(any(is.na(hits[[score_column]]))){
    hits <- hits[which(!is.na(hits[[score_column]])),]
  }
  if(any(is.na(hits[[gene_column]]))){
    hits <- hits[which(!is.na(hits[[gene_column]])),]
  }

  if(database != "CETSA"){
    ### load database
    ensembl <- NULL
    while(is.null(ensembl)){
      ensembl <- tryCatch(biomaRt::useMart("ensembl", dataset = biomart_data, port = ""),
                          error = function(e) message("Timeout reached for getting gene ensemble,
                                                    fetching it again.")) # genes_id / gene symbols
    }

    # get genes id from gene symbol
    hits_gene_id <- biomaRt::getBM(attributes = c("uniprot_gn_symbol", "entrezgene_id"),
                                   filters = "uniprot_gn_symbol",
                                   values = unique(sort(hits[[gene_column]])),
                                   bmHeader = TRUE,
                                   mart = ensembl)
    colnames(hits_gene_id) <- c(gene_column, "Gene_id")

    hits <- dplyr::left_join(hits, hits_gene_id, by = gene_column, multiple = "all")
  }

  # start enrichment analysis
  if(database == "WikiPathway"){
    wp <- get_wikipath(FALSE, species = species) # wiki pathway
    hits_enrich <- clusterProfiler::enricher(hits$Gene_id,
                                             TERM2GENE = wp,
                                             pvalueCutoff = pval_cutoff,
                                             minGSSize = minGSSize)
  }
  else if(database == "KEGG"){
    if(species == "human"){
      hits_enrich <- clusterProfiler::enrichKEGG(hits$Gene_id,
                                                 organism = "hsa",
                                                 pvalueCutoff = pval_cutoff,
                                                 minGSSize = minGSSize)
    }
    else if(species == "mouse"){
      hits_enrich <- clusterProfiler::enrichKEGG(hits$Gene_id,
                                                 organism = "mmu",
                                                 pvalueCutoff = pval_cutoff,
                                                 minGSSize = minGSSize)
    }
    rm(.KEGG_clusterProfiler_Env, envir=sys.frame()) # hidden object from clusterprofiler prevent dbplyr to load when in the environment
  }
  else if(database == "GO"){
    if(species == "human"){
      if(!("org.Hs.eg.db" %in% installed.packages())){
        message("Installing org.Hs.eg.db package")
        BiocManager::install("org.Hs.eg.db")
      }
      hits_enrich <- clusterProfiler::enrichGO(hits$Gene_id,
                                               OrgDb = "org.Hs.eg.db",
                                               pvalueCutoff = pval_cutoff,
                                               minGSSize = minGSSize)
    }
    else if(species == "mouse"){
      if(!("org.Mm.eg.db" %in% installed.packages())){
        message("Installing org.Mm.eg.db package")
        BiocManager::install("org.Mm.eg.db")
      }
      hits_enrich <- clusterProfiler::enrichGO(hits$Gene_id,
                                               OrgDb = "org.Mm.eg.db",
                                               pvalueCutoff = pval_cutoff,
                                               minGSSize = minGSSize)
    }
    rm(.GO_clusterProfiler_Env, .GOTERM_Env, envir=sys.frame()) # hidden object from clusterprofiler prevent dbplyr to load when in the environment
  }
  else if(database == "CETSA"){
    print(hits)
    hits_enrich <- clusterProfiler::enricher(hits[[gene_column]],
                                             TERM2GENE = cetsa_gsea_database[,c("name", "gene")], # data.frame of 2 columns with term and corresponding gene
                                             pvalueCutoff = pval_cutoff,
                                             minGSSize = minGSSize)
    print(hits_enrich)
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

    return(graph)
  }
  else if(nrow(hits_enrich@result) == 0){
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

    return(graph)
  }
  else{
    n_toshow <- length(which(hits_enrich@result$p.adjust <= pval_cutoff))
    if(n_toshow == 0){
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

      return(graph)
    }

    if(database != "CETSA"){
      hits_enrich@result$geneSymbol <- unlist(lapply(strsplit(hits_enrich@result$geneID, "/"),
                                                     function(x){x <- as.numeric(x)
                                                     g <- hits[[gene_column]][which(!is.na(match(hits$Gene_id, x)))]
                                                     g <- paste(g, collapse = "/");
                                                     g
                                                     }
      ))
    }
    else{
      hits_enrich@result$geneSymbol <- hits_enrich@result$geneID
    }

    hits_enrich@result$Description <- gsub("%.{1,}", "", hits_enrich@result$Description)
    colnames(hits_enrich@result)[c(8,10)] <- c("geneNumber", "geneID")

    fold <- hits[,c(gene_column, score_column)]
    fold <- as.data.frame(fold)
    if(any(duplicated(fold[[gene_column]]))){
      fold <- fold[-which(duplicated(fold[[gene_column]])),]
    }
    rownames(fold) <- fold[[gene_column]]
    fold[[gene_column]] <- NULL
    fold <- unlist(as.list(as.data.frame(t(fold))))


    graph <- enrichplot::cnetplot(hits_enrich, showCategory = n_toshow,
                                  color_category = "#0059FE",
                                  foldChange = fold) +
      scale_color_gradient(low = "#87FE00", high = "#FE0000") +
      labs(title = paste("Gene-concept network", treatment, "hits"),
           subtitle = paste(database, "p-value cutoff:", pval_cutoff),
           color = score_column,
           size = "Count")

    return(graph)
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

  date <- strsplit(Sys.Date(), "-")[[1]]
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
  if(!RCurl::url.exists(url_wiki)){
    date <- strsplit(date, "")[[1]][-(nchar(date)-1)]
    date <- paste0(date, collapse = "")
    date <- paste0(date, 1)
    url_wiki <- paste0("https://wikipathways-data.wmcloud.org/", date, "/gmt/wikipathways-",
                       date, "-gmt-", species, ".gmt")
  }
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
