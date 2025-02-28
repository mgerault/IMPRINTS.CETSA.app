#' run_gsea
#'
#' Function to run gene set enrichment analysis and plot a GSEA plot.
#'
#' @param hits A data.frame containing the genes id, a score value and preferably a treatment column but not necessary.
#' @param gene_column The name of the column that contains the genes. Default is 'Gene'.
#' @param score_column The name of the column that contains the score values. Default is 'IS'.
#' @param treatment_column The name of the column that contains the treatments. Default is NULL.
#' @param treatment The name of the treatment you want to keep. Default is NULL.
#' @param species Specify the species. Currently, only 'human' and 'mouse' are available.
#' @param pos_enrichment Logical to tell if you want to only look at positive enrichment score.
#'                       If FALSE, show only negative.
#' @param pval_cutoff The p-value cutoff for the enrichment analysis.
#' @param minGSSize minimal size of each gene set for analyzing. default here is 3
#' @param database Specify the database. Currently, WikiPathway, KEGG, GO and CETSA are available.
#'
#' @return A list that contains the results and the plot.
#'
#' @export
#'
#' @seealso \code{\link{clusterProfiler}}

run_gsea <- function(hits, gene_column = "Gene", score_column = "IS",
                     treatment_column = NULL, treatment = NULL,
                     species = c("human", "mouse"), pos_enrichment = TRUE,
                     pval_cutoff = 0.01, minGSSize = 3,
                     database = c("WikiPathway", "KEGG", "GO", "CETSA")){
  require(clusterProfiler)
  species <- tolower(species)
  species <- match.arg(species)
  database <- match.arg(database)

  if(!("KEGGREST" %in% installed.packages())){
    message("Installing KEGGREST package")
    BiocManager::install("KEGGREST")
  }
  if(species == "human"){
    if(!("org.Hs.eg.db" %in% installed.packages())){
      message("Installing org.Hs.eg.db package")
      BiocManager::install("org.Hs.eg.db")
    }
  }
  else if(species == "mouse"){
    if(!("org.Mm.eg.db" %in% installed.packages())){
      message("Installing org.Mm.eg.db package")
      BiocManager::install("org.Mm.eg.db")
    }
  }

  org_data <- ifelse(species == "human", "org.Hs.eg.db", "org.Mm.eg.db")

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
    ### convert gene IDs
    hits_gene_id <- clusterProfiler::bitr(unique(sort(hits[[gene_column]])),
                                          fromType = "SYMBOL", toType = c("ENTREZID"),
                                          OrgDb = org_data, drop = FALSE)
    colnames(hits_gene_id) <- c(gene_column, "Gene_id")

    hits <- dplyr::left_join(hits, hits_gene_id, by = gene_column, multiple = "all")
    if(any(is.na(hits$Gene_id))){
      hits <- hits[which(!is.na(hits$Gene_id)),]
    }
  }
  else{
    hits$Gene_id <- hits[[gene_column]]
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
                                      scoreType = "pos", pvalueCutoff = pval_cutoff,
                                      minGSSize = minGSSize)
  }
  else if(database == "KEGG"){
    if(species == "human"){
      gsea_res <- clusterProfiler::gseKEGG(hits, organism = "hsa", pvalueCutoff = pval_cutoff,
                                           minGSSize = minGSSize)
    }
    else if(species == "mouse"){
      gsea_res <- clusterProfiler::gseKEGG(hits, organism = "mmu", pvalueCutoff = pval_cutoff,
                                           minGSSize = minGSSize)
    }
    rm(.KEGG_clusterProfiler_Env, envir=sys.frame()) # hidden object from clusterprofiler prevent dbplyr to load when in the environment
  }
  else if(database == "GO"){
    if(species == "human"){
      gsea_res <- clusterProfiler::gseGO(hits, OrgDb = "org.Hs.eg.db", pvalueCutoff = pval_cutoff,
                                         ont = "BP", minGSSize = minGSSize)
    }
    else if(species == "mouse"){
      gsea_res <- clusterProfiler::gseGO(hits, OrgDb = "org.Mm.eg.db", pvalueCutoff = pval_cutoff,
                                         ont = "BP", minGSSize = minGSSize)
    }
    rm(.GO_clusterProfiler_Env, .GOTERM_Env, envir=sys.frame()) # hidden object from clusterprofiler prevent dbplyr to load when in the environment
  }
  else if(database == "CETSA"){
    gsea_res <- clusterProfiler::GSEA(hits,
                                      TERM2GENE=cetsa_gsea_database[,c("cetsa.id", "gene")],
                                      TERM2NAME=cetsa_gsea_database[,c("cetsa.id", "name")],
                                      scoreType = "pos", pvalueCutoff = pval_cutoff,
                                      minGSSize = minGSSize)
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
    res <- gsea_res@result
    if(database != "CETSA"){
      res$geneSymbol <- unlist(lapply(strsplit(res$core_enrichment, "/"),
                                                  function(x){x <- as.numeric(x)
                                                  g <- hits_gene_id[[gene_column]][which(!is.na(match(hits_gene_id$Gene_id, x)))]
                                                  g <- paste(g, collapse = "/");
                                                  g
                                                  }
                                                  )
                                           )
    }
    else{
      extra_info <- unique(cetsa_gsea_database[,c("cetsa.id", "function", "functional.hypothesis")])
      colnames(extra_info)[1] <- "ID"
      res <- dplyr::left_join(res, extra_info, by = "ID")
    }

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
    graph[[1]] <- graph[[1]] +
      theme(legend.position = "inside",
            legend.position.inside = c(0.8,0.9),
            legend.text = element_text(face = "bold"))

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

  date <- strsplit(as.character(Sys.Date()), "-")[[1]]
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

  wikipath <- tryCatch({
    if(wp){
      clusterProfiler::read.gmt.wp(url_wiki)
    }
    else{
      clusterProfiler::read.gmt(url_wiki)
    }
  },
  # if wrong date
  error = function(e){
    date <- strsplit(date, "")[[1]][-(nchar(date)-1)]
    date <- paste0(date, collapse = "")
    date <- paste0(date, 1)
    url_wiki <- paste0("https://wikipathways-data.wmcloud.org/", date, "/gmt/wikipathways-",
                       date, "-gmt-", species, ".gmt")
    url_wiki <- url(url_wiki)

    wikipath2 <- tryCatch({
      if(wp){
        clusterProfiler::read.gmt.wp(url_wiki)
      }
      else{
        clusterProfiler::read.gmt(url_wiki)
      }
    },
    # if wrong date again, retrieve one month
    error = function(e){
      date <- strsplit(date, "")[[1]]
      date <- c(paste(date[1:4], collapse = ""), # year
                paste(date[5:6], collapse = ""), # month
                paste(date[7:8], collapse = "")  # day
      )
      date[3] <- "10"
      date[2] <- as.numeric(date[2]) - 1
      date[2] <- ifelse(date[2] == "0", "12", date[2])
      date[2] <- ifelse(nchar(date[2]) == 1, paste0("0", date[2]), date[2])
      date[1] <- ifelse(date[2] == "12", as.numeric(date[1]) - 1, date[1])
      date <- paste(date, collapse = "")

      url_wiki <- paste0("https://wikipathways-data.wmcloud.org/", date, "/gmt/wikipathways-",
                         date, "-gmt-", species, ".gmt")
      url_wiki <- url(url_wiki)

      wikipath3 <- tryCatch({
        if(wp){
          clusterProfiler::read.gmt.wp(url_wiki)
        }
        else{
          clusterProfiler::read.gmt(url_wiki)
        }
      },
      # if wrong date
      error = function(e) {
        date <- strsplit(date, "")[[1]][-(nchar(date)-1)]
        date <- paste0(date, collapse = "")
        date <- paste0(date, 1)
        url_wiki <- paste0("https://wikipathways-data.wmcloud.org/", date, "/gmt/wikipathways-",
                           date, "-gmt-", species, ".gmt")
        url_wiki <- url(url_wiki)

        if(wp){
          wikipath4 <- clusterProfiler::read.gmt.wp(url_wiki)
        }
        else{
          wikipath4 <-clusterProfiler::read.gmt(url_wiki)
        }

        return(wikipath4)
      })

      return(wikipath3)
    })

    return(wikipath2)
  })


  close(url_wiki)
  return(wikipath)
}

