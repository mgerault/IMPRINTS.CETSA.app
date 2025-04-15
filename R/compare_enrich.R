#' compare_enrich
#'
#' Function to run over representation analysis on your hits and return
#' a plot that compare the pathways found between your treatments.
#'
#' @param hits The hitlist; a data.frame containing the genes id and preferably a treatment column but not necessary.
#' @param gene_column The name of the coulumn that contains the genes. Default is 'Gene'.
#' @param treatment_column The name of the column that contains the treatments. Default is NULL.
#' @param species Specify the species. Currently, only 'human' and 'mouse' are available.
#' @param n_pathway Number of pathway to show on plot. Default is 5.
#'                  For more info you, see \code{\link{compareCluster}}.
#' @param pval_cutoff The p-value cutoff for the enrichment analysis.
#' @param minGSSize minimal size of each gene set for analyzing. Default here is 3
#' @param database Specify the database. Currently, WikiPathway, KEGG, GO and CETSA are available.
#'
#' @return A list that contains the results and the plot.
#'
#' @export
#'
#' @seealso \code{\link{clusterProfiler}}

compare_enrich <- function(hits, gene_column = "Gene", treatment_column = NULL,
                       species = c("human", "mouse"), n_pathway = 5,
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


  if(is.null(treatment_column)){
    hits$treatment <- "treatment"
  }
  else{
    colnames(hits)[grep(treatment_column, colnames(hits))] <- "treatment"
  }

  # cluster compare enrichment analysis
  if(database == "WikiPathway"){
    wp <- get_wikipath(species = species) # wiki pathway
    hits_enrich <- clusterProfiler::compareCluster(Gene_id~treatment,
                                                   data = hits, fun = "enricher",
                                                   TERM2GENE=wp[,c("wpid", "gene")],
                                                   TERM2NAME=wp[,c("wpid", "name")],
                                                   pvalueCutoff = pval_cutoff,
                                                   minGSSize = minGSSize)
  }
  else if(database == "KEGG"){
    if(species == "human"){
      hits_enrich <- clusterProfiler::compareCluster(Gene_id~treatment,
                                                     data = hits, fun = "enrichKEGG",
                                                     organism = "hsa", pvalueCutoff = pval_cutoff,
                                                     minGSSize = minGSSize)
    }
    else if(species == "mouse"){
      hits_enrich <- clusterProfiler::compareCluster(Gene_id~treatment,
                                                     data = hits, fun = "enrichKEGG",
                                                     organism = "mmu", pvalueCutoff = pval_cutoff,
                                                     minGSSize = minGSSize)
    }
    rm(.KEGG_clusterProfiler_Env, envir=sys.frame()) # hidden object from clusterprofiler prevent dbplyr to load when in the environment
  }
  else if(database == "GO"){
    if(species == "human"){
      hits_enrich <- clusterProfiler::compareCluster(Gene_id~treatment,
                                                     data = hits, fun = "enrichGO", ont = "BP",
                                                     OrgDb = "org.Hs.eg.db", pvalueCutoff = pval_cutoff,
                                                     minGSSize = minGSSize)
    }
    else if(species == "mouse"){
      hits_enrich <- clusterProfiler::compareCluster(Gene_id~treatment,
                                                     data = hits, fun = "enrichGO",ont = "BP",
                                                     OrgDb = "org.Mm.eg.db", pvalueCutoff = pval_cutoff,
                                                     minGSSize = minGSSize)
    }
    rm(.GO_clusterProfiler_Env, .GOTERM_Env, envir=sys.frame()) # hidden object from clusterprofiler prevent dbplyr to load when in the environment
  }
  else if(database == "CETSA"){
    if(species != "human"){
      stop("Only human is available for CETSA database.")
    }
    hits$Gene_id <- hits[[gene_column]]
    hits_enrich <- clusterProfiler::compareCluster(Gene_id~treatment,
                                                   data = hits, fun = "enricher",
                                                   TERM2GENE=cetsa_gsea_database[,c("cetsa.id", "gene")],
                                                   TERM2NAME=cetsa_gsea_database[,c("cetsa.id", "name")],
                                                   pvalueCutoff = pval_cutoff,
                                                   minGSSize = minGSSize)
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
    if(any(res$Count < minGSSize)){
      res <- res[-which(res$Count < minGSSize),]
    }
    if(nrow(res)){
      res_data <- res
      if(any(!(unique(hits$treatment) %in% res$treatment))){
        res$Cluster <- as.character(res$Cluster)
        res$Description <- as.character(res$Description)
        for(i in unique(hits$treatment)[!(unique(hits$treatment) %in% res$treatment)]){
          res <- as.data.frame(rbind(res[1,], res))
          res[1,which(sapply(as.list(res[1,]), is.numeric))] <- NA
          res[1,which(!sapply(as.list(res[1,]), is.numeric))] <- ""
          res$Cluster[1] <- i
          res$treatment[1] <- i
        }
      }

      if(is.null(levels(hits$treatment))){
        res$treatment <- factor(res$treatment, levels = unique(res$treatment))
      }
      else{
        res$treatment <- factor(res$treatment, levels = levels(hits$treatment))
      }

      ord <- res %>% dplyr::select(Description, treatment) %>%
        dplyr::group_by(Description) %>%
        dplyr::summarise(n = length(unique(treatment)),
                         t = paste(as.numeric(treatment), collapse = ""),
                         ord = paste0(n, "_", t))
      res$Description <- factor(res$Description, levels = rev(ord$Description[order(ord$ord)]))
      if(any(res$Description == "")){
        res$Description[which(res$Description == "")] <- res$Description[nrow(res)]
      }

      graph <- ggplot(res, aes(treatment, Description,
                               fill = p.adjust,
                               size = GeneRatio)) +
        geom_point(color = "black", shape = 21, na.rm = TRUE) +
        scale_y_discrete(labels = enrichplot:::default_labeller(30)) +
        scale_fill_continuous(low = "#01DD05", high = "#B30000", name = "p.adjust",
                              guide = guide_colorbar(reverse = TRUE)) +
        scale_size(range = c(3,8)) +
        DOSE::theme_dose(12) +
        labs(subtitle = paste(database, " pvalueCutoff:", pval_cutoff)) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 30, hjust = 1, size = rel(1.6), face = "bold"),
              axis.text.y = element_text(size = rel(1.3)))

      if(database != "CETSA"){
        res_data$geneSymbol <- unlist(lapply(strsplit(res_data$geneID, "/"), function(x){
          if(length(x)){
            x <- as.numeric(x)
            g <- hits_gene_id[[gene_column]][which(!is.na(match(hits_gene_id$Gene_id, x)))]
            g <- paste(g, collapse = "/")
          }
          else{
            g <- ""
          };
          g})
          )
      }
      else{
        extra_info <- unique(cetsa_gsea_database[,c("cetsa.id", "function", "functional.hypothesis")])
        colnames(extra_info)[1] <- "ID"
        res_data <- dplyr::left_join(res_data, extra_info, by = "ID")
      }

      return(list("res" = res_data, "graph" = graph))
    }
    else{
      #no term enriched with enough genes
      graph <- ggplot(data.frame(x = c(0,1), y = c(0,1)), aes(x,y, label = "s")) +
        geom_text(x=0.5, y=0.5, label = paste("No term enriched \nwith more than", minGSSize, "genes"), size = 10) +
        cowplot::theme_cowplot() +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank())

      return(list("res" = NULL, "graph" = graph))
    }
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

