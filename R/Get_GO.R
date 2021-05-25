#' Get_GO
#'
#' Function to get the result list from the get_annotations or get_enrichment functions from STRINGdb
#' package.
#'
#' @param anno The output data from  get_annotations or get_enrichment functions from STRINGdb package.
#' @param enrich A logical to tell wether your data is the out put from
#'               get_enrichment function from STRINGdb package or not
#' @param all_cat A logical to tell if you want the results from all the database or not
#' @param sing_cat A character; if all_cat is set to FALSE, will select only this database
#'
#' @return A list of the results
#'
#' @export
#'
#' @seealso \code{\link{STRINGdb}}

Get_GO <- function(anno, enrich = FALSE, all_cat = TRUE, sing_cat = "Component"){
  final <-list()
  if(all_cat){
    catego <- unique(anno$category)
  }
  else{
    catego <- sing_cat
  }

  for (l in catego){
    print(l)
    anno_l <- anno[which(anno$category == l),]

    gene_l <- list()
    check_g <- c()
    for (i in 1:nrow(anno_l)){
      prf_g <- str_which(colnames(anno_l), "preferred")
      lg <- anno_l[,prf_g][i]
      lg_b <- str_split(lg, ",")[[1]]

      stringid <- str_which(colnames(anno_l), "inputGenes|string_ids")
      lst <- anno_l[,stringid][i]
      lst_b <- str_split(lst, ",")[[1]]
      for(k in 1:length(lg_b)){
        if (!(lg_b[k] %in% check_g)){
          idx <- which(lapply(anno_l[,prf_g], function(x) str_detect(x, lg_b[k])) == TRUE)
          if(enrich){
            gene_l[[paste(lg_b[k], lst_b[k], sep = ",")]] <- data.frame(p_value = anno_l$p_value[idx],
                                                                        fdr = anno_l$fdr[idx],
                                                                        description = anno_l$description[idx])
          }
          else{
            gene_l[[paste(lg_b[k], lst_b[k], sep = ",")]] <- anno_l$description[idx]
          }
          check_g <- append(check_g, lg_b[k])
        }
      }
    }
    final[[l]] <- gene_l
  }
  return(final)
}
