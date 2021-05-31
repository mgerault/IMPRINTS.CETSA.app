#' hit_for_cell
#'
#' Function to get subcellular location based on the database from protein Atlas, and coordinates
#' of each protein on the pplot cell (img2) according to their location
#'
#' @param data Usually, the output from hitlist function; but just has to be a data frame with the
#'             columns id, Condition and category
#'
#' @return Dataframe (format for hit_plotcell)
#'
#' @export
#'
#' @seealso \code{\link{hitplot_cell}}

hit_for_cell <- function(HIT){
  loca_hit <- pr_atlas[which(!is.na(match(pr_atlas$Uniprot, unique(HIT$id)))),]

  message("Assign proteins to their main location from protein Atlas")
  HIT[,c("gene.name", "main.location", "additional.location")] <- rep(NA, nrow(HIT))
  for(i in loca_hit$Uniprot){
    idx <- which(!is.na(match(HIT$id, i)))
    HIT[idx, c("gene.name", "main.location", "additional.location")] <-
      loca_hit[which(loca_hit$Uniprot == i), c("Gene", "Subcellular.main.location", "Subcellular.additional.location")]
  }
  HIT <- HIT[,!(colnames(HIT) %in% "additional.location")] #for now, not take additional locations

  HIT$main.location.cell <- rep(NA, nrow(HIT))
  for (i in 1:nrow(HIT)){
    ORG <- str_split(HIT$main.location[i], ", ")[[1]]
    corres <- c()
    for(k in 1:length(ORG)){
      corres <-  append(corres, orgatlas_match[ORG[k],2])
    }
    HIT$main.location.cell[i] <- paste(corres, collapse = ", ")
  }

  HIT$nb_location <- as.numeric(
    unlist(
      lapply(as.list(HIT$main.location.cell), function(x) {ifelse(x != "NA", length(str_split(x, ", ")[[1]]), NA)})))

  nwl <- nrow(HIT)
  HIT <- HIT[which(!is.na(HIT$nb_location)),] #remove prot without location
  nwl <- nwl - nrow(HIT)
  message(paste(nwl, "proteins were removed because no subcellular location has been found"))

  #duplicate row when more than 1 location
  message("Duplicate row when more than one location and separate them")
  #uncount : Based on a given column's value (x), repeat each row x times, thereby elongating the df.
  HIT <- HIT %>%
    tidyr::uncount(nb_location, .remove = FALSE)
  #then separate them
  mem <- c()
  for(i in 1:nrow(HIT)){
    if(HIT$nb_location[i] > 1){
      if(!(i %in% mem)){
        l <- HIT$nb_location[i]
        for(k in 1:l){
          HIT[(i+k-1), "main.location"] <- str_split(HIT$main.location[(i+k-1)], ", ")[[1]][k]
          HIT[(i+k-1), "main.location.cell"] <- str_split(HIT$main.location.cell[(i+k-1)], ", ")[[1]][k]
        }
        mem <- i:(i + l - 1)
      }
    }
  }

  message("Assign cooordinates to proteins according to their subcellular locations")
  idx_malo <- grep("main.location.cell", colnames(HIT))
  HIT[, c("x", "y")] <- t(apply(HIT, 1, function(x) point_in_orga(x[idx_malo])))

  HIT$main.location.cell <- as.factor(HIT$main.location.cell)
  HIT$nb_location <- as.factor(HIT$nb_location)
  HIT$category <- as.factor(HIT$category)

  HIT$txt <- paste("Uniprot ID :", HIT$id,
                    "<br> Gene name :", HIT$gene.name,
                    "<br> Category :", HIT$category,
                    "<br> Condition :", HIT$Condition,
                    "<br> Organelle :", HIT$main.location.cell,
                    "<br> Located in", HIT$nb_location, "organelles" )

  message("Your data are ready !")

  return(HIT)
}


