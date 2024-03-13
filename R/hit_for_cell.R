#' hit_for_cell
#'
#' Function to get subcellular location based on the database from protein Atlas, and coordinates
#' of each protein on the plot cell (img2) according to their location
#'
#' @param HIT Usually, the output from hitlist function; but just has to be a data frame with the
#'             columns id, treatment and category
#' @param organism The organism on which you performed your experiment. For now, only 'HUMAN' and 'MOUSE' are available.
#'
#' @return Dataframe (format for hit_plotcell)
#'
#' @export
#'
#' @seealso \code{\link{hitplot_cell}}

hit_for_cell <- function(HIT, organism = c("HUMAN", "MOUSE")){
  if(organism == "HUMAN"){
    loca_hit <- IMPRINTS.CETSA.app::pr_atlas[which(!is.na(match(IMPRINTS.CETSA.app::pr_atlas$Uniprot, unique(HIT$id)))),]
  }
  else if(organism == "MOUSE"){
    loca_hit <- IMPRINTS.CETSA.app::pr_atlas_mouse[which(!is.na(match(IMPRINTS.CETSA.app::pr_atlas_mouse$Uniprot, unique(HIT$id)))),]
  }
  else{
    stop("Please, choose a valid organism name : 'HUMAN' or 'MOUSE'")
  }

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
    ORG <- strsplit(HIT$main.location[i], ", ")[[1]]
    corres <- c()
    for(k in 1:length(ORG)){
      corres <-  append(corres, orgatlas_match[ORG[k],2])
    }
    HIT$main.location.cell[i] <- paste(corres, collapse = ", ")
  }

  HIT$nb_location <- as.numeric(
    unlist(
      lapply(as.list(HIT$main.location.cell), function(x) {ifelse(x != "NA", length(strsplit(x, ", ")[[1]]), NA)})))

  nwl <- length(unique(HIT$id))
  HIT <- HIT[which(!is.na(HIT$nb_location)),] #remove prot without location
  nwl <- nwl - length(unique(HIT$id))
  message(paste(nwl, "proteins were removed because no subcellular location has been found"))

  if(nrow(HIT) == 0){
    message("No proteins have been found in the Protein Atlas database")
    return(HIT)
  }

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
          HIT[(i+k-1), "main.location"] <- strsplit(HIT$main.location[(i+k-1)], ", ")[[1]][k]
          HIT[(i+k-1), "main.location.cell"] <- strsplit(HIT$main.location.cell[(i+k-1)], ", ")[[1]][k]
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

  message("Your data are ready !")

  return(HIT)
}





