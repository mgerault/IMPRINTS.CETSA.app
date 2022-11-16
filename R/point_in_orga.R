#' point_in_orga
#'
#' Function to generate random coordinate according to a subcellular location from img2.
#'
#' @param Organelle One of the valid organelle; To see the valid organelles type
#'                  \code{unique(mineCETSAapp::orgatlas_match$location.corres)}
#'
#' @return A double containing coordinates
#'
#' @export
#'
#' @seealso \code{\link{is_in_zone}}

point_in_orga <- function(Organelle = "Golgi"){
  if (str_detect(Organelle, "Mitochondrian|Vacuole|Vesicle")){
    if(sample(c(1,2), 1) == 1){
      Organelle <- paste0(Organelle, 2)
    }
  }
  if (str_detect(Organelle, "^endoplasmic_reticulum")){
    if(sample(c(1,2), 1) == 1){
      Organelle <- paste0("Rough_", Organelle)
    }
    else{
      Organelle <- paste0("Smooth_", Organelle)
    }
  }
  rg_orgx <- rg_list[[Organelle]]$x
  rg_orgy <- rg_list[[Organelle]]$y

  rnd_x <- runif(1, rg_orgx[1], rg_orgx[2])
  rnd_y <- runif(1, rg_orgy[1], rg_orgy[2])

  data <- loca_orga[which(loca_orga$organelle == Organelle), 1:2]
  while(!is_in_zone(data, c(rnd_x, rnd_y))){
    rnd_x <- runif(1, rg_orgx[1], rg_orgx[2])
    rnd_y <- runif(1, rg_orgy[1], rg_orgy[2])
  }

  return(c(rnd_x, rnd_y))
}
