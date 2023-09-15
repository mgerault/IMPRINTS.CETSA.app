library(IMPRINTS.CETSA.app)
library(IMPRINTS.CETSA)
library(stringr)

#add data from cell cyle paper from Dai Lingyun 2018
#elutriation file
elutriation <- ms_fileread("210827_1735_elutriation_imprints_caldiff.txt")
hitlist_elutriation <- openxlsx::read.xlsx("elutriation_Summary.xlsx")
colnames(hitlist_elutriation)[2] <- "treatment"
NN_elutriation <- openxlsx::read.xlsx("elutriation_NN.xlsx")
colnames(NN_elutriation)[3] <- "treatment"

#the chemarrest file
chemarrest <- ms_fileread("210827_1739_chemarrest_imprints_caldiff.txt")
hitlist_chemarrest <- openxlsx::read.xlsx("chemarrest_Summary.xlsx")
colnames(hitlist_chemarrest)[2] <- "treatment"
NN_chemarrest <- openxlsx::read.xlsx("chemarrest_NN.xlsx")
colnames(NN_chemarrest)[3] <- "treatment"

elutriation_ave <- imprints_average_sh(elutriation, FALSE)
chemarrest_ave <- imprints_average_sh(chemarrest, FALSE)

#the list named drug_data
drug_data <- list("data" = list("elutriation" = elutriation, "chemarrest" = chemarrest),
                  "data_ave" = list("elutriation" = elutriation_ave, "chemarrest" = chemarrest_ave),
                  "treat_level" = list("elutriation" = get_treat_level(elutriation), "chemarrest" = get_treat_level(chemarrest)),
                  "hitlist" = list("elutriation" = hitlist_elutriation, "chemarrest" = hitlist_chemarrest),
                  "NN" = list("elutriation" = NN_elutriation, "chemarrest" = NN_chemarrest))

#a graph to print on the tab interactive cell of the app
library(ggplot2)
ev_null_print <-
  ggplot(data.frame(x = c(0,1), y = c(0,1)), aes(x,y, label = "s")) +
  geom_text(x=0.5, y=0.5, label = "Click on a protein \nto see its bar plot !", size = 10) +
  theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())


#the interactive cell
img2 <- png::readPNG("./data-raw/Animal_Cell.png")


###get location
### HUMAN
#https://www.proteinatlas.org/api/search_download.php?search=Human&format=json&columns=g,up,relih,scl,scml,scal&compress=no
pr_atlas <- rjson::fromJSON(file = "./data-raw/Human.json")
pr_atlas <- lapply(pr_atlas, function(y) lapply(y, function(x) {
  if(length(x) != 0){
    if(length(x) > 1){
      x <- paste(x, collapse = ", ")
    }
  }
  else{
    x <- NA
  }; x})
)
pr_atlas <- lapply(pr_atlas, function(x) {x <- as.data.frame(x); x})
pr_atlas <- do.call(rbind, pr_atlas)


### Mouse
#https://www.proteinatlas.org/api/search_download.php?search=Mouse&format=json&columns=g,up,relih,scl,scml,scal&compress=no
pr_atlas_mouse <- rjson::fromJSON(file = "./data-raw/Mouse.json")
pr_atlas_mouse <- lapply(pr_atlas_mouse, function(y) lapply(y, function(x) {
  if(length(x) != 0){
    if(length(x) > 1){
      x <- paste(x, collapse = ", ")
    }
  }
  else{
    x <- NA
  }; x})
)
pr_atlas_mouse <- lapply(pr_atlas_mouse, function(x) {x <- as.data.frame(x); x})
pr_atlas_mouse <- do.call(rbind, pr_atlas_mouse)

### no subcellular location that we don't already have with Human
#setdiff(unique(unlist(str_split(unique(pr_atlas_mouse$Subcellular.main.location), ", "))),
 #       orgatlas_match$location.from.pratlas)


#get correspondances between organellar names
orgatlas_match <- openxlsx::read.xlsx("./data-raw/organellar list.xlsx")
orgatlas_match$location.corres <- mapply(function(x,y){
  if(purrr::is_empty(x)){
    z <- y
  }
  else{
    z <- x
  }
  ;z
},
stringr::str_extract_all(orgatlas_match$location.from.pratlas,  "(?<=\\().+?(?=\\))"),
as.list(orgatlas_match$location.from.pratlas))

orgatlas_match$location.from.pratlas <- stringr::str_replace_all(orgatlas_match$location.from.pratlas,  "\\(|(?<=\\().+?(?=\\))|\\)", "")
orgatlas_match <- orgatlas_match[,-c(1:2)]
rownames(orgatlas_match) <- orgatlas_match$location.from.pratlas


loca_orga

rg_list <- list()
for(i in levels(loca_orga$organelle)){
  idx <- which(loca_orga$organelle == i)
  rg_x <- range(loca_orga$x[idx])
  rg_y <- range(loca_orga$y[idx])
  rg_list[[i]] <- data.frame(x = rg_x, y = rg_y)
}
rg_list


### adding CETSA cluster database to run gene set enrichment analysis base on this database
cetsa_gsea_database <- openxlsx::read.xlsx("./data-raw/Inital_CETSA_clusters.xlsx")
cetsa_gsea_database$cetsa.id <- paste0("CETSA", 1:nrow(cetsa_gsea_database))
cetsa_gsea_database <- cetsa_gsea_database %>%
  group_by(Name.of.cluster, cetsa.id) %>%
  group_modify(~ {
    p <- .x$Proteins.in.cluster
    p <- strsplit(p, ", | |,")[[1]]
    p <- gsub(" ", "", p)
    p <- p[nzchar(p)]

    p <- data.frame(Function = .x$Function,
                    Proteins.in.cluster = p,
                    Functional.hypothesis.for.shifts = .x$Functional.hypothesis.for.shifts)

    return(p)
  })

colnames(cetsa_gsea_database) <- c("name", "cetsa.id", "function", "gene", "functional.hypothesis")
cetsa_gsea_database$species <- "Homo sapiens"


### save data created
usethis::use_data(elutriation, overwrite = TRUE)
usethis::use_data(elutriation_ave, overwrite = TRUE)
usethis::use_data(hitlist_elutriation, overwrite = TRUE)
usethis::use_data(NN_elutriation, overwrite = TRUE)
usethis::use_data(chemarrest, overwrite = TRUE)
usethis::use_data(chemarrest_ave, overwrite = TRUE)
usethis::use_data(hitlist_chemarrest, overwrite = TRUE)
usethis::use_data(NN_chemarrest, overwrite = TRUE)
usethis::use_data(drug_data, overwrite = TRUE)
usethis::use_data(ev_null_print, overwrite = TRUE)
usethis::use_data(img2, overwrite = TRUE)
usethis::use_data(pr_atlas, overwrite = TRUE)
usethis::use_data(pr_atlas_mouse, overwrite = TRUE)
usethis::use_data(orgatlas_match, overwrite = TRUE)
usethis::use_data(loca_orga, overwrite = TRUE)
usethis::use_data(rg_list, overwrite = TRUE)
usethis::use_data(cetsa_gsea_database, overwrite = TRUE)

