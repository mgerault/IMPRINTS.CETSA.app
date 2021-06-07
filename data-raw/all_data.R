library(mineCETSA)
library(stringr)

#the PI3K file
PI3K1h6h_file <- ms_fileread("210503_1729_PI3K1h6h.txt")
names(PI3K1h6h_file) <- str_replace(names(PI3K1h6h_file), ".x", "1h") #change condition name to fit the treatmentlevel
names(PI3K1h6h_file) <- str_replace(names(PI3K1h6h_file), ".y", "6h")

hitlist_PI3K1h6h <- read.csv("./Data/PI3K/PI3K1h6h_ForHitGeneration_1722_30-03-21_Summary.csv",
                             row.names = 1)
NN_PI3K1h6h <- read.csv("./Data/PI3K/PI3K1h6h_ForHitGeneration_1722_30-03-21_NN.csv",
                        row.names = 1)

#the TNF file
TNF_MOLM1316 <- ms_fileread("210519_1416_TNF_MOLM1316.txt")

hitlist_TNF16 <- read.csv("./Data/TNF/210409_1222_MOLM16_s1_1327_9-04-21_Summary.csv",
                          row.names = 1)
NN_TNF16 <- read.csv("./Data/TNF/210409_1222_MOLM16_s1_1327_9-04-21_NN.csv",
                     row.names = 1)
hitlist_TNF13 <- read.csv("./Data/TNF/210409_1223_MOLM13_s1_1328_9-04-21_Summary.csv",
                          row.names = 1)
NN_TNF13 <- read.csv("./Data/TNF/210409_1223_MOLM13_s1_1328_9-04-21_NN.csv",
                     row.names = 1)
hitlist_TNF <- rbind(hitlist_TNF13, hitlist_TNF16)
NN_TNF <- rbind(NN_TNF13, NN_TNF16)



#the list named drug_data
drug_data <- list("data" = list("PI3K" = PI3K1h6h_file, "TNF" = TNF_MOLM1316),
                  "treat_level" = list("PI3K" = get_treat_level(PI3K1h6h_file), "TNF" = get_treat_level(TNF_MOLM1316)),
                  "hitlist" = list("PI3K" = hitlist_PI3K1h6h, "TNF" = hitlist_TNF),
                  "NN" = list("PI3K" = NN_PI3K1h6h, "TNF" = NN_TNF))

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
library(EBImage)
img2 <- readImage("./data-raw/Animal_Cell.png")

img2@.Data <- aperm(img2@.Data, c(2,1,3))




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

usethis::use_data(PI3K1h6h_file, overwrite = TRUE)
usethis::use_data(hitlist_PI3K1h6h, overwrite = TRUE)
usethis::use_data(NN_PI3K1h6h, overwrite = TRUE)
usethis::use_data(TNF_MOLM1316, overwrite = TRUE)
usethis::use_data(hitlist_TNF, overwrite = TRUE)
usethis::use_data(NN_TNF, overwrite = TRUE)
usethis::use_data(drug_data, overwrite = TRUE)
usethis::use_data(ev_null_print, overwrite = TRUE)
usethis::use_data(img2, overwrite = TRUE)
usethis::use_data(pr_atlas, overwrite = TRUE)
usethis::use_data(pr_atlas_mouse, overwrite = TRUE)
usethis::use_data(orgatlas_match, overwrite = TRUE)
usethis::use_data(loca_orga, overwrite = TRUE)
usethis::use_data(rg_list, overwrite = TRUE)


