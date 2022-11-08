## ---- eval=FALSE--------------------------------------------------------------
#  if(!requireNamespace("BiocManager", quietly = TRUE)){ #check if you already have the BiocManager package
#   install.packages("BiocManager")  #if not, install it
#  }
#  BiocManager::install("STRINGdb")
#  BiocManager::install("EBImage")

## ---- eval=FALSE--------------------------------------------------------------
#  if(!requireNamespace("devtools", quietly = TRUE)){ #check if you already have the devtools package
#   install.packages("devtools")  #if not, install it
#  }
#  devtools::install_github("mgerault/mineCETSAapp")

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  library("mineCETSAapp")

## ---- eval=FALSE--------------------------------------------------------------
#  runCETSAapp()      #this function will directly start the app

## ---- eval=FALSE--------------------------------------------------------------
#  mydata_caldiff <- drug_data$data$elutriation  #save your data
#  mydata_ave <- drug_data$data_ave$elutriation
#  
#  View(mydata_caldiff)  #take a look at the data
#  View(mydata_ave)

## ---- eval=FALSE--------------------------------------------------------------
#  IMPRINTS_barplotting_sh(mydata_caldiff, ret_plot = FALSE, save_pdf = TRUE) #will save all the barplots from the data
#  
#  mydata_caldiff_chemarrest <- drug_data$data$chemarrest
#  all_drug <- list("elutriation" = mydata_caldiff, "chemarrest" = mydata_caldiff_chemarrest)
#  IMPRINTS_barplotting_sh(all_drug, ret_plot = FALSE, save_pdf = TRUE) #will save all the barplots from the data
#  #here, it will "join" the two pdf in one

## ---- eval=FALSE--------------------------------------------------------------
#  mydata_hit <- drug_data$hitlist$elutriation
#  #will map the proteins categorized as CN, CC and NC in S condition to some protein complex (core Corum database)
#  map_compl <- IMPRINTS_complex_mapping_sh(mydata_ave, mydata_hit, treatment = "S",
#                                        targetcategory = c("CN", "CC", "NC"))
#  View(map_compl) #if you want to check the data
#  
#  map_compl <- map_compl[, c("ComplexName", "subunitsNum", "subunitsIdentifiedNum",  #keep columns of interest
#                                   "id", "description", "gene", "category")]
#  
#  if(nrow(map_compl) !=0){
#    map_compl$description <- mineCETSAapp:::getProteinName(map_compl$description)  #keep only protein names in description
#  }

## ---- eval=FALSE--------------------------------------------------------------
#  data_l <- list() #initialize list
#  some_complex <- unique(map_compl$ComplexName)
#  some_complex <- some_complex[sample(1:length(some_complex), 4)]
#  for(i in some_complex){
#    cate_ <- map_compl[which(map_compl$ComplexName == i), ]
#    pr_comp <- cate_$id
#  
#    data_l[[i]] <- ms_subsetting(data, isfile = F, hitidlist = c(pr_comp)) #filter the proteins from those complexes
#  
#    data_l[[i]]$category <- cate_$category[which(!is.na(match(cate_$id, data_l[[i]]$id)))] #keep category
#  }
#  data <- data_l
#  
#  #each element of the list contains data from caldiff output for each protein complex selected
#  
#  IMPRINTS_barplotting_sh(data_l,save_pdf = TRUE, ret_plot = FALSE,
#                       toplabel = "IMPRINTS-CETSA bar plotting \nProtein complex :",  #title on each page (will add the complex name, cf. the names of the list)
#                       pdfname = "complex_barplot"
#                       )

## ---- eval=FALSE--------------------------------------------------------------
#  IMPRINTS_barplotting_simprof(mydata_caldiff, mydata_ave, #here we have the average data, but if it wasn't the case we would have let it NULL and let the function calculate it for us
#                            treatmentlevel = "S", protein_profile = "O43776",
#                            pdfname = "similar_O43776",
#                            use_score = "euclidean",
#                            score_threshold = 0.65, max_na_prow = 0
#                            )
#  #Here we choose the Euclidean distance score with a threshold of 0.65 and no missing values

## ---- eval=FALSE--------------------------------------------------------------
#  # According their category
#  #here we need the average data
#  IMPRINTS_heatmap(mydata_ave, hit_summary = mydata_hit,
#                treatment = "S", max_na = 0,
#                response = "both",
#                select_cat = c("CC", "CN", "NC"),
#                saveHeat = TRUE, file_type = "png")
#  #Will save a heatmap in a png file, of the proteins categorized as CN, CC or NC under the condition S
#  
#  # According to their protein complex
#  some_complex <- unique(map_compl$ComplexName)
#  some_complex <- some_complex[sample(1:length(some_complex), 4)]
#  cate_ <- map_compl[which(map_compl$ComplexName == i), ]
#  pr_comp <- cate_$id
#  
#  data <- ms_subsetting(mydata_ave, isfile = F, hitidlist = c(pr_comp)) #filter the proteins from those complexes
#  
#  IMPRINTS_heatmap(data, PRcomplex_data = cate_,
#                treatment = "S", max_na = 0,
#                response = "both",
#                saveHeat = TRUE, file_type = "png")
#  #Will save a heatmap in a png file, of the proteins of the complex from some_complex, under the condition S

## ---- eval=FALSE--------------------------------------------------------------
#  new_hits <- hitlist_outliers(mydata_caldiff, control = "G1", #the control
#                               basetemp = "37C",        #categorization based on the lowest tempearature
#                               format = "xlsx") #will save in xlsx format
#  #will ask you if you want to save or not the results

## ---- eval=FALSE--------------------------------------------------------------
#  library(STRINGdb)
#  string_db <- STRINGdb$new(version="11", species=9606,               #ID 9606 correspond to human
#                             score_threshold=200,
#                             input_directory=  file.path(getwd(), "STRING_data")) #will save the data in a folder named STRING_data

## ---- eval=FALSE--------------------------------------------------------------
#  data <- mydata_hit %>% dplyr::filter(Condition == "S")
#  network_hit <- string_db$map(data, "id", removeUnmappedRows = TRUE)         #will map the proteins to the string ID
#  
#  My_net(network_hit$STRING_id , inter = FALSE) #will plot the network

## ---- eval=FALSE--------------------------------------------------------------
#  enrich_hit <- string_db$get_enrichment(network_hit$STRING_id)
#  View(enrich_hit) #take a look at the data
#  #Each protein has been mapped to a category (Component, function, ...) with a description of this category.
#  #For example, with the category component we can have the description ribosome
#  
#  #sum up the result, will take only the category "Component"
#  df <- Get_GO(enrich_hit, enrich = TRUE, all_cat = FALSE, sing_cat = "Component") #a list of data frame
#  d_n <- as.list(lapply(df, names)[["Component"]]) #get names of the list : the string ID and the gene name
#  #as.list() to keep the class list and use it in mapply after
#  
#  #add the column id in each data frame (the names we got at the last step)
#  df <- mapply(function(x,y) {x$id <- rep(y, nrow(x)); x},
#                       df[["Component"]],
#                       d_n, SIMPLIFY = FALSE)
#  df <- do.call(rbind, df)                    #put all the data frame in one
#  rownames(df) <- 1:nrow(df)
#  View(df) #take a look at the data
#  
#  #keep the id and separate the gene name from the string ID
#  id_string <- do.call(rbind, str_split(df$id, ","))
#  colnames(id_string) <- c("gene.names", "STRING_id")
#  
#  #bind the two data frames
#  df <- cbind(id_string, df[,-ncol(df)])
#  View(df) #take a look at the final data
#  
#  #Let's see the network from some proteins mapped to a specific description
#  #for example the ribosome and the cytosolic ribosome
#  ribosome_net <- df$STRING_id[str_which(df$description, paste0("^", c("ribosome", "cytosolic ribosome"), "$"))]
#  My_net(ribosome_net , inter = FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  # Let's test it on the S condition
#  data <- mydata_hit %>% dplyr::filter(Condition == "S")
#  data <- data %>% dplyr::filter(!is.na(match(category, c("CC", "CN", "NC")))) #filter the category, but you could have take them all, even add also the NN
#  
#  cell_result <- hit_for_cell(data, organism = "HUMAN")
#  View(cell_result) #take a look at the data

## ---- eval=FALSE--------------------------------------------------------------
#  hit_plotcell(cell_result,
#               cat_col_list = list("CC" = "#FB4F0B",
#                                   "CN" = "#0FAEB9",
#                                   "NC" = "#E7B700"))
#  #it takes some time to run, since the plot contains a lot of information

## ---- eval=FALSE--------------------------------------------------------------
#  find_in_pubmed(mydata_hit, feat = "cell", imp_by_hitlist = TRUE, condition = "S",
#                 language = "english", year_rg = "2000:2021", your_API = NULL,
#                 newfolder_name = "elutriation_pubmed_search")

## ---- eval=FALSE--------------------------------------------------------------
#  get_treat_level(mydata_caldiff)
#  get_treat_level(mydata_ave)

## ---- eval=FALSE--------------------------------------------------------------
#  dl <- list(elutriation, chemarrest)
#  new_data <- join_cetsa(dl, new_names = c("elutriation", "chemarrest"))
#  new_data_nosuffix <- join_cetsa(dl, new_names = c("", "")) #without new names
#  
#  #just rename conditions
#  chemarrest_new_name <- join_cetsa(chemarrest, "new")

## ---- eval=FALSE--------------------------------------------------------------
#  join_drugdata(drug_data$data, by = c("id", "description"))

## ---- eval=FALSE--------------------------------------------------------------
#  #check between two data frames
#  pr_elutriation <- unique(elutriation$id)
#  pr_chemarrest <- unique(chemarrest$id)
#  pr_list <- list("elutriation" = pr_elutriation, "chemarrest" = pr_chemarrest)
#  
#  result_common <- com_protein_loop(pr_list)
#  result_common #A list which contains the proteins only identified in the elutriation experiment, the chemarrest experiment and both
#  
#  #check hits between conditions
#  all_hits <- rbind(hitlist_elutriation, hitlist_chemarrest)
#  pr_list <- list()
#  for(i in unique(all_hits$Condition)){
#    pr_list[[i]] <- (all_hits %>% dplyr::filter(Condition == i))$id  #get hits for each condition in a hit
#  }
#  pr_list
#  pr_list <- com_protein_loop(pr_list)  #A list which contains all the common and unique hits between all the condtions

## ---- eval=FALSE--------------------------------------------------------------
#  all_hits$drug <- rep("c", nrow(all_hits))
#  for (i in names(pr_list)){
#    all_hits$drug[which(!is.na(match(all_hits$id, pr_list[[i]])))] <- i
#  }
#  View(all_hits) #take a look at the data

## ---- fig.show='hold', echo=FALSE, fig.width=6, fig.height=6------------------
library(ggplot2)
ggplot(data.frame(x = c(1,1,2,2,1), y = c(1,2,2,1,1)), aes(x,y)) + geom_point() + geom_path() +
  xlim(c(0.75,2.25)) + ylim(c(0.75,2.25)) + 
  geom_point(data = data.frame(x = 1.5, y = 1.5), color = "red", size = 5) +
  geom_segment(data = data.frame(x1 = c(1,1), x2 = c(2,2), y1 = c(1,2), y2 = c(2,1)), 
               aes(x = x1, y = y1, xend = x2, yend = y2),
               color = "blue") +
 geom_curve(
  aes(x = x1, y = y1, xend = x2, yend = y2),
  data = data.frame(x1 = c(1.5, 1.75, 1.5, 1.25), x2 = c(1.75, 1.5, 1.25, 1.5), 
                    y1 = c(1.25, 1.5, 1.75, 1.5), y2 = c(1.5, 1.75, 1.5, 1.25)),
  arrow = arrow(length = unit(0.03, "npc"))
  ) + 
  geom_text(data = data.frame(x = c(1.15,1.5,1.85,1.5), y = c(1.5,1.15,1.5,1.85),
                              text = rep("90°", 4)),
            aes(x, y, label = text))

## ---- eval=FALSE--------------------------------------------------------------
#  square <- data.frame(x = c(1,1,2,2), y = c(1,2,2,1))
#  point_in <- c(1.5,1.5)
#  point_out <- c(2.1,2)
#  
#  is_in_zone(square, point_in)  #TRUE
#  is_in_zone(square, point_out) #FALSE
#  
#  #You can create way more complicated border. You'll just need a data frame of 2 columns named x and y.

