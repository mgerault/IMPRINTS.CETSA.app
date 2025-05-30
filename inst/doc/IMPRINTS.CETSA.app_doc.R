## ----eval=FALSE---------------------------------------------------------------
# if(!requireNamespace("BiocManager", quietly = TRUE)){ #check if you already have the BiocManager package
#  install.packages("BiocManager")  #if not, install it
# }
# BiocManager::install("STRINGdb")
# BiocManager::install("clusterProfiler")
# BiocManager::install("biomaRt")
# BiocManager::install("enrichplot")
# BiocManager::install("multtest")

## ----eval=FALSE---------------------------------------------------------------
# if(!requireNamespace("devtools", quietly = TRUE)){ #check if you already have the devtools package
#  install.packages("devtools")  #if not, install it
# }
# devtools::install_github("mgerault/IMPRINTS.CETSA.app")

## ----eval=FALSE---------------------------------------------------------------
# devtools::install_github("nkdailingyun/IMPRINTS.CETSA")

## ----message=FALSE, eval=FALSE------------------------------------------------
# library("IMPRINTS.CETSA.app")

## ----eval=FALSE---------------------------------------------------------------
# runCETSAapp()      #this function will directly start the app

## ----eval=FALSE---------------------------------------------------------------
# mydata_caldiff <- drug_data$data$elutriation  #save your data
# mydata_ave <- drug_data$data_ave$elutriation
# 
# View(mydata_caldiff)  #take a look at the data
# View(mydata_ave)

## ----eval=FALSE---------------------------------------------------------------
# imprints_barplotting_app(mydata_caldiff, ret_plot = FALSE, save_pdf = TRUE) #will save all the barplots from the data
# 
# mydata_caldiff_chemarrest <- drug_data$data$chemarrest
# all_drug <- list("elutriation" = mydata_caldiff, "chemarrest" = mydata_caldiff_chemarrest)
# imprints_barplotting_app(all_drug, ret_plot = FALSE, save_pdf = TRUE) #will save all the barplots from the data
# #here, it will "join" the two pdf in one

## ----eval=FALSE---------------------------------------------------------------
# library(IMPRINTS.CETSA)
# mydata_hit <- drug_data$hitlist$elutriation
# #will map the proteins categorized as CN, CC and NC in S treatment to some protein complex (core Corum database)
# map_compl <- imprints_complex_mapping(mydata_ave, mydata_hit, treatment = "S",
#                                       targetcategory = c("CN", "CC", "NC"))
# View(map_compl) #if you want to check the data
# 
# map_compl <- map_compl[, c("ComplexName", "subunitsNum", "subunitsIdentifiedNum",  #keep columns of interest
#                                  "id", "description", "gene", "category")]
# 
# if(nrow(map_compl) !=0){
#   map_compl$description <- unname(unlist(sapply(map_compl$description, IMPRINTS.CETSA .app:::getProteinName))) #keep only protein names in description
# }

## ----eval=FALSE---------------------------------------------------------------
# data_l <- list() #initialize list
# some_complex <- unique(map_compl$ComplexName)
# some_complex <- some_complex[sample(1:length(some_complex), 4)]
# for(i in some_complex){
#   cate_ <- map_compl[which(map_compl$ComplexName == i), ]
#   pr_comp <- unique(cate_$id)
# 
#   data_l[[i]] <- ms_subsetting(mydata_caldiff, isfile = F, hitidlist = c(pr_comp)) #filter the proteins from those complexes
# 
#   data_l[[i]]$category <- cate_$category[which(!is.na(match(cate_$id, data_l[[i]]$id)))] #keep category
# }
# data <- data_l
# 
# #each element of the list contains data from caldiff output for each protein complex selected
# 
# imprints_barplotting_app(data_l,save_pdf = TRUE, ret_plot = FALSE,
#                      toplabel = "IMPRINTS-CETSA bar plotting \nProtein complex :",  #title on each page (will add the complex name, cf. the names of the list)
#                      pdfname = "complex_barplot"
#                      )

## ----eval=FALSE---------------------------------------------------------------
# imprints_barplotting_simprof(mydata_caldiff, mydata_ave, #here we have the average data, but if it wasn't the case we would have let it NULL and let the function calculate it for us
#                           treatmentlevel = "S", protein_profile = "O43776",
#                           pdfname = "similar_O43776",
#                           use_score = "euclidean",
#                           score_threshold = 0.65, max_na_prow = 0
#                           )
# #Here we choose the Euclidean distance score with a threshold of 0.65 and no missing values

## ----eval=FALSE---------------------------------------------------------------
# # According their category
# #here we need the average data
# imprints_heatmap(mydata_ave, hit_summary = mydata_hit,
#               treatment = "S", max_na = 0,
#               response = "both",
#               select_cat = c("CC", "CN", "NC"),
#               saveHeat = TRUE, file_type = "png")
# #Will save a heatmap in a png file, of the proteins categorized as CN, CC or NC under the treatment S
# 
# # According to their protein complex
# some_complex <- unique(map_compl$ComplexName)
# some_complex <- some_complex[sample(1:length(some_complex), 4)]
# cate_ <- map_compl[which(map_compl$ComplexName == i), ]
# pr_comp <- cate_$id
# 
# data <- ms_subsetting(mydata_ave, isfile = F, hitidlist = c(pr_comp)) #filter the proteins from those complexes
# 
# imprints_heatmap(data, PRcomplex_data = cate_,
#               treatment = "S", max_na = 0,
#               response = "both",
#               saveHeat = TRUE, file_type = "png")
# #Will save a heatmap in a png file, of the proteins of the complex from some_complex, under the treatment S

## ----eval=FALSE---------------------------------------------------------------
# norm_data <- readr::read_tsv("path_to_your_normalized_data.txt")
# new_hits <- imprints_IS(norm_data, ctrl = "Vehicle")

## ----eval=FALSE---------------------------------------------------------------
# library(STRINGdb)
# dir.create("STRING_data")   # create folder to stor STRING data
# string_db <- STRINGdb$new(version="11", species=9606,               #ID 9606 correspond to human
#                            score_threshold=200,
#                            input_directory=  file.path(getwd(), "STRING_data")) #will save the data in a folder named STRING_data

## ----eval=FALSE---------------------------------------------------------------
# data <- mydata_hit %>% dplyr::filter(treatment == "S")
# network_hit <- string_db$map(data, "id", removeUnmappedRows = TRUE)         #will map the proteins to the string ID
# 
# get_net_app(network_hit$STRING_id , inter = FALSE) #will plot the network

## ----eval=FALSE---------------------------------------------------------------
# enrich_hit <- string_db$get_enrichment(network_hit$STRING_id)
# View(enrich_hit) #take a look at the data
# #Each protein has been mapped to a category (Component, function, ...) with a description of this category.
# #For example, with the category component we can have the description ribosome
# 
# #sum up the result, will take only the category "Component"
# df <- get_GO_app(enrich_hit, enrich = TRUE, all_cat = FALSE, sing_cat = "Component") #a list of data frame
# d_n <- as.list(lapply(df, names)[["Component"]]) #get names of the list : the string ID and the gene name
# #as.list() to keep the class list and use it in mapply after
# 
# #add the column id in each data frame (the names we got at the last step)
# df <- mapply(function(x,y) {x$id <- rep(y, nrow(x)); x},
#                      df[["Component"]],
#                      d_n, SIMPLIFY = FALSE)
# df <- do.call(rbind, df)                    #put all the data frame in one
# rownames(df) <- 1:nrow(df)
# View(df) #take a look at the data
# 
# #keep the id and separate the gene name from the string ID
# id_string <- do.call(rbind, stringr::str_split(df$id, ","))
# colnames(id_string) <- c("gene.names", "STRING_id")
# 
# #bind the two data frames
# df <- cbind(id_string, df[,-ncol(df)])
# View(df) #take a look at the final data
# 
# #Let's see the network from some proteins mapped to a specific description
# #for example the ribosome and the cytosolic ribosome
# ribosome_net <- df$STRING_id[stringr::str_which(df$description, "^Ribosome$")]
# get_net_app(ribosome_net , inter = FALSE)

## ----eval=FALSE---------------------------------------------------------------
# hits_G2 <- mydata_hit %>% dplyr::filter(treatment == "G2")
# network <- imprints_network(mydata_caldiff, hits_G2$id,
#                             GOterm = "Component",      # perform enrichment analysis from 'Component' database
#                             required_score = 900)
# 
# network  # to see it
# # it returns a visNetwork object so you can modify it with visNetwork functions

## ----eval=FALSE---------------------------------------------------------------
# # Let's test it on the S treatment
# data <- mydata_hit %>% dplyr::filter(treatment == "S")
# data <- data %>% dplyr::filter(!is.na(match(category, c("CC", "CN", "NC")))) #filter the category, but you could have take them all, even add also the NN
# 
# cell_result <- hit_for_cell(data, organism = "HUMAN")
# View(cell_result) #take a look at the data

## ----eval=FALSE---------------------------------------------------------------
# hit_plotcell(cell_result,
#              cat_col_list = list("CC" = "#FB4F0B",
#                                  "CN" = "#0FAEB9",
#                                  "NC" = "#E7B700"))
# #it takes some time to run, since the plot contains a lot of information

## ----eval=FALSE---------------------------------------------------------------
# # to use find_in_pubmed with parameter imp_by_hitlist set to TRUE, need the description column, i.e. the protein name
# mydata_hit$description <- mydata_caldiff$description[which(!is.na(match(mydata_hit$id, mydata_caldiff$id)))]
# mydata_hit <- mydata_hit[1:20,] # take a sample of the proteins
# find_in_pubmed(mydata_hit, feat = "cell", imp_by_hitlist = TRUE, treatment = "S",
#                language = "english", year_rg = "2021:2022", your_API = NULL,
#                newfolder_name = "elutriation_pubmed_search")

## ----eval=FALSE---------------------------------------------------------------
# get_treat_level(mydata_caldiff)
# get_treat_level(mydata_ave)

## ----eval=FALSE---------------------------------------------------------------
# dl <- list(elutriation, chemarrest)
# new_data <- join_cetsa(dl, new_names = c("elutriation", "chemarrest"))
# new_data_nosuffix <- join_cetsa(dl, new_names = c("", "")) #without new names
# 
# #just rename treatments
# chemarrest_new_name <- join_cetsa(chemarrest, "new")

## ----eval=FALSE---------------------------------------------------------------
# join_drugdata(drug_data$data, by = c("id", "description"))

## ----eval=FALSE---------------------------------------------------------------
# #check between two data frames
# pr_elutriation <- unique(elutriation$id)
# pr_chemarrest <- unique(chemarrest$id)
# pr_list <- list("elutriation" = pr_elutriation, "chemarrest" = pr_chemarrest)
# 
# result_common <- com_protein_loop(pr_list)
# result_common #A list which contains the proteins only identified in the elutriation experiment, the chemarrest experiment and both
# 
# #check hits between treatments
# all_hits <- rbind(hitlist_elutriation, hitlist_chemarrest)
# pr_list <- list()
# for(i in unique(all_hits$treatment)){
#   pr_list[[i]] <- (all_hits %>% dplyr::filter(treatment == i))$id  #get hits for each treatment in a hit
# }
# pr_list
# pr_list <- com_protein_loop(pr_list)  #A list which contains all the common and unique hits between all the condtions

## ----eval=FALSE---------------------------------------------------------------
# all_hits$drug <- rep("c", nrow(all_hits))
# for (i in names(pr_list)){
#   all_hits$drug[which(!is.na(match(all_hits$id, pr_list[[i]])))] <- i
# }
# View(all_hits) #take a look at the data

## ----fig.show='hold', echo=FALSE, fig.width=6, fig.height=6-------------------
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

## ----eval=FALSE---------------------------------------------------------------
# square <- data.frame(x = c(1,1,2,2), y = c(1,2,2,1))
# point_in <- c(1.5,1.5)
# point_out <- c(2.1,2)
# 
# is_in_zone(square, point_in)  #TRUE
# is_in_zone(square, point_out) #FALSE
# 
# #You can create way more complicated border. You'll just need a data frame of 2 columns named x and y.

## ----eval=FALSE---------------------------------------------------------------
# library(IMPRINTS.CETSA .app)
# 
# peptides <- imprints_read_peptides(peptides_files = list.files("Analysis_files",     # your files
#                                                                pattern = "Peptide",
#                                                                full.names = T),
#                                    treatment = c("B1_Vehicle", "B1_DrugA", "B1_DrugB",  # the treatments corresponding to your channel
#                                                  "B2_Vehicle","B2_DrugA", "B2_DrugB",
#                                                  "B3_Vehicle", "B3_DrugA", "B3_DrugB",
#                                                  "Mix"),
#                                    temperatures = c("37C", "47C", "50C", "52C", "54C", "57C"), # temperatures from your corresponding files
#                                    proteins = "Analysis_files/proteins_imprints_caldiff.txt") # The proteins from which you want to analyze the peptides, can be a dataframe

## ----eval=FALSE---------------------------------------------------------------
# peptides_norm <- imprints_normalize_peptides(peptides)

## ----eval=FALSE---------------------------------------------------------------
# # one of your hitlist
# proteins <- hitlist$id
# peptides_norm_diff <-  imprints_sequence_peptides(peptides_norm,
#                                                   proteins =  proteins,
#                                                   barplot = TRUE,
#                                                   control = "Vehicle") # needs to specify control from your experiment

## ----eval=FALSE---------------------------------------------------------------
# # no barplotting this time as it would be quite long
# peptidesall_norm_diff <-  imprints_sequence_peptides(peptides_norm, control = "Vehicle",
#                                                      barplot = FALSE)
# 
# potential_cleaved <- imprints_cleaved_peptides(peptides_norm, peptidesall_norm_diff,
#                                                control = "Vehicle", FDR = 0.01)

## ----eval=FALSE---------------------------------------------------------------
# # compute potential cleaved hit list
# potential_cleaved_notcat <- imprints_cleaved_peptides(peptides_norm, peptidesall_norm_diff,
#                                                       control = "Vehicle", FDR = 0.01,
#                                                       categorize = FALSE)
# 
# # perform categorization
# potential_cleaved_cat <- imprints_categorize_peptides(peptides_norm, potential_cleaved_notcat,
#                                                       control = "Vehicle")

## ----eval=FALSE---------------------------------------------------------------
# imprints_barplotting_categorize_peptides(peptides_norm, potential_cleaved_cat,
#                                          treatment = "DrugB", control = "Vehicle",
#                                          color = "red")

## ----eval=FALSE---------------------------------------------------------------
# ptms_ambiguity <- imprints_ptms_peptides(peptides_norm, peptidesall_norm_diff, control = "Vehicle",
#                                          PTM_direction = "both")

## ----eval=FALSE---------------------------------------------------------------
# isoform_ambiguity <- imprints_isoform_peptides(peptides_norm, peptidesall_norm_diff,
#                                                control = "Vehicle", fasta = "path/to/your/FASTA_file.fasta")

## ----eval=FALSE---------------------------------------------------------------
# imprints_plotting_isoform_peptides(peptides_norm, isoform_ambiguity,
#                                    control = "Vehicle", treatment = "DrugB",
#                                    ret_plot = FALSE, save_pdf = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# ## Drug A
# # keeping only drug A
# peptides_norm_drugA <- peptides_norm[,-grep("_DrugB$", colnames(peptides_norm))]
# potential_cleaved_drugA <- potential_cleaved %>%
#   dplyr::filter(treatment == "DrugA") %>%
#   select(id, cleaved_site)
# # computing fold-changes
# peptides_cleaved_drugA <-  imprints_sequence_peptides(peptides_norm_drugA,
#                                                       proteins =  potential_cleaved_drugA$id,
#                                                       sequence = potential_cleaved_drugA$cleaved_site, # give sequence = the potential cleavage sites
#                                                       control = "Vehicle") # needs to specify control from your experiment
# # removing peptides considred as cleavage sites
# peptides_cleaved_drugA <- imprints_remove_peptides(peptides_cleaved_drugA,
#                                                    proteins =  potential_cleaved_drugA$protein,  # the protein from which you want to remove the specific sequence
#                                                    sequence = potential_cleaved_drugA$cleaved_site) # the sequence you want to remove
# 
# ## Drug B
# peptides_norm_drugB <- peptides_norm[,-grep("_DrugA$", colnames(peptides_norm))]
# potential_cleaved_drugB <- potential_cleaved %>%
#   dplyr::filter(treatment == "DrugB") %>%
#   select(id, cleaved_site)
# # computing fold-changes
# peptides_cleaved_drugB <-  imprints_sequence_peptides(peptides_norm_drugB,
#                                                       proteins =  potential_cleaved_drugB$id,
#                                                       sequence = potential_cleaved_drugB$cleaved_site, # give sequence = the potential cleavage sites
#                                                       control = "Vehicle") # needs to specify control from your experiment
# # removing peptides considred as cleavage sites
# peptides_cleaved_drugB <- imprints_remove_peptides(peptides_cleaved_drugB,
#                                                    proteins =  potential_cleaved_drugB$protein,  # the protein from which you want to remove the sepecific sequence
#                                                    sequence = potential_cleaved_drugB$cleaved_site) # the sequence you want to remove

## ----eval=FALSE---------------------------------------------------------------
# peptides_cleaved <- list(peptides_cleaved_drugA, peptides_cleaved_drugB) # you could add more datasets
# 
# peptides_cleaved <- imprints_join_peptides(peptides_cleaved, mode = "partial")

## ----eval=FALSE---------------------------------------------------------------
# # save the imprints in a pdf
# imprints_barplotting_peptides(peptides_cleaved, format = "RESP_peptide",
#                               save_pdf = T, ret_plot = F)

