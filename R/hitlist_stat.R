#' hitlist_stat
#'
#' Function to categorize protein according to their Expression/Stability change
#' based on two treatment conditions and characteristic conditions for fold change.
#'
#' @param data The input data set from which categorization is performed on and hitlist is produced from, ie. caldiff output.
#'             You can also type a file path. The file must be in txt, xlsx or csv file. If in csv, check
#'             if in the description column there are no commas. Otherwise, it could lead to a misreading of the data.
#'             (The function will remind you to do so)
#' @param control The name(s) of the control(s). If several, separate them by '|'.
#' @param rem_proteomic The name of the proteomic channel (example : '36C'). Can be NULL.
#' @param basetemp The name of the base temperature (example : '37C').
#' @param RSS_cutoff The residual sum of squares cutoff between the base temperature and the others.
#'                   This cutoff is for differentitate the CN and the CC.
#'                   The RSS is reduced to a score between 0 and 1 (1/(1+RSS)).
#'                   The more it's close to one, the more the RSS is close to 0. So the more you choose
#'                   a cutoff close to one, the more you will consider small variations between the base temperatures and the
#'                   others as thermal stability change (categorized as CC then). Default is 0.985.
#' @param ErrScore The error score threshold, a score between 0 and 1; it is calculated thanks to the  standard error.
#'                 From the SE, you will get the error bar of each mean value for each protein for each temperature and 
#'                 conditions. Then, from this error value, the score is calculated as follow : 
#'                 \code{1/(1+abs((error/mean)))}. So the more it's close to one, the you're confident about your measurment.
#'                 You can choose a value or not. If not, the algorithm will plot the histogramm of the error  score of the hits
#'                 and take the error score that gives the highest frequency. If several conditions, will take the
#'                 maximum threshold value found between the conditions.
#' @param use_prompt Logical to tell if you want to use the prompt function (added for the shiny app)
#' @param exported Logical to tell if you want to export the results. Use it when use_prompt is FALSE.
#' @param format The format of the saved files. Format available : xlsx, txt or csv.
#'
#'
#' @return A list which contains the hitlist, the summary and the NN
#'
#' @examples \dontrun{
#' library(mineCETSAapp)
#' PI3K_data <- mineCETSAapp::PI3K1h6h_file
#' HITS <- hitlist_stat(PI3K_data, control = "Vehicle1h|Vehicle6h",
#'                      rem_proteomic = "36C", format = "xlsx")
#' }
#' 
#' @export
#'
#

hitlist_stat <- function (data = NULL, control = NULL,  
                          rem_proteomic = NULL, 
                          basetemp = "37C", RSS_cutoff = 0.985,
                          ErrScore = NULL,
                          use_prompt = TRUE,
                          exported = FALSE, format = c("xlsx", "txt", "csv")){
  
  dname <- deparse(substitute(data))
  dname <- str_remove_all(dname, ":")
  
  if(length(str_subset(class(data), "^character$")) >= 1){
    if(file.exists(data)){
      extension <- str_split(data, "\\.")[[1]]
      if(length(extension) == 1){
        stop("Don't forget to type the format of your file !")
      }
      
      dname <- extension[length(extension)-1]
      dname <- str_split(dname, "/")[[1]]
      dname <- dname[length(dname)]
      
      extension <- extension[length(extension)]
      if(extension == "txt"){
        data <- mineCETSA::ms_fileread(data)
      }
      else if(extension == "xlsx"){
        data <- dplyr::as_tibble(import(data))
      }
      else if(extension == "csv"){
        no_comma <- ''
        while(!(no_comma %in% c('YES','NO','Y','N')) ){
          no_comma <- toupper(readline(prompt = "Are you sure there is no comma in your description column? Otherwise, it can lead to a misreading of your data (Yes/No): "))
          if (no_comma %in% c('YES','Y')){
            message("Let's get this hits then !")
            data <- dplyr::as_tibble(import(data))
          }
          else if (no_comma %in% c('NO','N')) {
            message("Go check your file then and re-run the function after. I recommend you to use txt or xlsx files so.")
            return()
          }
          else {
            print("Invalid choice")
          }
        }
      }
      else{
        stop("This file format is not supported yet. Please, use txt, xlsx or eventually csv files.")
      }
    }
    else{
      stop("This file doesn't exists ! Check your file name")
    }
  }
  else if(sum(class(data) == "data.frame") == 0){
    stop("Please, provide either a data.frame or a file name")
  }
  
  if(use_prompt){
    exporting <- ''
    while(!(exporting %in% c('YES','NO','Y','N')) ){
      exporting <- toupper(readline(prompt = "Do you want to export the results (the hitlist, the summary and the NN) ? (Yes/No): "))
      if (exporting %in% c('YES','Y')){
        exported <- TRUE
      }
      else if (exporting %in% c('NO','N')) {
        exported <- FALSE
      }
      else {
        print("Invalid choice.")
      }
    }
  }
  
  
  if (length(grep("countNum", names(data)))) {
    countinfo1 <- data[,str_which(names(data), "^id$|^description$|^sumPSM|^countNum|^sumUniPeps")]
    data <- data[, -str_which(names(data), "^sumPSM|^countNum|^sumUniPeps")]
  }
  subset <- str_which(names(data), paste0("_", control, "$"))
  if (length(subset) > 0) {
    data1 <- data[, -subset]
  }
  else {
    stop("Please provide the right treatment keyword character")
  }
  
  if (length(grep("description", names(data1)))) {
    proteininfo <- unique(data1[, c("id", "description")])
    data1$description <- NULL
  }
  
  if(!is.null(rem_proteomic) & is.character(rem_proteomic)){
    data1 <- data1[,-str_which(names(data1), paste0("^", rem_proteomic, "_"))]
  }

  message("Getting mean, standard deviation, error bar, etc. This may take a while.")
  data1 <- tidyr::gather(data1, condition, reading, -id)
  a <- data1$condition[1]
  if (length(unlist(strsplit(a, "_"))) == 4) {
    data1 <- tidyr::separate(data1, condition, into = c("set",
                                                        "temperature", "replicate", "treatment"), sep = "_")
    if(length(str_subset(basetemp, unique(data1$temperature))) != 1){
      stop(paste("Please provide only one of this temperature :", paste(unique(data1$temperature), collapse = ", ")))
    }
    
    data1$id <- factor(data1$id, levels = unique(data1$id), ordered = TRUE) #preserve order
    cdata <- plyr::ddply(data1, c("id", "set", "temperature",
                                  "treatment"),
                         summarise, N = length(na.omit(reading)),
                         mean = mean(reading, na.rm = T), sd = sd(reading, na.rm = T), se = sd/sqrt(N))
    cdata$id <- as.character(cdata$id)
    
  }
  else if (length(unlist(strsplit(a, "_"))) == 3) {
    data1 <- tidyr::separate(data1, condition, into = c("temperature",
                                                        "replicate", "treatment"), sep = "_")
    
    if(length(str_subset(basetemp, unique(data1$temperature))) != 1){
      stop(paste("Please provide only one of this temperature :", paste(unique(data1$temperature), collapse = ", ")))
    }
    
    cdata <- plyr::ddply(data1, c("id", "temperature", "treatment"),
                         summarise, N = length(na.omit(reading)), mean = mean(reading,na.rm = T),
                         sd = sd(reading, na.rm = T), se = sd/sqrt(N))
  }
  else {
    stop("make sure the namings of the columns of the dasaset are correct.")
  }
  
  cdata$error <- abs((cdata$mean - cdata$se) - (cdata$mean + cdata$se))
  cdata$error.score <- 1/(1+abs((cdata$error/cdata$mean)))
  
  message("Getting outlier")
  
  hits <- list() #hitlist
  
  for(i in unique(cdata$treatment)){
    message(i)
    tr <- cdata %>% 
      dplyr::filter(treatment == i)
    for(k in unique(cdata$temperature)){ ### get hits for each temperature
      message(k)
      meh <- tr %>% 
        dplyr::filter(temperature == k)
      meh$temperature <- NULL
      meh <- na.omit(meh)           ### na will be removed anyway 
      meh <- meh[order(meh$mean),]  ### order proteins (not necessary)
      
      meh$idx <- 1:nrow(meh)
      is.outlier <- function (x) {
        x < quantile(x, .25) - 1.5 * IQR(x) |
          x > quantile(x, .75) + 1.5 * IQR(x)
      }
      
      meh$is_out <- is.outlier(meh$mean)
      hits[[k]] <- meh$id[meh$is_out]
    }
  }
  
  message("Prepare data for categorization")
  hits_list <- list()
  for(i in unique(cdata$temperature)){
    y <- cdata %>% dplyr::filter(temperature == i)
    y$is.outlier <- rep(FALSE, nrow(y))
    
    hits_list[[i]] <- y
  }
  hits_list <- mapply(function(x, y) {x$is.outlier <- x$id %in% y; x}, hits_list, hits, SIMPLIFY = FALSE)
  hits_list <- lapply(hits_list, function(y) full_join(y, countinfo1, by = "id"))
  hits_list <- do.call(rbind, hits_list)
  hits_list <- hits_list[order(hits_list$id),]
  rownames(hits_list) <- 1:nrow(hits_list)
  
  reorder_hits_list <- (1:ncol(hits_list))[!((1:ncol(hits_list)) %in% c(1,str_which(names(hits_list), "^description$")))]
  hits_list <- hits_list[,c(str_which(names(hits_list), "^id$|^description$"), reorder_hits_list)]
  
  hits_list <- hits_list %>% dplyr::group_by(id, treatment) %>% 
    dplyr::mutate(n_out = sum(is.outlier))
  
  NN <- hits_list %>% dplyr::filter(n_out == 0)
  NN$category <- "NN"
  NN <- NN %>% dplyr::select(id, description, treatment, temperature, mean, error, error.score, category)
  colnames(NN)[which(colnames(NN) == "treatment")] <- "Condition"
  
  
  hits_list <- hits_list %>% dplyr::filter(n_out >= 1)
  hits_list$ErrThreshold <- rep(1, nrow(hits_list))
  
  if(is.null(ErrScore)){
    err_thresh <- c()
    for(i in unique(hits_list$treatment)){
      err <- (hits_list %>% filter(treatment == i))$error.score
      s <- hist(err, breaks = 20, main = paste("Error score histogram \nCondition :", i))
      err <- (which.max(s$counts) - 1)/20
      err_thresh <- append(err_thresh, err)
      hits_list$ErrThreshold[which(hits_list$treatment == i)] <- err
    }
  }
  else{
    err_thresh <- ErrScore
    hits_list$ErrThreshold <- rep(ErrScore, nrow(hits_list))
  }
  hits_list <- hits_list %>% dplyr::group_by(treatment) %>% dplyr::mutate(ErrOK = error.score >= ErrThreshold)
  
  
  ### reorganize data
  hits_list <- hits_list[,c("id", "description", "treatment", "temperature", "mean", "is.outlier", "n_out", "error", "error.score", "ErrOK")]
  hits_list <- hits_list[order(hits_list$n_out, decreasing = TRUE),]
  
  n_chan <- length(unique(hits_list$temperature))
  
  hits_list <- hits_list %>% dplyr::group_by(id, treatment) %>%
    dplyr::mutate(nb_conf_out = ErrOK*is.outlier)
  hits_list <- hits_list %>% dplyr::group_by(id, treatment) %>% 
    dplyr::mutate(nb_conf_out = sum(nb_conf_out, na.rm = TRUE))
  
  hits_list <- hits_list %>% dplyr::group_by(id, treatment) %>%
    dplyr::mutate(nb_conf_notout = (!is.outlier)*ErrOK)
  hits_list <- hits_list %>% dplyr::group_by(id, treatment) %>% 
    dplyr::mutate(nb_conf_notout = sum(nb_conf_notout, na.rm = TRUE))
  
  hits_list <- hits_list %>% dplyr::group_by(id, treatment) %>% 
    dplyr::mutate(nb_conf_notout = nb_conf_notout/(n_chan - n_out), 
                  nb_conf_out = nb_conf_out/n_out)
  
  
  # Definition for hit
  colnames(hits_list)[which(colnames(hits_list) == "treatment")] <- "Condition"
  hits_definition <- hits_list
  
  #Unique id,cond pair for hits
  keys_hits <- hits_definition %>% dplyr::ungroup() %>%  dplyr::select(id,Condition) %>%  dplyr::distinct()
  
  # Reference hitlist
  hitlist <- hits_list %>% dplyr::right_join(keys_hits, by = c('id','Condition'))
  referencelist <- hitlist
  
  ## Extract Not determinable -  ND
  
  message("Getting ND")
  # Proteins with noisy 37C measurement
  ND_condition <- hitlist %>%  dplyr::filter(temperature == basetemp,
                                             ErrOK == FALSE)
  
  keys_ND <- ND_condition %>% dplyr::ungroup() %>% dplyr::select(id,Condition) %>%  dplyr::distinct()
  
  ND <- hitlist %>%  dplyr::right_join(keys_ND,by = c('id','Condition'))
  
  # Proteins without 37C measurements
  NDwo37 <- (hitlist %>% dplyr::filter(temperature == basetemp,
                                       is.na(mean) == TRUE))$id
  NDwo37 <- hitlist[which(!is.na(match(hitlist$id, NDwo37))),]
  
  # Combine all of ND
  ND <- ND %>%
    dplyr::full_join(NDwo37, by = names(ND)) %>%
    dplyr::group_by(category = 'ND')
  
  
  # Remove ND from hitlist
  hitlist <- hitlist %>% dplyr::anti_join(ND, by = c('id','Condition'))
  
  ## Stability change w/out expression change  - NC
  message("Getting NC")
  # High Temp. hits
  highT <- hitlist %>% dplyr::filter(temperature != basetemp,
                                     is.outlier == TRUE)
  keys_highT <- highT %>% dplyr::ungroup() %>%  dplyr::select(id, Condition) %>%  dplyr::distinct()
  
  # Low Temp., small mean and well measured
  lowT <- hitlist %>% dplyr::filter(temperature == basetemp,
                                    is.outlier == FALSE,
                                    ErrOK == TRUE)
  keys_lowT <- lowT %>%  dplyr::ungroup() %>%  dplyr::select(id,Condition) %>%  dplyr::distinct()
  
  # Match id,condition pair fulfilling NC condition
  NC_keys <- dplyr::inner_join(keys_highT,keys_lowT, by = c('id','Condition'))
  
  
  NC <- dplyr::right_join(hitlist,NC_keys,by = c('id','Condition')) %>%
    dplyr::group_by(category = 'NC')
  
  # Remove NC from hitlist
  hitlist <- hitlist %>%  dplyr::anti_join(NC, by = c('id','Condition'))
  
  
  if(nrow(hitlist) > 0){
    message("Getting CN and CC")
    ## Expression change w/out stability change - CN
    
    ### DEVIATION FROM REFERENCE MEAN CONDITION
    
    # 37C with high mean
    CN_referencekeys <- hitlist %>%
      dplyr::filter(temperature == basetemp, is.outlier == TRUE) %>%
      dplyr::ungroup () %>%
      dplyr::select(id,Condition) %>%
      dplyr::distinct()
    
    CN_keep_meandev <- c()
    CN_discard_meandev <- c() # Will be part of CC
    
    for (idx in 1:nrow(CN_referencekeys)) {
      iterTibble <- hitlist %>% dplyr::filter(id == CN_referencekeys$id[idx],
                                              Condition == CN_referencekeys$Condition[idx])
      
      refmean <- (iterTibble %>% dplyr::filter(temperature == basetemp))$mean
      deviateTibble <- iterTibble %>% dplyr::filter(temperature != basetemp)
      
      ### R squared on h line 
      rss <- sum((deviateTibble$mean - refmean)**2, na.rm = TRUE)
      # if rss == 0, means all other temperature are NAs
      rss <- 1/(1+rss)
      
      #ref_over_mean close to one means we might have CN or a CC with base_temperature value close to the mean 
      #so even with a clear CC we would have a R2 close to one 
      #to prevent from this, use only the rss, score it (0 to 1)
      #keep if at least one errok and deviation threshold
      #maybe change any in percentage
      if((any(deviateTibble$ErrOK, na.rm = TRUE) & rss >= RSS_cutoff) | rss == 1){
        CN_keep_meandev <- rbind(CN_keep_meandev,CN_referencekeys[idx,])
      }
      else{
        CN_discard_meandev <- rbind(CN_discard_meandev,CN_referencekeys[idx,])
      }
    }
    
    # Proteins that do not fulfill either of CN conditions are CC - Expression and stability change
    CCkeys <- CN_discard_meandev %>% dplyr::select(id, Condition)
    CC <- dplyr::right_join(hitlist, CCkeys,by = c('id','Condition')) %>%
      dplyr::group_by(category = 'CC')
    
    
    # Remove all CC and reveal CNs
    CN <- dplyr::anti_join(hitlist,CC, by = c('id','Condition')) %>%
      dplyr::group_by(category = 'CN')
    
  }
  else{
    CC <- as.data.frame(matrix(ncol = ncol(ND), nrow = 0))
    colnames(CC) <- colnames(ND)
    CC <- dplyr::as_tibble(CC)
    #match the class of each column
    a <- sapply(ND, function(x) class(x))
    a <- paste0("as.", a)
    for(i in 1:ncol(CC)){
      CC[,i] <- do.call(a[i], list(CC[,i]))
    }
    CN <- CC
  }
  
  ###    Data Wrangling: Compiling and exporting
  
  # Complete set of all categorized proteins
  categorized_full <- CC %>%
    dplyr::full_join(CN, by = colnames(CN)) %>%
    dplyr::full_join(NC, by = colnames(NC)) %>%
    dplyr::full_join(ND, by = colnames(ND))
  
  hitlist <- categorized_full 
  
  hitlist$ord_outlier <- rep(NA, nrow(hitlist))
  for(k in unique(hitlist$Condition)){
    for(i in unique(hitlist$temperature)){
      idx <- which(hitlist$temperature == i & hitlist$Condition == k)
      hitlist$ord_outlier[idx] <- rank(abs(hitlist$mean[idx]), na.last = FALSE)
    }
  }
  
  hitlist <- hitlist %>% dplyr::group_by(id, Condition) %>%
    dplyr::mutate(score_outlier = sum(ord_outlier*is.outlier)/(length(unique(hitlist$id))*n_out))
  hitlist <- hitlist %>% dplyr::group_by(id, Condition) %>%
    dplyr::mutate(cumul_error_out = sum(is.outlier*error.score, na.rm = TRUE))
  hitlist <- hitlist %>% dplyr::group_by(id, Condition) %>%
    dplyr::mutate(cumul_error_out = cumul_error_out/n_out)
  hitlist$score <- hitlist$score_outlier*hitlist$cumul_error_out
  
  hit_summary <- hitlist %>% dplyr::group_by(id,Condition,category, score) %>%  dplyr::summarise()
  
  results <- list()
  results$hitlist <- hitlist
  results$summary <- hit_summary
  results$NN <- NN
  
  g <- ggplot(hitlist %>% filter(temperature == "37C"), 
              aes(mean, error.score, color = category)) + geom_point() +
    labs(x = "Abundance difference at 37°C", y = "Error score",
         color = "category") +
    ylim(c(0,1)) +
    scale_color_manual(values = c("CN" = "#0FAEB9", "NC" = "#E7B700", "CC" = "#FB4F0B", 
                                  "ND" = "#8F3A8461")) +
    facet_grid(~Condition)
  
  if(exported){
    message("Start saving your results")
    if(!file.exists(paste0(dname, "_Hits"))){
      dir.create(paste0(dname, "_Hits"))
    }
    filename <- paste(dname,
                      format(Sys.time(), "%H%M_%e-%m-%y"),
                      sep = "_")
    filename <- file.path(paste0(dname, "_Hits"), filename)
    
    ggsave(paste(filename, "_summaryPlots.png"), g)
    if(format == "xlsx"){
      openxlsx::write.xlsx(hitlist, paste(filename, "_Hitlist.xlsx"))
      openxlsx::write.xlsx(hit_summary, paste(filename, "_Summary.xlsx"))
      openxlsx::write.xlsx(NN, paste(filename, "_NN.xlsx"))
      message("Results saved succesfully !")
    }
    else if(format == "txt"){
      write.table(hitlist, paste(filename, "_Hitlist.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
      write.table(hit_summary, paste(filename, "_Summary.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
      write.table(NN, paste(filename, "_NN.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
      message("Results saved succesfully !")
    }
    else if(format == "csv"){
      write.csv(hitlist, paste(filename, "_Hitlist.csv"), quote = FALSE, row.names = FALSE)
      write.csv(hit_summary, paste(filename, "_Summary.csv"), quote = FALSE, row.names = FALSE)
      write.csv(NN, paste(filename, "_NN.csv"), quote = FALSE, row.names = FALSE)
      message("Results saved succesfully !")
    }
    else{
      message("The results couldn't have been saved ! Please, provide a valid file format (xlsx, txt or csv)")
    }
  }
  
  TR <- unique(hit_summary$Condition)
  nb_hits <- list()
  for(i in TR){
    nb_hits[[i]] <- nrow(hit_summary %>% dplyr::filter(Condition == i))
  }
  if(length(err_thresh) == 1){
    message(paste(paste(nb_hits, "hits have been found under the condition", names(nb_hits), collapse = ", "), 
                  ", with an error score threshold of", err_thresh))
  }
  else if(length(err_thresh) > 1){
    message(paste(nb_hits, "hits have been found under the condition", names(nb_hits), 
                  "with an error score threshold of", err_thresh, collapse = ", "))
  }
  
  
  return(results) 
}

