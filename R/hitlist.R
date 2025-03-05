#' hitlist
#'
#' Function to categorize protein according to their Expression/Stability change
#' based on two treatment conditions and characteristic conditions for fold change.
#'
#' @param inputdata The input data set from which categorization is performed on and hitlist is produced from.
#' @param meancutoff The threshold value for the mean of a triplicate(First filter of hits).
#' @param boundedness Factor of SD/SEM, the mean needs to be within.
#' @param qualitycutoff The threshold value for the SD/SEM to separate well-measured and illmeasured measurements.
#' @param meandev The deviation of individual mean from the group mean in CN categorization.
#' @param use_prompt Logical to tell if you want to use the prompt function (added for the shiny app)
#' @param exported Logical to tell if you want to export the results. Use it when use_prompt is FALSE.
#'
#' @details
#' This function was a first draft of hit selection from which I, Marc-Antoine Gerault, wasn't the author.
#' After corrections and mosifications of tghe scripts, it has been kept in IMPRINTS.CETSA.app for comparison
#' purpose but is not advised for hit selection.
#'
#'
#' @return Dataframe of the export.
#'
#' @export
#'
#'

hitlist <- function(inputdata, meancutoff = 0.25, boundedness = 4,
                    qualitycutoff = 0.15, meandev = 0.15,
                    use_prompt = TRUE, exported = FALSE){
  ##################################################
  #                                                #
  #     Loading data and precursory settings       #
  #                                                #
  ##################################################

  # Handling input data
  if(typeof(inputdata) == 'character'){
    extension <- sub(".*\\.", "", inputdata)
    treatment <- sub(paste0("\\.", extension, "$"), "", inputdata)

    # Import data
    if(extension == 'txt'){
      data <- readr::read_tsv(inputdata, show_col_types = FALSE)
    }
    else if(extension == 'csv'){
      data <- readr::read_csv(inputdata, show_col_types = FALSE)
    }
    else if(extension == 'xlsx'){
      data <- openxlsx::read.xlsx(inputdata)
    }
    else {
      stop('This file format is not supported.')
    }
  }
  else if("data.frame" %in% class(inputdata)){
    treatment <- deparse(substitute(inputdata))
    data <- inputdata
  }
  else{
    stop("Data can either be path to a file or a data.frame")
  }

  if(use_prompt){
    exporting <- ''
    while(!(exporting %in% c('YES','NO','Y','N')) ){
      exporting <- toupper(readline(prompt = "Do you want to export the datasets? (Yes/No): "))
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

  ##################################################
  #                                                #
  #         Data cleaning and preparation          #
  #                                                #
  ##################################################

  # remove quantitative proteomics columns
  if(any(grepl("^36C_", colnames(data)))){
    data <- data[,-grep("^36C_", colnames(data))]
  }

  ###  Saving protein information and count data  ###
  if (length(grep('^id$', colnames(data)))) {
    proteininfo <- (unique(data[,c('id','description')]))
    data <- data[,-grep('^description$', colnames(data))]
  }
  if (length(grep('countNum', colnames(data)))) {
    countinfo <- unique(data[,grep("^id$|^sumPSM|^countNum|^sumUniPeps", colnames(data))])
    data <- data[,-grep("^sumPSM|^countNum|^sumUniPeps", colnames(data))]
  }

  control <- data[1,grep("^\\d{2}C_", colnames(data))]
  control <- colnames(control)[which(control[1,grep("^\\d{2}C_", colnames(control))] == 0)]
  if(length(control)){
    control <- unique(sub(".*_", "", control))
    data <- data[,-grep(paste0("_", control, "$"), colnames(data))]
  }

  ###  Cleaning data and calculating statistics  ###
  data <- data %>%
    tidyr::gather(treatment,value,-id)

  # Cleaning and calculating according to treatment format
  if (length(strsplit(data$treatment[1], '_')[[1]]) == 4){
    data <- data %>%
      tidyr::separate(treatment, into = c('set', 'Temperature', 'Replicate', 'treatment'), sep = '_')
    data_clean <- data %>%
      dplyr::group_by(id, set, Temperature, treatment) %>%
      dplyr::reframe(Mean = mean(value,na.rm=T),
              SD = sd(value,na.rm=T),
              SEM = SD/sqrt(length(na.omit(value))))
  }
  else if (length(strsplit(data$treatment[1],'_')[[1]]) == 3){
    data <- data %>%
        tidyr::separate(treatment, into = c('Temperature', 'Replicate', 'treatment'), sep = '_')

    data_clean <- data %>%
        dplyr::group_by(id,Temperature, treatment) %>%
        dplyr::reframe(Mean = mean(value,na.rm=T),
                SD = sd(value,na.rm=T),
                SEM = SD/sqrt(length(na.omit(value))))

  }
  else if (length(strsplit(data$treatment[1],'_')[[1]]) < 3){
    stop('Make sure the format/label is correct.')
  }


  ##################################################
  #                                                #
  #   Data Wrangling: 'Rule' based categorization  #
  #                                                #
  ##################################################

  ### Categorization according to expression/stability change  ###
  ## Initial hitlist (Reference)
  selection_metrics <- data_clean %>%
      dplyr::group_by(mean_threshold = (abs(Mean) > meancutoff),
                      bounded = (abs(Mean)- boundedness*SEM>0),
                      wellmeasured = (SEM < qualitycutoff))
  hits_definition <- selection_metrics %>%
    dplyr::filter(mean_threshold & bounded)

  #Unique id,cond pair for hits
  keys_hits <- hits_definition %>% dplyr::ungroup() %>%
    dplyr::select(id,treatment) %>%
    dplyr::distinct()

  # Reference hitlist
  hitlist <- selection_metrics %>% dplyr::right_join(keys_hits, by = c('id','treatment'))
  referencelist <- hitlist

  # Separate hits and NN from data
  NN <- selection_metrics %>%
    dplyr::anti_join(hitlist, by = c('id','treatment')) %>%
    dplyr::group_by(category = 'NN')

  ## Extract Not determinable -  ND
  # Proteins with noisy 37C measurement
  ND_condition <- hitlist %>%
    dplyr::filter(Temperature == '37C' & !bounded & !wellmeasured == FALSE)

  keys_ND <- ND_condition %>% dplyr::ungroup() %>%
    dplyr::select(id,treatment) %>%  dplyr::distinct()

  ND <- hitlist %>%
    dplyr::right_join(keys_ND, by = c('id','treatment'))

  # Proteins without 37C measurements
  NDwo37 <- hitlist %>%
    dplyr::filter(Temperature == '37C' & is.na(Mean))

  # Combine all of ND
  ND <- ND %>%
    dplyr::full_join(NDwo37, by = names(ND)) %>%
    dplyr::group_by(category = 'ND')

  # Remove ND from hitlist
  hitlist <- hitlist %>% dplyr::anti_join(ND, by = c('id','treatment'))

  ## Stability change w/out expression change  - NC
  # High Temp. hits
  highT <- hitlist %>%
    dplyr::filter(Temperature != '37C' & mean_threshold & bounded)
  keys_highT <- highT %>% dplyr::ungroup() %>%
    dplyr::select(id, treatment) %>%  dplyr::distinct()

  # Low Temp., small mean and well measured
  lowT <- hitlist %>%
    dplyr::filter(Temperature == '37C' & !mean_threshold == FALSE & wellmeasured == TRUE)
  keys_lowT <- lowT %>% dplyr::ungroup() %>%
    dplyr::select(id, treatment) %>%  dplyr::distinct()

  # Match id,treatment pair fulfilling NC condition
  NC_keys <- dplyr::inner_join(keys_highT, keys_lowT, by = c('id','treatment'))

  NC <- dplyr::right_join(hitlist, NC_keys, by = c('id','treatment')) %>%
    dplyr::group_by(category = 'NC')

  # Remove NC from hitlist
  hitlist <- hitlist %>%
    dplyr::anti_join(NC, by = c('id','treatment'))

  if(nrow(hitlist) > 0){
    ## Expression change w/out stability change - CN
    ### INTERVAL CONDITION FOR CN
    CN_intervals <- hitlist %>%
      dplyr::group_by(upperbound = Mean + SEM, lowerbound = Mean - SEM)
    CN_keys <- CN_intervals %>%
      dplyr::ungroup() %>%
      dplyr::select(id, treatment) %>% dplyr::distinct()

    CN_keep_interval <- c()
    CN_discard_interval <- c() # Will be part of CC

    for (idx in 1:nrow(CN_keys)){
      iterTibble <- CN_intervals %>%
        dplyr::filter(id == CN_keys$id[idx], treatment == CN_keys$treatment[idx]) %>%
        dplyr::select(upperbound, lowerbound)
      overlap <- c()
      for (i in 1:nrow(iterTibble)) {
        refinterval <- iterTibble[i,]
        testinterval <- iterTibble[-i,]
        overlap <- c(overlap, sum(refinterval$lowerbound <= testinterval$upperbound & refinterval$upperbound >= testinterval$lowerbound) > 0)
      }
      if(all(overlap, na.rm = TRUE)){
        CN_keep_interval <- rbind(CN_keep_interval, CN_keys[idx,])
      }
      else {
        CN_discard_interval <- rbind(CN_discard_interval, CN_keys[idx,])
      }
    }


    ### DEVIATION FROM REFERENCE MEAN treatment
    # 37C with high mean
    CN_referencekeys <- hitlist %>%
      dplyr::filter(Temperature =='37C' & mean_threshold) %>%
      dplyr::ungroup () %>%
      dplyr::select(id,treatment) %>% dplyr::distinct()

    CN_keep_meandev <- c()
    CN_discard_meandev <- c() # Will be part of CC

    for (idx in 1:nrow(CN_referencekeys)) {
      iterTibble <- hitlist %>% dplyr::filter(id == CN_referencekeys$id[idx],
                                       treatment == CN_referencekeys$treatment[idx])
      if (nrow(iterTibble %>% dplyr::filter(Temperature !='37C' & bounded)) == 0){
        next # Not bounded high temps are disregarded
      }
      else {
        refmean <- (iterTibble %>% dplyr::filter(Temperature == '37C'))$Mean
        boundediterTibble <- iterTibble %>% dplyr::filter(Temperature != '37C' & bounded)
        deviateTibble <- iterTibble %>% dplyr::filter(Temperature != '37C' & wellmeasured)
        # Bounded high temp. deviate at most meandev AND
        # wellmeasured high temp don't deviate more than meandev
        if ((sum((abs(boundediterTibble$Mean - refmean)< meandev*abs(refmean))) == nrow(boundediterTibble)) &&
             (!any(abs(deviateTibble$Mean - refmean)> meandev*abs(refmean)))  ) {
          CN_keep_meandev <- rbind(CN_keep_meandev,CN_referencekeys[idx,])
        }
        else{
          CN_discard_meandev <- rbind(CN_discard_meandev,CN_referencekeys[idx,])
        }
      }
    }

    # Proteins that do not fulfill either of CN conditions are CC - Expression and stability change
    CCkeys <- dplyr::inner_join(CN_discard_interval, CN_discard_meandev, by = c('id','treatment'))
    CC <- dplyr::right_join(hitlist,CCkeys,by = c('id','treatment')) %>%
      dplyr::group_by(category = 'CC')

    # Remove all CC and reveal CNs
    CN <- dplyr::anti_join(hitlist, CC, by = c('id','treatment')) %>%
      dplyr::group_by(category = 'CN')
  }
  else{
    CC <- as.data.frame(matrix(ncol = ncol(ND), nrow = 0))
    colnames(CC) <- colnames(ND)
    CC <- as_tibble(CC)
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

  ###  Export : Combine with protein information and write to csv  ###
  if (exported){
    # Output content
    out_reference <- dplyr::inner_join(proteininfo,referencelist, by = 'id')
    out_categorized <- dplyr::inner_join(proteininfo,categorized_full, by = 'id')
    out_summary <- categorized_full %>% dplyr::group_by(id,treatment,category) %>%  dplyr::summarize()
    out_NN <- dplyr::inner_join(proteininfo,NN, by = 'id')

    # Output name formatting
    filename_reference <- paste(treatment, format(Sys.time(), "%H%M_%e-%m-%y"), 'Hitlist', sep = '_')
    filename_categorized <- paste(treatment, format(Sys.time(), "%H%M_%e-%m-%y"), 'Categorized', sep = '_')
    filename_NN <- paste(treatment, format(Sys.time(), "%H%M_%e-%m-%y"), 'NN', sep = '_')
    filename_summary <- paste(treatment, format(Sys.time(), "%H%M_%e-%m-%y"), 'Summary', sep = '_')

    # Output csv
    write.csv(out_reference, file = paste0(filename_reference,'.csv'))
    write.csv(out_categorized, file = paste0(filename_categorized,'.csv'))
    write.csv(out_summary, file = paste0(filename_summary,'.csv'))
    write.csv(out_NN, file = paste0(filename_NN,'.csv'))
  }

  ### Return : Summary of run and dataframe with categories for further handling. ###
  summary_categorized <- categorized_full %>%
    dplyr::group_by(id,treatment,category) %>% dplyr::reframe()
  summary_NN <- NN %>%
    dplyr::group_by(id,treatment,category) %>% dplyr::reframe()
  tablecate <- table(summary_categorized[,2:3])
  tableNN <- table(summary_NN[,2:3])

  print('This dataset has the following count for the categories:')
  print(cbind(tablecate,tableNN))

  print('The categories were obtained with the following values for the parameters:')
  print(sprintf(fmt = 'meancutoff: %s',meancutoff))
  print(sprintf(fmt = 'boundedness: %s', boundedness))
  print(sprintf(fmt = 'qualitycutoff: %s',qualitycutoff))
  print(sprintf(fmt = 'meandev: %s',meandev))

  results <- list()
  results$hitlist <- referencelist
  results$ND <- ND
  results$NC <- NC
  results$CN <- CN
  results$CC <- CC
  results$NN <- NN
  return(results)
}
