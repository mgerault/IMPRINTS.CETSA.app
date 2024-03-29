#' hitlist
#'
#' Function to categorize protein according to their Expression/Stability change
#' based on two treatment conditions and characteristic conditions for fold change.
#'
#' @param inputdata The input data set from which categorization is performed on and hitlist is produced from.
#' @param meancutoff The threshold value for the mean of a triplicate(First filter of hits).
#' @param boundedness Number of SD/SEM, the mean needs to be within.
#' @param qualitycutoff The threshold value for the SD/SEM to separate well-measured and illmeasured measurements.
#' @param meandev The deviation of individual mean from the group mean in CN categorization.
#' @param expstr List of strings with treatments/drugs and control group
#' @param use_prompt Logical to tell if you want to use the prompt function (added for the shiny app)
#' @param exported Logical to tell if you want to export the results. Use it when use_prompt is FALSE.
#'
#'
#' @return Dataframe of the export.
#'
#' @export
#'
#'

hitlist <- function(inputdata,
                    meancutoff = 0.25,
                    boundedness = 4,
                    qualitycutoff = 0.15,
                    meandev = 0.15,
                    expstr = c(),
                    use_prompt = TRUE,
                    exported = FALSE){

  ### Package handling
  # Wanted packages
  packages <- c("dplyr", "tidyr", "stringr")

  # Install packages not yet installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
    }

  # Load all wanted packages
  invisible(lapply(packages, library, character.only = TRUE))



  ##################################################
  #                                                #
  #     Loading data and precursory settings       #
  #                                                #
  ##################################################

  # Handling input data
  if(typeof(inputdata) == 'character'){

    inputdata_strsplit <- unlist(strsplit(inputdata,'\\W+'))
    treatment <- inputdata_strsplit[1]
    extension <- inputdata_strsplit[2]

    # Import data
    if(extension == 'txt'){
      data <- read.delim(inputdata)
      }
    else if(extension == 'csv'){
      data <- read.csv2(file = inputdata, sep = ',')
      data <- data[,-1]
      }
    else {
      stop('This file format is not supported.')
      }
    }
  else {
    treatment <- deparse(substitute(inputdata))
    data <- inputdata # Loading data directly from environment

    if (sum(str_detect(names(data), "^\\d{1}")) >= 1){
      names(data)[str_which(names(data), "^\\d{1}")] <- paste0("X", str_subset(names(data), "^\\d{1}"))
    }
    }

  ########## Will ppotentially be removed as SEM will always be used ##########
  dispmeas <- 'SEM'
  while( !(dispmeas %in% c('SEM','SD')) ){
    dispmeas <- toupper(readline(prompt = "SD or SEM ?"))
    if ( !(dispmeas %in% c('SEM','SD')) ){
      print("Make sure to choose a valid measure of dispersion.")
    }
    }
  #############################################################################


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
  if(sum(stringr::str_detect(colnames(data), "36C_"))){
    data <- data[,-stringr::str_which(colnames(data), "36C_")]
  }

  ###  Saving protein information and count data  ###
  if (length(grep('id',colnames(data)))) {
    proteininfo <- (unique(data[,c('id','description')]))
    data <- data[,!(names(data) %in% 'description')]
    }
  if (length(grep('countNum',colnames(data)))) {
    countinfo <- unique(data[,str_which(names(data), "^id$|^sumPSM|^countNum|^sumUniPeps")])
    data <- data[,-str_which(names(data), "^sumPSM|^countNum|^sumUniPeps")]
    }


  ###  Cleaning data and calculating statistics  ###
  data <- data %>%
    gather(treatment,value,-id)

  # Cleaning and calculating according to treatment format
  if (length(unlist(strsplit(data$treatment[1],'_')))==4) {
    data <- data %>%
      separate(treatment, into = c('set',
                                 'Temperature',
                                 'Replicate',
                                 'treatment'),sep = '_')
    data_clean <- data %>%
      dplyr::group_by(id, set,Temperature, treatment) %>%
      dplyr::reframe(Mean = mean(value,na.rm=T),
              SD = sd(value,na.rm=T),
              SEM = SD/sqrt(length(na.omit(value))))
    }
  else if (length(unlist(strsplit(data$treatment[1],'_')))==3) {
    if(any(grepl('[a-zA-Z][0-9]{2}[a-zA-Z]',unlist(strsplit(data$treatment[1],'_'))))){
      data[,2] <- gsub("X","",data[,2],perl = TRUE)
      data <- data %>%
        separate(treatment, into = c('Temperature',
                                   'Replicate',
                                   'treatment'),sep = '_')
      data_clean <- data %>%
        dplyr::group_by(id,Temperature, treatment) %>%
        dplyr::reframe(Mean = mean(value,na.rm=T),
                SD = sd(value,na.rm=T),
                SEM = SD/sqrt(length(na.omit(value))))
      }
    else {
      identifiers <- data.frame(str_split(data$treatment,
                                        '\\_',
                                        simplify=T))

      identifiers[,2] <- gsub("(...)(?!$)", "\\1_", identifiers[,2], perl=TRUE)
      identifiers <- separate(data = identifiers,
                            col = 2,
                            into = c('temp','repl'),
                            sep = '\\_')

      data[,2] <- unite(identifiers,'treatment',sep = '_',remove=T)
      data <- data %>%
        separate(treatment, into = c('set',
                                   'Temperature',
                                   'Replicate',
                                   'treatment'),sep = '_')

      data_clean <- data %>%
        dplyr::group_by(id,set,Temperature, treatment) %>%
        dplyr::reframe(Mean = mean(value,na.rm=T),
                SD = sd(value,na.rm=T),
                SEM = SD/sqrt(length(na.omit(value))))
      }
    }
  else if (length(unlist(strsplit(data$treatment[1],'_'))) < 3  ) {
    print('Make sure the format/label is correct.')
    }


  # Handling user inputs string of treatments/drugs and control group
  veh_loc <- (data_clean$Mean==0)
  default_expstr <- na.omit(c(unique(data_clean$treatment[!veh_loc]),unique(data_clean$treatment[veh_loc])))

  # Placing vehicle string in the last position
  if(any(length(expstr))){
    for (i in seq_along(expstr)){
      data_clean$treatment <- gsub(default_expstr[i],
                                   expstr[i],
                                   data_clean$treatment,
                                   perl = TRUE)
      }
    }
  else {
    expstr <- default_expstr
    }



  ##################################################
  #                                                #
  #   Data Wrangling: 'Rule' based categorization  #
  #                                                #
  ##################################################

  ### Categorization according to expression/stability change  ###


  ######### Partially changed since SEM is used more frequently #########
  ## Impose metric based conditions


  if (dispmeas == 'SD') {
    selection_metrics <- data_clean %>%
      dplyr::group_by(mean_threshold = (abs(Mean) > meancutoff),
                      bounded = (abs(Mean)- boundedness*SD>0),
                      wellmeasured = (SD < qualitycutoff))
    }
  else {
    selection_metrics <- data_clean %>%
      dplyr::group_by(mean_threshold = (abs(Mean) > meancutoff),
                      bounded = (abs(Mean)- boundedness*SEM>0),
                      wellmeasured = (SEM < qualitycutoff))
    }
  ########################################################################


  ## Initial hitlist (Reference)

  # Definition for hit
  hits_definition <- selection_metrics %>%
    dplyr::filter(mean_threshold == T, bounded == T)

  #Unique id,cond pair for hits
  keys_hits <- hits_definition %>% dplyr::ungroup() %>%  dplyr::select(id,treatment) %>%  dplyr::distinct()

  # Reference hitlist
  hitlist <- selection_metrics %>% dplyr::right_join(keys_hits, by = c('id','treatment'))
  referencelist <- hitlist



  # Separate hits and NN from data
  NN <- selection_metrics %>%
    dplyr::anti_join(hitlist, by = c('id','treatment')) %>%
    dplyr::filter(treatment != expstr[length(expstr)]) %>%
    dplyr::group_by(category = 'NN')


  ## Extract Not determinable -  ND

  # Proteins with noisy 37C measurement
  ND_condition <- hitlist %>%  dplyr::filter(Temperature == '37C',
                                    bounded == FALSE,
                                    wellmeasured == FALSE)

  keys_ND <- ND_condition %>% dplyr::ungroup() %>% dplyr::select(id,treatment) %>%  dplyr::distinct()

  ND <- hitlist %>%  dplyr::right_join(keys_ND,by = c('id','treatment'))

  # Proteins without 37C measurements
  NDwo37 <- hitlist %>% dplyr::filter(Temperature == '37C',
                             is.na(Mean) == T)

  # Combine all of ND
  ND <- ND %>%
    dplyr::full_join(NDwo37, by = names(ND)) %>%
    dplyr::group_by(category = 'ND')

  # Remove ND from hitlist
  hitlist <- hitlist %>% dplyr::anti_join(ND, by = c('id','treatment'))


  ## Stability change w/out expression change  - NC

  # High Temp. hits
  highT <- hitlist %>% dplyr::filter(Temperature != '37C',
                             mean_threshold == TRUE,
                             bounded == TRUE)
  keys_highT <- highT %>% dplyr::ungroup() %>%  dplyr::select(id, treatment) %>%  dplyr::distinct()

  # Low Temp., small mean and well measured
  lowT <- hitlist %>% dplyr::filter(Temperature == '37C',
                            mean_threshold == FALSE,
                            wellmeasured == TRUE)
  keys_lowT <- lowT %>%  dplyr::ungroup() %>%  dplyr::select(id,treatment) %>%  dplyr::distinct()

  # Match id,treatment pair fulfilling NC condition
  NC_keys <- dplyr::inner_join(keys_highT,keys_lowT, by = c('id','treatment'))


  NC <- dplyr::right_join(hitlist,NC_keys,by = c('id','treatment')) %>%
    dplyr::group_by(category = 'NC')

  # Remove NC from hitlist
  hitlist <- hitlist %>%  dplyr::anti_join(NC, by = c('id','treatment'))

  if(nrow(hitlist) > 0){

    ## Expression change w/out stability change - CN

    ### INTERVAL CONDITION FOR CN
    CN_intervals <- hitlist %>%  dplyr::group_by(upperbound = Mean + SEM, lowerbound = Mean - SEM)
    CN_keys <- CN_intervals %>%  dplyr::ungroup() %>%  dplyr::select(id, treatment) %>%  dplyr::distinct()

    CN_keep_interval <- c()
    CN_discard_interval <- c() # Will be part of CC

    for (idx in 1:nrow(CN_keys)){
      iterTibble <- CN_intervals %>%
        dplyr::filter(id == CN_keys$id[idx],treatment == CN_keys$treatment[idx]) %>%
        dplyr::select(upperbound,lowerbound)
      overlap <- c()
      for (i in 1:nrow(iterTibble)) {
        refinterval <- iterTibble[i,]
        testinterval <- iterTibble[-i,]
        overlap <- c(overlap,sum(refinterval$lowerbound <= testinterval$upperbound & refinterval$upperbound >= testinterval$lowerbound)>0)
      }
      if(all(overlap, na.rm = TRUE)){
        CN_keep_interval <- rbind(CN_keep_interval,CN_keys[idx,])
      }
      else {
        CN_discard_interval <- rbind(CN_discard_interval,CN_keys[idx,])
      }
    }


    ### DEVIATION FROM REFERENCE MEAN treatment

    # 37C with high mean
    CN_referencekeys <- hitlist %>%
      dplyr::filter(Temperature =='37C', mean_threshold == TRUE) %>%
      dplyr::ungroup () %>%
      dplyr::select(id,treatment) %>%
      dplyr::distinct()

    CN_keep_meandev <- c()
    CN_discard_meandev <- c() # Will be part of CC

    for (idx in 1:nrow(CN_referencekeys)) {
      iterTibble <- hitlist %>% dplyr::filter(id == CN_referencekeys$id[idx],
                                       treatment == CN_referencekeys$treatment[idx])
      if (nrow(iterTibble %>% dplyr::filter(Temperature !='37C',bounded == TRUE)) == 0){
        next # Not bounded high temps are disregarded
      } else {
        refmean <- (iterTibble %>% dplyr::filter(Temperature == '37C'))$Mean
        boundediterTibble <- iterTibble %>% dplyr::filter(Temperature != '37C', bounded == TRUE)
        deviateTibble <- iterTibble %>% dplyr::filter(Temperature != '37C', wellmeasured == TRUE)
        # Bounded high temp. deviate at most meandev AND
        # wellmeasured high temp don't deviate more than meandev
        if ( (sum((abs(boundediterTibble$Mean - refmean)< meandev*abs(refmean))) == nrow(boundediterTibble)) &&
             (!any(abs(deviateTibble$Mean - refmean)> meandev*abs(refmean)))  ) {
          CN_keep_meandev <- rbind(CN_keep_meandev,CN_referencekeys[idx,])
        }
        else{
          CN_discard_meandev <- rbind(CN_discard_meandev,CN_referencekeys[idx,])
        }
      }
    }



    # Proteins that do not fulfill either of CN conditions are CC - Expression and stability change
    CCkeys <- dplyr::inner_join(CN_discard_interval,CN_discard_meandev, by = c('id','treatment'))
    CC <- dplyr::right_join(hitlist,CCkeys,by = c('id','treatment')) %>%
      dplyr::group_by(category = 'CC')


    # Remove all CC and reveal CNs
    CN <- dplyr::anti_join(hitlist,CC, by = c('id','treatment')) %>%
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
    filename_reference <- paste(treatment,
                                format(Sys.time(), "%H%M_%e-%m-%y"),
                                'Hitlist',
                                sep = '_')
    filename_categorized <- paste(treatment,
                                  format(Sys.time(), "%H%M_%e-%m-%y"),
                                  'Categorized',
                                  sep = '_')
    filename_NN <- paste(treatment,
                         format(Sys.time(), "%H%M_%e-%m-%y"),
                         'NN',
                         sep = '_')
    filename_summary <- paste(treatment,
                              format(Sys.time(), "%H%M_%e-%m-%y"),
                              'Summary',
                              sep = '_')
    # Output csv
    write.csv(out_reference, file = str_c(filename_reference,'.csv'))
    write.csv(out_categorized, file = str_c(filename_categorized,'.csv'))
    write.csv(out_summary, file = str_c(filename_summary,'.csv'))
    write.csv(out_NN, file = str_c(filename_NN,'.csv'))
    }


  ### Summary: Return to user

  ### Return : Summary of run and dataframe with categories for further handling. ###
  summary_categorized <- categorized_full %>% dplyr::group_by(id,treatment,category) %>%  dplyr::reframe()
  summary_NN <- NN %>% dplyr::group_by(id,treatment,category) %>%  dplyr::reframe()
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
