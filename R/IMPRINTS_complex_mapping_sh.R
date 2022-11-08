#' IMPRINTS_complex_mapping_sh
#'
#' Function to generate 2D bar plot and pdf file with multipanel bar plots for 2D-CETSA data of
#' proteins which have similar profile from a selected protein.
#' It is totally based on the function ms_2D_barplotting from the mineCETSA package.
#'
#' @param data dataset after calculating the relative protein abundance differences and
#'             deriving the average of the profile, ie IMPRINTS_average_sh,
#'             make sure the columns with readings are named in the format like "37C_Treatment".
#' @param categorytable Dataset containing the category of each proteins from data, ie hitlist.
#'                      The column that contains the conditions must be named treatment.
#' @param set A single character to indicate the sample name if any
#' @param treatment A single character to indicate the sample name
#' @param complexdatabase By default is core Corum database, 03.09.2018 Corum 3.0 current release
#' @param organism Choose from Human/Mouse/Rat, default is Human
#' @param complexID A vector of complexID to retrieve from the data, if specified
#' @param targetcategory A character vector to select the category you want. Default value is c("NC","CN","CC")
#' @param minsubunitsIdentified The minimal number of subunits should be identified, default value is 3
#' @param cutisoformbycount When there are several isoforms mapped in a complex,
#'                          by default to filter away the isoforms with comparatively fewer count number
#'
#' @return A dataframe with the assigned protein complex
#'
#' @export
#'

IMPRINTS_complex_mapping_sh <- function (data, categorytable = NULL, set = NULL, treatment = NULL,
                                      complexdatabase = "Corum", organism = "Human", complexID = NULL,
                                      targetcategory = c("NC", "CN", "CC"), minsubunitsIdentified = 3,
                                      cutisoformbycount = TRUE)
{
  if (complexdatabase != "Corum") {
    stop("Only core Corum database is supported for now")
  }
  if (!organism %in% c("Human", "Mouse", "Rat")) {
    stop("Provide the right organism")
  }
  complexdb <- subset(mineCETSA::coreCorum, Organism == organism)
  if (length(complexID)) {
    if (class(complexID) == "character") {
      complexID <- as.numeric(complexID)
    }
    complexdb <- subset(complexdb, ComplexID %in% complexID)
    if (nrow(complexdb) == 0) {
      stop("Make sure the provided complexID is present in core Corum database")
    }
  }

  dataname <- deparse(substitute(data))
  outdir <- mineCETSA:::ms_directory(data, dataname)$outdir
  data <- mineCETSA:::ms_directory(data, dataname)$data

  if(length(grep("Condition", names(categorytable)))){
    names(categorytable)[grep("Condition", names(categorytable))] <- "treatment"
  }

  if (length(treatment) == 1) {
    treat = treatment
    message(paste0("To check the data under the treatment of ", treat))
    data <- subset(data, treatment == treat)
    sel <- grep(treatment, names(data))
    if (length(set) == 1) {
      message(paste0("To check the data within the set of ", set))
      sel1 <- grep(set, names(data))
      sel <- dplyr::intersect(sel, sel1)
    }
    if (length(sel)) {
      data <- data[, c(1, 2, sel, (ncol(data) - 2):ncol(data))] #take the id descr, next treat, next sumUni, etc.
    }
  }
  else {
    stop("Need to specify a single character to indicate the condition name")
  }
  if (length(categorytable) > 0 & length(grep("category", names(categorytable)))) {
    if (length(targetcategory) > 0) {
      categorytable <- subset(categorytable, category %in%  targetcategory)             #filter the prot with catego selected
    }
    categorytable <- subset(categorytable, treatment == treat)
    if (length(set) & length(grep("category", names(categorytable)))) {
      setname = set
      categorytable <- subset(categorytable, set == setname)
    }
    data <- merge(data, categorytable[, c("id", "category")], by = "id")
  }
  if (nrow(data) != 0) {
    message(paste0(nrow(data), " protein entries to be assigned to Corum database..."))
  }
  else {
    stop("Make sure the right categories were specified.")
  }
  comps <- data.frame(ComplexID = character(), ComplexName = character(),
                      subunitsUniprot_IDs = character(), subunitsNum = integer())
  data2 <- data[F, ]
  for (h in 1:nrow(complexdb)) {
    subunitNames <- complexdb[h, "subunitsUniProt_IDs"]
    subunitname <- unique(unlist(str_split(as.character(subunitNames),";")))
    for (s in subunitname) {
      pos <- grep(s, data$id, value = FALSE)
      if (length(pos) > 0) {
        data2 <- rbind(data2, data[pos, ])
        for (i in 1:length(pos)) {
          comps <- rbind(comps, complexdb[h, c(1, 2, 6, 21)])
        }
      }
    }
  }
  if (nrow(data2) == 0) {
    tab_name <- c("ComplexName", "subunitsNum", "subunitsIdentifiedNum",
                  "id", "description", "gene", "category")
    data2 <- as.data.frame(matrix(nrow = 0, ncol = length(tab_name)))
    colnames(data2) <- tab_name
    return(data2)
  }
  data2 <- data2 %>% dplyr::rowwise() %>% dplyr::mutate(gene = mineCETSAapp:::getGeneName(description))
  comps3 <- cbind(comps, data2)
  comp_table <- comps3 %>% dplyr::add_count(ComplexID, gene)
  cleariso <- FALSE
  if (sum(comp_table$n > 1)) {
    message("Be cautious about the possible isoforms in the complex mapping...")
    message("Double check the potential_isoform_in_complex table in the working directory.")
    comp_table <- subset(comp_table, n > 1)
    comp_table$n <- NULL
    ms_filewrite(comp_table, paste0(dataname, "_potential_isoform_in_complex.txt"),
                 outdir = outdir, withdescription = F)
    cleariso <- TRUE
  }
  if (cutisoformbycount) {
    if (cleariso) {
      message("Would remove the isoforms with comparatively fewer count number.")
    }

    names(comps3)[str_which(names(comps3), "^countNum")] <- "countNum"
    names(data)[str_which(names(data), "^countNum")] <- "countNum"
    comps3 <- comps3 %>% dplyr::left_join(data[, c("id", "countNum")]) %>%
      dplyr::group_by(ComplexID, gene) %>% dplyr::top_n(1, countNum) %>%
      dplyr::ungroup() %>% dplyr::group_by(ComplexID) %>% dplyr::mutate(subunitsIdentifiedNum = n()) %>%
      dplyr::rowwise() %>% dplyr::mutate(subunitsIdentifiedPerc = subunitsIdentifiedNum/subunitsNum) %>%
      dplyr::filter(subunitsIdentifiedNum >= minsubunitsIdentified) %>%
      dplyr::ungroup()
  }
  message(paste0("The number of identified multi-protein complexes with at least ",
                 minsubunitsIdentified, " subunits is: ", length(unique(comps3$ComplexID))))
  return(comps3)
}



