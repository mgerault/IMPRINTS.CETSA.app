#' imprints_network
#'
#' Function to plot an interactive STRING network with the IMPRINTS barplots in each nodes
#'
#' @param data Dataset after imprints_caldiff
#' @param hits The Uniprots ids from which you want to plot the network.
#'   If NULL, it will take all the proteins from the data.
#' @param treatment The treatment from which you want to plot to plot the barplots.
#'   If NULL, it will plot all the treatments present in the data.
#' @param GOterm Either a data frame containing the columns 'id', 'GOterm' and eventually 'color';
#'   a character corresponding to a data base or a NULL object.
#'   If you choose a data base, it will perform an enrichment analysis and assign a GO term
#'   with the lowest FDR to each protein.
#'   The data bases now supported are: COMPARTMENTS, Process, Component, Function, TISSUES,
#'   Keyword, KEGG, SMART, PMID, RCTM, WikiPathways and NetworkNeighborAL
#' @param required_score A numeric between 0 and 1000 to set the required score
#'   for the interectation in the STRING network.
#' @param witherrorbar Logical to tell if you want to print the error bar or not on the bar plots.
#' @param FC_border Logical to tell if you want to color the nodes borders according to the
#'   maximum Fold-Change from each protein.
#' @param colorbar A vector of colors corresponding to each treatment.
#' @param colorFC If FC_border is set to TRUE, a vector of three color corresponding
#'   to the min, mid and max value from the maximum fold change.
#' @param species The species; either human, mouse or rat.
#' @param physics_type A string corresponding to the solver of the network; default is forceAtlas2Based.
#'   Possible options: 'barnesHut', 'repulsion', 'hierarchicalRepulsion', 'forceAtlas2Based'.
#' @param physics_enabled Logical to toggle the physics system on or off.
#' @param ... Options for network visualization. Options available are:
#'   borderWidth, imagePadding, node_size, edge_length, edge_color
#'   font_background, font_color, legend_width, file_type
#'
#'
#' @return visNetwork object
#'
#' @export

imprints_network <- function(data, hits = NULL, treatment = NULL, GOterm = NULL,
                             required_score = 400, witherrorbar = TRUE, FC_border = TRUE,
                             colorbar = NULL, colorFC = c("#0041FF", "#FFFFFF", "#FF0000"),
                             species = c("human", "mouse", "rat"),
                             physics_type = c("forceAtlas2Based", "barnesHut",
                                              "repulsion", "hierarchicalRepulsion"),
                             physics_enabled = TRUE, ...){
  # load if necessary packages
  if(!("STRINGdb" %in% names(sessionInfo()$otherPkgs))){
    library(STRINGdb)
  }
  if(!("visNetwork" %in% names(sessionInfo()$otherPkgs))){
    library(visNetwork)
  }

  # define options for network
  opts <- list(borderWidth = 5, imagePadding = 8,
               node_size = 40, edge_length = 60,
               edge_color = "#00000075", legend_width = 0.2,
               font_background = "#A6A6A669", font_color = "#343434",
               file_type = "png")
  dot_opts <- match.call(expand.dots = FALSE)$...
  if(!is.null(dot_opts)){
    if(is.null(names(dot_opts))){
      warning("The options parameters should be named")
    }
    else{
      opts_to_change <- names(opts)[names(opts) %in% names(dot_opts)]
      if(length(opts_to_change)){
        opts[opts_to_change] <- dot_opts[opts_to_change]
      }
      else{
        warning(paste("The options parameters can only be", paste(names(opts), collapse = ", ")))
      }
    }
  }

  if(is.null(treatment)){
    treatment <- get_treat_level(data)
  }
  else{
    treatment_in <- treatment %in% get_treat_level(data)
    if(!all(treatment_in)){
      message(paste("Error: ", treatment[!treatment_in], "wasn't found in your data."))
      return()
    }
  }

  if(!is.null(colorbar)){
    if(length(colorbar) != length(treatment)){
      message("Error: The number of colors you gave doesn't match the number of conditons.")
      return()
    }
  }

  if(length(colorFC) != 3){
    message("Error: You need to provide 3 colors for colorFC.")
    return()
  }

  species <- match.arg(species)
  if(species == "human"){
    species <- 9606
  }
  else if(species == "mouse"){
    species <- 10090
  }
  else if(species == "rat"){
    species <- 10116
  }

  if(!is.null(hits)){
    hits_in <- hits %in% data$id

    if(all(!hits_in)){
      message("Error: None of the proteins you selected are in your data")
      return()
    }
    if(any(!hits_in)){
      message(paste(paste(hits[!hits_in], collapse = ", "), "were not found in your data and hence, have been removed."))
      hits <- hits[hits_in]
    }

    data <- data[which(!is.na(match(data$id, hits))),]
  }

  # extracting Gene info
  data_genes <- data[,c("id", "description")] %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(description = stringr::word(stringr::str_extract(description, "(?<=GN=).*(?=)"), 1))
  colnames(data_genes)[2] <- "Gene"
  data_genes <- as.data.frame(data_genes) # STRINGdb cannot handle tibble

  run_enrich <- FALSE
  if(inherits(GOterm, "data.frame")){
    good_col <- c("id", "GOterm") %in% colnames(GOterm)
    if(any(!good_col)){
      message(paste("Error:", c("id", "GOterm")[!good_col], "wasn't found in your GOterm data.frame"))
      return()
    }
    if(!("color" %in% colnames(GOterm))){
      GOterm$color <- scales::col_factor("Spectral", domain = NULL)(GOterm$GOterm)
    }
  }
  else if(inherits(GOterm, "character")){
    if(length(GOterm) > 1){
      GOterm <- GOterm[1]
      message(paste("Only", GOterm, "has been kept; you should provide only one data base."))
    }
    GOterm_in <- GOterm %in% c("COMPARTMENTS", "Process", "Component",
                               "Function", "TISSUES", "Keyword",
                               "KEGG", "SMART", "PMID",
                               "RCTM", "WikiPathways", "NetworkNeighborAL")
    if(GOterm_in){
      run_enrich <- TRUE
    }
    else{
      message(paste(GO_term, "is not in the database available yet. Please, choose one from\n
                    COMPARTMENTS, Process, Component, Function, TISSUES, Keyword, KEGG,\n
                    SMART, PMID, RCTM, WikiPathways or NetworkNeighborAL"))
      return()
    }
  }

  message("Start to map genes...")
  if(!file.exists("STRING_data")){
    dir.create("STRING_data")
  }
  string_db <- STRINGdb$new(version="11.5", species=species,               #ID 9606 correspond to human
                            score_threshold=200,
                            input_directory=file.path(getwd(), "STRING_data"))
  string_id <- string_db$map(data_genes, "id", removeUnmappedRows = TRUE)

  message("Get STRING network...")
  interact <- string_db$get_interactions(string_id$STRING_id)
  if(nrow(interact) == 0){
    message("Warning: No interactions were found between the proteins you selected !")
    return(NULL)
  }
  interact$from <- sapply(interact$from, function(x) string_id$id[which(string_id$STRING_id == x)][1])
  interact$to <- sapply(interact$to, function(x) string_id$id[which(string_id$STRING_id == x)][1])
  interact <- interact[-which(duplicated(interact)),] # remove potential duplicated rows
  if(length(which(interact$combined_score >= required_score)) == 0){
    message("Warning: No interactions passed the required interaction score ! Try to decrease it")
    return(NULL)
  }

  message("Generates barplots...")
  if(is.null(colorbar)){
    bar <- imprints_barplotting_app(data, treatmentlevel = treatment,
                                    witherrorbar = witherrorbar,
                                    printBothName = FALSE)
  }
  else{
    bar <- imprints_barplotting_app(data, treatmentlevel = treatment,
                                    colorpanel = colorbar, witherrorbar = witherrorbar,
                                    printBothName = FALSE)
  }

  bar <- lapply(bar, function(g){
    g <- g +
      theme(legend.position = "none",
            plot.title = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_line(color = NA),
            panel.background = element_rect(fill = "transparent", color = NA),
            plot.background = element_rect(fill = "transparent", color = NA)
      );
    list(g)
  })
  names(bar) <- data$id

  # preparing legend if needed
  lnodes <- data.frame(matrix(nrow = 0, ncol = 5))
  colnames(lnodes) <- c("label", "color.border", "color.background", "shape", "borderWidth")
  if(FC_border){
    to_rm <- get_treat_level(data)[!(get_treat_level(data) %in% treatment)]
    if(length(to_rm)){
      data <- data[,-stringr::str_which(colnames(data), paste0("_", to_rm, "$", collapse = "|"))]
    }

    if(length(treatment) > 1){
      message("More than one treatment has been selected. The biggest maximum FC in absolute value will be taken.")
    }

    FC <- data %>%
      dplyr::select(-countNum, -sumUniPeps, -sumPSMs) %>%
      tidyr::gather("treatment", "value", -id, -description) %>%
      tidyr::separate(treatment, into = c("temperature", "biorep", "treatment"), sep = "_") %>%
      dplyr::group_by(id, description, treatment) %>%
      dplyr::reframe(FC = ifelse(max(value, na.rm = TRUE) < abs(min(value, na.rm = TRUE)),
                                 min(value, na.rm = TRUE), max(value, na.rm = TRUE))) %>%
      dplyr::group_by(id, description) %>%
      dplyr::reframe(FC = FC[which.max(abs(FC))]) %>%
      dplyr::mutate(description = stringr::word(stringr::str_extract(description, "(?<=GN=).*(?=)"), 1)) %>%
      dplyr::rename(name = id)

    FC$description <- NULL

    rng_FC <- round(range(FC$FC, na.rm = TRUE), 1)
    rng_FC <- append(rng_FC, 0, 1)

    FC$colorFC <- scales::gradient_n_pal(colorFC, rng_FC, "Lab")(FC$FC)
    if(length(which(is.na(FC$colorFC)))){
      FC$colorFC[which(is.na(FC$colorFC))] <- "#666666"
    }

    # updating legend
    lnodes <- as.data.frame(rbind(lnodes,
                                  data.frame(label = paste(rng_FC),
                                             color.border = colorFC,
                                             color.background = "#CACACA",
                                             shape = "circle",
                                             borderWidth = opts$borderWidth))
                            )

  }

  if(run_enrich){
    message("Start enrichment analysis...")
    enrich <- string_db$get_enrichment(string_id$STRING_id, category = GOterm)

    if(nrow(enrich) == 0){
      message("Warning: No enrichment were found !")
      enrich <- NULL
    }
    else{
      enrich <- enrich[,c("preferredNames", "fdr", "description")] %>%
        dplyr::group_by(description, fdr) %>%
        dplyr::reframe(Gene = stringr::str_split(preferredNames, ",")[[1]])
      colnames(enrich)[1] <- "GOterm"

      enrich <- enrich %>%
        dplyr::ungroup() %>% dplyr::group_by(Gene) %>%
        dplyr::group_modify(~ {
          if(nrow(.x) > 1){
            y <- .x[which.min(.x$fdr),]
            return(y)
          }
          else{
            return(.x)
          }
        }) %>%
        dplyr::select(-fdr)

      enrich <- enrich %>%
        dplyr::right_join(data_genes, by = "Gene") %>%
        dplyr::rename(name = id)
      enrich$Gene <- NULL

      enrich$GOterm[which(is.na(enrich$GOterm))] <- "none"
      enrich$color <- scales::col_factor("Spectral", domain = NULL)(enrich$GOterm)
      enrich$color[which(enrich$GOterm == "none")] <- "#CACACA"
    }
  }
  else if(!is.null(GOterm)){
    enrich <- GOterm  %>%
      dplyr::right_join(data_genes, by = "id") %>%
      dplyr::rename(name = id)
    enrich$Gene <- NULL

    enrich$GOterm[which(is.na(enrich$GOterm))] <- "none"
    enrich$color[which(enrich$GOterm == "none")] <- "#CACACA"
  }
  else{
    enrich <- NULL
  }

  ### creating graph
  graph <- interact %>%
    dplyr::filter(combined_score >= required_score) %>%
    dplyr::mutate(combined_score = combined_score/1000)
  colnames(graph)[3] <- "value"

  edges <- data.frame(graph)

  # adding info
  data_genes <- data_genes %>%
    tibble::column_to_rownames(var = "id")
  graph <-  graph %>% dplyr::select(from, to) %>%
    tidyr::gather("key", "name") %>%
    dplyr::select(name) %>%
    unique() %>%
    dplyr::group_by(name) %>%
    dplyr::mutate(bar = bar[[name]],
                  Gene = data_genes[name,"Gene"]
                  )
  if(FC_border){
    graph <- graph %>%
      dplyr::left_join(FC, by = "name")
  }
  if(!is.null(enrich)){
    graph <- graph %>%
      dplyr::left_join(enrich, by = "name")

    # updating legend
    enrich_leg <- unique(graph[,c("GOterm", "color")])
    lnodes <- as.data.frame(rbind(lnodes,
                                  data.frame(label = enrich_leg$GOterm,
                                             color.border = enrich_leg$color,
                                             color.background = enrich_leg$color,
                                             shape = "ellipse",
                                             borderWidth = opts$borderWidth))
    )
  }

  # generates img for each barplot
  f_rnd <- paste0("img_", paste(sample(c(1:9, letters), 12, TRUE), collapse = ""))
  dir.create(f_rnd)
  lapply(graph$bar, function(img){
    n <- img$data$id[1]
    n <- stringr::word(n, 1, sep = "\n")
    ggsave(plot = img, paste0(f_rnd, "/", n , ".png"),
           width = 1, height = 1)
  })
  # convert it in base 64 so js can read them
  bar_base64 <- sapply(graph$name, function(x){
    f <- paste0(f_rnd, "/", x, ".png")
    res <- RCurl::base64Encode(readBin(f, "raw", file.info(f)[1, "size"]), "txt")
    res <- paste0("data:image/png;base64,", res[1])
  })
  unlink(f_rnd, recursive = TRUE) # remove directory created


  nodes <- data.frame(id = graph$name,
                      label = graph$Gene,
                      shape = "circularImage",
                      image = bar_base64
                      )
  if(FC_border){
    nodes$color.border <- graph$colorFC
  }
  if(!is.null(enrich)){
    nodes$group <- graph$GOterm
    nodes$color.background <- graph$color
  }

  physics_type <- match.arg(physics_type)
  net <- visNetwork(nodes, edges, width = "100%", height = "800px") %>%
    visNodes(borderWidth = opts$borderWidth, shapeProperties = list(useBorderWithImage = TRUE),
             imagePadding = opts$imagePadding,
             font = list(background = opts$font_background, color = opts$font_color),
             size = opts$node_size) %>%
    visEdges(color = opts$edge_color, length = opts$edge_length) %>%
    visPhysics(solver = physics_type, enabled = physics_enabled) %>%
    visInteraction(navigationButtons = TRUE) %>%
    visEvents(selectNode = "function(properties) {
                             if(properties.event.srcEvent.shiftKey){
                                var win = window.open('https://www.uniprot.org/uniprot/' + properties.nodes);
                                win.focus();
                              }
                            }") %>%
    visExport(type = opts$file_type, name = paste0(format(Sys.time(), "%y%m%d_%H%M_"), "network"),
              float = "left", label = paste("Save network as", opts$file_type),
              background = "#FFFFFF00", style= "") %>%
    visOptions(manipulation = TRUE)

  if(nrow(lnodes) != 0){
    net <- net %>%
      visLegend(addNodes = lnodes, stepY = 50,
                width = opts$legend_width, useGroups = FALSE)
  }

  message("Done !")

  return(net)
}

