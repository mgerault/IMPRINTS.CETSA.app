library(IMPRINTS.CETSA)
library(stringr)
library(rio)
library(DT)
library(cowplot)
library(plotly)

library(pubmedR)
library(bibliometrix)
library(officer)
library(magrittr)
library(STRINGdb)
library(visNetwork)
library(clusterProfiler)

library(shiny)
library(shinyjs)
library(shinyMatrix)
library(shinydashboard)
library(shinycssloaders)
library(colourpicker)

#increase the max request size for uploading files
options(shiny.maxRequestSize = 1000*1024^2)
#set options for the spinner when things are loading
options(spinner.color = "#518CE2", spinner.color.background = "000000", spinner.size = 2)

# javascript code to collapse box
jscode <- "
shinyjs.collapse = function(boxid) {
$('#' + boxid).closest('.box').find('[data-widget=collapse]').click();
}
"

ui <-  navbarPage(title = img(src="logo.png", height = "28px"),
                 id = "navBar",
                 theme = "paper.css", # file in www
                 collapsible = TRUE, # usefull when viewing on smaller screen
                 inverse = FALSE, # true: use a dark background and light text for the navigation bar
                 windowTitle = "IMPRINTS.CETSA.app", # just name of tab
                 position = "fixed-top",
                 footer = includeHTML("./www/include_footer.html"), # bottom of the page/site
                 header = tagList(
                   shinyWidgets::useShinydashboard(),      # allow to render the boxes from shinydashboard
                   tags$style("body {padding-top: 75px;}")
                 ),

                 tabPanel("Home", value = "home",
                          tags$head(tags$script(HTML('
                                                       var fakeClick = function(tabName) {
                                                       var dropdownList = document.getElementsByTagName("a");
                                                       for (var i = 0; i < dropdownList.length; i++) {
                                                       var link = dropdownList[i];
                                                       if(link.getAttribute("data-value") == tabName) {
                                                       link.click();
                                                       };
                                                       }
                                                       };
                                                       '))
                                    ),

                          fluidRow(style = "height:20px;"),
                          fluidRow(
                            column(12,
                                   shiny::HTML("<h1>About CETSA</h1><br>"),
                                   shiny::HTML("<h5>In 2013, we published the first paper describing the transformative Cellular Thermal Shift Assay (CETSA, Martinez Molina et al, 2013, Science).
                                                    CETSA is the first widely applicable method for assessing direct drug binding in cells and tissues.
                                                    This method has had a major impact on drug discovery and is broadly applied in academia and
                                                    industry to establish the correct ligand-protein relationship, improve the efficacy and quality of drug
                                                    candidates, and has the potential to serve as an important clinical diagnostic of drug efficacy in the future.
                                                    Recently we have published papers demonstrating that a new generation of multi-dimensional CETSA
                                                    provides a highly resolved means to study the interactions of proteins with other cellular components
                                                    (including metabolites, nucleic acids and other proteins) in intact cells and tissues at the proteome level.
                                                    This approach provides a completely new perspective on how cellular processes are executed and we expect
                                                    it to have a large impact on understanding disease processes and drug action in the future.
                                                    It is also a new way to discover new functional protein targets and biomarkers for intervention and therapy.

                                                    <br><br> To learn more about CETSA and our lab, click <a href='https://www.cetsa.org/'>here</a></h5>"
                                               )
                                   )
                            ),

                          tags$hr(),

                          fluidRow(
                            column(12,
                                   shiny::HTML("<h1>IMPRINTS.CETSA.app</h1><br>"),
                                   shiny::HTML("<h5>IMPRINTS.CETSA.app is an R package that includes a shiny app that can help you
                                                    easily analyze your IMPRINTS-CETSA data. In this app, you will be able
                                                    to process the TMT multiplexing-based quantification files, narrow down protein targets,
                                                    visualize the results in bar plot or heat map formats, run gene ontology enrichment analysis and more. <br>
                                                    The app also includes 2 sample datasets comparing different phases of the cell cycle
                                                    originally published in Dai et al. 2018, Cell, for easy visualization and comparison.
                                                    You can add new datasets at your ease from your local machine. <br>
                                                    To learn how to use the app, you are encouraged to watch the tutorial video by clicking
                                                    on the question mark icon in the upper right corner.
                                                    <br><br><br>
                                                    <center> Otherwise, start your analysis here ! </center></h5>")
                                   )
                            ),
                          fluidRow(
                            column(3),
                            column(6,
                                   tags$div(align = "center",
                                            tags$a("Start",
                                                   onclick="fakeClick('proteins')", # take you to other tab with value 'analysis'
                                                   class="btn btn-primary btn-lg")
                                   )
                            ),
                            column(3)
                          ),

                          tags$hr()
                          ),

                 navbarMenu("Analysis",
                            tabPanel("Peptides", value = "peptides",
                                     shinyjs::useShinyjs(),
                                     tags$style(type = 'text/css',
                                                '#modal1_pep .modal-dialog { width: fit-content !important; overflow-x: initial !important}
                                                 #modal1_pep .modal-body { width: 150vh; overflow-x: auto;}'
                                                ),
                                     tags$style(type = 'text/css',
                                                '#modal2_pep .modal-dialog { width: fit-content !important; overflow-x: initial !important}
                                                 #modal2_pep .modal-body { width: 150vh; overflow-x: auto;}'
                                                ),
                                     tags$style(type = 'text/css',
                                                '#modal3_pep .modal-dialog { width: fit-content !important; overflow-x: initial !important}
                                                 #modal3_pep .modal-body { width: 150vh; overflow-x: auto;}'
                                                ),

                                     fluidRow(style = "height:20px;"),
                                     fluidRow(column(1),
                                              column(8, checkboxInput("step_peptides", tags$strong("Did you already performed normalization ?"), FALSE, width = "100%"))
                                              ),

                                     shinyjs:::useShinyjs(),
                                     shinyjs:::extendShinyjs(text = jscode, functions = "collapse"),
                                     fluidRow(box(id="upload_peptide", title = "Upload and concatenate your data", status = "primary",
                                                  solidHeader = TRUE, collapsible = TRUE, width = 12,
                                                  tags$u(h3("Upload your data")),
                                                  checkboxInput("got_data_pep", "Do you already have the concatenate peptide file ?", FALSE),

                                                  conditionalPanel(condition = "input.got_data_pep",
                                                                   fileInput("file_data_pep", "Select your concatenate peptide file",
                                                                             accept = ".txt")
                                                                   ),
                                                  conditionalPanel(condition = "!input.got_data_pep",
                                                                   fileInput("PD_data_pep", "Select PD txt files for your analysis",
                                                                             accept = ".txt", multiple = TRUE),
                                                                   tags$hr(),

                                                                   conditionalPanel(condition = "output.pep_fileup",
                                                                                    fluidRow(column(4, shiny::HTML("<br><h5>On the table on your right, you can type the tempeature
                                                                                                                   you used for each of your files. It needs to start with a digit
                                                                                                                   and also it's better to finish by a letter like this:
                                                                                                                   '37C' or '99F'. Remember to not use '_'.
                                                                                                                   <br>Also, if you have a quantitative proteomic file, it is
                                                                                                                   advised to name its 'temperature' as '36C' for easier handle
                                                                                                                   in the other functions from IMPRINTS.CETSA.app.</h5>")),
                                                                                             column(8, uiOutput("temp_nameui_pep"))
                                                                                             ),
                                                                                    tags$hr(),
                                                                                    fluidRow(column(4, shiny::HTML("<br><h5>On the table on your right, you can type the name
                                                                                                                    of the sample to the corresponding channel. The underscore '_'
                                                                                                                    will be used as a separator between temperatures, bioreplicates and
                                                                                                                    treatments in all further functions, so make sure of your spelling.
                                                                                                                    <br><br>Here, you'll need to type first the bioreplicate and then
                                                                                                                    the treatment, like this : 'B1_Vehicle', 'B1_treatment', etc.
                                                                                                                    <br>Also, if you have a 'Mix' channel; it needs to be named explicitely
                                                                                                                    as 'Mix'.</h5>")),
                                                                                             column(8, uiOutput("treat_nameui_pep"))
                                                                                             ),
                                                                                    tags$hr(),
                                                                                    fluidRow(column(4, shiny::HTML("<br><h5>On your right, you can upload a data frame that contains
                                                                                                                   some proteins of interest, the one you want to analyze their peptides.
                                                                                                                   This data frame should contain the column 'id' and the column 'description',
                                                                                                                   the Uniprot IDs and the protein description respectively.<br>
                                                                                                                   So you can use the caldiff file output you have from your protein analysis.
                                                                                                                   <br><br>If you upload nothing, all proteins from the pepetides files
                                                                                                                   will be kept but the protein description will be missing.</h5>")),
                                                                                             column(8, fileInput("prot_data_pep", "Select a protein file",
                                                                                                                 accept = ".txt", multiple = TRUE))
                                                                                             ),
                                                                                    tags$hr(),
                                                                                    fluidRow(column(4, shiny::HTML("<br><h5>For now, you cannot choose the peptide modification
                                                                                                                   you want to keep or remove. By default, only the TMT modifications
                                                                                                                   are kept.</h5>")),
                                                                                             column(4, textInput("dname_pep", "Type a name for your dataset so the saved
                                                                                                                               file name can contain it. Can be NULL", value = ""))
                                                                                             ),
                                                                                    tags$hr(),
                                                                                    actionButton("read_pep", "Read your peptides data", class = "btn-primary")
                                                                                    )
                                                                   ),
                                                  tags$hr(),
                                                  actionButton("see1_pep", "View data uploaded")
                                                  )
                                              ),

                                     conditionalPanel(condition = "output.pep_dataup | input.step_peptides",
                                                      fluidRow(box(title = "Normalize your peptides data", status = "primary",
                                                                   solidHeader = TRUE, collapsible = TRUE, width = 12,
                                                                   tags$u(h3("Normalize your data")),
                                                                   fluidRow(column(6, checkboxInput("got_norm_pep", "Do you already have the peptide file NormPeptides ?")),
                                                                            column(6, conditionalPanel(condition = "!input.got_norm_pep",
                                                                                                       actionButton("NORM_pep", "Start Normalization", class = "btn-primary")
                                                                                                       ),
                                                                                   conditionalPanel(condition = "input.got_norm_pep",
                                                                                                    fileInput("normfile_pep", "Select the NormPeptides file", accept = ".txt")
                                                                                                    )
                                                                                   )
                                                                            ),
                                                                   tags$hr(),
                                                                   actionButton("see2_pep", "View normalized data")
                                                                   )
                                                               )
                                                      ),
                                     conditionalPanel(condition = "output.norm_pep_dataup",
                                                      fluidRow(box(title = "Compute fold change and save bar plots", status = "primary",
                                                                   solidHeader = TRUE, collapsible = TRUE, width = 12,
                                                                   tags$u(h3("Fold-change calculation")),
                                                                   checkboxInput("sequence_file", "Import a file with proteins and sequences"),
                                                                   tags$hr(),
                                                                   conditionalPanel(condition = "input.sequence_file",
                                                                                    fluidRow(column(4, shiny::HTML("<br><h5>This file needs to contains at leats one column named
                                                                                                                    'protein' and eventually another one named 'sequence'.
                                                                                                                    <br>The 'protein' column contains the Uniprot ID from the
                                                                                                                    protein you want to compute fold changes at the peptide level.
                                                                                                                    <br>The 'sequence' column contains the peptide position you want to highlight.
                                                                                                                    Every peptides before this sequence will be summed, same for the ones after.
                                                                                                                    The position needs to be in this format precisely: a number followed by a
                                                                                                                    dash and another number; like this for example '208-221'.
                                                                                                                    <br>If you import nothing, it will select all proteins from your
                                                                                                                   peptides data and will not select specific sequences. Which means
                                                                                                                   it will compute and plot fold change for all the peptides in your data.</h5>")),
                                                                                             column(8, fileInput("protseq_file_pep", "Import the file"))
                                                                                             )
                                                                                    ),
                                                                   conditionalPanel(condition = "!input.sequence_file",
                                                                                    shiny::HTML("<h5>Here, you can select the protein from which you want to compute fold
                                                                                                     changes at the peptide level.
                                                                                                     <br>The sequence you can select correspond to the peptide position you want to highlight.
                                                                                                     Every peptides before this sequence will be summed, same for the ones after.
                                                                                                     <br>The position needs to be in this format precisely: a number followed by a
                                                                                                     dash and another number; like this for example '208-221'.
                                                                                                     <br>If you select nothing, it will select all proteins from your
                                                                                                     peptides data and will not select specific sequences. Which means
                                                                                                     it will compute and plot fold change for all the peptides in your data.</h5>"),
                                                                                    fluidRow(column(6, selectizeInput("protseq_pep", "Select some proteins (if NULL, will select all)", choices = NULL, multiple = TRUE)),
                                                                                             column(6, uiOutput("selectSequenceui_pep"))
                                                                                             )
                                                                                    ),
                                                                   fluidRow(column(6, selectInput("control_pep", "Select the control from your experiment", choices = NULL)),
                                                                            column(6, textInput("dnamediff_pep", "Type a name for your dataset so the saved
                                                                                                                  file name can contain it. Can be NULL", value = ""))),
                                                                   tags$hr(),
                                                                   fluidRow(column(6, actionButton("SEQU_pep", "Compute and plot fold change", class = "btn-primary")),
                                                                            column(6, textOutput("diag_pep_sequence"))
                                                                            ),
                                                                   tags$hr(),
                                                                   fluidRow(style = "height:10px;"),

                                                                   tags$u(h3("Find potential cleaved sites")),
                                                                   fluidRow(column(6, checkboxInput("got_FCfile_pep", "Import your own fold-change file.
                                                                                                     If not, will use the one obtained in previous step.", value = FALSE)),
                                                                            conditionalPanel(condition = "input.got_FCfile_pep",
                                                                                             column(6, fileInput("FCfile_pep", "Import a peptide fold change file (txt)", accept = ".txt")))
                                                                           ),
                                                                   tags$hr(),

                                                                   conditionalPanel(condition = "output.sequence_pep_dataup",
                                                                                    fluidRow(column(3, numericInput("R2cleaved_pep", "Choose the maximum R-squared from the
                                                                                                                                      cumulative sum of the fold change.",
                                                                                                                    value = 0.9, min = 0, max = 1, step = 0.01),
                                                                                                    shiny::HTML("<br><h5>The higher it is, the less stringent your are. It means that you
                                                                                                                can accept more 'cumulative sum profile' with a linear evolution, i.e. it is
                                                                                                                less likely that there is a cleaved site.</h5>")),
                                                                                             column(3, selectInput("controlcleaved_pep", "Select the control from your experiment", choices = NULL)),
                                                                                             column(3, numericInput("propValcleaved_pep", "Choose the minimum proportion of valid values per peptide
                                                                                                                                           per treatment; i.e. if 6 temperatures and 0.5, it can't
                                                                                                                                           have more than 3 missing values.",
                                                                                                                    value = 0.4, min = 0, max = 1, step = 0.01)),
                                                                                             column(3, actionButton("CLEAVED_pep", "Search for potential cleaved site", class = "btn-primary"))
                                                                                             )
                                                                                    )
                                                                   )
                                                               )
                                                      ),
                                     fluidRow(box(title = "Join peptides datasets, remove peptides and barplots", status = "primary",
                                                  solidHeader = TRUE, collapsible = TRUE, width = 12, collapsed = TRUE,
                                                  tags$u(h3("Filter your dataset")),

                                                  shiny::HTML("<br><h5>Here, you can filter out some treatments from your peptide fold-change data
                                                              and remove some specific sequences like the cleaved sites found.
                                                              <br>For selecting the sequence you want to remove,
                                                              you can import the same file you used previously. A file that contains the protein from
                                                              which you want to remove the specific sequence (column named 'protein') and the sequence
                                                              you want to remove (column named 'sequence'). The sequence needs to be in the this format:
                                                              a number followed by a dash followed by a number, like this for example '208-221'.
                                                              <br>Once the filtration is done a txt file is saved.</h5>"),
                                                  fileInput("filter_joinpep", "Import a fold-change peptide file (txt)", accept = ".txt"),
                                                  tags$hr(),
                                                  checkboxInput("sequence_file_joinpep", "Import a file with proteins and sequences"),
                                                  tags$hr(),
                                                  conditionalPanel(condition = "input.sequence_file_joinpep",
                                                                   fluidRow(column(6, fileInput("protseq_file_joinpep", "Import the file"))
                                                                            )
                                                                   ),
                                                  conditionalPanel(condition = "!input.sequence_file_joinpep",
                                                                    fluidRow(column(6, selectizeInput("protseq_joinpep", "Select some proteins (if NULL, will select all)", choices = NULL, multiple = TRUE)),
                                                                             column(6, uiOutput("selectSequenceui_joinpep"))
                                                                             )
                                                                   ),
                                                  fluidRow(column(6, selectInput("remcond_joinpep", "Select some treatments you want to remove. If NULL, nothing is removed.",
                                                                                 choices = NULL, multiple = TRUE)),
                                                           column(6, actionButton("gofilter_joinpep", "Filter your dataset", class = "btn-primary"),
                                                                  textOutput("diag_pep_filter"))
                                                           ),

                                                  tags$hr(),
                                                  tags$u(h3("Join your datasets")),
                                                  shiny::HTML("<br><h5>Here you can import as much peptides dataset as you want and join them.
                                                              <br>This feature has been mainly made to join dataset after you checked for cleaved sites
                                                              and computed fold changes. For example if you checked for the same potential cleaved site for
                                                              one protein for several drugs and you want now to compare their effect.
                                                              <br><br>Once you joined your data, you can plot their bar plots if you want</h5>"),
                                                  tags$hr(),
                                                  fluidRow(column(6, fileInput("joinFC_file_pep", "Import the peptides files you want to join", multiple = TRUE)),
                                                           conditionalPanel(condition = "output.tojoin_pep_dataup",
                                                                            column(3, actionButton("JOIN_pep", "Start joining datasets", class = "btn-primary"),
                                                                                      textOutput("diag_pep_join")),
                                                                            column(3, actionButton("see3_pep", "View joined data"))
                                                                            )
                                                           ),

                                                  tags$hr(),
                                                  fluidRow(style = "height:10px;"),

                                                  tags$u(h3("Plot your data")),
                                                  fluidRow(column(6, radioButtons("join_data_pep", label = "", choices = c("Use the joined data" = "join_app",
                                                                                                                           "Use a file" = "join_file"))
                                                                  ),
                                                           conditionalPanel(condition = "input.join_data_pep == 'join_file'",
                                                                            column(6, fileInput("joined_file_pep", "Import your joined data file (txt)", accept = ".txt"))
                                                                            )
                                                           ),
                                                  conditionalPanel(condition = "output.toplot_pep_dataup",
                                                                   fluidRow(column(4, selectInput("condition_plotjoinpep", "Select one or more treatments",
                                                                                                  choices = NULL,
                                                                                                  multiple = TRUE)
                                                                                   ),
                                                                            column(4, numericInput("lay_bar1_plotjoinpep", "Type the number of plot per row", min = 1, max = 10, step = 1, value = 2),
                                                                                      numericInput("lay_bar2_plotjoinpep", "Type the number of plot per column", min = 1, max = 10, step = 1, value = 2)),
                                                                            column(4, textInput("pdftit_plotjoinpep", "Choose a name for your pdf file", "barplot"))
                                                                            ),
                                                                   fluidRow(column(4, checkboxInput("ch_own_col_plotjoinpep", "Choose your own color", FALSE)),
                                                                            conditionalPanel(condition = "input.ch_own_col_plotjoinpep",
                                                                                             column(4, textOutput("n_cond_sel_plotjoinpep"),
                                                                                                       colourpicker::colourInput("own_color_pick_plotjoinpep", NULL, "#FF2B00",
                                                                                                                                 allowTransparent = TRUE, closeOnClick = TRUE),
                                                                                                       textOutput("own_color_plotjoinpep")
                                                                                                    ),
                                                                                             column(4, actionButton("add_col_plotjoinpep", "Add the color"),
                                                                                                       actionButton("rem_col_plotjoinpep", "Remove the last color")
                                                                                                    )
                                                                                             )
                                                                            ),
                                                                   tags$hr(),
                                                                   fluidRow(column(6, actionButton("getbar_plotjoinpep", "Save bar plots", class = "btn-primary")),
                                                                            column(6, textOutput("diag_bar_plotjoinpep"))
                                                                            )
                                                                   )

                                                  )
                                              )
                                     ),
                            tabPanel("Proteins", value = "proteins",
                                      shinyjs::useShinyjs(),
                                      tags$style(type = 'text/css',
                                                 '#modal1 .modal-dialog { width: fit-content !important; overflow-x: initial !important}
                                                 #modal1 .modal-body { width: 150vh; overflow-x: auto;}'
                                                 ),
                                      tags$style(type = 'text/css',
                                                 '#modal2 .modal-dialog { width: fit-content !important; overflow-x: initial !important}
                                                 #modal2 .modal-body { width: 150vh; overflow-x: auto;}'
                                                 ),
                                      tags$style(type = 'text/css',
                                                 '#modal3 .modal-dialog { width: fit-content !important; overflow-x: initial !important}
                                                 #modal3 .modal-body { width: 150vh; overflow-x: auto;}'
                                                 ),
                                      tags$style(type = 'text/css',
                                                 '#modal4 .modal-dialog { width: fit-content !important; overflow-x: initial !important}
                                                 #modal4 .modal-body { width: 150vh; overflow-x: auto;}'
                                                 ),
                                      tags$style(type = 'text/css',
                                                 '#modal5 .modal-dialog { width: fit-content !important; overflow-x: initial !important}
                                                 #modal5 .modal-body { width: 150vh; overflow-x: auto;}'
                                                 ),
                                      tags$style(type = 'text/css',
                                                 '#modal6 .modal-dialog { width: fit-content !important; overflow-x: initial !important}
                                                 #modal6 .modal-body { width: 150vh; overflow-x: auto;}'
                                                 ),
                                      tags$style(type = 'text/css',
                                                 '#modal7 .modal-dialog { width: fit-content !important; overflow-x: initial !important}
                                                 #modal7 .modal-body { width: 150vh; overflow-x: auto;}'
                                                 ),

                                      tags$style(type = 'text/css',
                                                 '#modal8_FC .modal-dialog { width: fit-content !important; overflow-x: initial !important}
                                                  #modal8_FC .modal-body { width: 150vh; overflow-x: auto;}'
                                                 ),
                                      tags$style(type = 'text/css',
                                                 '#modal9_IS .modal-dialog { width: fit-content !important; overflow-x: initial !important}
                                                  #modal9_IS .modal-body { width: 150vh; overflow-x: auto;}'
                                                 ),

                                      fluidRow(style = "height:20px;"),

                                      h1(tags$u(class = "main-1", "The IMPRINTS.CETSA analysis")),

                                      tags$br(),

                                      radioButtons("step_cetsa", "At which step do you want to start your analysis ?",
                                                   choices = c("From the beginning" = "1begin",
                                                               "Consolidate isoforms and rearrange your data" = "2conso_ISO",
                                                               "Normalize your data" = "3NORM",
                                                               "Get the protein abundance difference" = "4DIFF",
                                                               "Get your hitlist" = "5HIT"),
                                                   inline = TRUE),

                                      fluidRow(box(title = "Upload and clean your data", status = "primary",
                                                   solidHeader = TRUE, collapsible = TRUE, width = 12,
                                                   tags$u(h3("Upload your data")),
                                                   fluidRow(column(4, selectInput("n_chan", "Select the number of channels", choices = c(10,11,16,18), selected = 10),
                                                                      shiny::HTML("<br><h5>On the table on your right, you can type the name
                                                                                  of the sample to the corresponding channel. The underscore '_'
                                                                                  will be used as a separator between temperatures, bioreplicates and
                                                                                  treatments in all further functions, so make sure of your spelling.
                                                                                  <br><br>Here, you'll need to type first the bioreplicate and then
                                                                                  the treatment, like this : 'B1_Vehicle', 'B1_treatment', etc. Make sure that
                                                                                  all names are different !
                                                                                  <br>If you have a 'Mix' channel; it needs to explicitely start with
                                                                                      'Mix'.
                                                                                  <br>If you have any 'Empty' channels; it needs to explicitely start with
                                                                                      'Empty'.</h5>")
                                                                   ),

                                                            column(8, uiOutput("treat_nameui"))
                                                            ),


                                                   fileInput("PD_data", "Select txt files for your analysis",
                                                             accept = ".txt", multiple = TRUE),


                                                   conditionalPanel(condition = "output.cetsa_fileup",
                                                                    textOutput("diag_rawread"),
                                                                    actionButton("see1_cetsa", "View data uploaded"),
                                                                    tags$hr(),
                                                                    tags$u(h3("Rename your treatments and clean your data")),
                                                                    tags$hr(),

                                                                    fluidRow(column(2, shiny::HTML("<br><h5>On the table on your right, you can rename your temperatures.
                                                                                                   For example, like this: '37C', '47C', etc. <br>Remember to not use '_'.
                                                                                                   <br>Also, if you have a quantitative proteomic file, it is
                                                                                                   advised to name its 'temperature' as '36C' for easier handle
                                                                                                   in the other functions from IMPRINTS.CETSA.app.</h5>")),
                                                                             column(4, uiOutput("temp_nameui")),
                                                                             column(2, checkboxInput("rem_mix", "Remove the 'Mix' channel if any", TRUE),
                                                                                       checkboxInput("rem_empty", "Remove the 'Empty' channels if any", TRUE),
                                                                                       checkboxInput("clean_data", "Remove proteins without quantitative information", TRUE)),
                                                                             column(2, actionButton("str_ren", "Rename the treatments", class = "btn-primary")),
                                                                             column(2, actionButton("see2_cetsa", "View data renamed"))
                                                                             )
                                                                    )
                                                   )
                                               ),
                                      tags$hr(),

                                      conditionalPanel(condition = "output.cetsa_cleanup | input.step_cetsa > '2' ",
                                                       fluidRow(box(title = "Isoform ambiguity cleanup, rearrange and normalization", status = "primary",
                                                                    solidHeader = TRUE, collapsible = TRUE, width = 12,
                                                                    tags$u(h3("Isoform ambiguity cleanup and rearrange")),
                                                                    tags$hr(),

                                                                    checkboxInput("got_ISO_cetsa", "Do you already have the file isoform_resolved ?", FALSE),
                                                                    conditionalPanel(condition = "!input.got_ISO_cetsa",
                                                                                     actionButton("ISO", "Resolve isoform", class = "btn-primary")
                                                                    ),
                                                                    conditionalPanel(condition = "input.got_ISO_cetsa",
                                                                                     fileInput("ISOresfile_cetsa", "Select the file named isoform_resolved", accept = ".txt")
                                                                    ),


                                                                    conditionalPanel(condition = "output.cetsa_isoup | input.step_cetsa > '3' ",
                                                                                     tags$hr(),
                                                                                     actionButton("see3_cetsa", "View data with isoform resolved"),
                                                                                     checkboxInput("got_rearr_cetsa", "Do you already have the file data_pre_normalization
                                                                                                                       (output after rarranging your data) ?", FALSE),
                                                                                     conditionalPanel(condition = "!input.got_rearr_cetsa",
                                                                                                      fluidRow(column(6, checkboxInput("iso_conso", "Perform isoform consolidate", TRUE),
                                                                                                                      tags$hr(),
                                                                                                                      conditionalPanel(condition = "input.iso_conso",
                                                                                                                                       numericInput("n_chan2", "Type the number of reading channels", value = 9, min = 1),
                                                                                                                                       fileInput("tab_conso", "Upload the txt file containing an isoform substitution matching table",
                                                                                                                                                 accept = ".txt")
                                                                                                                                       ),
                                                                                                                      conditionalPanel(condition = "!input.iso_conso",
                                                                                                                                       tags$br()
                                                                                                                                       )
                                                                                                                      ),
                                                                                                               column(6, checkboxInput("iso_rearr", "Rearrange data", TRUE), tags$hr()),
                                                                                                               conditionalPanel(condition = "input.iso_rearr",
                                                                                                                                column(3, numericInput("n_chan3", "Type the number of reading channels", value = 9, min = 1),
                                                                                                                                          numericInput("count_thr", "Type the minimal threshold number
                                                                                                                                                        of associated abundance count of proteins",
                                                                                                                                                       value = 2, min = 1, step = 1),
                                                                                                                                       checkboxInput("wit_37", "Whether the kept proteins should have readings at 37C", FALSE)
                                                                                                                                       ),
                                                                                                                                column(3, numericInput("rep_thr", "Type the minimal percentage threshold of
                                                                                                                                                                   protein being sampled from multiple runs",
                                                                                                                                                       value = 0.1, min = 0, max = 1, step = 0.01),
                                                                                                                                          checkboxInput("avgcount_abd", "Take the median average of abundance count
                                                                                                                                                      across temperature", TRUE)
                                                                                                                                       )
                                                                                                                                )
                                                                                                               ),
                                                                                                      actionButton("ISO2", "Consolidate isoform and/or rearrange", class = "btn-primary"),
                                                                                                      textOutput("diag_rearrange")
                                                                                     ),
                                                                                     conditionalPanel(condition = "input.got_rearr_cetsa",
                                                                                                      fileInput("rearrfile_cetsa", "Select the file named data_pre_normalization", accept = ".txt")
                                                                                     ),
                                                                                     tags$hr(),
                                                                                     fluidRow(column(3, actionButton("see4_cetsa", "View consolidated data")),
                                                                                              column(3, actionButton("see5_cetsa", "View rearranged data"))),

                                                                                     tags$hr(),

                                                                                     tags$u(h3("Normalize your data")),
                                                                                     tags$hr(),

                                                                                     fluidRow(column(6, checkboxInput("got_norm_cetsa", "Do you already have the file data_post_normalization ?")),
                                                                                              column(6, conditionalPanel(condition = "!input.got_norm_cetsa",
                                                                                                                         actionButton("NORM", "Start Normalization", class = "btn-primary")
                                                                                              ),
                                                                                              conditionalPanel(condition = "input.got_norm_cetsa",
                                                                                                               fileInput("normfile_cetsa", "Select the file named data_post_normalization", accept = ".txt")
                                                                                              )
                                                                                              )
                                                                                     )
                                                                    )
                                                       )
                                                       ),

                                                       tags$hr(),
                                                       conditionalPanel(condition = "output.cetsa_normup | input.step_cetsa > '4' ",
                                                                        fluidRow(box(title = "Abundance difference calculation and hitlist", status = "primary",
                                                                                     solidHeader = TRUE, collapsible = TRUE, width = 12,
                                                                                     actionButton("see6_cetsa", "View normalized data"),
                                                                                     tags$hr(),
                                                                                     tags$u(h3("Calculate the pair-wise protein abundance differences")),
                                                                                     tags$hr(),

                                                                                     checkboxInput("got_diff_cetsa", "Do you already have the file imprints_caldiff ?", FALSE),
                                                                                     conditionalPanel(condition = "!input.got_diff_cetsa",
                                                                                                      fluidRow(column(4, selectInput("ctrl_name2", "Select the treatment that corresponds to your control",
                                                                                                                                     choices = NULL)
                                                                                                                      ),
                                                                                                               column(4, checkboxInput("wit_rep", "Whether the calculation of the relative protein
                                                                                                                                       abundance difference should be within the same biorep", TRUE)),
                                                                                                               column(4, actionButton("CAL_DIF", "Start difference calculation", class = "btn-primary"))
                                                                                                      )
                                                                                     ),
                                                                                     conditionalPanel(condition = "input.got_diff_cetsa",
                                                                                                      fileInput("difffile_cetsa", "Select the file named imprints_caldiff", accept = ".txt")),

                                                                                     tags$hr(),

                                                                                     conditionalPanel(condition = "output.cetsa_difup | input.step_cetsa > '5' ",
                                                                                                      actionButton("see7_cetsa", "View caldiff output"),
                                                                                                      tags$hr(),
                                                                                                      tags$u(h3("Get the protein hitlist")),
                                                                                                      tags$hr(),

                                                                                                      conditionalPanel(condition = "!input.calc_diff",
                                                                                                                       radioButtons("hitmethod_cetsa", "Choose a method to get your hitlist",
                                                                                                                                    choices = c("Intercept Score" = "IS",
                                                                                                                                                "IMPRINTS score" = "ImpS",
                                                                                                                                                "Fold Change cutoff" = "FC"
                                                                                                                                                ),
                                                                                                                                    selected = "IS",
                                                                                                                                    inline = TRUE),
                                                                                                                       conditionalPanel(condition = "input.hitmethod_cetsa == 'ImpS'",
                                                                                                                                        tags$hr(),
                                                                                                                                        fluidRow(column(4, selectInput("formatImpS_cetsa", "Choose how many categories to segregate the results", choices = c("4", "9"), selected = "4")),
                                                                                                                                                 column(4, numericInput("FDRImpS_cetsa", "Choose the FDR", value = 0.01, min = 0, max = 1, step = 0.01)),
                                                                                                                                                 column(4, numericInput("cvcutoff_cetsa", "Choose the significance level of threshold used for CV quality control", value = 2.5, min = 0, max = 10, step = 0.25))
                                                                                                                                                 ),
                                                                                                                                        tags$hr(),
                                                                                                                                        textOutput("diag_ImpS"),
                                                                                                                                        tags$hr()
                                                                                                                                        ),
                                                                                                                       conditionalPanel(condition = "input.hitmethod_cetsa == 'IS'",
                                                                                                                                        actionButton("see9_cetsa", "See more information"),
                                                                                                                                        tags$hr(),
                                                                                                                                        fluidRow(column(4, numericInput("IScut_cetsa", "Choose a Intercept Score cutoff", value = 1.5, min = 0, step = 0.1)),
                                                                                                                                                 column(4, numericInput("FDR_cetsa", "Choose the FDR", value = 0.01, min = 0, max = 1, step = 0.01)),
                                                                                                                                                 column(4, numericInput("validval_cetsa", "Choose the minimum proportion of valid values", value = 0, min = 0, max = 1, step = 0.05))
                                                                                                                                                 ),
                                                                                                                                        tags$hr(),
                                                                                                                                        textOutput("diag_IS"),
                                                                                                                                        tags$hr()
                                                                                                                                        ),
                                                                                                                       conditionalPanel(condition = "input.hitmethod_cetsa == 'FC'",
                                                                                                                                        actionButton("see8_cetsa", "See more information"),
                                                                                                                                        tags$hr(),
                                                                                                                                        fluidRow(column(4, numericInput("meancut_cetsa", "Choose a mean cutoff", value = 0.25, min = 0, step = 0.01)),
                                                                                                                                                 column(4, numericInput("bound_cetsa", "Choose the boundedness", value = 4)),
                                                                                                                                                 column(4, checkboxInput("save_hit", "Save the hitlist", TRUE))
                                                                                                                                                 )
                                                                                                                                        ),

                                                                                                                       actionButton("str_calchitlist", "Start calculation", class = "btn-primary"),
                                                                                                                       tags$hr(),

                                                                                                                       radioButtons("HIT", h3("Choose a result to print"),
                                                                                                                                    choices = c("hitlist", "CC", "CN", "NC", "ND"),
                                                                                                                                    selected = "hitlist", inline = TRUE),

                                                                                                                       DT::dataTableOutput("hit_out")
                                                                                                      )
                                                                                     )
                                                                        )
                                                                        ),

                                                                        h3("Go check your file in your working directory, all the results from your analysis should be saved !"),
                                                                        tags$hr()
                                                                        )
                                                       )
                                     )
                            ),

                 tabPanel("Database", value = "database",
                          shinyjs::useShinyjs(),
                          fluidRow(style = "height:20px;"),

                          h1(tags$u(class = "main-1", "Add new dataset and remove old ones")),
                          tags$hr(),

                          htmlOutput("info_daba"),
                          tags$hr(),

                          HTML("<p><h5>In order to add new dataset, you need to import two files.<br>
                              This files are : <br>
                              - The output from the imprints_caldiff function from the IMPRINTS.CETSA package <br>
                              - The file named 'Summary', from the hitlist function output OR,
                                if you have it, the analysis tab from IS function that contains the IS for all proteins<br>
                              Once you uploaded these two files, choose a name for your dataset (like 'elutriation' for example),
                              click on the button 'Add dataset', and you're good to go !</h5></p>
                              <p><h5>If you want to remove a dataset, select one of the dataset available from the database,
                              and click on the button 'Remove dataset', in the box below. Beware, this operation cannot be undone !</h5></p>
                              <br><h5>Once you made your changements, don't forget to click on the button 'Reload the database' to use it directly.</h5>"
                          ),

                          fluidRow(box(title = "Add new dataset", status = "success", solidHeader = TRUE, collapsible = TRUE, width = 12,
                                       fluidRow(column(6, fileInput("caldif_daba", "Import the output from imprints_caldiff"),
                                                       checkboxInput("gave_daba", "Don't have the imprints_average output
                                                                                   (will calculate and save it)", TRUE),
                                                       conditionalPanel(condition = "!input.gave_daba",
                                                                        fileInput("AVE_dabafile", "Import the output from imprints_average")
                                                       )
                                       ),
                                       column(6, fileInput("hitsum_daba", "Import the summary file OR the analysis tab from the hitlist outputs")),
                                       ),
                                       conditionalPanel(condition = "output.DIFdaba_fileup & output.AVEdaba_fileup & output.HITdaba_fileup",
                                                        fluidRow(column(6, textInput("name_daba", "Type a name for your new dataset")),
                                                                 column(6, actionButton("add_daba", "Add dataset", class = "btn-success btn-lg"))
                                                        )
                                       )

                          )),

                          fluidRow(box(title = "Rename your treatments", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
                                       fluidRow(column(6, uiOutput("davai2_daba_ui")),
                                                column(6, htmlOutput("condfrom_daba"),
                                                       textInput("condnew_daba", "Type the new names of the treatments
                                                                   (same order; separated by a comma; if empty, no changement)"))
                                       ),
                                       actionButton("changename_daba", "Change the name of your treatments", class = "btn-primary btn-lg")
                          )
                          ),

                          fluidRow(box(title = "Remove dataset", status = "danger", solidHeader = TRUE, collapsible = TRUE, width = 12,
                                       fluidRow(column(6, uiOutput("davai_daba_ui")),
                                                column(6, actionButton("rem_daba", "Remove dataset", class = "btn-danger btn-lg"))
                                       )
                          )),

                          actionButton("up_daba", "Reload the database", class = "btn-primary btn-lg"),

                          tags$hr()
                 ),

                 navbarMenu("Visualization",
                            tabPanel("Barplot", value = "visu",
                                     shinyjs::useShinyjs(),
                                     tags$style(HTML(".tabbable > .nav > li > a                  {background-color: #A1BAC8;  color:#FFFFFF}
                                                                  .tabbable > .nav > li[class=active]    > a {background-color: #3C8DBC; color:#FFFFFF}
                                                                  .tabbable > .nav > li    > a:hover {background-color: #3BAAE6; color:#FFFFFF}
                                                                  .tabbable > .nav > li[class=active]    > a:hover {background-color: #3C8DBC; color:#FFFFFF}")
                                     ),

                                     fluidRow(style = "height:20px;"),

                                     tabsetPanel(type = "tabs", selected = "2D Bar plot",
                                                 tabPanel("2D Bar plot",
                                                          h1(tags$u(class = "main-1", "Get the 2D bar plot")),
                                                          tags$hr(),

                                                          fluidRow(box(title = "2D bar plot parameters", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,

                                                                       radioButtons("drug", "Choose a dataset",
                                                                                    choices = c("Database" = "base",
                                                                                                "Your data" = "dat"),
                                                                                    selected = "base", inline = TRUE),
                                                                       conditionalPanel(condition = "input.drug == 'base'",
                                                                                        uiOutput("drug2_ui")
                                                                       ),


                                                                       tags$hr(),

                                                                       conditionalPanel(condition = "input.drug == 'dat' ",

                                                                                        fluidRow(column(6, fileInput("data_barplot", "Upload the 'imprints_caldiff' file (log2 fold change)",
                                                                                                                     accept = c(".txt", ".csv", ".xlsx"))
                                                                                        ),
                                                                                        column(6, checkboxInput("calc_hitlist", "Find the hitlist from your data file", FALSE),
                                                                                               conditionalPanel(condition = "!input.calc_hitlist",
                                                                                                                fileInput("data_hitlist", "Upload your own hitlist, the summary file from the hitlist outputs",
                                                                                                                          accept = c(".txt", ".csv", ".xlsx"), multiple = TRUE)),
                                                                                               conditionalPanel(condition = "input.calc_hitlist",
                                                                                                                fluidRow(column(2, numericInput("meancut_bar", "Choose a mean cutoff", value = 0.25, min = 0, step = 0.01)),
                                                                                                                         column(2, numericInput("bound_bar", "Choose the boundedness", value = 4)),
                                                                                                                         column(2, checkboxInput("save_hit_bar", "Save the hitlist", FALSE))
                                                                                                                ),
                                                                                                                actionButton("str_calchit", "Start calculation"))
                                                                                        )
                                                                                        ),

                                                                                        tags$hr()
                                                                       ),

                                                                       fluidRow(
                                                                         column(4, checkboxInput("protlist_bar", "Import a protein list", FALSE),
                                                                                conditionalPanel(condition = "!input.protlist_bar",
                                                                                                 checkboxInput("hit", "Only take the hited proteins", FALSE),
                                                                                                 conditionalPanel(condition = "input.hit",
                                                                                                                  selectInput("cond_fhit", "Select hits treatment", choices = NULL, multiple = TRUE))
                                                                                ),
                                                                                conditionalPanel(condition = "input.protlist_bar",
                                                                                                 fileInput("prlist_file_bar", "Import your protein list (txt file)", accept = ".txt"),
                                                                                ),
                                                                                checkboxInput("ALL_prot", "Select all the proteins", FALSE),
                                                                                conditionalPanel(condition = "!input.ALL_prot",
                                                                                                 selectizeInput("prot", "Select a protein", choices = NULL, multiple = TRUE))

                                                                         ),
                                                                         column(4, conditionalPanel(condition = "input.cond_sel != 'cat' ",
                                                                                                    checkboxInput("rem_con", "Remove the controls", FALSE),
                                                                                                    conditionalPanel(condition = "input.rem_con",
                                                                                                                     textInput("con_name", "Type the name of your controls (if several names, separate them by |)", "G1")
                                                                                                                     )
                                                                                                    )
                                                                                ),
                                                                         column(4, radioButtons("cond_sel", "Selection type",
                                                                                                choices = c("Select the treatment level" = "treat",
                                                                                                            "Select treatment level by category" = "cat",
                                                                                                            "Select all the treatment level" = "all_cond"),
                                                                                                selected = "treat"),

                                                                                conditionalPanel(condition = "input.cond_sel != 'all_cond' ",
                                                                                                 selectInput("cond", "Select one or more treatments",
                                                                                                             choices = NULL,
                                                                                                             multiple = TRUE)
                                                                                                 )
                                                                                )
                                                                       ),

                                                                       tags$hr(),
                                                                       fluidRow(column(4, checkboxInput("ch_own_col", "Choose your own color", FALSE)),
                                                                                column(4, conditionalPanel(condition = "input.ch_own_col",
                                                                                                           textOutput("n_cond_sel"),
                                                                                                           colourpicker::colourInput("own_color_pick", NULL, "#FF2B00",
                                                                                                                       allowTransparent = TRUE, closeOnClick = TRUE),
                                                                                                           textOutput("own_color")
                                                                                )
                                                                                ),
                                                                                conditionalPanel(condition = "input.ch_own_col",
                                                                                                 column(4, actionButton("add_col", "Add the color"),
                                                                                                        actionButton("rem_col", "Remove the last color"))
                                                                                )
                                                                       ),

                                                                       tags$hr(),

                                                                       fluidRow(column(4, checkboxInput("save_bar", "Save the bar plots in a pdf file", FALSE),
                                                                                       conditionalPanel(condition = "input.save_bar",
                                                                                                        textInput("pdftit", "Choose a name for your pdf file", "barplot")
                                                                                                        )
                                                                                       ),
                                                                                conditionalPanel(condition = "input.save_bar",
                                                                                                 column(4, numericInput("lay_bar1", "Type the number of plot per row",
                                                                                                                        min = 1, max = 10, step = 1, value = 4),
                                                                                                          numericInput("lay_bar2", "Type the number of plot per column",
                                                                                                                       min = 1, max = 10, step = 1, value = 3)),
                                                                                                 column(4, numericInput("pdfw", "Type the width of the pdf page",
                                                                                                                        min = 1, step = 1, value = 12),
                                                                                                        numericInput("pdfh", "Type the height of the pdf page",
                                                                                                                     min = 1, step = 1, value = 12))
                                                                                )
                                                                       ),

                                                                       tags$hr(),

                                                                       fluidRow(column(4, checkboxInput("werb", "Print error bar", TRUE),
                                                                                          checkboxInput("wpts", "Print point of each replicate", FALSE)),
                                                                                column(4, checkboxInput("grad", "Use color gradient", FALSE)),
                                                                                column(4, checkboxInput("line", "Use line instead of bar", FALSE))
                                                                       )
                                                          )),

                                                          DT::dataTableOutput("pr_info"),
                                                          conditionalPanel(condition = "output.identifcomp_barup",
                                                                           downloadButton("downtabidentif_barplot", "Download the identification comparison tab")),

                                                          tags$hr(),
                                                          actionButton("barp", "See bar plot", class = "btn-primary btn-lg"),
                                                          tags$hr(),
                                                          textOutput("diag_bar"),
                                                          tags$hr(),

                                                          withSpinner(plotOutput("bar_plot", height = "800px"), type = 6),
                                                          fluidRow(column(2, tags$div(style="line-height:175%;",
                                                                                      tags$br()
                                                                                      ),
                                                                             downloadButton("downbar", "Download plot")),
                                                                   column(2, selectInput("downbar_format", "Download as", choices = c("png", "pdf")))
                                                                   ),

                                                          tags$hr()
                                                 ),

                                                 tabPanel("Protein complex",
                                                          h1(tags$u(class = "main-1", "Protein complex and 2D bar plot")),
                                                          tags$hr(),

                                                          fluidRow(box(title = "Map proteins to known protein complex", status = "primary",
                                                                       solidHeader = TRUE, collapsible = TRUE, width = 12,

                                                                       radioButtons("drug_compl", "Choose a dataset",
                                                                                    choices = c("Database" = "base",
                                                                                                "Your data" = "dat"),
                                                                                    selected = "base", inline = TRUE),
                                                                       conditionalPanel(condition = "input.drug_compl == 'base'",
                                                                                        uiOutput("drug2ui_compl")
                                                                       ),

                                                                       conditionalPanel(condition = "input.drug_compl == 'dat' ",
                                                                                        fluidRow(column(6, fileInput("caldif_compl", "Import the output from imprints_caldiff")),
                                                                                                 column(6, fileInput("hitsum_compl", "Import the summary file from the hitlist outputs"))
                                                                                                 ),
                                                                                        fluidRow(column(6, checkboxInput("gave_compl", "Don't have the imprints_average output
                                                                                                                                        (will calculate and save it)", TRUE)),
                                                                                                 conditionalPanel(condition = "!input.gave_compl",
                                                                                                                  column(6, fileInput("avef_compl", "Import the output from imprints_average"))
                                                                                                 )
                                                                                        )
                                                                       ),
                                                                       tags$hr(),

                                                                       conditionalPanel(condition = "output.DIFcompl_fileup & output.HITcompl_fileup & output.NNcompl_fileup & output.AVEcompl_fileup",
                                                                                        fluidRow(column(4, selectInput("condsel_compl", "Select a treatment", choices = NULL)),
                                                                                                 column(4, selectInput("catego_compl", "Select some categories", choices = NULL, multiple = TRUE)),
                                                                                                 column(4, selectInput("organism_compl", "Choose an organism", choices = c("Human", "Mouse", "Rat"), selected = "Human"))
                                                                                        ),

                                                                                        actionButton("ave_map_compl", "Map proteins to known protein complex", class = "btn-primary btn-lg"),
                                                                                        textOutput("diagmapping_compl"),

                                                                                        tags$hr(),

                                                                                        conditionalPanel(condition = "output.resmappingcompl_fileup",
                                                                                                         DT::dataTableOutput("tabmap_compl"),
                                                                                                         downloadButton("downrestab_compl")
                                                                                        )
                                                                       )

                                                          )
                                                          ),

                                                          conditionalPanel(condition = "output.resmappingcompl_fileup",
                                                                           fluidRow(box(title = "2D bar plot parameter", status = "primary",
                                                                                        solidHeader = TRUE, collapsible = TRUE, width = 12,

                                                                                        fluidRow(column(4,selectInput("allcomplex_compl", "Select some protein complex", choices = NULL, multiple = TRUE)),
                                                                                                 column(4, checkboxInput("ALL_prot_compl", "Select all the proteins", FALSE)),
                                                                                                 conditionalPanel(condition = "!input.ALL_prot_compl",
                                                                                                                  column(4,selectizeInput("prot_compl", "Select a protein", choices = NULL, multiple = TRUE))
                                                                                                 )
                                                                                        ),

                                                                                        tags$hr(),
                                                                                        fluidRow(column(4, checkboxInput("ch_own_col_compl", "Choose your own color", FALSE)),
                                                                                                 column(4, conditionalPanel(condition = "input.ch_own_col_compl",
                                                                                                                            colourpicker::colourInput("own_color_pick_compl", NULL, "#FF2B00",
                                                                                                                                        allowTransparent = TRUE, closeOnClick = TRUE)
                                                                                                                            )
                                                                                                        )
                                                                                                 ),

                                                                                        tags$hr(),

                                                                                        fluidRow(column(4, checkboxInput("save_bar_compl", "Save the bar plots in a pdf file", FALSE),
                                                                                                        conditionalPanel(condition = "input.save_bar_compl",
                                                                                                                         textInput("pdftit_compl", "Choose a name for your pdf file", "barplot")
                                                                                                                         )
                                                                                                        ),
                                                                                                 conditionalPanel(condition = "input.save_bar_compl",
                                                                                                                  column(4, numericInput("lay_bar1_compl", "Type the number of plot per row",
                                                                                                                                         min = 1, max = 10, step = 1, value = 4),
                                                                                                                         numericInput("lay_bar2_compl", "Type the number of plot per column",
                                                                                                                                      min = 1, max = 10, step = 1, value = 3)),
                                                                                                                  column(4, numericInput("pdfw_compl", "Type the width of the pdf page",
                                                                                                                                         min = 1, step = 1, value = 12),
                                                                                                                         numericInput("pdfh_compl", "Type the height of the pdf page",
                                                                                                                                      min = 1, step = 1, value = 12))
                                                                                                 )
                                                                                        ),

                                                                                        tags$hr(),

                                                                                        fluidRow(column(4, checkboxInput("werb_compl", "Print error bar", TRUE),
                                                                                                           checkboxInput("wpts_compl", "Print point of each replicate", FALSE)),
                                                                                                 column(4, checkboxInput("grad_compl", "Use color gradient", FALSE)),
                                                                                                 column(4, checkboxInput("line_compl", "Use line instead of bar", FALSE))
                                                                                        ),

                                                                                        tags$hr(),
                                                                                        actionButton("barp_compl", "See bar plot", class = "btn-primary btn-lg"),
                                                                                        tags$hr(),
                                                                                        textOutput("diag_bar_compl"),
                                                                                        tags$hr(),

                                                                                        withSpinner(plotOutput("bar_plot_compl", height = "800px"), type = 6),
                                                                                        fluidRow(column(2, tags$div(style="line-height:175%;",
                                                                                                                    tags$br()
                                                                                                                    ),
                                                                                                        downloadButton("downbar_compl", "Download plot")),
                                                                                                 column(2, selectInput("downbar_compl_format", "Download as", choices = c("png", "pdf")))
                                                                                        )
                                                                           )
                                                                           ),

                                                          ),

                                                          tags$hr()

                                                 ),

                                                 tabPanel("Similar profiles",
                                                          h1(tags$u(class = "main-1", "Find similar profiles")),
                                                          tags$hr(),

                                                          fluidRow(box(title = "2D bar plot parameters", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,

                                                                       radioButtons("drug_simpf", "Choose a dataset",
                                                                                    choices = c("Database" = "base",
                                                                                                "Your data" = "dat"),
                                                                                    selected = "base", inline = TRUE),
                                                                       conditionalPanel(condition = "input.drug_simpf == 'base'",
                                                                                        uiOutput("drug2ui_simpf")
                                                                       ),

                                                                       conditionalPanel(condition = "input.drug_simpf == 'dat' ",
                                                                                        fluidRow(column(4, fileInput("cdiff_simpf", "Import the output from imprints_caldiff")),
                                                                                                 column(4, checkboxInput("gave_simpf", "Don't have the imprints_average output
                                                                                                                                       (will calculate and save it)", TRUE)),
                                                                                                 conditionalPanel(condition = "!input.gave_simpf",
                                                                                                                  column(4, fileInput("avef_simpf", "Import the output from imprints_average"))
                                                                                                 )
                                                                                        )
                                                                       )
                                                                       ,

                                                                       conditionalPanel(condition = "output.AVEsimpf_fileup & output.DIFsimpf_fileup",
                                                                                        fluidRow(column(3, selectInput("treat_simpf", "Select a treatment", choices = NULL)),
                                                                                                 column(3, selectizeInput("prot_simpf", "Select a protein from which you want to get
                                                                                                                                         the similar profiles", choices = NULL)),
                                                                                                 column(3, sliderInput("maxna_simpf", "Choose a maximum number of
                                                                                                                                       missing values per rows", value = 0, min = 0, max = 10, step = 1)),
                                                                                                 column(3, selectInput("scoremeth_simpf", "Select a method for calculating the similarity score",
                                                                                                                       choices = c("Euclidean distance score" = "euclidean",
                                                                                                                                   "Pearson correlation" = "pearson"), selected = "euclidean"),
                                                                                                        numericInput("scothr_simpf", "Choose a threshold for the similarity score",
                                                                                                                     value = 0.9, min = 0, max = 1, step = 0.01))
                                                                                        ),
                                                                                        tags$hr(),

                                                                                        checkboxInput("infoscmeth_simpf", h4("See some informations about the method for calculating the similarity score"), FALSE),
                                                                                        conditionalPanel(condition = "input.infoscmeth_simpf",
                                                                                                         HTML("<p><h5>You have actually two methods for calculating the similarity score : <br>
                                                                                                               - The euclidean distance score <br>
                                                                                                               - The Pearson correlation <br>
                                                                                                               This score will determine which proteins got a similar profile from the one you selected. <br>
                                                                                                               The euclidean distance score : <br>
                                                                                                               With this method, every euclidean distance between the value from each protein profile and the
                                                                                                               selected one will be calculated. Then for each distance a score between 0 and 1 is calculated
                                                                                                               by dividing 1 by 1 + d, where d is the euclidean distance. <br>
                                                                                                               This score means that you will search for protein profile with similar values from the the one you selected.
                                                                                                               So the profile with a similar shape but with lower or higher values will not have a good score.
                                                                                                               It also means that with a high score (~0.9) you're not very likely to find a lot of proteins.<br><br>
                                                                                                               This is not the case with Pearson correlation. For this score, each covariance and standard deviation
                                                                                                               between the protein you selected and all the other proteins will be calculated. Then, the covariance is divided
                                                                                                               by the product of the two standard devation. It gives you score between -1 and 1. -1 means the data are negatively
                                                                                                               correlated, 1 positively correlated and 0 not correlated. <br>
                                                                                                               Because you calculate a correlation score, you will search for all proteins profile with a similar shape from the one
                                                                                                               you selected, no matter their values. It also means that with a high score (~0.95) you may find a lot of proteins.
                                                                                                              </h5></p>")),


                                                                                        tags$hr(),

                                                                                        fluidRow(column(4, checkboxInput("ch_own_col_simpf", "Choose your own color", FALSE)),
                                                                                                 conditionalPanel(condition = "input.ch_own_col_simpf",
                                                                                                                  column(4, colourpicker::colourInput("own_color_pick_simpf", NULL, "#FF2B00",
                                                                                                                                        allowTransparent = TRUE, closeOnClick = TRUE)
                                                                                                                  )
                                                                                                 )
                                                                                        ),

                                                                                        tags$hr(),

                                                                                        fluidRow(column(4, checkboxInput("save_bar_simpf", "Save the bar plots in a pdf file", TRUE),
                                                                                                        checkboxInput("save_prot_simpf", "Save the list of proteins ID with similar profiles (will save in a xlsx file)", TRUE),
                                                                                                        conditionalPanel(condition = "input.save_bar_simpf",
                                                                                                                         textInput("pdftit_simpf", "Choose a name for your pdf file", "barplot"))
                                                                                                        ),
                                                                                                 conditionalPanel(condition = "input.save_bar_simpf",
                                                                                                                  column(4, numericInput("lay_bar1_simpf", "Type the number of plot per row",
                                                                                                                                         min = 1, max = 10, step = 1, value = 4),
                                                                                                                         numericInput("lay_bar2_simpf", "Type the number of plot per column",
                                                                                                                                      min = 1, max = 10, step = 1, value = 3)),
                                                                                                                  column(4, numericInput("pdfw_simpf", "Type the width of the pdf page",
                                                                                                                                         min = 1, step = 1, value = 12),
                                                                                                                         numericInput("pdfh_simpf", "Type the height of the pdf page",
                                                                                                                                      min = 1, step = 1, value = 12))
                                                                                                                  )
                                                                                        ),

                                                                                        checkboxInput("seeprsel_simpf", "See barplot from the selected protein", FALSE),

                                                                                        tags$hr(),

                                                                                        fluidRow(column(4, checkboxInput("werb_simpf", "Print error bar", TRUE)),
                                                                                                 column(4, checkboxInput("grad_simpf", "Use color gradient", FALSE)),
                                                                                                 column(4, checkboxInput("line_simpf", "Use line instead of bar", FALSE))
                                                                                        ),

                                                                                        tags$hr(),

                                                                                        actionButton("getsimi_simpf", "Get similar profile !", class = "btn-primary btn-lg"),
                                                                                        tags$hr(),
                                                                                        textOutput("diag_bar_simpf"),
                                                                                        tags$hr(),




                                                                                        withSpinner(plotOutput("bar_plot_simpf", height = "800px"), type = 6),
                                                                                        fluidRow(column(2, tags$div(style="line-height:175%;",
                                                                                                                    tags$br()
                                                                                                                    ),
                                                                                                        downloadButton("downbar_simpf", "Download plot")),
                                                                                                 column(2, selectInput("downbar_simpf_format", "Download as", choices = c("png", "pdf")))
                                                                                        )
                                                                       )

                                                          )
                                                          ),

                                                          tags$hr()
                                                 )
                                     )
                            ),
                            tabPanel("Heatmap", value = "visu",
                                     shinyjs::useShinyjs(),
                                     tags$style(HTML(".tabbable > .nav > li > a                  {background-color: #A1BAC8;  color:#FFFFFF}
                                                                  .tabbable > .nav > li[class=active]    > a {background-color: #3C8DBC; color:#FFFFFF}
                                                                  .tabbable > .nav > li    > a:hover {background-color: #3BAAE6; color:#FFFFFF}
                                                                  .tabbable > .nav > li[class=active]    > a:hover {background-color: #3C8DBC; color:#FFFFFF}")
                                     ),
                                     fluidRow(style = "height:20px;"),

                                     tabsetPanel(type = "tabs",

                                                 tabPanel("Heatmap",
                                                          h2(tags$u(class = "main-1", "Get heatmaps from your data")),
                                                          tags$hr(),

                                                          fluidRow(box(title = "Import your data and Heatmap parameter", status = "primary",
                                                                       solidHeader = TRUE, collapsible = TRUE, width = 12,

                                                                       radioButtons("drug_heat", "Choose a dataset",
                                                                                    choices = c("Database" = "base",
                                                                                                "Your data" = "dat"),
                                                                                    selected = "base", inline = TRUE),
                                                                       conditionalPanel(condition = "input.drug_heat == 'base'",
                                                                                        uiOutput("drug2ui_heat")
                                                                       ),

                                                                       conditionalPanel(condition = "input.drug_heat == 'dat' ",
                                                                                        fluidRow(column(6, conditionalPanel(condition = "input.gave_heat",
                                                                                                                         fileInput("filedif_heat", "Choose an imprints_caldiff output")
                                                                                                                         ),
                                                                                                        conditionalPanel(condition = "!input.gave_heat",
                                                                                                                         fileInput("fileave_heat", "Choose an imprints_average output")
                                                                                                                         ),
                                                                                                        checkboxInput("gave_heat", "Don't have the imprints_average output
                                                                                                                                       (will calculate and save it)", TRUE)
                                                                                                        ),
                                                                                                 column(6, fileInput("summary_heat", "Choose the summary file from the hitlist output"))
                                                                                                 )
                                                                                        ),

                                                                       tags$hr(),

                                                                       conditionalPanel(condition = "output.heat_fileup & output.HITheat_fileup & output.NNheat_fileup",
                                                                                        fluidRow(column(3, selectInput("cond_heat", "Select a treatment", choices = NULL)),
                                                                                                 column(3, selectInput("resp_heat", "Select a response to the drug",
                                                                                                                       choices = c("Stabilization" = "S",
                                                                                                                                   "Destabilization" = "D",
                                                                                                                                   "Both" = "both"), selected = "both")),
                                                                                                 column(3, selectInput("catego_heat", "Select some categories", choices = NULL, multiple = TRUE)),
                                                                                                 column(3, sliderInput("maxna_heat", "Choose a maximum number of
                                                                                                                                      missing values per rows", value = 0, min = 0, max = 7, step = 1))
                                                                                                 ),

                                                                                        fluidRow(column(3, textInput("titleH_heat", "Type a title for your heatmap", "Heatmap")),
                                                                                                 column(3, colourpicker::colourInput("backcol_heat", "Choose a background color", "#FFFFFF",
                                                                                                                       allowTransparent = TRUE, closeOnClick = TRUE)),
                                                                                                 column(3, colourpicker::colourInput("bordercol_heat", "Choose a border color (can be NULL)", NULL,
                                                                                                                       allowTransparent = TRUE, closeOnClick = TRUE)),
                                                                                                 column(3,
                                                                                                        colourpicker::colourInput("grad1col_heat", "Choose the low gradient color", "#005EFF",
                                                                                                                    allowTransparent = TRUE, closeOnClick = TRUE),
                                                                                                        colourpicker::colourInput("grad2col_heat", "Choose the middle gradient color", "#FFFFFF",
                                                                                                                    allowTransparent = TRUE, closeOnClick = TRUE),
                                                                                                        colourpicker::colourInput("grad3col_heat", "Choose the high gradient color", "#FF0000",
                                                                                                                    allowTransparent = TRUE, closeOnClick = TRUE))
                                                                                        ),

                                                                                        fluidRow(column(4, checkboxInput("saveH_heat", "Save the heatmap", TRUE)),
                                                                                                 conditionalPanel(condition = "input.saveH_heat",
                                                                                                                  column(4, textInput("fnameH_heat", "Type your file name", "My_heatmap")),
                                                                                                                  column(4, selectInput("formatH_heat", "Choose a format for your file",
                                                                                                                                        choices = c("png", "pdf"), selected = "png"))
                                                                                                 )
                                                                                        )
                                                                       )
                                                          )
                                                          ),


                                                          conditionalPanel(condition = "output.heat_fileup & output.HITheat_fileup & output.NNheat_fileup",
                                                                           actionButton("getH_heat", "See heatmap", class = "btn-primary btn-lg"),
                                                                           tags$hr(),

                                                                           textOutput("diagl_heat"),
                                                                           tags$hr(),

                                                                           withSpinner(plotOutput("H_heat", height = "800px"), type = 6)
                                                          )
                                                 ),

                                                 tabPanel("Protein complex",
                                                          h2(tags$u(class = "main-1", "Protein complex and heatmap")),
                                                          tags$hr(),

                                                          fluidRow(box(title = "Map proteins to known protein complex", status = "primary",
                                                                       solidHeader = TRUE, collapsible = TRUE, width = 12,

                                                                       radioButtons("drug_heatcom", "Choose a dataset",
                                                                                    choices = c("Database" = "base",
                                                                                                "Your data" = "dat"),
                                                                                    selected = "base", inline = TRUE),
                                                                       conditionalPanel(condition = "input.drug_heatcom == 'base'",
                                                                                        uiOutput("drug2ui_heatcom")
                                                                       ),

                                                                       conditionalPanel(condition = "input.drug_heatcom == 'dat' ",
                                                                                        fluidRow(column(4, checkboxInput("gave_heatcom", "Don't have the imprints_average output
                                                                                                                                          (will calculate and save it)", TRUE)),
                                                                                                 column(4,
                                                                                                        conditionalPanel(condition = "input.gave_heatcom",
                                                                                                                         fileInput("filedif_heatcom", "Choose an imprints_caldiff output")
                                                                                                        ),
                                                                                                        conditionalPanel(condition = "!input.gave_heatcom",
                                                                                                                         fileInput("fileave_heatcom", "Choose an imprints_average output")
                                                                                                        )
                                                                                                 )
                                                                                        )
                                                                       ),

                                                                       tags$hr(),

                                                                       conditionalPanel(condition = "output.heatcom_fileup",
                                                                                        fluidRow(column(4, selectInput("cond_heatcom", "Select a treatment", choices = NULL)),
                                                                                                 column(4, selectInput("organism_heatcom", "Choose an organism", choices = c("Human", "Mouse", "Rat"), selected = "Human"))
                                                                                        ),

                                                                                        actionButton("ave_map_heatcom", "Map proteins to known protein complex", class = "btn-primary btn-lg"),
                                                                                        textOutput("diagmapping_heatcom"),

                                                                                        tags$hr(),

                                                                                        conditionalPanel(condition = "output.resmappingheatcom_fileup",
                                                                                                         DT::dataTableOutput("tabmap_heatcom"),
                                                                                                         downloadButton("downrestab_heatcom")
                                                                                                         )
                                                                                        )
                                                          )
                                                          ),

                                                          conditionalPanel(condition = "output.resmappingheatcom_fileup",
                                                                           fluidRow(box(title = "Heatmap parameter", status = "primary",
                                                                                        solidHeader = TRUE, collapsible = TRUE, width = 12,

                                                                                        fluidRow(column(4,selectInput("allcomplex_heatcom", "Select some protein complex", choices = NULL, multiple = TRUE)),
                                                                                                 column(4, selectInput("resp_heatcom", "Select a response to the drug",
                                                                                                                       choices = c("Stabilization" = "S",
                                                                                                                                   "Destabilization" = "D",
                                                                                                                                   "Both" = "both"), selected = "both")),
                                                                                                 column(4, sliderInput("maxna_heatcom", "Choose a maximum number of
                                                                                                                                         missing values per rows", value = 0, min = 0, max = 7, step = 1))
                                                                                                 ),

                                                                                        fluidRow(column(3, textInput("titleH_heatcom", "Type a title for your heatmap", "Heatmap")),
                                                                                                 column(3, colourpicker::colourInput("backcol_heatcom", "Choose a background color", "#FFFFFF",
                                                                                                                       allowTransparent = TRUE, closeOnClick = TRUE)),
                                                                                                 column(3, colourpicker::colourInput("bordercol_heatcom", "Choose a border color (can be NULL)", NULL,
                                                                                                                       allowTransparent = TRUE, closeOnClick = TRUE)),
                                                                                                 column(3,
                                                                                                        colourpicker::colourInput("grad1col_heatcom", "Choose the low gradient color", "#005EFF",
                                                                                                                    allowTransparent = TRUE, closeOnClick = TRUE),
                                                                                                        colourpicker::colourInput("grad2col_heatcom", "Choose the middle gradient color", "#FFFFFF",
                                                                                                                    allowTransparent = TRUE, closeOnClick = TRUE),
                                                                                                        colourpicker::colourInput("grad3col_heatcom", "Choose the high gradient color", "#FF0000",
                                                                                                                    allowTransparent = TRUE, closeOnClick = TRUE))
                                                                                                 ),

                                                                                        fluidRow(column(4, checkboxInput("saveH_heatcom", "Save the heatmap", TRUE)),
                                                                                                 conditionalPanel(condition = "input.saveH_heatcom",
                                                                                                                  column(4, textInput("fnameH_heatcom", "Type your file name", "My_heatmap")),
                                                                                                                  column(4, selectInput("formatH_heatcom", "Choose a format for your file",
                                                                                                                                        choices = c("png", "pdf"), selected = "png"))
                                                                                                                  )
                                                                                                 )
                                                                           )
                                                                           ),

                                                                           actionButton("getH_heatcom", "See heatmap", class = "btn-primary btn-lg"),
                                                                           tags$hr(),

                                                                           textOutput("diagl_heatcom"),
                                                                           tags$hr(),

                                                                           withSpinner(plotOutput("H_heatcom", height = "800px"), type = 6)
                                                          )
                                                 )
                                     )
                                     )
                            ),

                 navbarMenu("Gene Ontology",
                            tabPanel("STRING", value = "string",
                                     shinyjs::useShinyjs(),
                                     tags$style(HTML(".tabbable > .nav > li > a                  {background-color: #A1BAC8;  color:#FFFFFF}
                                                      .tabbable > .nav > li[class=active]    > a {background-color: #3C8DBC; color:#FFFFFF}
                                                      .tabbable > .nav > li    > a:hover {background-color: #3BAAE6; color:#FFFFFF}
                                                      .tabbable > .nav > li[class=active]    > a:hover {background-color: #3C8DBC; color:#FFFFFF}")
                                     ),
                                     fluidRow(style = "height:20px;"),

                                     tabsetPanel(type = "tabs",

                                                 tabPanel("STRING network",
                                                          h1(tags$u(class = "main-1", "Network and enrichment analysis from STRING")),
                                                          tags$hr(),

                                                          fluidRow(box(title = "Import your data and start the analysis", status = "primary",
                                                                       solidHeader = TRUE, collapsible = TRUE, width = 12,

                                                                       radioButtons("drug_stri", "Choose a dataset",
                                                                                    choices = c("Database" = "base",
                                                                                                "Your data" = "dat"),
                                                                                    selected = "base", inline = TRUE),
                                                                       conditionalPanel(condition = "input.drug_stri == 'base'",
                                                                                        uiOutput("drug2ui_stri"),
                                                                                        fluidRow(column(6, selectInput("cond_fhitB_stri", "Select some treatments to filter your proteins",
                                                                                                                       choices = NULL, multiple = TRUE)),
                                                                                                 column(6, selectInput("cat_fhitB_stri", "Select some categories to filter your proteins (If NULL, will select all)",
                                                                                                                       choices = c("CN", "NC", "CC", "ND", "NN"), multiple = TRUE))
                                                                                                 )
                                                                                        ),

                                                                       conditionalPanel(condition = "input.drug_stri == 'dat' ",
                                                                                        checkboxInput("impfile_stri", "Import a file", TRUE),
                                                                                        conditionalPanel(condition = "input.impfile_stri",
                                                                                                         fluidRow(column(4, fileInput("file_stri", "Choose a file")),
                                                                                                                  column(4, checkboxInput("ishit_stri", "Do you import a hitlist ? (needs a column named 'treatment')", TRUE),
                                                                                                                         conditionalPanel(condition = "!input.ishit_stri",
                                                                                                                                          textInput("idfile_stri", "What is the name of the column of
                                                                                                                                                    your file which contains the protein IDs ?")
                                                                                                                                          )
                                                                                                                         ),
                                                                                                                  conditionalPanel(condition = "input.ishit_stri",
                                                                                                                                   column(4, selectInput("cond_fhit_stri", "Select some treatments to filter your hits",
                                                                                                                                                         choices = NULL, multiple = TRUE))
                                                                                                                                   )
                                                                                                                  )
                                                                                                         ),
                                                                                        conditionalPanel(condition = "!input.impfile_stri",
                                                                                                         textInput("txt_stri", "Type some protein ID separated by a comma")
                                                                                                         )
                                                                                        ),

                                                                       conditionalPanel(condition = "output.file_stri_up | !input.impfile_stri",
                                                                                        fluidRow(column(6, selectInput("species_string", "Choose an organism",
                                                                                                                       choices = c("Human" = 9606,
                                                                                                                                   "Mouse" = 10090,
                                                                                                                                   "Rat" = 10116), selected = 9606)),
                                                                                                 column(6,  actionButton("start_string", "Start to map genes", class = "btn-primary btn-lg"))
                                                                                                 )
                                                                                        ),

                                                                       conditionalPanel(condition = "output.data_stri_up",
                                                                                        tags$hr(),
                                                                                        fluidRow(column(3, checkboxInput("intnet_stri", "Interactive network", FALSE)),
                                                                                                 column(3, checkboxInput("hidnet1_stri", "Hide network", FALSE)),
                                                                                                 column(3, selectInput("edgetype1_stri", "Meaning of network edges",
                                                                                                                       choices = c("evidence", "confidence", "actions"),
                                                                                                                       selected = "evidence")),
                                                                                                 column(3, numericInput("intscore1_stri", "Required interaction score",
                                                                                                                        min = 0, max = 1000, step = 100, value = 400))
                                                                                                 ),
                                                                                        fluidRow(column(3, actionButton("netbase_stri", "See network", class = "btn-primary btn-lg"))
                                                                                                 ),
                                                                                        tags$hr(),
                                                                                        actionButton("go_enrich", "Start the enrichment analysis", class = "btn-primary btn-lg"),
                                                                                        tags$hr(),

                                                                                        conditionalPanel(condition = "!input.hidnet1_stri",
                                                                                                         conditionalPanel(condition = "input.intnet_stri",
                                                                                                                          withSpinner(plotlyOutput("netInt_stri", height = "800px"), type = 6)
                                                                                                                          ),
                                                                                                         conditionalPanel(condition = "!input.intnet_stri",
                                                                                                                          withSpinner(plotOutput("net_stri", height = "800px"), type = 6),
                                                                                                                          fluidRow(column(2, tags$div(style="line-height:175%;",
                                                                                                                                                      tags$br()
                                                                                                                                                      ),
                                                                                                                                          downloadButton("downnet_stri", "Download plot")),
                                                                                                                                   column(2, selectInput("downnet_stri_format", "Download as", choices = c("png", "pdf")))
                                                                                                                                   )
                                                                                                                          )
                                                                                                         )
                                                                                        )
                                                                       )
                                                                   ),


                                                          tags$hr(),

                                                          conditionalPanel(condition = "output.enrich_stri_up",
                                                                           fluidRow(box(title = "Results from enrichment analysis", status = "primary",
                                                                                        solidHeader = TRUE, collapsible = TRUE, width = 12,

                                                                                        fluidRow(
                                                                                          column(6, selectInput("catego_stri", "Select a category for enrichment", choices = NULL)),
                                                                                          column(6, checkboxInput("hidtab_stri", "Hide the enrichment tab", FALSE))
                                                                                          ),

                                                                                        conditionalPanel(condition = "!input.hidtab_stri",
                                                                                                         DT::dataTableOutput("enrich_table_stri"),
                                                                                                         tags$hr(),
                                                                                                         downloadButton("downenrich_stri", "Download the tab as xlsx file")
                                                                                                         ),
                                                                                        conditionalPanel(condition = "output.enrich_res_tab_up",
                                                                                                         tags$hr(),
                                                                                                         fluidRow(
                                                                                                           column(3, selectInput("descri_stri", "Select a description to filter proteins", choices = NULL)),
                                                                                                           column(3, checkboxInput("hidnet2_stri", "Hide new network", FALSE)),
                                                                                                           column(3, selectInput("edgetype2_stri", "Meaning of network edges",
                                                                                                                                 choices = c("evidence", "confidence", "actions"),
                                                                                                                                 selected = "evidence")),
                                                                                                           column(3, numericInput("intscore2_stri", "Minimum interaction score",
                                                                                                                                  min = 0, max = 1000, step = 100, value = 400))
                                                                                                           ),
                                                                                                         fluidRow(column(6, actionButton("netfilt_stri", "See new network", class = "btn-primary btn-lg"))),
                                                                                                         tags$hr(),

                                                                                                         conditionalPanel(condition = "output.enrich_res_tab_up",
                                                                                                                          conditionalPanel(condition = "!input.hidnet2_stri",
                                                                                                                                           conditionalPanel(condition = "input.intnet_stri",
                                                                                                                                                            withSpinner(plotlyOutput("netInt2_stri", height = "800px"), type = 6)
                                                                                                                                                            ),
                                                                                                                                           conditionalPanel(condition = "!input.intnet_stri",
                                                                                                                                                            withSpinner(plotOutput("net2_stri", height = "800px"), type = 6),
                                                                                                                                                            fluidRow(column(2, tags$div(style="line-height:175%;",
                                                                                                                                                                                        tags$br()
                                                                                                                                                                                        ),
                                                                                                                                                                            downloadButton("downnetfilt_stri", "Download plot")),
                                                                                                                                                                     column(2, selectInput("downnetfilt_stri_format", "Download as", choices = c("png", "pdf")))
                                                                                                                                                                     )
                                                                                                                                                            )
                                                                                                                                           )
                                                                                                                          )
                                                                                                         )
                                                                                        )
                                                                                    )
                                                                           )
                                                          ),

                                                 tabPanel("Barplots network",
                                                          h1(tags$u(class = "main-1", "Interactive barplots network")),
                                                          tags$hr(),
                                                          HTML("<h5>In this tab, you can select some proteins and plot their STRING network.<br>
                                                               The network will be interactive and inside each node, their corresponding barplots
                                                               with the treatments you selected will be plotted. You can also color the node
                                                               according GO term from an enrichment analysis or any other category and color the
                                                               nodes border accoring their corresponding maximum fold change.</h5>"),
                                                          tags$hr(),

                                                          fluidRow(box(title = "Networks data parameters", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,

                                                                       radioButtons("drug_barnet", "Choose a dataset",
                                                                                    choices = c("Database" = "base",
                                                                                                "Your data" = "dat"),
                                                                                    selected = "base", inline = TRUE),
                                                                       conditionalPanel(condition = "input.drug_barnet == 'base'",
                                                                                        uiOutput("drug2_ui_barnet")
                                                                                        ),

                                                                       conditionalPanel(condition = "input.drug_barnet == 'dat' ",
                                                                                        fluidRow(column(6, fileInput("caldiff_barnet", "Upload your imprints caldiff file",
                                                                                                                     accept = c(".txt", ".csv", ".xlsx"))
                                                                                                        ),
                                                                                                 column(6, fileInput("hits_barnet", "Upload your hitlist, the summary file from the hitlist outputs",
                                                                                                                     accept = c(".txt", ".csv", ".xlsx"))
                                                                                                        )
                                                                                                 ),
                                                                                        ),
                                                                       tags$hr(),

                                                                       fluidRow(
                                                                         column(4, checkboxInput("importprot_barnet", "Import a protein list", FALSE),
                                                                                conditionalPanel(condition = "!input.importprot_barnet",
                                                                                                 checkboxInput("onlyhit_barnet", "Only take the hits proteins", FALSE),
                                                                                                 conditionalPanel(condition = "input.onlyhit_barnet",
                                                                                                                  selectInput("cond_fhit_barnet", "Select hits treatment", choices = NULL, multiple = TRUE))
                                                                                                 ),
                                                                                conditionalPanel(condition = "input.importprot_barnet",
                                                                                                 fileInput("prlist_file_barnet", "Import your protein list (txt file)", accept = ".txt")
                                                                                                 ),
                                                                                selectizeInput("prot_barnet", "Select some proteins (if NULL, will select all)", choices = NULL, multiple = TRUE)
                                                                                ),

                                                                         column(4, selectInput("condition_barnet", "Select one or more treatments (if NULL, selet all)",
                                                                                               choices = NULL,  multiple = TRUE)
                                                                                ),
                                                                         column(4, checkboxInput("importGO_barnet", "Import a file with GO terms", FALSE),
                                                                                conditionalPanel(condition = "input.importGO_barnet",
                                                                                                 HTML("<h5>This file must contain the column 'id' and the column 'GOterm'
                                                                                                      corresponding to the Uniprot IDs and some GO terms or any other terms
                                                                                                      respectively. Optionnally, it can contain a column 'color' to specify the
                                                                                                      node colors.</h5>"),
                                                                                                 fileInput("GOtermfile_barnet", "")
                                                                                                 ),
                                                                                conditionalPanel(condition = "!input.importGO_barnet",
                                                                                                 selectInput("GOtype_barnet", "Perform an enrichment analysis from",
                                                                                                             choices = c("none", "COMPARTMENTS", "Process", "Component",
                                                                                                                         "Function", "TISSUES", "Keyword", "KEGG", "SMART",
                                                                                                                         "PMID", "RCTM", "WikiPathways", "NetworkNeighborAL")
                                                                                                             )
                                                                                                 )
                                                                                )
                                                                         ),

                                                                       tags$hr(),
                                                                       fluidRow(column(4, checkboxInput("ch_own_col_barnet", "Choose your own color for the bar plots", FALSE)),
                                                                                column(4, conditionalPanel(condition = "input.ch_own_col_barnet",
                                                                                                           textOutput("n_cond_sel_barnet"),
                                                                                                           colourpicker::colourInput("own_color_pick_barnet", NULL, "#FF2B00",
                                                                                                                                     allowTransparent = TRUE, closeOnClick = TRUE),
                                                                                                           textOutput("own_color_barnet")
                                                                                                           )
                                                                                       ),
                                                                                conditionalPanel(condition = "input.ch_own_col_barnet",
                                                                                                 column(4, actionButton("add_col_barnet", "Add the color"),
                                                                                                        actionButton("rem_col_barnet", "Remove the last color"))
                                                                                                 )
                                                                                ),

                                                                       tags$hr(),

                                                                       fluidRow(column(4, selectInput("species_barnet", "Select a species",
                                                                                                      choices = c("human", "mouse", "rat"), selected = "human")),
                                                                                column(4, numericInput("reqscore_barnet", "Required interaction score",
                                                                                                       min = 200, max = 1000, value = 900, step = 10)),
                                                                                column(4, checkboxInput("werb_barnet", "Print error bar", TRUE))
                                                                                ),

                                                                       tags$hr(),

                                                                       fluidRow(column(3, checkboxInput("FCborder_barnet", "Color node border according maximum Fold-Change", TRUE)),
                                                                                conditionalPanel(condition = "input.FCborder_barnet",
                                                                                                 column(3, colourpicker::colourInput("FCbordercolorlow_barnet", NULL, "#0041FF",
                                                                                                                                     allowTransparent = TRUE, closeOnClick = TRUE)
                                                                                                        ),
                                                                                                 column(3, colourpicker::colourInput("FCbordercolormid_barnet", NULL, "#FFFFFF",
                                                                                                                                     allowTransparent = TRUE, closeOnClick = TRUE)
                                                                                                        ),
                                                                                                 column(3, colourpicker::colourInput("FCbordercolorhigh_barnet", NULL, "#FF0000",
                                                                                                                                     allowTransparent = TRUE, closeOnClick = TRUE)
                                                                                                        )
                                                                                                 )
                                                                                )
                                                                       )
                                                                   ),

                                                          tags$hr(),
                                                          actionButton("plotnet_barnet", "Plot network", class = "btn-primary btn-lg"),
                                                          tags$hr(),
                                                          shinyjs::useShinyjs(),
                                                          textOutput("diag_barnet"),
                                                          tags$br(),

                                                          conditionalPanel(condition = "output.netready_barnet",
                                                                           tags$hr(),

                                                                           fluidRow(
                                                                             column(2,
                                                                                    selectizeInput("select_groups_barnet", "Change color from", choices = NULL, size = 5),
                                                                                    colourpicker::colourInput("node_bordercolor_barnet", label = "Node border color", "#2B7CE9",
                                                                                                              allowTransparent = TRUE, closeOnClick = TRUE),
                                                                                    numericInput("node_size_barnet", "Node size", min = 1, max = 200,
                                                                                                 value = 40, step = 1),
                                                                                    colourpicker::colourInput("edge_color_barnet", label = "Edge color", "#00000075",
                                                                                                              allowTransparent = TRUE, closeOnClick = TRUE),
                                                                                    colourpicker::colourInput("font_color_barnet", label = "Label color", "#343434",
                                                                                                              allowTransparent = TRUE, closeOnClick = TRUE),
                                                                                    numericInput("label_size_barnet", "Label size", min = 1, max = 100, step = 0.5, value = 14),
                                                                                    selectInput("physics_type_barnet", "Physics",
                                                                                                choices = c('forceAtlas2Based', 'barnesHut', 'repulsion', 'hierarchicalRepulsion'),
                                                                                                selected = 'forceAtlas2Based'),
                                                                                    checkboxInput("button_enable_barnet", "Interaction button", value = TRUE),
                                                                                    tags$hr(),
                                                                                    downloadButton("down_barnet", "Save as html")
                                                                                    ),
                                                                             column(2,
                                                                                    colourpicker::colourInput("node_color_barnet", label = "Node color", "#D2E5FF",
                                                                                                              allowTransparent = TRUE, closeOnClick = TRUE),
                                                                                    column(12, style = "height:6px;"),
                                                                                    numericInput("node_borderwidth_barnet", "Node border width", min = 1, max = 50,
                                                                                                 value = 5, step = 0.5),
                                                                                    numericInput("node_imgpadding_barnet", "Node padding", min = 1, max = 50,
                                                                                                 value = 8, step = 0.5),
                                                                                    numericInput("edge_length_barnet", "Edge length", min = 1, max = 500,
                                                                                                 value = 60, step = 1),
                                                                                    colourpicker::colourInput("font_backcolor_barnet", label = "Label background color", "#A6A6A669",
                                                                                                              allowTransparent = TRUE, closeOnClick = TRUE),
                                                                                    column(12, style = "height:105px;"),
                                                                                    checkboxInput("physics_enable_barnet", "Enable physics", value = TRUE)
                                                                                    ),
                                                                             column(8,
                                                                                    HTML("<h5>Click+shift on a node to go to its Uniprot page</h5>"),
                                                                                    withSpinner(visNetworkOutput("network_barnet", height = "800px"), type = 6)
                                                                                    )
                                                                             ),

                                                                           tags$br(), tags$br(), tags$br(),
                                                                           )
                                                          )
                                                 )
                                     ),

                            tabPanel("ClusterProfiler", value = "clusprof",
                                     shinyjs::useShinyjs(),
                                     tags$style(HTML(".tabbable > .nav > li > a                  {background-color: #A1BAC8;  color:#FFFFFF}
                                                                  .tabbable > .nav > li[class=active]    > a {background-color: #3C8DBC; color:#FFFFFF}
                                                                  .tabbable > .nav > li    > a:hover {background-color: #3BAAE6; color:#FFFFFF}
                                                                  .tabbable > .nav > li[class=active]    > a:hover {background-color: #3C8DBC; color:#FFFFFF}")
                                     ),

                                     fluidRow(style = "height:20px;"),

                                     h1(tags$u(class = "main-1", "Enrichment analysis and visualization")),
                                     tags$hr(),

                                     fluidRow(box(title = "Import your data and start the analysis", status = "primary",
                                                  solidHeader = TRUE, collapsible = TRUE, width = 12,

                                                  radioButtons("drug_clus", h3("Choose a dataset"),
                                                               choices = c("Database" = "base",
                                                                           "Your data" = "dat"),
                                                               selected = "base", inline = TRUE),
                                                  conditionalPanel(condition = "input.drug_clus == 'base'",
                                                                   uiOutput("drug2ui_clus"),
                                                                   fluidRow(column(6, selectInput("cond_fhitB_clus", "Select some treatments to filter your proteins",
                                                                                                  choices = NULL, multiple = TRUE)),
                                                                            column(6, selectInput("cat_fhitB_clus", "Select some categories to filter your proteins (If NULL, will select all)",
                                                                                                  choices = c("CN", "NC", "CC", "ND", "NN"), multiple = TRUE))
                                                                            )
                                                                   ),

                                                  conditionalPanel(condition = "input.drug_clus == 'dat' ",
                                                                   fluidRow(column(4, fileInput("file_clus", "Choose a file")),
                                                                            column(4, checkboxInput("ishit_clus", "Do you import a hitlist ? (needs a column named 'treatment')", TRUE),
                                                                                   conditionalPanel(condition = "!input.ishit_clus",
                                                                                                    textInput("idfile_clus", "What is the name of the column of
                                                                                                              your file which contains the genes ?")
                                                                                                    )
                                                                                   ),
                                                                            conditionalPanel(condition = "input.ishit_clus",
                                                                                             column(4, selectInput("cond_fhit_clus", "Select some treatments to filter your hits",
                                                                                                                   choices = NULL, multiple = TRUE))
                                                                                             )
                                                                            )
                                                                   ),
                                                  fluidRow(column(3, selectInput("species_clus", "Specify the species from your data",
                                                                                 choices = c("Human", "Mouse"), selected = "Human")),
                                                           column(3, selectInput("database_clus", "Choose a database to perform the enrichment analysis",
                                                                                 choices = c("WikiPathway", "KEGG", "GO", "CETSA"), selected = "WikiPathway")),
                                                           column(3, numericInput("pvcut_clus", "Choose a p-value cutoff for gene set enrichment analysis",
                                                                                  value = 0.01, min = 0, max = 1, step = 0.01)),
                                                           column(3, numericInput("minGNsize_clus", "Choose a minimal size of genes annotated for testing",
                                                                                  value = 3, min = 1, max = 100, step = 1))
                                                           )

                                                  )
                                              ),

                                     tabsetPanel(type = "tabs",

                                                 tabPanel("Compare cluster",
                                                          fluidRow(box(title = "Compare the enrichment results of your data between treatments", status = "primary",
                                                                       solidHeader = TRUE, collapsible = TRUE, width = 12,
                                                                       fluidRow(column(6, sliderInput("npath_clus", "Choose the maximum number of pathway found to show",
                                                                                                      min = 1, max = 100, step = 1, value = 5)
                                                                                       ),
                                                                                column(6, actionButton("gocomp_clus", "Compare pathways", class = "btn-primary btn-lg"))
                                                                                ),
                                                                       DT::dataTableOutput("comptab_clus"),
                                                                       downloadButton("downcomptab_clus"),
                                                                       tags$hr(),
                                                                       withSpinner(plotOutput("compplot_clus", height = "800px"), type = 6),
                                                                       fluidRow(column(2, tags$div(style="line-height:175%;",
                                                                                                   tags$br()
                                                                                                   ),
                                                                                       downloadButton("downcomplot_clus", "Download plot")),
                                                                                column(2, selectInput("downcomplot_clus_format", "Download as", choices = c("png", "pdf")))
                                                                                )
                                                                       )
                                                                   )
                                                          ),
                                                 tabPanel("GSEA",
                                                          fluidRow(box(title = "Perform GSEA on your data", status = "primary",
                                                                       solidHeader = TRUE, collapsible = TRUE, width = 12,
                                                                       fluidRow(column(4, uiOutput("scorenameui_clus")),
                                                                                column(4, checkboxInput("onlypos_clus", "Show only enrcihment set with positive enrichment score", TRUE)),
                                                                                column(4, actionButton("gogsea_clus", "Start GSEA", class = "btn-primary btn-lg"))
                                                                                ),
                                                                       DT::dataTableOutput("gseatab_clus"),
                                                                       downloadButton("downgseatab_clus"),
                                                                       tags$hr(),
                                                                       withSpinner(plotOutput("gseaplot_clus", height = "800px"), type = 6),
                                                                       fluidRow(column(2, tags$div(style="line-height:175%;",
                                                                                                   tags$br()
                                                                                                   ),
                                                                                       downloadButton("downgsealot_clus", "Download plot")),
                                                                                column(2, selectInput("downgsealot_clus_format", "Download as", choices = c("png", "pdf")))
                                                                                )
                                                                       )
                                                                   )
                                                          ),
                                                 tabPanel("Gene concept network",
                                                          fluidRow(box(title = "Plot a gene concept network from your data", status = "primary",
                                                                       solidHeader = TRUE, collapsible = TRUE, width = 12,
                                                                       fluidRow(column(6, uiOutput("scorename2ui_clus")),
                                                                                column(6, actionButton("gogeneconc_clus", "See gene concept network", class = "btn-primary btn-lg"))
                                                                                ),
                                                                       withSpinner(plotOutput("geneplot_clus", height = "800px"), type = 6),
                                                                       fluidRow(column(2, tags$div(style="line-height:175%;",
                                                                                                   tags$br()
                                                                                                   ),
                                                                                       downloadButton("downgenelot_clus", "Download plot")),
                                                                                column(2, selectInput("downgenelot_clus_format", "Download as", choices = c("png", "pdf")))
                                                                                )
                                                                       )
                                                                   )
                                                          )
                                                 )
                                     )

                            ),

                 tabPanel("Cell", value = "cell",
                          fluidRow(style = "height:20px;"),
                          shinyjs::useShinyjs(),

                          h1(tags$u(class = "main-1", "Proteins localization")),
                          tags$hr(),
                          HTML("<h5>In this tab, you can assign protins to their subcellular location
                               from Protein Atlas database and then plot them on a cell.<br>
                               By clicking on the points on the plot, you  select a protein and then
                               can plot its IMPRINTS profile.</h5>"),
                          tags$hr(),
                          fluidRow(
                            box(title = "Get subcellular location from your hitlist", status = "primary",
                                solidHeader = TRUE, collapsible = TRUE, width = 12,

                                radioButtons("drug_cell", "Choose a dataset",
                                             choices = c("Database" = "base",
                                                         "Your data" = "dat"),
                                             selected = "base", inline = TRUE),
                                conditionalPanel(condition = "input.drug_cell == 'base'",
                                                 uiOutput("drug2ui_cell")
                                ),
                                fluidRow(column(6, selectInput("organism_cell", "Choose an organism",
                                                               choices = c("Human" = "HUMAN",
                                                                           "Mouse" = "MOUSE"), selected = "HUMAN")),
                                         conditionalPanel(condition = "input.drug_cell == 'dat' ",
                                                          column(6, fileInput("hitl_cell", "Import the summary file from the hitlist output"))
                                         )
                                ),

                                conditionalPanel(condition = "output.hitdata_cell_up",
                                                 fluidRow(column(4, selectInput("condhit_cell", "Select a treatment", choices = NULL)),
                                                          column(4, selectInput("cathit_cell", "Select some categories (if NULL, will select all)", choices = NULL, multiple = TRUE)),
                                                          column(4, actionButton("goloca_cell", "Get subcellular location", class = "btn-primary btn-lg"))
                                                 )
                                ),

                                textOutput("diagl_cell"),

                                tags$hr(),
                                conditionalPanel(condition = "output.resdata_cell_up",
                                                 DT::dataTableOutput("locatab_cell"),
                                                 downloadButton("down_prl_cell")
                                )

                            )
                          ),

                          conditionalPanel(condition = "output.resdata_cell_up",
                                           fluidRow(
                                             box(title = "The cell", status = "primary",
                                                 solidHeader = TRUE, collapsible = TRUE, width = 12,
                                                 fluidRow(column(4, textInput("titp_cell", "Type a title for the plot", "elutriation data in the cell")),
                                                          column(4, selectInput("condp_cell", "Select some treatments", multiple = TRUE, choices = NULL)),
                                                          column(4, actionButton("gop_cell", "See the plot", class = "btn-primary btn-lg"))
                                                 ),
                                                 tags$hr(),

                                                 withSpinner(plotlyOutput("cell_p", height = "700px"), type = 6),
                                                 downloadButton("downthe_cell", "Download the plot as html file")
                                             )
                                           )
                          ),

                          tags$hr(),

                          fluidRow(
                            box(title = "Bar plot", status = "primary",
                                solidHeader = TRUE, collapsible = TRUE, width = 12,
                                conditionalPanel(condition = "input.drug_cell == 'dat'",
                                                 fileInput("filebarp_cell", "If you want to see the bar plot from the protein you clicked on,
                                                                             please import the imprints_caldiff output file which correspond to your hitlist.")
                                ),
                                conditionalPanel(condition = "!input.selpr_loca_cell",
                                                 htmlOutput("prsel_p_cell")
                                ),

                                conditionalPanel(condition = "output.barpdata_cell_up | input.drug_cell == 'base'",
                                                 fluidRow(column(4, checkboxInput("selpr_loca_cell", "Select proteins according to their subcellular location", FALSE)),
                                                          conditionalPanel(condition = "input.selpr_loca_cell",
                                                                           column(4, selectInput("selorga_cell", "Select some organelles", multiple = TRUE, choices = NULL)),
                                                                           column(4, checkboxInput("allpr_cell", "Select all the proteins from the organelles you selected", FALSE),
                                                                                  conditionalPanel(condition = "!input.allpr_cell",
                                                                                                   selectizeInput("selectpr_cell", "Select some proteins", multiple = TRUE, choices = NULL)
                                                                                  )
                                                                           )
                                                          )
                                                 ),

                                                 radioButtons("cond_sel_cell", "Selection type",
                                                              choices = c("Select the treatment level" = "treat",
                                                                          "Select treatment level by category" = "cat",
                                                                          "Select all the treatment level" = "all_cond"),
                                                              selected = "treat", inline = TRUE
                                                 ),

                                                 conditionalPanel(condition = "input.cond_sel_cell != 'all_cond' ",
                                                                  selectInput("cond_cell", "Select one or more treatments",
                                                                              choices = NULL,
                                                                              multiple = TRUE)
                                                 ),

                                                 tags$hr(),

                                                 fluidRow(column(4, checkboxInput("save_bar_cell", "Save the bar plots in a pdf file", FALSE),
                                                                 conditionalPanel(condition = "input.save_bar_cell",
                                                                                  textInput("pdftit_cell", "Choose a name for your pdf file", "barplot")
                                                                                  )
                                                                 ),
                                                          conditionalPanel(condition = "input.save_bar_cell",
                                                                           column(4, numericInput("lay_bar1_cell", "Type the number of plot per row",
                                                                                                  min = 1, max = 10, step = 1, value = 4),
                                                                                  numericInput("lay_bar2_cell", "Type the number of plot per column",
                                                                                               min = 1, max = 10, step = 1, value = 3)),
                                                                           column(4, numericInput("pdfw_cell", "Type the width of the pdf page",
                                                                                                  min = 1, step = 1, value = 12),
                                                                                  numericInput("pdfh_cell", "Type the height of the pdf page",
                                                                                               min = 1, step = 1, value = 12))
                                                                           )
                                                 ),

                                                 tags$hr(),
                                                 fluidRow(column(4, checkboxInput("ch_own_col_cell", "Choose your own color", FALSE)),
                                                          column(4, conditionalPanel(condition = "input.ch_own_col_cell",
                                                                                     textOutput("n_cond_sel_cell"),
                                                                                     colourpicker::colourInput("own_color_pick_cell", NULL, "#FF2B00",
                                                                                                 allowTransparent = TRUE, closeOnClick = TRUE),
                                                                                     textOutput("own_color_cell")
                                                          )
                                                          ),
                                                          conditionalPanel(condition = "input.ch_own_col_cell",
                                                                           column(4, actionButton("add_col_cell", "Add the color"),
                                                                                  actionButton("rem_col_cell", "Remove the last color"))
                                                          )
                                                 ),
                                                 tags$hr(),

                                                 fluidRow(column(4, checkboxInput("werb_cell", "Print error bar", TRUE),
                                                                    checkboxInput("wpts_cell", "Print point of each replicate", FALSE)),
                                                          column(4, checkboxInput("grad_cell", "Use color gradient", FALSE)),
                                                          column(4, checkboxInput("line_cell", "Use line instead of bar", FALSE))
                                                 ),

                                                 tags$hr(),

                                                 actionButton("barp_cell", "See bar plot", class = "btn-primary btn-lg"),
                                                 tags$hr(),
                                                 textOutput("diag_bar_cell"),
                                                 tags$hr(),

                                                 withSpinner(plotOutput("bar_pr_cell", height = "800px"), type = 6),
                                                 fluidRow(column(2, tags$div(style="line-height:175%;",
                                                                             tags$br()
                                                                             ),
                                                                 downloadButton("downbar_cell", "Download plot")),
                                                          column(2, selectInput("downbar_cell_format", "Download as", choices = c("png", "pdf")))
                                                 )
                                )
                            )
                          )
                          ),

                 tabPanel("PubMed search", value = "pubmed",
                          fluidRow(style = "height:20px;"),
                          shinyjs::useShinyjs(),

                          h1(tags$u(class = "main-1", "Search publications in PubMed")),
                          tags$hr(),
                          HTML("<h5>In this tab, you can look for potential pubmed publication
                               related to the keywords you selected. <br>
                               If a/some publicaitons are found, their title, author and abstract are saved in
                               one word file in the folder you named.</h5>"),
                          tags$hr(),

                          fluidRow(box(title = "Search parameters", status = "primary",
                                       solidHeader = TRUE, collapsible = TRUE, width = 12,
                                       fluidRow(column(3, checkboxInput("impc_pubmed", "Import a file", TRUE),
                                                       conditionalPanel(condition = "input.impc_pubmed",
                                                                        fileInput("data_pubmed", "Import your data")),
                                                       conditionalPanel(condition = "!input.impc_pubmed",
                                                                        textInput("dtext_pubmed", "Type some protein names or any other words, separated by a comma", "proteomics"))
                                                       ),
                                       column(3, textInput("feat_pubmed", "Type your second research word (can be null)", "")),
                                       column(3, textInput("LA_pubmed", "Type a language to match (can be null)", "english")),
                                       column(3, textInput("Y_pubmed", "Type a year range to match (can be null, format is Y1:Y2)", "2022:2023"))
                                       ),

                                       fluidRow(column(3, textInput("api_pubmed", "Type your NCBI API if you have an account")),
                                                column(3, textInput("fname_pubmed", "Type the name of the folder that will be created", "pubmed_search")),
                                                column(3, conditionalPanel(condition = "input.impc_pubmed",
                                                                           checkboxInput("hit_pubmed", "Do you import a hitlist ? (need description column)", TRUE))
                                                       ),
                                                conditionalPanel(condition = "input.hit_pubmed",
                                                                 column(3, textInput("cond_pubmed", "Type a treatment from you hitlist (if null, will take all the treatments)"))
                                                                 )
                                                ),
                                       fluidRow(column(3, checkboxInput("save_in_word", "Save publication found for each query in
                                                                        a word file (title, authors and abstract)", TRUE))
                                                ),
                                       tags$hr(),

                                       conditionalPanel(condition = "output.pubmed_fileup | !input.impc_pubmed",
                                                        actionButton("go_pub", "Start searching", class = "btn-primary btn-lg"),
                                                        tags$hr()
                                                        ),

                                       textOutput("diag")
                                       )
                                   ),

                          DT::dataTableOutput("pubmed_out"),

                          downloadButton("down_pubmed"),

                          tags$hr()
                          ),

                 bslib::nav_item(tags$a(href = "https://github.com/mgerault/IMPRINTS.CETSA.app",
                        icon("github"),
                        title = "See source code to the github repository"), class = "icon1"),
                 bslib::nav_item(tags$a(href = "https://youtu.be/djpP8nc_JUE",
                        icon("question-circle"),
                        title = "See the tutorial video of the app"), class = "icon2"),
                 bslib::nav_item(tags$a(href = "mailto:marco.gerault@gmail.com",
                        icon("envelope"),
                        title = "Any questions, suggestions or bug report ? Feel free to send me an e-mail !"), class = "icon3")


)

server <- function(input, output, session){
  setwd(WD)

  ### analysis tab - Peptides

  # PD peptides files
  pep_file_data <- reactive({
    File <- input$PD_data_pep
    if (is.null(File)){
      return(NULL)
    }
    File
  })
  #check if a files are uploaded
  output$pep_fileup <- reactive({
    return(!is.null(pep_file_data()))
  })
  outputOptions(output, "pep_fileup", suspendWhenHidden = FALSE)

  # temperatures
  output$temp_nameui_pep <- renderUI({
    if(!is.null(pep_file_data())){
      files_temp <- pep_file_data()$name
    }
    else{
      files_temp <- NULL
    }
    m <- matrix("", length(files_temp), 1,
                dimnames = list(files_temp, "Temperatures"))

    matrixInput("temp_name_pep", "Type the name you want for your temperatures",
                value = m,
                rows = list(names = TRUE),
                cols = list(names = TRUE)
                )
  })

  # treatments
  output$treat_nameui_pep <- renderUI({
    if(!is.null(pep_file_data())){
      TMT <- colnames(readr::read_tsv(pep_file_data()$datapath[1], n_max = 0, progress = F, show_col_types = F)) # only read header
      TMT <- unique(unlist(stringr::str_extract_all(TMT, "(?<=: )\\d{3}[C|N](?=,)"))) # extract TMT channels --> not 126
      TMT <- c("126", TMT)
      if(length(TMT) == 9){
        TMT <- c(TMT, "131")
      }
    }
    else{
      TMT <- NULL
    }

    m <- matrix("", length(TMT), 1,
                dimnames = list(TMT, "Treatment"))

    matrixInput("treat_name_pep", "Type the name of your channels",
                value = m,
                rows = list(names = TRUE),
                cols = list(names = TRUE)
    )
  })

  # protein file
  prot_data_pep <- reactive({
    File <- input$prot_data_pep
    if (is.null(File)){
      return(NULL)
    }
    readr::read_tsv(File$datapath)
  })

  # read the data
  pep_data <- reactiveValues(
    x = NULL
  )
  observeEvent(input$file_data_pep,{
    if(input$got_data_pep){
      File <- input$file_data_pep
      if(!is.null(File)){
        pep_data$x <- readr::read_tsv(File$datapath)
      }
    }
  })
  observeEvent(input$read_pep, {
    df <- NULL
    treat <- as.character(input$treat_name_pep[,1])
    temp <- as.character(input$temp_name_pep[,1])
    if(any(stringr::str_length(treat) == 0)){
      showNotification("Type the treatment names !", type = "error")
      return(NULL)
    }
    else if(any(stringr::str_length(temp) == 0)){
      showNotification("Type the temperature names !", type = "error")
      return(NULL)
    }
    else{
      showNotification("Reading files...", type = "message")
      df <- imprints_read_peptides(pep_file_data()$datapath, treatment = treat,
                                   temperatures = temp,
                                   proteins = prot_data_pep(),
                                   dataset_name = input$dname_pep)
    }
    pep_data$x <- df
  })

  # check if pep data available
  output$pep_dataup <- reactive({
    return(!is.null(pep_data$x))
  })
  outputOptions(output, "pep_dataup", suspendWhenHidden = FALSE)

  # see the peptides data
  observeEvent(input$see1_pep,{
    if(!is.null(pep_data$x)){
      showModal(tags$div(id="modal1_pep", modalDialog(
        DT::renderDataTable({DT::datatable(pep_data$x,
                                           caption = htmltools::tags$caption(
                                             style = 'caption-side: top; text-align: left;',
                                             htmltools::strong("Peptides data")
                                           ),
                                           rownames = FALSE,
                                           options = list(lengthMenu = c(10,20,30), pageLength = 10,
                                                          scrollX = TRUE))
        }),
        footer = NULL,
        easyClose = TRUE
      )))
    }
    else{
      showNotification("The data are currently NULL, try to refresh.", type = "error")
    }
  })


  # normalization
  observeEvent(input$step_peptides, {
    js$collapse("upload_peptide")
    updateCheckboxInput(session, "got_norm_pep", value = input$step_peptides)
  }, ignoreInit = TRUE)

  norm_pep_data <- reactiveValues(
    x = NULL
  )
  observeEvent(input$normfile_pep,{
    if(input$got_norm_pep){
      File <- input$normfile_pep
      if(!is.null(File)){
        norm_pep_data$x <- readr::read_tsv(File$datapath)
      }
    }
  })
  observeEvent(input$NORM_pep, {
    df <- NULL

    showNotification("Starting normalization, this may take a while", type = "message")
    df <- imprints_normalize_peptides(pep_data$x)

    norm_pep_data$x <- df
    showNotification("Normalized data saved !",  type = "message")
  })

  # check if norm pep data available
  output$norm_pep_dataup <- reactive({
    return(!is.null(norm_pep_data$x))
  })
  outputOptions(output, "norm_pep_dataup", suspendWhenHidden = FALSE)

  # see the norm peptides data
  observeEvent(input$see2_pep,{
    if(!is.null(norm_pep_data$x)){
      showModal(tags$div(id="modal2_pep", modalDialog(
        DT::renderDataTable({DT::datatable(norm_pep_data$x,
                                           caption = htmltools::tags$caption(
                                             style = 'caption-side: top; text-align: left;',
                                             htmltools::strong("Normalized peptides data")
                                           ),
                                           rownames = FALSE,
                                           options = list(lengthMenu = c(10,20,30), pageLength = 10,
                                                          scrollX = TRUE))
        }),
        footer = NULL,
        easyClose = TRUE
      )))
    }
    else{
      showNotification("The data are currently NULL, try to refresh.", type = "error")
    }
  })


  ## sequence function
  observe({
    updateSelectInput(session, "control_pep", choices = get_treat_level(norm_pep_data$x))
  })

  protseq_file_pep <- reactive({
    File <- input$protseq_file_pep
    if (is.null(File)){
      return(NULL)
    }
    df <- rio::import(File$datapath, header = TRUE)
    if(!("protein" %in% colnames(df))){
      showNotification("Your file needs to contain the protein column !", type = "error")
      return(NULL)
    }
    else if(!("sequence" %in% colnames(df))){
      showNotification("Your file doesn't contain the sequence column, only protein are kept.", type = "warning")
      df <- df[,"protein", drop = FALSE]
    }
    else{
      df <- df[,c("protein", "sequence")]
    }
    df
  })
  observe({
    updateSelectizeInput(session, "protseq_pep", choices = norm_pep_data$x$`Master Protein Accessions`, server = TRUE)
  })

  # handling sequence selection
  output$selectSequenceui_pep <- renderUI({
    if(!is.null(input$protseq_pep)){
      if(length(input$protseq_pep) <= 20){
        prot <- input$protseq_pep
        m <- matrix("", length(prot), 1,
                    dimnames = list(prot, "sequence"))

        matrixInput("selectSequence_pep", "Type the sequences (i.e. peptide position) you want to highlight.
                                           Press enter or click outside the table when you're done.",
                    value = m,
                    rows = list(names = TRUE),
                    cols = list(names = TRUE)
                    )
      }
      else{
        shiny::HTML("<h5>If you want to select a specific sequence for more than 20 proteins,
                         you need to import a file (check box on the top).</h5>")
      }
    }
    else{
      textInput("selectSequence_pep", "Type the sequences (i.e. peptide position) you want to highlight.
                                       It will be applied for all proteins.")
    }
  })

  sequence_pep_data <- reactiveValues(
    x = NULL
  )
  observeEvent(input$FCfile_pep,{
    if(input$got_FCfile_pep){
      File <- input$FCfile_pep
      if(!is.null(File)){
        sequence_pep_data$x <- readr::read_tsv(File$datapath)
      }
    }
  })
  observeEvent(input$SEQU_pep, {
    prot <- NULL
    sequ <- NULL
    if(input$sequence_file){
      if(!is.null(protseq_file_pep())){
        prot <- protseq_file_pep()$protein
        sequ <- protseq_file_pep()$sequence
      }
      else{
        return(NULL)
      }
    }
    else{
      prot <- input$protseq_pep
      if(!is.null(input$selectSequence_pep)){
        if(inherits(input$selectSequence_pep, "matrix")){
          sequ <- as.character(input$selectSequence_pep[,1])
        }
        else{
          if(stringr::str_detect(input$selectSequence_pep, "^\\d{1,}-\\d{1,}$") & stringr::str_length(input$selectSequence_pep) == 0){
            sequ <- input$selectSequence_pep
          }
          else{
            showNotification("The sequence you wrote isn't in the right format.
                              No sequence has been selected.", duration = 8, type = "warning")
          }
        }
      }
    }

    withCallingHandlers({
      shinyjs::html("diag_pep_sequence", "")
      sequence_pep_data$x <- imprints_sequence_peptides(norm_pep_data$x,
                                                        proteins = prot, sequence = sequ,
                                                        control = input$control_pep,
                                                        dataset_name = input$dnamediff_pep)
      showNotification("Fold change and bar plot saved !",  type = "message")
      },
      message = function(m) {
        shinyjs::html(id = "diag_pep_sequence", html = paste(m$message, "<br>", sep = ""), add = FALSE)
        }
      )
  })

  # check if FC pep data available
  output$sequence_pep_dataup <- reactive({
    return(!is.null(sequence_pep_data$x))
  })
  outputOptions(output, "sequence_pep_dataup", suspendWhenHidden = FALSE)


  # potential cleaved
  observe({
    updateSelectInput(session, "controlcleaved_pep", choices = get_treat_level(sequence_pep_data$x))
  })

  cleaved_pep_data <- reactiveValues(
    x = NULL
  )
  observeEvent(input$CLEAVED_pep, {
    showNotification("Searching for cleaved sites", type = "message")
    n_prot <- length(unique(sequence_pep_data$x$`Master Protein Accessions`))

    cleaved_pep_data$x <- imprints_cleaved_peptides(sequence_pep_data$x,
                                                    R2 = input$R2cleaved_pep,
                                                    control = input$controlcleaved_pep,
                                                    min_ValidValue = input$propValcleaved_pep)
    cleaved_pep_data$x <- cleaved_pep_data$x %>% dplyr::ungroup()

    showModal(
      modalDialog(
        selectInput("conditioncleaved_pep", "Choose a treatment",
                    choices = unique(cleaved_pep_data$x$treatment),
                    selected = unique(cleaved_pep_data$x$treatment)[1]),
        DT::renderDataTable({DT::datatable(cleaved_pep_data$x[cleaved_pep_data$x$treatment %in% input$conditioncleaved_pep,],
                                           caption = htmltools::tags$caption(
                                             style = 'caption-side: top; text-align: left;',
                                             htmltools::strong(paste("Potentially cleaved -", input$conditioncleaved_pep))
                                           ),
                                           rownames = FALSE,
                                           options = list(lengthMenu = c(10,20,30), pageLength = 10,
                                                          scrollX = TRUE))
          }),
        tags$hr(),
        textOutput("diag_pep_cleaved"),

        title = tags$strong("Do you want to compute and plot fold change from the potentially cleaved proteins ?"),
        footer = tagList(actionButton("confirm_cleavedplot", tags$strong("Yes"), class = "btn-lg btn-success"),
                         modalButton(tags$strong("No"))
        )
      )
    )
  })
  observeEvent(input$confirm_cleavedplot, {
    showNotification("Start computing and plotting fold change", type = "message")
    withCallingHandlers({
      shinyjs::html("diag_pep_cleaved", "")
      cleaved_pepTab <- cleaved_pep_data$x %>% dplyr::filter(treatment == input$conditioncleaved_pep)
      foo <- imprints_sequence_peptides(norm_pep_data$x,
                                                        proteins = cleaved_pepTab$protein,
                                                        sequence = cleaved_pepTab$cleaved_site,
                                                        control = input$controlcleaved_pep,
                                                        dataset_name = "potentially_cleaved")
      },
      message = function(m) {
        shinyjs::html(id = "diag_pep_cleaved", html = paste(m$message, "<br>", sep = ""), add = FALSE)
        }
    )

    removeModal()
    showNotification("Fold change computed !", type = "message")
  })


  # filter peptides data
  info_filterpep <- reactiveValues(
    name = NULL
  )
  tofilter_pep_data <- reactive({
    File <- input$filter_joinpep
    if (is.null(File)){
      return(NULL)
    }
    info_filterpep$name <- File$name
    readr::read_tsv(File$datapath)
  })


  protseq_file_joinpep <- reactive({
    File <- input$protseq_file_joinpep
    if (is.null(File)){
      return(NULL)
    }
    df <- rio::import(File$datapath, header = TRUE)
    if(!("protein" %in% colnames(df)) | !("sequence" %in% colnames(df))){
      showNotification("Your file needs to contain the protein column and the sequence column !", type = "error")
      return(NULL)
    }
    else{
      df <- df[,c("protein", "sequence")]
    }
    df
  })
  observe({
    updateSelectizeInput(session, "protseq_joinpep",
                         choices = unique(stringr::str_remove_all(tofilter_pep_data()$`Master Protein Accessions`, "\\s.*")),
                         server = TRUE)
  })
  observe({
    updateSelectInput(session, "remcond_joinpep", choices = get_treat_level(tofilter_pep_data()))
  })

  # handling sequence selection
  output$selectSequenceui_joinpep <- renderUI({
    if(!is.null(input$protseq_joinpep)){
      if(length(input$protseq_joinpep) <= 20){
        prot <- input$protseq_joinpep
        m <- matrix("", length(prot), 1,
                    dimnames = list(prot, "sequence"))

        matrixInput("selectSequence_joinpep", "Type the sequences (i.e. peptide position) you want to highlight.
                                               Press enter or click outside the table when you're done.",
                    value = m,
                    rows = list(names = TRUE),
                    cols = list(names = TRUE)
        )
      }
      else{
        shiny::HTML("<h5>If you want to select a specific sequence for more than 20 proteins,
                         you need to import a file (check box on the top).</h5>")
      }
    }
    else{
      textInput("selectSequence_joinpep", "Type the sequences (i.e. peptide position) you want to highlight.
                                       It will be applied for all proteins.")
    }
  })

  observeEvent(input$gofilter_joinpep, {
    prot <- NULL
    sequ <- NULL
    if(input$sequence_file_joinpep){
      if(!is.null(protseq_file_joinpep())){
        prot <- protseq_file_joinpep()$protein
        sequ <- protseq_file_joinpep()$sequence
      }
      else{
        return(NULL)
      }
    }
    else{
      prot <- input$protseq_joinpep
      if(!is.null(input$selectSequence_joinpep)){
        if(inherits(input$selectSequence_joinpep, "matrix")){
          sequ <- as.character(input$selectSequence_joinpep[,1])
        }
        else{
          if(stringr::str_detect(input$selectSequence_joinpep, "^\\d{1,}-\\d{1,}$") & stringr::str_length(input$selectSequence_joinpep) == 0){
            sequ <- input$selectSequence_joinpep
          }
          else{
            showNotification("The sequence you wrote isn't in the right format.
                              No sequence has been selected.", duration = 8, type = "warning")
          }
        }
      }
    }

    df_filtered <- tofilter_pep_data()
    withCallingHandlers({
      shinyjs::html("diag_pep_filter", "")
      if(!is.null(prot) & !is.null(sequence)){
        message("Rmoving specific peptides")
        df_filtered <- imprints_remove_peptides(tofilter_pep_data(),
                                                proteins = prot,
                                                sequence = sequ)
      }
      if(!is.null(input$remcond_joinpep)){
        message("Removing treatments")
        df_filtered <- df_filtered[,-stringr::str_which(colnames(df_filtered), paste0("_", input$remcond_joinpep,
                                                                                      "$", collapse = "|")
                                                        )
                                   ]
      }
      message("Saving filtered data")
      f_name <- stringr::str_replace(info_filterpep$name, "\\.txt", "_filtered.txt")
      f_name <- str_replace_all(f_name, "\\d{6}_\\d{4}_", format(Sys.time(), "%y%m%d_%H%M_"))
      readr::write_tsv(df_filtered, file = f_name)
      message("Filtered data saved !")
      showNotification("Filtered data saved !",  type = "message")
      },
      message = function(m) {
        shinyjs::html(id = "diag_pep_filter", html = paste(m$message, "<br>", sep = ""), add = FALSE)
        }
    )
  })

  # join peptides datasets
  tojoin_pep_data <- reactive({
    File <- input$joinFC_file_pep
    if (is.null(File)){
      return(NULL)
    }
    File$datapath
  })
  # check if all files are uploaded
  output$tojoin_pep_dataup <- reactive({
    return(!is.null(tojoin_pep_data()))
  })
  outputOptions(output, "tojoin_pep_dataup", suspendWhenHidden = FALSE)

  joined_pep_data <- reactiveValues(
    x = NULL
  )
  observeEvent(input$JOIN_pep, {
    withCallingHandlers({
      shinyjs::html("diag_pep_join", "")
      message("Reading and joining data")
      tojoin <- lapply(tojoin_pep_data(), readr::read_tsv)
      joined_pep_data$x <- imprints_join_peptides(tojoin)

      message("Saving joined dataset")
      readr::write_tsv(joined_pep_data$x,
                       file = paste0(format(Sys.time(), "%y%m%d_%H%M_"),
                                     "JoinedPeptides.txt")
                       )
      message("Joined data saved !")
      showNotification("Joined data saved !",  type = "message")
    },
    message = function(m) {
      shinyjs::html(id = "diag_pep_join", html = paste(m$message, "<br>", sep = ""), add = FALSE)
    }
    )
  })

  # see the joined peptides data
  observeEvent(input$see3_pep,{
    if(!is.null(joined_pep_data$x)){
      showModal(tags$div(id="modal3_pep", modalDialog(
        DT::renderDataTable({DT::datatable(joined_pep_data$x,
                                           caption = htmltools::tags$caption(
                                             style = 'caption-side: top; text-align: left;',
                                             htmltools::strong("Joined peptides data")
                                           ),
                                           rownames = FALSE,
                                           options = list(lengthMenu = c(10,20,30), pageLength = 10,
                                                          scrollX = TRUE))
        }),
        footer = NULL,
        easyClose = TRUE
      )))
    }
    else{
      showNotification("The data are currently NULL, try to refresh.", type = "error")
    }
  })


  # plot joined data
  toplot_pep_data <- reactive({
    if(input$join_data_pep == "join_app"){
      return(joined_pep_data$x)
    }
    else if(input$join_data_pep == "join_file"){
      File <- input$joined_file_pep
      if (is.null(File)){
        return(NULL)
      }
      else{
        df <- readr::read_tsv(File$datapath)
        return(df)
      }
    }
  })
  # check if data available
  output$toplot_pep_dataup <- reactive({
    return(!is.null(toplot_pep_data()))
  })
  outputOptions(output, "toplot_pep_dataup", suspendWhenHidden = FALSE)

  observe({
    updateSelectInput(session, "condition_plotjoinpep", choices = get_treat_level(toplot_pep_data()))
  })

  # handling color selection
  output$n_cond_sel_plotjoinpep <- renderText({
    if(input$ch_own_col_plotjoinpep){
      paste("You selected", length(input$condition_plotjoinpep), "treatments, please enter the same number of colors")
    }
    else{
      NULL
    }
  })

  OWN_color_plotjoinpep <- reactiveValues(
    ch = c()
  )
  observeEvent(input$add_col_plotjoinpep, {
    OWN_color_plotjoinpep$ch <- append(OWN_color_plotjoinpep$ch, input$own_color_pick_plotjoinpep)
  })
  observeEvent(input$rem_col_plotjoinpep, {
    if(length(OWN_color_plotjoinpep$ch) <= 1){
      OWN_color_plotjoinpep$ch <- c()
    }
    else{
      OWN_color_plotjoinpep$ch <- OWN_color_plotjoinpep$ch[1:(length(OWN_color_plotjoinpep$ch)-1)]
    }
  })
  output$own_color_plotjoinpep <- renderText({
    paste("You selected this colors :", paste(OWN_color_plotjoinpep$ch, collapse = ", "))
  })


  # plot peptides bar plot !
  observeEvent(input$getbar_plotjoinpep, {
    data_toplot <- toplot_pep_data()
    data_toplot$`Master Protein Accessions` <- paste(data_toplot$`Positions in Master Proteins`, "\n", "\n")
    colnames(data_toplot)[1:5] <- c("id", "description", "sumUniPeps", "sumPSMs", "countNum")

    withCallingHandlers({
      shinyjs::html("diag_bar_plotjoinpep", "")
      if(length(input$condition_plotjoinpep)){
        if(input$ch_own_col_plotjoinpep){
          nbc <- length(input$condition_plotjoinpep)
          COL <- OWN_color_plotjoinpep$ch
          if(nbc == length(COL)){
            imprints_barplotting_app(data_toplot, save_pdf = TRUE, ret_plot = FALSE,
                                    colorpanel = COL, treatmentlevel = input$condition_plotjoinpep,
                                    layout = c(input$lay_bar1_plotjoinpep, input$lay_bar2_plotjoinpep),
                                    pdfname = input$pdftit_plotjoinpep)
            showNotification("Bar plot saved !",  type = "message")
          }
          else{
            showNotification("The number of colors given doesn't match the number of treatment selected !", type = "error")
          }
        }
        else{
          imprints_barplotting_app(data_toplot, save_pdf = TRUE, ret_plot = FALSE,
                                  treatmentlevel = input$condition_plotjoinpep,
                                  layout = c(input$lay_bar1_plotjoinpep, input$lay_bar2_plotjoinpep),
                                  pdfname = input$pdftit_plotjoinpep)
          showNotification("Bar plot saved !",  type = "message")
        }
      }
      else{
        showNotification("Don't forget to select some treatments !", type = "error")
      }
    },
    message = function(m) {
      shinyjs::html(id = "diag_bar_plotjoinpep", html = paste(m$message, "<br>", sep = ""), add = FALSE)
      }
    )
  })




  ### analysis tab - Proteins
  output$treat_nameui <- renderUI({
    TMT <- list("10" = c("126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131"),
                "11" = c("126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C"),
                "16" = c("126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N"),
                "18" = c("126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N", "134C", "135N")
    )
    m <- matrix("", as.numeric(input$n_chan), 1,
                dimnames = list(c(TMT[[input$n_chan]]), "Treatment"))

    matrixInput("treat_name", "Type the name of your channels",
                value = m,
                rows = list(names = TRUE),
                cols = list(names = TRUE)
    )
  })

  cetsa_data <- reactive({
    File <- input$PD_data
    if (is.null(File) | is.null(input$treat_name)){
      return(NULL)
    }
    else if(any(apply(input$treat_name, 1, function(x) stringr::str_length(x) == 0))){
      return(NULL)
    }
    else if(length(unique(as.character(input$treat_name[,1]))) != nrow(input$treat_name)){
      return(NULL)
    }

    withCallingHandlers({
      shinyjs::html("diag_rawread", "")
      df <- imprints_rawread(File$datapath,
                             #the name of each treatment
                             treatment = as.character(input$treat_name[,1]),
                             #number of channel and the name of it
                             nread = as.numeric(input$n_chan),
                             channels = rownames(input$treat_name)
                             )
    },
    message = function(m) {
      m <- m$message
      if(grepl('protein', m)){
        shinyjs::html(id = "diag_rawread",
                      html = paste0("<span style='color:red;'>", m, "</span><br>"),
                      add = TRUE)
      }
    }
    )

    df
  })

  #check if a file is upload
  output$cetsa_fileup <- reactive({
    return(!is.null(cetsa_data()))
  })
  outputOptions(output, "cetsa_fileup", suspendWhenHidden = FALSE)

  observeEvent(input$see1_cetsa,{
    if(!is.null(cetsa_data())){
      showModal(tags$div(id="modal1", modalDialog(
        DT::renderDataTable({DT::datatable(cetsa_data(),
                                           caption = htmltools::tags$caption(
                                             style = 'caption-side: top; text-align: left;',
                                             htmltools::strong("Base data")
                                           ),
                                           rownames = FALSE,
                                           options = list(lengthMenu = c(10,20,30), pageLength = 10,
                                                          scrollX = TRUE))
          }),
        footer = NULL,
        easyClose = TRUE
      )))
    }
    else{
      showNotification("The data are currently NULL, try to refresh.", type = "error")
    }
  })


  output$temp_nameui <- renderUI({
    if(!is.null(cetsa_data())){
      temp <- unique(cetsa_data()$condition)
    }
    else{
      temp <- NULL
    }
    m <- matrix("", length(temp), 1,
                dimnames = list(temp, "Temperatures"))

    matrixInput("temp_name", "Type the name you want for your temperatures",
                value = m,
                rows = list(names = TRUE),
                cols = list(names = TRUE)
    )
  })
  cetsa_data_clean <- eventReactive(input$str_ren, {
    d1 <- NULL
    temp <- as.character(input$temp_name[,1])
    if(any(sapply(temp, stringr::str_length) == 0)){
      showNotification("Some temperatures names are empty !", type = "error", duration = 5)
    }
    else if(length(unique(temp)) != length(temp)){
      showNotification("Some temperatures names are equal !", type = "error", duration = 5)
    }
    else{
      d1 <- ms_conditionrename(cetsa_data(),
                               incondition = unique(cetsa_data()$condition),
                               outcondition = temp
                               )
      showNotification("Renaming done !", type = "message", duration = 3)

      if(input$rem_mix){
        if(length(grep("^Mix", names(d1)))){
          d1 <- d1[, -grep("^Mix", names(d1))]
        }
        else{
          showNotification("The column 'Mix' was not found in your data !", type = "warning", duration = 5)
        }
      }

      if(input$rem_empty){
        if(length(grep("^Empty", names(d1)))){
          d1 <- d1[, -grep("^Empty", names(d1))]
        }
        else{
          showNotification("The column 'Empty' was not found in your data !", type = "warning", duration = 5)
        }
      }

      if(input$clean_data){
        d1 <- ms_clean(d1, nread = as.numeric(input$n_chan))
        showNotification("Cleaning done !", type = "message", duration = 3)
      }
    }
    d1
  })
  #check if a file is upload
  output$cetsa_cleanup <- reactive({
    return(!is.null(cetsa_data_clean()))
  })
  outputOptions(output, "cetsa_cleanup", suspendWhenHidden = FALSE)

  observeEvent(input$see2_cetsa,{
    if(!is.null(cetsa_data_clean())){
      showModal(tags$div(id="modal2", modalDialog(
        DT::renderDataTable({DT::datatable(cetsa_data_clean(),
                                           caption = htmltools::tags$caption(
                                             style = 'caption-side: top; text-align: left;',
                                             htmltools::strong("Base data")
                                           ),
                                           rownames = FALSE,
                                           options = list(lengthMenu = c(10,20,30), pageLength = 10,
                                                          scrollX = TRUE))
          }),
        footer = NULL,
        easyClose = TRUE
      )))
    }
    else{
      showNotification("The data are currently NULL, try to refresh.", type = "error")
    }
  })

  observe({
    if(input$step_cetsa > 2){
      updateCheckboxInput(session, "got_ISO_cetsa", value = TRUE)
    }
    if(input$step_cetsa > 3){
      updateCheckboxInput(session, "got_rearr_cetsa", value = TRUE)
    }
    if(input$step_cetsa > 4){
      updateCheckboxInput(session, "got_norm_cetsa", value = TRUE)
    }
    if(input$step_cetsa > 5){
      updateCheckboxInput(session, "got_diff_cetsa", value = TRUE)
    }
    if(input$step_cetsa < 2){
      updateCheckboxInput(session, "got_ISO_cetsa", value = FALSE)
    }
    if(input$step_cetsa < 3){
      updateCheckboxInput(session, "got_rearr_cetsa", value = FALSE)
    }
    if(input$step_cetsa < 4){
      updateCheckboxInput(session, "got_norm_cetsa", value = FALSE)
    }
    if(input$step_cetsa < 5){
      updateCheckboxInput(session, "got_diff_cetsa", value = FALSE)
    }
  })

  cetsa_isoform <- reactiveValues(
    x = NULL,
    y = NULL,
    norm = NULL,
    dif = NULL,
    conso = NULL,
    rearr = NULL
  )
  observeEvent(input$ISO, {
    showNotification("Start resolving isoform, this may take a while. Please wait a few minutes",
                     type = "message", duration = 5)
    x <- ms_isoform_resolve(cetsa_data_clean())

    cetsa_isoform$x <- x
  })
  ISOresdata_cetsa <- reactive({
    File <- input$ISOresfile_cetsa
    if (is.null(File) | !input$got_ISO_cetsa)
      return(NULL)

    ms_fileread(File$datapath)
  })
  observe({
    if(input$got_ISO_cetsa){
      cetsa_isoform$x <- ISOresdata_cetsa()
    }
  })
  #check if a file is upload
  output$cetsa_isoup <- reactive({
    return(!is.null(cetsa_isoform$x))
  })
  outputOptions(output, "cetsa_isoup", suspendWhenHidden = FALSE)

  observeEvent(input$see3_cetsa,{
    if(!is.null(cetsa_isoform$x)){
      showModal(tags$div(id="modal3", modalDialog(
        DT::renderDataTable({DT::datatable(cetsa_isoform$x,
                                           caption = htmltools::tags$caption(
                                             style = 'caption-side: top; text-align: left;',
                                             htmltools::strong("Base data")
                                           ),
                                           rownames = FALSE,
                                           options = list(lengthMenu = c(10,20,30), pageLength = 10,
                                                          scrollX = TRUE))
          }),
        footer = NULL,
        easyClose = TRUE
      )))
    }
    else{
      showNotification("The data are currently NULL, try to refresh.", type = "error")
    }
  })

  observeEvent(input$ISO2, {
    d2 <- cetsa_isoform$x
    if(!is.null(d2)){
      if(input$iso_conso){
        showNotification("Start consolidating isoform, this may take a while. Please wait a few minutes",
                         type = "message", duration = 5)
        File <- input$tab_conso
        if (is.null(File))
          return(NULL)

        d2 <- ms_isoform_consolidate(d2,
                                     nread = input$n_chan2,
                                     matchtable = File$datapath)
        cetsa_isoform$conso <- d2
        showNotification("Consolidation succeed !", type = "message", duration = 5)
      }

      if(input$iso_rearr){
        showNotification("Start rearranging data", type = "message", duration = 3)
        withCallingHandlers({
          shinyjs::html("diag_rearrange", "")
          d2 <- imprints_rearrange(d2, nread = input$n_chan3,
                                   repthreshold = input$rep_thr,
                                   averagecount = input$avgcount_abd,
                                   countthreshold = input$count_thr,
                                   withabdreading = input$wit_37)
          message("Done rearranging !")
        },
        message = function(m) {
          shinyjs::html(id = "diag_rearrange", html = paste(m$message, "<br>", sep = ""), add = TRUE)
        }
        )
        cetsa_isoform$rearr <- d2
        showNotification("Rearrangement succeed !", type = "message", duration = 5)
      }

      cetsa_isoform$y <- d2
    }
    else{
      showNotification("Don't forget to import a file or start the analysis", type = "error")
    }
  })
  rearrdata_cetsa <- reactive({
    File <- input$rearrfile_cetsa
    if (is.null(File) | !input$got_rearr_cetsa)
      return(NULL)

    ms_fileread(File$datapath)
  })
  observe({
    if(input$got_rearr_cetsa){
      cetsa_isoform$y <- rearrdata_cetsa()
    }
  })

  observeEvent(input$see4_cetsa,{
    if(!is.null(cetsa_isoform$conso)){
      showModal(tags$div(id="modal4", modalDialog(
        DT::renderDataTable({DT::datatable(cetsa_isoform$conso,
                                           caption = htmltools::tags$caption(
                                             style = 'caption-side: top; text-align: left;',
                                             htmltools::strong("Base data")
                                           ),
                                           rownames = FALSE,
                                           options = list(lengthMenu = c(10,20,30), pageLength = 10,
                                                          scrollX = TRUE))
          }),
        footer = NULL,
        easyClose = TRUE
      )))
    }
    else{
      showNotification("The data are currently NULL, try to refresh.", type = "error")
    }
  })

  observeEvent(input$see5_cetsa,{
    if(!is.null(cetsa_isoform$rearr)){
      showModal(tags$div(id="modal5", modalDialog(
        DT::renderDataTable({DT::datatable(cetsa_isoform$rearr,
                                           caption = htmltools::tags$caption(
                                             style = 'caption-side: top; text-align: left;',
                                             htmltools::strong("Base data")
                                           ),
                                           rownames = FALSE,
                                           options = list(lengthMenu = c(10,20,30), pageLength = 10,
                                                          scrollX = TRUE))
          }),
        footer = NULL,
        easyClose = TRUE
      )))
    }
    else{
      showNotification("The data are currently NULL, try to refresh.", type = "error")
    }
  })
  observeEvent(input$NORM, {
    if(is.null(cetsa_isoform$y)){
      d <- cetsa_isoform$x
    }
    else{
      d <- cetsa_isoform$y
    }

    if(!is.null(d)){
      showNotification("Start normalizing, this may take a while. Please wait a few minutes",
                       type = "message", duration = 5)
      d <- imprints_normalization(d)
      cetsa_isoform$norm <- d
      showNotification("Normalization succeed !", type = "message", duration = 5)
    }
    else{
      showNotification("Don't forget to import a file or start the analysis", type = "error")
    }

  })
  normdata_cetsa <- reactive({
    File <- input$normfile_cetsa
    if (is.null(File) | !input$got_norm_cetsa)
      return(NULL)

    ms_fileread(File$datapath)
  })
  observe({
    if(input$got_norm_cetsa){
      cetsa_isoform$norm <- normdata_cetsa()
    }
  })
  #check if a file is upload
  output$cetsa_normup <- reactive({
    return(!is.null(cetsa_isoform$norm))
  })
  outputOptions(output, "cetsa_normup", suspendWhenHidden = FALSE)

  observeEvent(input$see6_cetsa,{
    if(!is.null(cetsa_isoform$norm)){
      showModal(tags$div(id="modal6", modalDialog(
        DT::renderDataTable({DT::datatable(cetsa_isoform$norm,
                                           caption = htmltools::tags$caption(
                                             style = 'caption-side: top; text-align: left;',
                                             htmltools::strong("Base data")
                                           ),
                                           rownames = FALSE,
                                           options = list(lengthMenu = c(10,20,30), pageLength = 10,
                                                          scrollX = TRUE))
          }),
        footer = NULL,
        easyClose = TRUE
      )))
    }
    else{
      showNotification("The data are currently NULL, try to refresh.", type = "error")
    }
  })

  observe({
    if(!is.null(cetsa_isoform$norm)){
      tr_level <- get_treat_level(cetsa_isoform$norm)
      updateSelectInput(session, "ctrl_name2", choices = tr_level, selected = tr_level[1])
    }
  })
  observeEvent(input$CAL_DIF, {
    if(!is.null(cetsa_isoform$norm)){
      showNotification("Start fold-change calculation, this may take a while. Please wait a few minutes",
                       type = "message", duration = 5)
      tr_level <- get_treat_level(cetsa_isoform$norm)
      tr_level <- tr_level[-which(tr_level == input$ctrl_name2)]
      tr_level <- c(input$ctrl_name2, tr_level)
      d <- imprints_caldiff_f(cetsa_isoform$norm,
                           reftreatment = tr_level,
                           withinrep = input$wit_rep
                         )

      cetsa_isoform$dif <- d
      message("Done to calculate the pair-wise (per replicate and temperature)
               protein abundance differences")
      showNotification("Difference calculation succeed !", type = "message", duration = 5)

    }
    else{
      showNotification("Don't forget to import a file or start the analysis", type = "error")
    }
  })
  diffdata_cetsa <- reactive({
    File <- input$difffile_cetsa
    if (is.null(File) | !input$got_diff_cetsa)
      return(NULL)

    read.delim(File$datapath, check.names = FALSE)
  })
  observe({
    if(input$got_diff_cetsa){
      cetsa_isoform$dif <- diffdata_cetsa()
    }
  })
  #check if a file is upload
  output$cetsa_difup <- reactive({
    return(!is.null(cetsa_isoform$dif))
  })
  outputOptions(output, "cetsa_difup", suspendWhenHidden = FALSE)

  observeEvent(input$see7_cetsa,{
    if(!is.null(cetsa_isoform$dif)){
      showModal(tags$div(id="modal7",modalDialog(
        DT::renderDataTable({DT::datatable(cetsa_isoform$dif,
                                           caption = htmltools::tags$caption(
                                             style = 'caption-side: top; text-align: left;',
                                             htmltools::strong("Base data")
                                           ),
                                           rownames = FALSE,
                                           options = list(lengthMenu = c(10,20,30), pageLength = 10,
                                                          scrollX = TRUE))
          }),
        footer = NULL,
        easyClose = TRUE
      )))
    }
    else{
      showNotification("The data are currently NULL, try to refresh.", type = "error")
    }
  })


  # hitlist doc
  observeEvent(input$see8_cetsa,{
    showModal(tags$div(id="modal8_FC",modalDialog(
        shiny::HTML("<h1>Hitlist generation: using fold-change cutoff</h1><br>
                    <br>This method is simple and only needs the caldiff output, i.e. the fold-changes data.
                    <br>It defines a protein as a hit if in at least one of the temperatures, a fold change
                    passes some criterias. In the shiny app, you can choose two different cutoffs: a cutoff
                    on the mean and the acceptable boundedness.
                    <br> By default, the values are 0.25 and 4. It means that if a protein has one mean fold change
                    superior than 0.25 (in absolute value) and if this absolute value is superior than 4*SEM
                    (its Standard Error of the Mean), it will be considered as a hits.
                    <br><br>For the categorization, if a hits is an ND, it means that its 37°C value is not well
                    measured (i.e. has a missing value or has an SEM > 0.15). Indeed the categorization totally depends
                    on the behaviour of the protein at 37°C.
                    <br><br><h6><em>This function hasn't been written by me</em></h6>"),
        footer = NULL,
        easyClose = TRUE
      )))
  })
  observeEvent(input$see9_cetsa,{
    showModal(tags$div(id="modal9_IS",modalDialog(
        shiny::HTML("<h1>Hitlist generation: Intercept Score</h1><br>
                    <br>This method compute a score for each protein and returns a volcano plot.
                    Since it will compute p-value, it needs the post-normalization file. It will
                    also compute the fold-changes but if you already have the caldiff output, you
                    will gain some time.
                    <br>For each protein, we want to extract a p-value and a score, the Intercept Score (IS).
                    <br>Let's consider our data with <var>n</var> proteins, <var>T</var> temperatures,
                    1 control and <var>C</var> treatments, <var>B</var> bio-replicates where <var>B</var> >= 3.
                    <br><br><u>1. The p-value</u><br><br>
                    For each treatment <var>c</var> = 1, ..., <var>C</var>, we compute a moderated t-test for each temperature
                    <var>t</var> = 1, ..., <var>T</var>. We now have 𝑇 p-values for each protein <var>p</var> = 1, ..., <var>n</var>.
                    We then only keep the two p-values from the two biggest mean fold changes (in absolute value) from
                    each protein <var>p</var>. We now compute a Fisher’s t-test on these two p-values which return
                    one final p-value. (Figure 1. A.)
                    <br>Finally, we have <var>n</var> p-values for each treatment <var>c</var> = 1, ..., <var>C</var>.
                    <br><br><u>2. IS</u><br><br>
                    For each protein <var>p</var> = 1, ..., <var>n</var> and for each treatment <var>c</var> = 1, ..., <var>C</var>,
                    we extract the <var>T</var> mean fold changes. We take the absolute value from these and then order them in
                    the decreasing order. We now fit a weighted linear regression to these data with
                    weights five times higher for the two biggest fold changes, i.e. the two first data
                    points. Which means more importance is given to the two biggest fold changes but
                    we still take into account the whole profile. (Figure 1. A.)
                    <br>From this regression we extract the intercept with the y-axis. This intercept will
                    always be positive and to keep track if the protein is destabilized or stabilized, we
                    multiply this value by the sign of the mean of all the fold changes (either -1 or 1).
                    <br>Finally, we apply a z-score normalization on all IS for each treatment <var>c</var> = 1, ..., <var>C</var>.
                    <br><br>In the end, we plot the -log10(p-value) vs IS which gives us a volcano plot. (Figure 1. B.)
                    <br><br><img src='IS_figure1.png' alt='IS figure', width='1180' height='660'>
                    <br><br><u>Set the cutoffs</u><br><br>
                    We have two cutoffs to set in order to choose what are the significant hits: an IS cutoff and
                    a p-value cutoff.
                    <br>For the p-value, we apply the Benjamini and Hochberg correction on the p-values with a
                    chosen FDR (1% by default) which gives us the p-value cutoff to set for the corresponding FDR.
                    <br>To set the IS cutoff, we apply the default method used in a classical volcano plot. We
                    compute the median of the p-values <var>p<sub>m</sub></var> and then we compute the median of
                    the IS, <var>IS<sub>m</sub></var> with a corresponding p-value <var>p</var> < <var>p<sub>m</sub></var>.
                    Finally, to the default value chosen (typically 2 but here 1.5 is the default value), we add and remove
                    <var>IS<sub>m</sub></var>; so we have the cutoffs − 1. 5 − <var>IS<sub>m</sub></var> and
                    1. 5 + <var>IS<sub>m</sub></var>.
                    <br>In the end, a curve is computed from these cutoffs with by default a curvature of 0.5 which
                    gives us the final ‘cutoff curve’. Every protein that is above this curve in the volcano plot
                    is then considered as a significant hit. (Figure 1. B.)"),
        footer = NULL,
        easyClose = TRUE
      )))
  })
  # hitlist calculation
  hit_pr <- reactiveValues(
    hitlist = NULL,
    ND = NULL,
    NC = NULL,
    CN = NULL,
    CC = NULL
  )
  observeEvent(input$str_calchitlist, {
    if(!is.null(cetsa_isoform$dif)){
      if(input$hitmethod_cetsa == "ImpS"){
        showNotification("Calculation started, this may take a while. Please wait a few minutes !",
                         type = "message", duration = 5)

        Dif <- cetsa_isoform$dif
        if(any(stringr::str_detect(colnames(Dif), "^X\\d{2}C_"))){
          bad_col <- stringr::str_which(colnames(Dif), "^X\\d{2}C_")
          colnames(Dif)[bad_col] <- stringr::str_remove_all(colnames(Dif)[bad_col], "^X")
        }
        if(length(grep("^36C_", colnames(Dif)))){   # remove QP if in data
          Dif <- Dif[,-grep("^36C_", colnames(Dif))]
        }
        withCallingHandlers({
          shinyjs::html("diag_ImpS", "")
          h <- imprints_score(Dif, format = input$formatImpS_cetsa,
                              fdrthreshold = input$FDRImpS_cetsa,
                              cvcutoffthreshold = input$cvcutoff_cetsa)
          message("Done !")
        },
        message = function(m) {
          shinyjs::html(id = "diag_ImpS", html = paste(m$message, "<br>", sep = ""), add = FALSE)

        }
        )

        hit_pr$hitlist <- h[which(h$category != "NN"),]
        hit_pr$ND <- h[which(h$category == "ND"),]
        hit_pr$NC <- h[grep("NC", h$category),]
        hit_pr$CN <- h[grep("CN", h$category),]
        hit_pr$CC <- h[grep("CC", h$category),]
      }
      else if(input$hitmethod_cetsa == "IS"){
        if(!is.null(cetsa_isoform$norm)){
          showNotification("Calculation started, this may take a while. Please wait a few minutes !",
                           type = "message", duration = 5)

          Dif <- cetsa_isoform$dif
          if(any(stringr::str_detect(colnames(Dif), "^X\\d{2}C_"))){
            bad_col <- stringr::str_which(colnames(Dif), "^X\\d{2}C_")
            colnames(Dif)[bad_col] <- stringr::str_remove_all(colnames(Dif)[bad_col], "^X")
          }
          ctrl <- Dif[,stringr::str_which(colnames(Dif), "^\\d{1,}")]
          idx_ctrl <- which(apply(ctrl, 1, function(x) all(!is.na(x))))[1]
          ctrl <- ctrl[idx_ctrl,]
          ctrl <- ctrl %>% tidyr::gather("key", "value") %>%
            tidyr::separate(key, into = c("t", "b", "cond"), sep = "_") %>%
            dplyr::group_by(cond) %>%
            dplyr::reframe(ctrl = all(value == 0))
          ctrl <- ctrl$cond[ctrl$ctrl]

          withCallingHandlers({
            shinyjs::html("diag_IS", "")
            h <- imprints_IS(cetsa_isoform$norm, Dif, ctrl = ctrl,
                             valid_val = input$validval_cetsa,
                             IS_cutoff = input$IScut_cetsa,
                             FDR = input$FDR_cetsa)
          },
          message = function(m) {
            shinyjs::html(id = "diag_IS", html = paste(m$message, "<br>", sep = ""), add = FALSE)

          }
          )

          hit_pr$hitlist <- h
          hit_pr$ND <- h %>% filter(category == "ND")
          hit_pr$NC <- h %>% filter(category == "NC")
          hit_pr$CN <- h %>% filter(category == "CN")
          hit_pr$CC <- h %>% filter(category == "CC")
        }
        else{
          showNotification("You also need the post normalization data for this method.
                           Import a file or start the analysis.", type = "error", duration = 8)
        }
      }
      else if(input$hitmethod_cetsa == "FC"){
        showNotification("Calculation started, this may take a while. Please wait a few minutes !",
                         type = "message", duration = 5)

        h <- hitlist(cetsa_isoform$dif, meancutoff = input$meancut_cetsa, boundedness = input$bound_cetsa,
                     use_prompt = FALSE, exported = input$save_hit)

        hit_pr$hitlist <- h$hitlist
        hit_pr$ND <- h$ND
        hit_pr$NC <- h$NC
        hit_pr$CN <- h$CN
        hit_pr$CC <- h$CC
      }
    }
    else{
      showNotification("Don't forget to import a file or start the analysis", type = "error")
    }
  })

  output$hit_out <- DT::renderDataTable({
    DT::datatable(hit_pr[[input$HIT]],
                  caption = htmltools::tags$caption(
                    style = 'caption-side: top; text-align: left;',
                    htmltools::strong(input$HIT)
                  ),
                  rownames = FALSE,
                  options = list(lengthMenu = c(10,20,30), pageLength = 10,
                                 scrollX = TRUE))
  })



  ### Data base ###
  drug_data_sh <- reactiveValues()
  if(exists("drug_data2")){
    drug_data_sh$y <- drug_data2
  }
  else{
    if(file.exists("drug_data")){
      drug_data_sh$y <- load_data()
    }
    else{
      drug_data_sh$y <- drug_data
    }
  }

  output$davai_daba_ui <- renderUI({
    selectInput("davai_daba", "Choose a dataset to remove", choices = names(drug_data_sh$y$data),
                selected = "elutriation")
  })
  output$davai2_daba_ui <- renderUI({
    selectInput("davai2_daba", "Choose a dataset in which you want to rename the treatments", choices = names(drug_data_sh$y$data),
                selected = "elutriation")
  })
  output$drug2_ui <- renderUI({
    selectInput("drug2", "Choose a dataset", choices = names(drug_data_sh$y$data), multiple = TRUE, selected = "elutriation")
  })

  observeEvent(input$up_daba,{
    showNotification("Start loading the data, this may take a while", type = "message")
    a <- load_data()
    if(!is.null(a)){
      drug_data_sh$y <- a
      drug_data2 <<- drug_data_sh$y

      updateSelectInput(session, "davai_daba", choices = names(a$data), selected = names(a$data)[1])
      updateSelectInput(session, "drug2", choices = names(a$data), selected = names(a$data)[1])

      showNotification("Data loaded !", type = "message")
    }
    else{
      showNotification("No folder named 'drug_data' has been created; there is nothing to update.", type = "error")
    }
  })

  output$info_daba <- renderText({
    dn <- length(drug_data_sh$y$data)
    d <- names(drug_data_sh$y$data)

    if(dn > 1){
      d <- paste(c(paste(d[1:(dn-1)], collapse = ", "), d[dn]), collapse = " and ")
    }

    HTML(paste("<p><h4>Your database contains at the moment", dn, "datasets :", d, ".</h4></p>"))
  })


  DIF_daba <- reactive({
    File <- input$caldif_daba
    if (is.null(File))
      return(NULL)

    ms_fileread(File$datapath)
  })
  #check if a file is upload
  output$DIFdaba_fileup <- reactive({
    return(!is.null(DIF_daba()))
  })
  outputOptions(output, "DIFdaba_fileup", suspendWhenHidden = FALSE)

  AVE_daba <- reactive({
    if(!input$gave_daba){
      File <- input$AVE_dabafile
      if (is.null(File))
        return(NULL)

      ms_fileread(File$datapath)
    }
    else{
      1  #simplify treatment is.null
    }
  })
  #check if a file is upload
  output$AVEdaba_fileup <- reactive({
    return(!is.null(AVE_daba()))
  })
  outputOptions(output, "AVEdaba_fileup", suspendWhenHidden = FALSE)

  NN_daba <- reactiveValues(
    x = NULL
  )
  HIT_daba <- reactive({
    File <- input$hitsum_daba
    if (is.null(File))
      return(NULL)

    dat <- import(File$datapath, header = TRUE)
    nv_nam <- str_subset(names(dat), "^V\\d{1}$")
    if(length(nv_nam)){
      dat <- dat[, !(names(dat) %in% nv_nam)]
    }
    if(!("treatment" %in% colnames(dat))){ # means that the analysis tab was imported
      dat <- dat[,stringr::str_which(colnames(dat), "^id$|^Fisher_|^IS_|^GlobalScore_|^category_")]
      dat <- dat %>% tidyr::gather("key", "value", -id) %>%
        tidyr::separate(key, into = c("key", "treatment"), sep = "_") %>%
        tidyr::spread(key, value)

      nn <- dat %>% dplyr::filter(category == "NN")
      dif <- DIF_daba()[,1:2]
      nn <- dplyr::left_join(nn, dif, by = "id")
      nn <- nn[,c("id", "description", "treatment", "category", "Fisher", "IS", "GlobalScore")]
      NN_daba$x <- nn

      dat <- dat %>% dplyr::filter(category != "NN")
    }
    else{
      dif <- DIF_daba()[,1:2]
      dat <- dat[,c("id", "treatment", "category")]
      nn <- lapply(unique(dat$treatment), function(x){
        x <- dat %>% dplyr::filter(treatment == x) %>%
          dplyr::right_join(dif, by = "id") %>%
          dplyr::filter(is.na(category)) %>%
          dplyr::mutate(category = "NN",
                              treatment = x);
        x
      })
      nn <- as.data.frame(Reduce(rbind, nn))
      nn <- nn[,c("id", "description", "treatment", "category")]
      NN_daba$x <- nn
    }
    dat
  })
  #check if a file is upload
  output$HITdaba_fileup <- reactive({
    return(!is.null(HIT_daba()))
  })
  outputOptions(output, "HITdaba_fileup", suspendWhenHidden = FALSE)


  observeEvent(input$add_daba, {
    if(input$name_daba %in% names(drug_data_sh$y$data)){
      showNotification("This name is already taken ! Please, choose anoter one.", type = "error")
    }
    else{
      ave_data <- NULL
      if(!input$gave_daba & !is.null(AVE_daba())){
        ave_data <- AVE_daba()
      }
      else{
        showNotification("Getting average dataset, this may take a while.", type = "message")
        ave_data <- imprints_average(DIF_daba(), savefile = TRUE)
        showNotification("Average calculation succeed !", type = "message")
      }
      showNotification("Start saving dataset, this may take a while.", type = "message")
      save_data(drug_data_sh$y, new_add = list("data" = DIF_daba(),
                                              "data_ave" = ave_data,
                                              "hitlist" = HIT_daba(),
                                              "NN" = NN_daba$x,
                                              "treat_level" = get_treat_level(DIF_daba())),
               input$name_daba)


      showNotification("New dataset added !", type = "message")
    }
  })

  output$condfrom_daba <- renderUI({
    cd_info <- NULL
    df <- NULL
    if(!is.null(input$davai2_daba)){
      df <- drug_data_sh$y$data[[input$davai2_daba]]
    }

    if(!is.null(df) & length(df)){
      cd <- get_treat_level(df)
      cd_1 <- cd[-length(cd)]
      cd_e <- cd[length(cd)]

      cd_info <- paste(paste(cd_1, collapse = ", "), cd_e, sep = " and ")
    }
    HTML(paste("<p>Your current condition names for the drug", input$davai2_daba, "are :", paste0("<b>", cd_info, "</b>"), "</p>"))
  })

  observeEvent(input$changename_daba, {
    showNotification("Checking new names", type = "message", duration = 2)
    nm <- input$condnew_daba
    nm <- str_split(nm, ",")[[1]]
    nm <- str_remove_all(nm, " ")
    if(sum(str_detect(nm, "_|/")) > 0){
      showNotification("The character '_' and '/' are not alllowed. Please, verify your new names.", type = "error")
    }
    else{
      if(!is.null(input$davai2_daba)){
        cd <- get_treat_level(drug_data_sh$y$data[[input$davai2_daba]])
      }

      for(i in 1:length(cd)){
        if(str_length(nm[i]) == 0){
          nm[i] <- cd[i]
        }
      }
      change <- cd[!(nm %in% cd)]
      if(length(change) & !is.null(change)){
        new <- nm[!(nm %in% cd)]
        showNotification(paste("You decided to change :", paste(change, collapse = ", "),
                               "In :", paste(new, collapse = ", ")), type = "message")

        df <- drug_data_sh$y$data[[input$davai2_daba]]
        df_ave <- drug_data_sh$y$data_ave[[input$davai2_daba]]
        dh <- drug_data_sh$y$hitlist[[input$davai2_daba]]
        dnn <- drug_data_sh$y$NN[[input$davai2_daba]]

        n_df <- names(df)[str_detect(names(df), paste(paste0("_", change, "$"), collapse = "|"))]
        n_df_ave <- names(df_ave)[str_detect(names(df_ave), paste(paste0("_", change, "$"), collapse = "|"))]
        n_dh <- dh$treatment
        n_dnn <- dnn$treatment

        for(i in 1:length(change)){
          n_df <- str_replace_all(n_df, paste0("_", change[i], "$"), paste0("_", new[i]))
          n_df_ave <- str_replace_all(n_df_ave, paste0("_", change[i], "$"), paste0("_", new[i]))
          n_dh <- str_replace_all(n_dh, paste0("^", change[i], "$"), new[i])
          n_dnn <- str_replace_all(n_dnn, paste0("^", change[i], "$"), new[i])
        }

        names(df)[str_detect(names(df), paste(paste0("_", change, "$"), collapse = "|"))] <- n_df
        names(df_ave)[str_detect(names(df_ave), paste(paste0("_", change, "$"), collapse = "|"))] <- n_df_ave
        dh$treatment <- n_dh
        dnn$treatment <- n_dnn
        dt <- get_treat_level(df)

        showNotification("Start saving changes, this may take a while.", type = "message")
        save_data(drug_data_sh$y, new_add = list("data" = df,
                                                "data_ave" = df_ave,
                                                "hitlist" = dh,
                                                "NN" = dnn,
                                                "treat_level" = dt),
                 input$davai2_daba)


        showNotification("Names changed !", type = "message")

      }
      else{
        showNotification("You didn't make any changement !", type = "error")
      }
    }
  })

  observeEvent(input$rem_daba, {
    showModal(
      modalDialog(
        title="Are you sure you want to remove this dataset ?",
        footer = tagList(actionButton("confirmRem", "Remove"),
                         modalButton("Cancel")
        )
      )
    )
  })
  observeEvent(input$confirmRem, {
    showNotification("Start removing dataset", type = "message")
    rem_data(drug_data_sh$y, input$davai_daba)

    removeModal()

    showNotification("Dataset removed !", type = "message")
  })


  ### BAR PLOT TAB ###

  barplot_data <- reactive({
    File <- input$data_barplot
    if (is.null(File))
      return(NULL)

    ms_fileread(File$datapath)
  })
  #check if a file is upload
  output$barplot_dataup <- reactive({
    return(!is.null(barplot_data()))
  })
  outputOptions(output, "barplot_dataup", suspendWhenHidden = FALSE)

  barhit_data <- reactive({
    File <- input$data_hitlist
    if (is.null(File))
      return(NULL)

    d <- import(File$datapath, header = TRUE)
    if(!all(c("id", "treatment", "category") %in% colnames(d))){
      missing_columns <- c("id", "treatment", "category")
      missing_columns <- missing_columns[!(c("id", "treatment", "category") %in% colnames(d))]
      missing_columns <- paste(missing_columns, collapse = ", ")
      verb <- ifelse(length(missing_columns) > 1, "are", "is")
      showNotification(paste(missing_columns, verb, "not in your summary file. Please check your columns names !"),
                       type = "error", duration = 8)
      d <- NULL
    }
    d
  })
  #check if a file is upload
  output$barhit_dataup <- reactive({
    return(!is.null(barhit_data()))
  })
  outputOptions(output, "barhit_dataup", suspendWhenHidden = FALSE)

  hit_bar <- reactiveValues(
    summa = NULL,
    NN = NULL
  )

  observeEvent(input$str_calchit, {
    showNotification("Calculation started, this may take a while. Please wait a few minutes !",
                     type = "message", duration = 5)

    h <- hitlist(barplot_data(), meancutoff = input$meancut_bar, boundedness = input$bound_bar,
                 use_prompt = FALSE, exported = input$save_hit_bar)
    h_s <- rbind(h$CC, h$CN, h$NC, h$ND)

    hit_bar$summa <- h_s %>% group_by(id,treatment,category) %>%  summarize()
    hit_bar$NN <- h$NN
  })


  Sel_cond_fhit_SUMMA <- reactiveValues(
    choice = NULL,
    hit = NULL
  )
  Sel_cond_fhit <- reactive({
    HIT <- NULL

    if(input$hit){
      if(input$drug == "base" & length(input$drug2) >= 1){
        HIT <- do.call(rbind, lapply(drug_data_sh$y$hitlist[input$drug2],
                                     function(x) x[,c("id", "treatment", "category")])
                       )
      }
      else if(input$drug == "dat"){
        if(is.null(hit_bar$summa)){
          HIT <- barhit_data()
        }
        else{
          HIT <- hit_bar$summa
        }
      }
      c_idx <- str_which(colnames(HIT), "treatment")
      if(length(c_idx)){
        HIT_summup <- list()
        for(i in unique(HIT[, c_idx])){
          HIT_summup[[i]] <- (HIT %>% dplyr::filter(treatment == i))$id
        }
        HIT_summup <- com_protein_loop(HIT_summup)

        for (i in names(HIT_summup)){
          HIT[which(!is.na(match(HIT$id, HIT_summup[[i]]))), c_idx] <- i
        }
        HIT <- unique(HIT)
      }
    }

    HIT
  })

  observe({
    if(!is.null(Sel_cond_fhit())){
      c_idx <- str_which(colnames(Sel_cond_fhit()), "treatment")
      Sel_cond_fhit_SUMMA$hit <- Sel_cond_fhit()
      if(length(c_idx)){
        Sel_cond_fhit_SUMMA$choice <- unique(Sel_cond_fhit()[,c_idx])
      }
      updateSelectInput(session, "cond_fhit", choices = Sel_cond_fhit_SUMMA$choice, selected = Sel_cond_fhit_SUMMA$choice[1])
    }
  })

  sel_prot <- reactive({
    pr <- NULL
    if(input$drug == "base"){
      if(input$protlist_bar){
        File <- input$prlist_file_bar
        if (is.null(File)){
          pr <- NULL
        }
        else{
          pr <- unique(read.delim(File$datapath, header = FALSE)[[1]])

          prcheck <- ""
          if(length(input$drug2) == 1){
            prcheck <- drug_data_sh$y$data[[input$drug2]][,c("id", "description")]
          }
          else if(length(input$drug2) > 1){
            prcheck <- plyr::join_all(drug_data_sh$y$data[input$drug2], by = c("id", "description"), type = "full")[,c("id", "description")]
          }

          a <- pr[!(pr %in% prcheck$id)]
          if(length(a)){
            pr <- pr[(pr %in% prcheck$id)]
            showNotification(paste(paste(a, collapse = ", "), "wasn't in the data and had to be removed."),
                             type = "error")
          }

          pr <- data.frame(id = pr)
          pr <- unique(pr)
          pr <- dplyr::left_join(pr, prcheck,
                                 by = "id")
          pr$description <- unname(sapply(pr$description, IMPRINTS.CETSA.app:::getGeneName))
          pr <- pr %>%
            group_by(description) %>%
            mutate(id = paste0(id, ":", description))
          pr <- pr$id
          pr <- unique(pr)
        }
      }
      else{
        if(length(input$drug2) == 1 & all(input$drug2 %in% names(drug_data_sh$y$data))){
          df <- drug_data_sh$y$data[[input$drug2]][,c("id", "description")]
          if(input$hit & !is.null(input$cond_fhit)){
            pr <- Sel_cond_fhit_SUMMA$hit
            pr <- pr %>% dplyr::filter(!is.na(match(treatment, c(input$cond_fhit))))
            pr <- pr[,"id", drop = FALSE]
            pr <- unique(pr)
            pr <- dplyr::left_join(pr, df,
                                   by = "id")
            pr$description <- unname(sapply(pr$description, IMPRINTS.CETSA.app:::getGeneName))
            pr <- pr %>%
              group_by(description) %>%
              mutate(id = paste0(id, ":", description))
            pr <- pr$id
            pr <- unique(pr)
          }
          else{
            pr <- df
            pr$description <- unname(sapply(pr$description, IMPRINTS.CETSA.app:::getGeneName))
            pr <- pr %>%
              group_by(description) %>%
              mutate(id = paste0(id, ":", description))
            pr <- pr$id
            pr <- unique(pr)
          }
        }
        else if(length(input$drug2) > 1 & all(input$drug2 %in% names(drug_data_sh$y$data))){
          df <- plyr::join_all(drug_data_sh$y$data[input$drug2], by = c("id", "description"), type = "full")
          if(input$hit & !is.null(input$cond_fhit)){
            pr <- Sel_cond_fhit_SUMMA$hit
            pr <- pr %>% dplyr::filter(!is.na(match(treatment, c(input$cond_fhit))))
            pr <- pr[,"id", drop = FALSE]
            pr <- unique(pr)
            pr <- dplyr::left_join(pr, df,
                                   by = "id")
            pr$description <- unname(sapply(pr$description, IMPRINTS.CETSA.app:::getGeneName))
            pr <- pr %>%
              group_by(description) %>%
              mutate(id = paste0(id, ":", description))
            pr <- pr$id
            pr <- unique(pr)
          }
          else{
            pr <- df
            pr$description <- unname(sapply(pr$description, IMPRINTS.CETSA.app:::getGeneName))
            pr <- pr %>%
              group_by(description) %>%
              mutate(id = paste0(id, ":", description))
            pr <- pr$id
            pr <- unique(pr)
          }
        }
      }
    }

    else if(input$drug == "dat"){
      if(input$protlist_bar){
        File <- input$prlist_file_bar
        if (is.null(File)){
          pr <- NULL
        }
        else{
          pr <- unique(read.delim(File$datapath, header = FALSE)[[1]])

          prcheck <- ""
          prcheck <- barplot_data()[,c("id", "description")]

          a <- pr[!(pr %in% prcheck$id)]
          if(length(a)){
            pr <- pr[(pr %in% prcheck$id)]
            showNotification(paste(paste(a, collapse = ", "), "wasn't in the data and had to be removed."),
                             type = "error")
          }

          pr <- data.frame(id = pr)
          pr <- unique(pr)
          pr <- dplyr::left_join(pr, prcheck,
                                 by = "id")
          pr$description <- unname(sapply(pr$description, IMPRINTS.CETSA.app:::getGeneName))
          pr <- pr %>%
            group_by(description) %>%
            mutate(id = paste0(id, ":", description))
          pr <- pr$id
          pr <- unique(pr)
        }
      }
      else{
        df <- barplot_data()[,c("id", "description")]
        if(!is.null(df)){
          if(input$hit  & !is.null(input$cond_fhit)){
            pr <- Sel_cond_fhit_SUMMA$hit
            pr <- pr %>% dplyr::filter(!is.na(match(treatment, c(input$cond_fhit))))
            pr <- pr[,"id", drop = FALSE]
            pr <- unique(pr)
            pr <- dplyr::left_join(pr, df,
                                   by = "id")
            pr$description <- unname(sapply(pr$description, IMPRINTS.CETSA.app:::getGeneName))
            pr <- pr %>%
              group_by(description) %>%
              mutate(id = paste0(id, ":", description))
            pr <- pr$id
            pr <- unique(pr)
          }
          else{
            pr <- df
            pr$description <- unname(sapply(pr$description, IMPRINTS.CETSA.app:::getGeneName))
            pr <- pr %>%
              group_by(description) %>%
              mutate(id = paste0(id, ":", description))
            pr <- pr$id
            pr <- unique(pr)
          }
        }
      }
    }
    pr
  })

  observe({
    updateSelectizeInput(session, "prot", choices = sel_prot(), server = TRUE)
  })

  Sel_cond <- reactive({
    if(input$ALL_prot){
      PROT <- sel_prot()
    }
    else{
      PROT <- input$prot
    }
    if(!is.null(PROT)){
      PROT <- unname(sapply(PROT, function(x) strsplit(x, ":")[[1]][1]))
    }

    if(input$drug == "base" & length(input$drug2) >= 1){
      HIT <- do.call(rbind, lapply(drug_data_sh$y$hitlist[input$drug2],
                                   function(x) x[,c("id", "treatment", "category")])
                     )
      NN <- do.call(rbind, lapply(drug_data_sh$y$NN[input$drug2],
                                  function(x) x[,c("id", "description", "treatment", "category")])
                    )
    }
    else if(input$drug == "dat"){
      if(is.null(hit_bar$summa)){

        HIT <- barhit_data()
        if(!is.null(HIT)){
          NN <- lapply(unique(HIT$treatment), function(z){
            z <- HIT %>% dplyr::filter(treatment == z) %>%
              dplyr::right_join(barplot_data()[,1:2], by = "id") %>%
              dplyr::filter(is.na(category)) %>%
              dplyr::mutate(category = "NN",
                            treatment = z);
            z
          })
          NN <- as.data.frame(Reduce(rbind, NN))
          NN <- NN[,c("id", "description", "treatment", "category")]
        }
      }
      else{
        HIT <- hit_bar$summa
        NN <- hit_bar$NN
      }
    }


    tr <- NULL
    if(input$cond_sel == "cat"){
      if(length(input$drug2) >= 1){
        trh <- HIT[which(!is.na(match(HIT$id, PROT))),c("treatment", "category")]
        tr <- NN[which(!is.na(match(NN$id, PROT))), c("treatment", "category")]
        tr <- tr[!duplicated(tr),]
        tr <- rbind(trh, tr)
      }
    }
    else if(input$cond_sel == "treat"){
      if(input$drug == "base" & length(input$drug2) >= 1){
        a <- join_drugdata(drug_data_sh$y$data[input$drug2], by = c("id", "description"))
        TRE <- get_treat_level(a)

        if(input$rem_con){
          tr <- TRE[-grep(input$con_name, TRE)]
        }
        else{
          tr <- TRE
        }
      }
      else if(input$drug == "dat"){
        if(input$rem_con){
          tr <- get_treat_level(barplot_data())[-grep(input$con_name, get_treat_level(barplot_data()))]
        }
        else{
          tr <- get_treat_level(barplot_data())
        }
      }
    }
    else{
      tr <- NULL
    }

    tr
  })

  observe({
    updateSelectInput(session, "cond", choices = Sel_cond())
  })

  observe({
    if(input$cond_sel == "cat"){
      updateSelectInput(session, "cond", choices = unique(Sel_cond()$category))
    }
    else{
      updateSelectInput(session, "cond", choices = Sel_cond())
    }
  })

  data <- reactive({
    if(input$ALL_prot){
      PROT <- sel_prot()
    }
    else{
      PROT <- input$prot
    }
    if(!is.null(PROT)){
      PROT <- unname(sapply(PROT, function(x) strsplit(x, ":")[[1]][1]))
    }

    if (input$drug == "base"){
      data <- join_drugdata(drug_data_sh$y$data[input$drug2], by = c("id", "description"))
      TREAT <- get_treat_level(data)
    }

    else if(input$drug == "dat"){
      data <- barplot_data()
      TREAT <- get_treat_level(barplot_data())
    }

    if(input$rem_con){
      data <- data[,-grep(input$con_name,  names(data))]
      TREAT <- get_treat_level(data)
    }

    data <- data[which(!is.na(match(data$id, PROT))),]


    if(input$cond_sel == "treat"){
      notsel_cond <- TREAT[!(TREAT %in% input$cond)]
      if(length(notsel_cond)){
        notsel_cond <- paste0("_", notsel_cond, "$")
        notsel_cond <- paste(notsel_cond, collapse = "|")
        data <- data[,-str_which(names(data), notsel_cond)]
      }

      id_sel <- str_which(names(data), paste(input$cond, collapse = "|"))
      w <- 1:ncol(data)
      w <- w[!(w %in% id_sel)]

      ord <- unlist(lapply(input$cond, function(x) str_which(names(data), paste0("_", x, "$"))))

      data <- data[,c(w,ord)]
    }
    else if(input$cond_sel == "cat"){
      sele_cond <- Sel_cond()$treatment[which(!is.na(match(Sel_cond()$category, input$cond)))]
      notsel_cond <- TREAT[!(TREAT %in% sele_cond)]
      notsel_cond <- paste(notsel_cond, collapse = "|")

      data <- data[,-str_which(names(data), notsel_cond)]
    }

    data
  })

  DAT_text <- reactive({
    DAT <- NULL
    if(input$drug == "base" & length(input$drug2) > 1){
      only_dat <- drug_data_sh$y$data[input$drug2]

      DAT <- join_drugdata(only_dat, by = c("id", "description"))
      DAT$drug <- rep(paste(input$drug2, collapse = " and "), nrow(DAT))
      DAT_id <- lapply(only_dat, function(x) x$id)

      com_pr <- com_protein_loop(DAT_id)

      for (i in names(com_pr)){
        DAT$drug[which(!is.na(match(DAT$id, com_pr[[i]])))] <- i
      }

      DAT <- DAT[,c("id", "drug")]
    }

  })
  tabident_bar <- reactiveValues(
    r = NULL
  )
  tabidentreac_bar <- reactive({
    if(input$ALL_prot){
      PROT <- sel_prot()
    }
    else{
      PROT <- input$prot
    }
    if(!is.null(PROT)){
      PROT <- unname(sapply(PROT, function(x) strsplit(x, ":")[[1]][1]))
    }
    DR <- NULL
    if(input$drug == "base" & length(PROT)){
      if(length(input$drug2) > 1){
        DR <- DAT_text()[which(!is.na(match(DAT_text()$id, PROT))),]
        DR$drug <- paste("has been identified in the experiment", DR$drug)
      }
      else{
        DR <- data.frame(id = PROT)
      }
    }
    else{
      NULL
    }
    if(!is.null(DR) & !is.null(Sel_cond_fhit_SUMMA$hit)){
      hit_info <- Sel_cond_fhit_SUMMA$hit
      hit_info <- hit_info[, !(names(hit_info) %in% "category")]
      names(hit_info)[!(names(hit_info) %in% "id")] <- "Hits_Info"

      DR <- left_join(DR, hit_info, by = "id")
    }
    if(!is.null(DR)){
      if(ncol(DR) <= 1){
        DR <- NULL
      }
    }
    unique(DR)
  })
  observe({
    tabident_bar$r <- tabidentreac_bar()
  })
  output$pr_info <- DT::renderDataTable({
    DT::datatable(tabident_bar$r,
                  caption = htmltools::tags$caption(
                    style = 'caption-side: top; text-align: left;',
                    htmltools::strong("Identification comparison")
                  ),
                  rownames = FALSE,
                  options = list(lengthMenu = c(10,20,30), pageLength = 10,
                                 scrollX = TRUE))
  })
  output$identifcomp_barup <- reactive({
    return(!is.null(tabident_bar$r))
  })
  outputOptions(output, "identifcomp_barup", suspendWhenHidden = FALSE)

  output$downtabidentif_barplot <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%y%m%d_%H%M_"), "Identification_comparison", ".xlsx")
    },
    content = function(file){
      openxlsx::write.xlsx(tabident_bar$r, file, row.names = FALSE)
    }
  )

  output$n_cond_sel <- renderText({
    if(input$ch_own_col){
      if (input$cond_sel  == "all_cond"){
        paste("You selected", length(get_treat_level(data())), "treatments, please enter the same number of colors")
      }
      else{
        paste("You selected", length(input$cond), "treatments, please enter the same number of colors")
      }
    }
    else{
      NULL
    }
  })

  OWN_color <- reactiveValues(
    ch = c()
  )
  observeEvent(input$add_col, {
    OWN_color$ch <- append(OWN_color$ch, input$own_color_pick)
  })
  observeEvent(input$rem_col, {
    if(length(OWN_color$ch) <= 1){
      OWN_color$ch <- c()
    }
    else{
      OWN_color$ch <- OWN_color$ch[1:(length(OWN_color$ch)-1)]
    }
  })
  output$own_color <- renderText({
    paste("You selected this colors :", paste(OWN_color$ch, collapse = ", "))
  })

  BAR <- reactiveValues(
    ch = NULL
  )

  Bar_one <- reactive({
    withCallingHandlers({
      shinyjs::html("diag_bar", "")
      if(input$ch_own_col){
        nbc <- ifelse(input$cond_sel == "all_cond", length(get_treat_level(data())), length(input$cond))
        COL <- OWN_color$ch
        if(nbc == length(COL)){
          imprints_barplotting_app(data(), witherrorbar = input$werb,
                                   withpoint = input$wpts,
                                   usegradient = input$grad, linegraph = input$line,
                                   save_pdf = input$save_bar, ret_plot = !input$save_bar,
                                   colorpanel = COL,
                                   layout = c(input$lay_bar1, input$lay_bar2),
                                   pdfname = input$pdftit,
                                   pdfwidth = input$pdfw, pdfheight = input$pdfh)
        }
        else{
          showNotification("The number of colors given doesn't match the number of treatment selected !", type = "error")
        }

      }
      else{
        imprints_barplotting_app(data(), witherrorbar = input$werb,
                                 withpoint = input$wpts,
                                 usegradient = input$grad, linegraph = input$line,
                                 save_pdf = input$save_bar, ret_plot = !input$save_bar,
                                 layout = c(input$lay_bar1, input$lay_bar2),
                                 pdfname = input$pdftit,
                                 pdfwidth = input$pdfw, pdfheight = input$pdfh)
      }

    },
    message = function(m) {
      shinyjs::html(id = "diag_bar", html = paste(m$message, "<br>", sep = ""), add = FALSE)

    }
    )

  })


  observeEvent(input$barp, {
    if(input$cond_sel == "cat"){
      if (length(unique(Sel_cond()$category)) > 1){
        che <- length(input$cond) == length(unique(Sel_cond()$category))
      }
      else{
        che <- FALSE
      }
    }
    else if(input$cond_sel == "all_cond"){
      che <- FALSE
    }

    che2 <- FALSE
    if (is.null(input$cond)){
      if(input$cond_sel != "all_cond"){
        che2 <- TRUE
      }
      else{
        che2 <- FALSE
      }
    }

    if(input$drug == "base" & length(input$drug2) < 1){
      showNotification("Don't forget to select a drug !", type = "error")
    }
    if(length(input$prot) == 0 & !input$ALL_prot){
      showNotification("Don't forget to select a protein !", type = "error")
    }
    else if (che2){
      showNotification("Don't forget to select a treatment !", type = "error")

    }
    else{
      BAR$ch <- Bar_one()
    }

  })

  output$bar_plot <- renderPlot({
    BAR$ch
  })

  output$downbar <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%y%m%d_%H%M_"), "2D_barplot_",  sub(":", "%", input$prot[length(input$prot)]), ".", input$downbar_format)
    },
    content = function(file){
      ggsave(file, plot = BAR$ch[[1]], device = input$downbar_format)
    }
  )




  ### PROTEIN COMPLEX
  output$drug2ui_compl <- renderUI({
    selectInput("drug2_compl", "Choose a drug", choices = names(drug_data_sh$y$data),
                multiple = TRUE, selected = "elutriation")
  })

  DIF_compl <- reactive({
    if(input$drug_compl == "dat"){
      File <- input$caldif_compl
      if (is.null(File))
        return(NULL)

      ms_fileread(File$datapath)
    }
    else if(input$drug_compl == "base" & length(input$drug2_compl) >= 1){
      join_drugdata(drug_data_sh$y$data[input$drug2_compl], by = c("id", "description"))
    }
    else{
      NULL
    }
  })
  #check if a file is upload
  output$DIFcompl_fileup <- reactive({
    return(!is.null(DIF_compl()))
  })
  outputOptions(output, "DIFcompl_fileup", suspendWhenHidden = FALSE)

  HIT_compl <- reactive({
    if(input$drug_compl == "dat"){
      File <- input$hitsum_compl
      if (is.null(File))
        return(NULL)

      dat <- import(File$datapath, header = TRUE)
      nv_nam <- str_subset(names(dat), "^V\\d{1}$")
      if(length(nv_nam)){
        dat <- dat[, !(names(dat) %in% nv_nam)]
      }
      if(!all(c("id", "treatment", "category") %in% colnames(dat))){
        missing_columns <- c("id", "treatment", "category")
        missing_columns <- missing_columns[!(c("id", "treatment", "category") %in% colnames(dat))]
        missing_columns <- paste(missing_columns, collapse = ", ")
        verb <- ifelse(length(missing_columns) > 1, "are", "is")
        showNotification(paste(missing_columns, verb, "not in your summary file. Please check your columns names !"),
                         type = "error", duration = 8)
        dat <- NULL
      }
      dat
    }
    else if(input$drug_compl == "base" & length(input$drug2_compl) >= 1){
      do.call(rbind, lapply(drug_data_sh$y$hitlist[input$drug2_compl],
                            function(x) x[,c("id", "treatment", "category")])
              )
    }
    else{
      NULL
    }
  })
  #check if a file is upload
  output$HITcompl_fileup <- reactive({
    return(!is.null(HIT_compl()))
  })
  outputOptions(output, "HITcompl_fileup", suspendWhenHidden = FALSE)

  NN_compl <- reactive({
    if(input$drug_compl == "dat"){
      nn <- NULL
      if(!is.null(HIT_compl()) & !is.null(DIF_compl())){
        dif <- DIF_compl()[,1:2]
        nn <- lapply(unique(HIT_compl()$treatment), function(x){
          x <- HIT_compl() %>% dplyr::filter(treatment == x) %>%
            dplyr::right_join(dif, by = "id") %>%
            dplyr::filter(is.na(category)) %>%
            dplyr::mutate(category = "NN",
                          treatment = x);
          x
        })
        nn <- as.data.frame(Reduce(rbind, nn))
        nn <- nn[,c("id", "description", "treatment", "category")]
        nn
      }
    }
    else if(input$drug_compl == "base" & length(input$drug2_compl) >= 1){
      do.call(rbind, lapply(drug_data_sh$y$NN[input$drug2_compl],
                            function(x) x[,c("id", "description", "treatment", "category")])
              )
    }
    else{
      NULL
    }
  })
  #check if a file is upload
  output$NNcompl_fileup <- reactive({
    return(!is.null(NN_compl()))
  })
  outputOptions(output, "NNcompl_fileup", suspendWhenHidden = FALSE)

  AVE_compl <- reactive({
    if(input$drug_compl == "dat"){
      File <- input$avef_compl
      if (is.null(File))
        return(NULL)

      ms_fileread(File$datapath)
    }
    else if(input$drug_compl == "base" & length(input$drug2_compl) >= 1){
      join_drugdata(drug_data_sh$y$data_ave[input$drug2_compl], by = c("id", "description"))
    }
    else{
      NULL
    }
  })
  #check if a file is upload
  output$AVEcompl_fileup <- reactive({
    if(input$gave_compl){
      return(TRUE)
    }
    else{
      return(!is.null(AVE_compl()))
    }
  })
  outputOptions(output, "AVEcompl_fileup", suspendWhenHidden = FALSE)

  observe({
    if(!is.null(HIT_compl()) & !is.null(NN_compl())){
      updateSelectInput(session, "condsel_compl", choices = unique(HIT_compl()$treatment))
      updateSelectInput(session, "catego_compl", choices = append(unique(HIT_compl()$category),  "NN"),
                        selected = unique(HIT_compl()$category)[1])
    }
  })

  resmapping_compl <- reactiveValues(
    ch = NULL
  )
  observeEvent(input$ave_map_compl, {
    if(length(input$catego_compl) == 0){
      showNotification("Don't forget to select a category !", type = "error")
    }
    else {
      showNotification("Start mapping proteins, this may take a while", type = "message")

      if(input$gave_compl & input$drug_compl == "dat"){
        data_ave <- imprints_average(DIF_compl(), savefile = TRUE)
        showNotification("Average calculation succeed !", type = "message")
      }
      else{
        data_ave <- AVE_compl()
      }

      cat_tab <- HIT_compl() %>% dplyr::group_by(id, treatment, category) %>% dplyr::reframe()

      cat_tabNN <- NN_compl()
      cat_tabNN <- cat_tabNN %>% dplyr::group_by(id, treatment, category) %>% dplyr::reframe()

      cat_tab <- rbind(cat_tab, cat_tabNN)

      withCallingHandlers({
        shinyjs::html("diagmapping_compl", "")
        map_compl <- imprints_complex_mapping(data_ave, cat_tab, treatment = input$condsel_compl,
                                              targetcategory = input$catego_compl,
                                              organism = input$organism_compl)
      },
      message = function(m) {
        shinyjs::html(id = "diagmapping_compl", html = paste(m$message, "<br>", sep = ""), add = TRUE)
      }
      )

      map_compl <- map_compl[, c("ComplexName", "subunitsNum", "subunitsIdentifiedNum",
                                 "id", "description", "gene", "category")]

      if(nrow(map_compl) !=0){
        map_compl$description <- unname(sapply(map_compl$description, IMPRINTS.CETSA.app:::getProteinName))
      }


      resmapping_compl$ch <- map_compl
      output$tabmap_compl <- DT::renderDataTable({
        DT::datatable(resmapping_compl$ch,
                      caption = htmltools::tags$caption(
                        style = 'caption-side: top; text-align: left;',
                        htmltools::strong("Mapping proteins results")
                      ),
                      rownames = FALSE,
                      options = list(lengthMenu = c(10,20,30), pageLength = 10,
                                     scrollX = TRUE))
      })

      if(nrow(map_compl) !=0){
        showNotification("Mapping protein succeed !", type = "message")
      }
      else{
        showNotification("No proteins could be mapped !
                         Try to add more category in order to have more proteins", type = "error")
      }


      updateSelectInput(session, "allcomplex_compl", choices = unique(resmapping_compl$ch$ComplexName))

    }

  })
  #check if a file is upload
  output$resmappingcompl_fileup <- reactive({
    return(!is.null(resmapping_compl$ch))
  })
  outputOptions(output, "resmappingcompl_fileup", suspendWhenHidden = FALSE)

  output$downrestab_compl <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%y%m%d_%H%M_"), "ProteinComplexMapping", ".xlsx")
    },
    content = function(file){
      openxlsx::write.xlsx(resmapping_compl$ch, file, row.names = FALSE)
    }
  )


  sel_prot_compl <- reactive({
    pr <- NULL
    if(!is.null(resmapping_compl$ch)){
      pr <- resmapping_compl$ch[which(!is.na(match(resmapping_compl$ch$ComplexName, input$allcomplex_compl))), c("id", "gene")]
      pr$id <- paste0(pr$id, ":", pr$gene)
      pr <- unique(pr$id)
    }

  })
  observe({
    updateSelectizeInput(session, "prot_compl", choices = sel_prot_compl(), server = TRUE)

  })

  data_compl <- reactive({
    if(input$ALL_prot_compl){
      PROT <- sel_prot_compl()
    }
    else{
      PROT <- input$prot_compl
    }
    if(!is.null(PROT)){
      PROT <- unname(sapply(PROT, function(x) strsplit(x, ":")[[1]][1]))
    }

    data <- DIF_compl()
    TREAT <- get_treat_level(data)

    cate <- resmapping_compl$ch[which(!is.na(match(resmapping_compl$ch$ComplexName, input$allcomplex_compl))),]
    notsel_cond <- TREAT[!(TREAT %in% input$condsel_compl)]
    notsel_cond <- paste(notsel_cond, collapse = "|")

    if(input$save_bar_compl){
      data_l <- list()
      for(i in input$allcomplex_compl){
        cate_ <- cate[which(cate$ComplexName == i), ]

        pr_comp <- cate_$id
        pr_comp <- pr_comp[which(!is.na(match(pr_comp, PROT)))]

        data_l[[i]] <- data[which(!is.na(match(data$id, pr_comp))),]

        data_l[[i]] <- data_l[[i]][,-str_which(names(data_l[[i]]), notsel_cond)]
        data_l[[i]] <- dplyr::left_join(data_l[[i]], unique(cate_[,c("id", "category")]),
                                        by = "id")
      }

      data <- data_l

    }
    else{
      data <- data[which(!is.na(match(data$id, PROT))),]

      data <- data[,-str_which(names(data), notsel_cond)]
      data <- dplyr::left_join(data, unique(cate[,c("id", "category")]),
                                        by = "id")
    }

    data
  })

  BAR_compl <- reactiveValues(
    ch = NULL
  )

  Bar_one_compl <- reactive({
    withCallingHandlers({
      shinyjs::html("diag_bar_compl", "")

      COL <- ifelse(input$ch_own_col_compl, input$own_color_pick_compl, "#18FF00")


      imprints_barplotting_app(data_compl(), witherrorbar = input$werb_compl,
                               withpoint = input$wpts_compl,
                               usegradient = input$grad_compl, linegraph = input$line_compl,
                               save_pdf = input$save_bar_compl, colorpanel = COL,
                               ret_plot = !input$save_bar_compl,
                               layout = c(input$lay_bar1_compl, input$lay_bar2_compl),
                               toplabel = "IMPRINTS-CETSA bar plotting \nProtein complex :",
                               pdfname = input$pdftit_compl,
                               pdfwidth = input$pdfw_compl, pdfheight = input$pdfh_compl
                               )

    },
    message = function(m) {
      shinyjs::html(id = "diag_bar_compl", html = paste(m$message, "<br>", sep = ""), add = FALSE)

    }
    )

  })


  observeEvent(input$barp_compl, {
    if(length(input$prot_compl) == 0 & !input$ALL_prot_compl){
      showNotification("Don't forget to select a protein !", type = "error")
    }
    else{
      BAR_compl$ch <- Bar_one_compl()
    }

  })

  output$bar_plot_compl <- renderPlot({
    BAR_compl$ch
  })

  output$downbar_compl <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%y%m%d_%H%M_"), "2D_barplot_", paste(str_remove_all(input$allcomplex_compl, " "), sep = "_"), ".", input$downbar_compl_format)
    },
    content = function(file){
      ggsave(file, plot = BAR_compl$ch[[1]], device = input$downbar_compl_format)
    }
  )




  ### SIMILAR PROFILE
  output$drug2ui_simpf <- renderUI({
    selectInput("drug2_simpf", "Choose a drug", choices = names(drug_data_sh$y$data),
                multiple = TRUE, selected = "elutriation")
  })

  DIF_simpf <- reactive({
    if(input$drug_simpf == "dat"){
      File <- input$cdiff_simpf
      if (is.null(File))
        return(NULL)

      ms_fileread(File$datapath)
    }
    else if(input$drug_simpf == "base" & length(input$drug2_simpf) >= 1){
      join_drugdata(drug_data_sh$y$data[input$drug2_simpf], by = c("id", "description"))
    }
    else{
      NULL
    }
  })
  #check if a file is upload
  output$DIFsimpf_fileup <- reactive({
    return(!is.null(DIF_simpf()))
  })
  outputOptions(output, "DIFsimpf_fileup", suspendWhenHidden = FALSE)

  AVE_simpf <- reactive({
    if(input$drug_simpf == "dat"){
      File <- input$avef_simpf
      if (is.null(File))
        return(NULL)

      ms_fileread(File$datapath)
    }
    else if(input$drug_simpf == "base" & length(input$drug2_simpf) >= 1){
      join_drugdata(drug_data_sh$y$data_ave[input$drug2_simpf], by = c("id", "description"))
    }
    else{
      NULL
    }
  })
  #check if a file is upload
  output$AVEsimpf_fileup <- reactive({
    if(input$gave_simpf){
      return(TRUE)
    }
    else{
      return(!is.null(AVE_simpf()))
    }
  })
  outputOptions(output, "AVEsimpf_fileup", suspendWhenHidden = FALSE)

  observe({
    if(!is.null(DIF_simpf())){
      updateSelectInput(session, "treat_simpf", choices = get_treat_level(DIF_simpf()))
      x <- DIF_simpf()[,c("id", "description")]
      x$description <- unname(sapply(x$description, IMPRINTS.CETSA.app:::getGeneName))
      x <- x %>%
        group_by(description) %>%
        mutate(id = paste0(id, ":", description))
      updateSelectizeInput(session, "prot_simpf", choices = unique(x$id), server = TRUE)
    }
  })
  observe({
    if(!is.null(DIF_simpf())){
      nc <- str_subset(names(DIF_simpf()), paste0("_", input$treat_simpf, "$"))
      nc <- str_split(nc, "B\\d{1}_")
      nc <- lapply(nc, function(x) paste(x, collapse = ""))
      nc <- length(unique(as.character(nc)))
      updateSliderInput(session, "maxna_simpf", max = nc)
    }
  })


  BAR_simpf <- reactiveValues(
    ch = NULL
  )

  Bar_one_simpf <- reactive({
    if(input$gave_simpf & input$drug_simpf == "dat"){
      average <- NULL
    }
    else{
      average <- AVE_simpf()
    }
    COL <- ifelse(input$ch_own_col_simpf, input$own_color_pick_simpf, "#18FF00")

    withCallingHandlers({
      shinyjs::html("diag_bar_simpf", "")
      imprints_barplotting_simprof(DIF_simpf(), average, witherrorbar = input$werb_simpf,
                                treatmentlevel = input$treat_simpf, protein_profile = strsplit(input$prot_simpf, ":")[[1]][1],
                                usegradient = input$grad_simpf, linegraph = input$line_simpf,
                                use_score = input$scoremeth_simpf, score_threshold = input$scothr_simpf,
                                max_na_prow = input$maxna_simpf,
                                ret_plot = input$seeprsel_simpf, save_pdf = input$save_bar_simpf,
                                colorpanel = COL, withprompt = FALSE, save_prlist = input$save_prot_simpf,
                                layout = c(input$lay_bar1_simpf, input$lay_bar2_simpf),
                                toplabel = paste0("IMPRINTS-CETSA bar plotting \nMethod :", input$scoremeth_simpf),
                                pdfname = input$pdftit_simpf)


    },
    message = function(m) {
      shinyjs::html(id = "diag_bar_simpf", html = paste(m$message, "<br>", sep = ""), add = FALSE)

    }
    )

  })


  geting_data_simpf <- reactiveValues(
    ch = NULL
  )

  observeEvent(input$getsimi_simpf, {
    showNotification("Getting similar profiles, this may take a while.", type = "message")

    if(input$gave_simpf & input$drug_simpf == "dat"){
      average <- NULL
    }
    else{
      average <- AVE_simpf()
    }
    COL <- ifelse(input$ch_own_col_simpf, input$own_color_pick_simpf, "#18FF00")

    withCallingHandlers({
      shinyjs::html("diag_bar_simpf", "")
      geting_data_simpf$ch <- imprints_barplotting_simprof(DIF_simpf(), average, witherrorbar = input$werb_simpf,
                                                        treatmentlevel = input$treat_simpf, protein_profile = strsplit(input$prot_simpf, ":")[[1]][1],
                                                        usegradient = input$grad_simpf, linegraph = input$line_simpf,
                                                        use_score = input$scoremeth_simpf, score_threshold = input$scothr_simpf,
                                                        max_na_prow = input$maxna_simpf,
                                                        ret_plot = input$seeprsel_simpf, save_pdf = FALSE,
                                                        withpopup = TRUE, continue = FALSE, modvar = "",
                                                        colorpanel = COL,  save_prlist = FALSE,
                                                        layout = c(input$lay_bar1_simpf, input$lay_bar2_simpf),
                                                        toplabel = paste0("IMPRINTS-CETSA bar plotting \nMethod :", input$scoremeth_simpf),
                                                        pdfname = input$pdftit_simpf)



    },
    message = function(m) {
      shinyjs::html(id = "diag_bar_simpf", html = paste(m$message, "<br>", sep = ""), add = FALSE)

    }
    )



  })

  observeEvent(input$ok, {
    removeModal()
    COL <- ifelse(input$ch_own_col_simpf, input$own_color_pick_simpf, "#18FF00")

    withCallingHandlers({
      shinyjs::html("diag_bar_simpf", "")
      BAR_simpf$ch <- imprints_barplotting_simprof(geting_data_simpf$ch, witherrorbar = input$werb_simpf,
                                                treatmentlevel = input$treat_simpf, protein_profile = strsplit(input$prot_simpf, ":")[[1]][1],
                                                usegradient = input$grad_simpf, linegraph = input$line_simpf,
                                                use_score = input$scoremeth_simpf, score_threshold = input$scothr_simpf,
                                                max_na_prow = input$maxna_simpf,
                                                ret_plot = input$seeprsel_simpf, save_pdf = input$save_bar_simpf,
                                                withpopup = TRUE, continue = FALSE, modvar = "Y", got_it = TRUE,
                                                colorpanel = COL, save_prlist = input$save_prot_simpf,
                                                layout = c(input$lay_bar1_simpf, input$lay_bar2_simpf),
                                                toplabel = paste0("IMPRINTS-CETSA bar plotting \nMethod :", input$scoremeth_simpf),
                                                pdfname = input$pdftit_simpf)
    },
    message = function(m) {
      shinyjs::html(id = "diag_bar_simpf", html = paste(m$message, "<br>", sep = ""), add = FALSE)

    }
    )
  })
  observeEvent(input$cancel, {
    removeModal()
    COL <- ifelse(input$ch_own_col_simpf, input$own_color_pick_simpf, "#18FF00")

    withCallingHandlers({
      shinyjs::html("diag_bar_simpf", "")
      BAR_simpf$ch <- imprints_barplotting_simprof(geting_data_simpf$ch, witherrorbar = input$werb_simpf,
                                                treatmentlevel = input$treat_simpf, protein_profile = strsplit(input$prot_simpf, ":")[[1]][1],
                                                usegradient = input$grad_simpf, linegraph = input$line_simpf,
                                                use_score = input$scoremeth_simpf, score_threshold = input$scothr_simpf,
                                                max_na_prow = input$maxna_simpf,
                                                ret_plot = input$seeprsel_simpf, save_pdf = FALSE,
                                                withpopup = TRUE, continue = FALSE, modvar = "N", got_it = TRUE,
                                                colorpanel = COL, save_prlist = FALSE,
                                                layout = c(input$lay_bar1_simpf, input$lay_bar2_simpf),
                                                toplabel = paste0("IMPRINTS-CETSA bar plotting \nMethod :", input$scoremeth_simpf),
                                                pdfname = input$pdftit_simpf)
    },
    message = function(m) {
      shinyjs::html(id = "diag_bar_simpf", html = paste(m$message, "<br>", sep = ""), add = FALSE)

    }
    )
  })

  output$bar_plot_simpf <- renderPlot({
    BAR_simpf$ch
  })

  output$downbar_simpf <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%y%m%d_%H%M_"), "2D_barplot_", paste0("similar_", sub(":", "%", input$prot_simpf)), ".", input$downbar_simpf_format)
    },
    content = function(file){
      ggsave(file, plot = BAR_simpf$ch[[1]], device = input$downbar_simpf_format)
    }
  )


  ### HEATMAP
  output$drug2ui_heat <- renderUI({
    selectInput("drug2_heat", "Choose a drug", choices = names(drug_data_sh$y$data),
                multiple = TRUE, selected = "elutriation")
  })

  DIF_heat <- reactive({
    if(input$drug_heat == "dat"){
      File <- input$filedif_heat
      if (is.null(File) | !input$gave_heat)
        return(NULL)

      ms_fileread(File$datapath)
    }
    else if(input$drug_heat == "base" & length(input$drug2_heat) >= 1){
      join_drugdata(drug_data_sh$y$data[input$drug2_heat], by = c("id", "description"))
    }
    else{
      NULL
    }
  })

  AVE_heat <- reactive({
    if(input$drug_heat == "dat"){
      File <- input$fileave_heat
      if (is.null(File)  | input$gave_heat)
        return(NULL)

      ms_fileread(File$datapath)
    }
    else if(input$drug_heat == "base" & length(input$drug2_heat) >= 1){
      join_drugdata(drug_data_sh$y$data_ave[input$drug2_heat], by = c("id", "description"))
    }
    else{
      NULL
    }
  })
  #check if a file is upload
  output$heat_fileup <- reactive({
    return(!is.null(AVE_heat()) | !is.null(DIF_heat()))
  })
  outputOptions(output, "heat_fileup", suspendWhenHidden = FALSE)

  HIT_heat <- reactive({
    if(input$drug_heat == "dat"){
      File <- input$summary_heat
      if (is.null(File))
        return(NULL)

      dat <- import(File$datapath, header = TRUE)
      nv_nam <- str_subset(names(dat), "^V\\d{1}$")
      if(length(nv_nam)){
        dat <- dat[, !(names(dat) %in% nv_nam)]
      }
      if(!all(c("id", "treatment", "category") %in% colnames(dat))){
        missing_columns <- c("id", "treatment", "category")
        missing_columns <- missing_columns[!(c("id", "treatment", "category") %in% colnames(dat))]
        missing_columns <- paste(missing_columns, collapse = ", ")
        verb <- ifelse(length(missing_columns) > 1, "are", "is")
        showNotification(paste(missing_columns, verb, "not in your summary file. Please check your columns names !"),
                         type = "error", duration = 8)
        dat <- NULL
      }
      dat
    }
    else if(input$drug_heat == "base" & length(input$drug2_heat) >= 1){
      do.call(rbind, lapply(drug_data_sh$y$hitlist[input$drug2_heat],
                            function(x) x[,c("id", "treatment", "category")])
              )
    }
    else{
      NULL
    }
  })
  #check if a file is upload
  output$HITheat_fileup <- reactive({
    return(!is.null(HIT_heat()))
  })
  outputOptions(output, "HITheat_fileup", suspendWhenHidden = FALSE)

  NN_heat <- reactive({
    if(input$drug_heat == "dat"){
      nn <- NULL
      if(!is.null(HIT_heat()) & !is.null(DIF_heat())){
        dif <- DIF_heat()[,1:2]
        nn <- lapply(unique(HIT_heat()$treatment), function(x){
          x <- HIT_heat() %>% dplyr::filter(treatment == x) %>%
            dplyr::right_join(dif, by = "id") %>%
            dplyr::filter(is.na(category)) %>%
            dplyr::mutate(category = "NN",
                          treatment = x);
          x
        })
        nn <- as.data.frame(Reduce(rbind, nn))
        nn <- nn[,c("id", "description", "treatment", "category")]
        nn
      }
    }
    else if(input$drug_heat == "base" & length(input$drug2_heat) >= 1){
      do.call(rbind, lapply(drug_data_sh$y$NN[input$drug2_heat],
                            function(x) x[,c("id", "description", "treatment", "category")])
              )
    }
    else{
      NULL
    }
  })
  #check if a file is upload
  output$NNheat_fileup <- reactive({
    if(input$drug_heat == "dat"){
      return(!is.null(NN_heat()) & !is.null(HIT_heat()))
    }
    else{
      return(!is.null(NN_heat()))
    }
  })
  outputOptions(output, "NNheat_fileup", suspendWhenHidden = FALSE)

  observe({
    if(!is.null(HIT_heat())){
      c_idx <- str_which(colnames(HIT_heat()), "treatment")
      cat_idx <- str_which(colnames(HIT_heat()), "^[C|c]ategory")

      tr <- NULL
      if(length(c_idx)){
        tr <- HIT_heat()[, c_idx]
        tr <- unique(tr)
      }
      cat <- c()
      if(length(cat_idx)){
        cat <- HIT_heat()[, cat_idx]
        cat <- unique(cat)
      }
      if(!is.null(NN_heat())){
        cat <- append(cat, "NN")
      }

      updateSelectInput(session, "cond_heat", choices = tr)
      updateSelectInput(session, "catego_heat", choices = cat, selected = cat[1])
    }
  })
  observe({
    if(!is.null(DIF_heat()) & input$drug_heat == "dat"){
      nc <- str_subset(names(DIF_heat()), paste0("_", input$cond_heat, "$"))
      nc <- str_split(nc, "B\\d{1}_")
      nc <- lapply(nc, function(x) paste(x, collapse = ""))
      nc <- length(unique(as.character(nc)))
      updateSliderInput(session, "maxna_heat", max = nc)
    }
    if(!is.null(AVE_heat())){
      nc <- str_subset(names(AVE_heat()), paste0("_", input$cond_heat, "$"))
      nc <- length(unique(nc))
      updateSliderInput(session, "maxna_heat", max = nc)
    }
  })

  pH_heat <- reactiveValues(
    g = NULL
  )

  plotH_heat <- reactive({
    dat <- NULL
    if(!is.null(AVE_heat())){
      dat <- AVE_heat()
    }
    else if(!is.null(DIF_heat()) & input$drug_heat == "dat"){
      showNotification("Start average calculation, this mays take a while.", type = "message")
      dat <- imprints_average(DIF_heat(), savefile = TRUE)
    }

    withCallingHandlers({
      shinyjs::html("diagl_heat", "")
      h <- imprints_heatmap(dat, HIT_heat(), NN_data = NN_heat(),
                         treatment = input$cond_heat, max_na = input$maxna_heat,
                         response = input$resp_heat, select_cat = input$catego_heat,
                         gradient_color = c(input$grad1col_heat, input$grad2col_heat, input$grad3col_heat),
                         titleH = input$titleH_heat,
                         saveHeat = input$saveH_heat, file_type = input$formatH_heat, file_name = input$fnameH_heat,
                         cat_color = NULL, back_color = input$backcol_heat, border_color = input$bordercol_heat)
    },
    message = function(m) {
      shinyjs::html(id = "diagl_heat", html = paste(m$message, "<br>", sep = ""), add = TRUE)

    }
    )

  })

  observeEvent(input$getH_heat, {
    if(is.null(input$catego_heat)){
      showNotification("Don't forget to select a category !", type = "error")
    }
    else{
      showNotification("Getting heatmap", type = "message")
      pH_heat$g <- plotH_heat()
    }
  })

  output$H_heat <- renderPlot({
    if(!is.null(pH_heat$g))
      return(plot(pH_heat$g))
    else
      NULL
  })



  ### HEATMAP PROTEIN COMPLEX
  output$drug2ui_heatcom <- renderUI({
    selectInput("drug2_heatcom", "Choose a drug", choices = names(drug_data_sh$y$data),
                multiple = TRUE, selected = "elutriation")
  })

  DIF_heatcom <- reactive({
    if(input$drug_heatcom == "dat"){
      File <- input$filedif_heatcom
      if (is.null(File) | !input$gave_heatcom)
        return(NULL)

      ms_fileread(File$datapath)
    }
    else if(input$drug_heatcom == "base" & length(input$drug2_heatcom) >= 1){
      join_drugdata(drug_data_sh$y$data[input$drug2_heatcom], by = c("id", "description"))
    }
    else{
      NULL
    }
  })

  AVE_heatcom <- reactive({
    if(input$drug_heatcom == "dat"){
      File <- input$fileave_heatcom
      if (is.null(File)  | input$gave_heatcom)
        return(NULL)

      ms_fileread(File$datapath)
    }
    else if(input$drug_heatcom == "base" & length(input$drug2_heatcom) >= 1){
      join_drugdata(drug_data_sh$y$data_ave[input$drug2_heatcom], by = c("id", "description"))
    }
    else{
      NULL
    }
  })
  #check if a file is upload
  output$heatcom_fileup <- reactive({
    return(!is.null(AVE_heatcom()) | !is.null(DIF_heatcom()))
  })
  outputOptions(output, "heatcom_fileup", suspendWhenHidden = FALSE)


  observe({
    if(!is.null(DIF_heatcom()) | !is.null(AVE_heatcom())){
      tr <- get_treat_level(DIF_heatcom())
      if(is.null(tr)){
        tr <- get_treat_level(AVE_heatcom())
      }

      updateSelectInput(session, "cond_heatcom", choices = tr)
    }
  })


  resmapping_heatcom <- reactiveValues(
    ch = NULL
  )
  resAVE_heatcom <- reactiveValues(
    d = NULL
  )
  observeEvent(input$ave_map_heatcom, {
    showNotification("Start mapping proteins, this may take a while", type = "message")

    if(input$gave_heatcom  & input$drug_heatcom == "dat"){
      data_ave <- imprints_average(DIF_heatcom(), savefile = TRUE)
      resAVE_heatcom$d <- data_ave
      showNotification("Average calculation succeed !", type = "message")
    }
    else{
      data_ave <- AVE_heatcom()
    }

    withCallingHandlers({
      shinyjs::html("diagmapping_heatcom", "")
      map_heatcom <- imprints_complex_mapping(data_ave, categorytable = NULL, treatment = input$cond_heatcom,
                                              targetcategory = NULL,
                                              organism = input$organism_heatcom)
    },
    message = function(m) {
      shinyjs::html(id = "diagmapping_heatcom", html = paste(m$message, "<br>", sep = ""), add = TRUE)
    }
    )

    map_heatcom <- map_heatcom[, c("ComplexName", "subunitsNum", "subunitsIdentifiedNum",
                                   "id", "description", "gene")]

    if(nrow(map_heatcom) !=0){
      map_heatcom$description <- IMPRINTS.CETSA.app:::getProteinName(map_heatcom$description)
    }


    resmapping_heatcom$ch <- map_heatcom
    output$tabmap_heatcom <- DT::renderDataTable({
      DT::datatable(resmapping_heatcom$ch,
                    caption = htmltools::tags$caption(
                      style = 'caption-side: top; text-align: left;',
                      htmltools::strong("Mapping proteins results")
                    ),
                    rownames = FALSE,
                    options = list(lengthMenu = c(10,20,30), pageLength = 10,
                                   scrollX = TRUE))
    })

    if(nrow(map_heatcom) != 0){
      showNotification("Mapping protein succeed !", type = "message")
    }
    else{
      showNotification("No proteins could be mapped !
                         Try to add more categories in order to have more proteins", type = "error")
    }


    updateSelectInput(session, "allcomplex_heatcom", choices = unique(resmapping_heatcom$ch$ComplexName))
  })
  #check if a file is upload
  output$resmappingheatcom_fileup <- reactive({
    return(!is.null(resmapping_heatcom$ch))
  })
  outputOptions(output, "resmappingheatcom_fileup", suspendWhenHidden = FALSE)

  output$downrestab_heatcom <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%y%m%d_%H%M_"), "ProteinComplexMapping", ".xlsx")
    },
    content = function(file){
      openxlsx::write.xlsx(resmapping_heatcom$ch, file, row.names = FALSE)
    }
  )


  observe({
    if(!is.null(DIF_heatcom()) & input$drug_heatcom == "dat"){
      nc <- str_subset(names(DIF_heatcom()), paste0("_", input$cond_heatcom, "$"))
      nc <- str_split(nc, "B\\d{1}_")
      nc <- lapply(nc, function(x) paste(x, collapse = ""))
      nc <- length(unique(as.character(nc)))
      updateSliderInput(session, "maxna_heatcom", max = nc)
    }
    if(!is.null(AVE_heatcom())){
      nc <- str_subset(names(AVE_heatcom()), paste0("_", input$cond_heatcom, "$"))
      nc <- length(unique(nc))
      updateSliderInput(session, "maxna_heatcom", max = nc)
    }
  })

  pH_heatcom <- reactiveValues(
    g = NULL
  )

  plotH_heatcom <- reactive({
    dat <- NULL
    if(!is.null(AVE_heatcom())){
      dat <- AVE_heatcom()
    }
    else if(!is.null(resAVE_heatcom$d)){
      dat <- resAVE_heatcom$d
    }

    pr <- NULL
    if(!is.null(resmapping_heatcom$ch)){
      pr <- resmapping_heatcom$ch$id[which(!is.na(match(resmapping_heatcom$ch$ComplexName, input$allcomplex_heatcom)))]
      pr <- unique(pr)
    }

    dat <- dat[which(!is.na(match(dat$id, pr))),]
    PRcompl <- resmapping_heatcom$ch[which(!is.na(match(resmapping_heatcom$ch$ComplexName, input$allcomplex_heatcom))),]

    withCallingHandlers({
      shinyjs::html("diagl_heatcom", "")
      h <- imprints_heatmap(dat, NULL, NN_data = NULL, PRcomplex_data = PRcompl,
                         treatment = input$cond_heatcom, max_na = input$maxna_heatcom,
                         response = input$resp_heatcom,
                         gradient_color = c(input$grad1col_heatcom, input$grad2col_heatcom, input$grad3col_heatcom),
                         titleH = input$titleH_heatcom,
                         saveHeat = input$saveH_heatcom, file_type = input$formatH_heatcom, file_name = input$fnameH_heatcom,
                         cat_color = NULL, back_color = input$backcol_heatcom, border_color = input$bordercol_heatcom)
    },
    message = function(m) {
      shinyjs::html(id = "diagl_heatcom", html = paste(m$message, "<br>", sep = ""), add = TRUE)

    }
    )

  })

  observeEvent(input$getH_heatcom, {
    if(is.null(input$allcomplex_heatcom)){
      showNotification("Don't forget to select a complex !", type = "error")
    }
    else{
      showNotification("Getting heatmap", type = "message")
      pH_heatcom$g <- plotH_heatcom()
    }
  })

  output$H_heatcom <- renderPlot({
    if(!is.null(pH_heatcom$g))
      return(plot(pH_heatcom$g))
    else
      NULL
  })



  ### STRINGdb
  output$drug2ui_stri <- renderUI({
    selectInput("drug2_stri", "Choose a drug", choices = names(drug_data_sh$y$data),
                multiple = TRUE, selected = "elutriation")
  })

  stri_data <- reactive({
    if(input$drug_stri == "dat"){
      if(input$impfile_stri){
        File <- input$file_stri
        if (is.null(File))
          return(NULL)

        import_list(File$datapath, header = TRUE)[[1]]
      }
      else{
        if (str_length(input$txt_stri) == 0)
          return(NULL)

        i <- str_remove_all(i, " ")
        i <- str_split(input$txt_stri, ",")[[1]]

        data.frame(id = i)
      }
    }
    else if(input$drug_stri == "base" & length(input$drug2_stri) >= 1){
      h <- do.call(rbind, lapply(drug_data_sh$y$hitlist[input$drug2_stri],
                                 function(x) x[,c("id", "treatment", "category")])
                   )
      n <- do.call(rbind, lapply(drug_data_sh$y$NN[input$drug2_stri],
                                 function(x) x[,c("id", "description", "treatment", "category")])
                   )
      n <- unique(n[,c("id", "treatment", "category")])

      rbind(h,n)
    }
    else{
      NULL
    }
  })
  #check if a file is upload
  output$file_stri_up <- reactive({
    return(!is.null(stri_data()))
  })
  outputOptions(output, "file_stri_up", suspendWhenHidden = FALSE)

  Sel_cond_fhit_stri <- reactive({
    tr <- NULL

    if((input$ishit_stri & input$impfile_stri) | input$drug_stri == "base"){
      HIT <- stri_data()

      c_idx <- str_which(colnames(HIT), "treatment")

      if(length(c_idx)){
        tr <- HIT[, c_idx]
        tr <- unique(tr)
      }
    }

    tr
  })

  observe({
    if(input$drug_stri =="dat"){
      updateSelectInput(session, "cond_fhit_stri", choices = Sel_cond_fhit_stri(), selected = Sel_cond_fhit_stri()[1])
    }
    else if(input$drug_stri =="base"){
      updateSelectInput(session, "cond_fhitB_stri", choices = Sel_cond_fhit_stri(), selected = Sel_cond_fhit_stri()[1])
    }
  })

  string_res <- reactiveValues(
    x = NULL
  )
  observeEvent(input$start_string, {
    showNotification("Getting the STRING ids, this may take a while", type = "message")

    if(!file.exists("STRING_data")){
      dir.create("STRING_data")
    }

    if(input$species_string == 9606){
      if(!exists("string_db_human")){
        string_db_human <<- STRINGdb$new(version="11.5", species=9606,               #ID 9606 correspond to human
                                         score_threshold=200,
                                         input_directory=  file.path(getwd(), "STRING_data"))
      }
      string_db <<- string_db_human
    }
    else if(input$species_string == 10090){
      if(!exists("string_db_mouse")){
        string_db_mouse <<- STRINGdb$new(version="11.5", species=10090,               #ID 10090 correspond to mouse
                                         score_threshold=200,
                                         input_directory=  file.path(getwd(), "STRING_data"))
      }
      string_db <<- string_db_mouse
    }
    else if(input$species_string == 10116){
      if(!exists("string_db_rat")){
        string_db_rat <<- STRINGdb$new(version="11.5", species=10116,               #ID 10116 correspond to rat
                                       score_threshold=200,
                                       input_directory=  file.path(getwd(), "STRING_data"))
      }
      string_db <<- string_db_rat
    }


    string_res$x <- NULL

    if (!is.null(stri_data())){
      if(input$drug_stri == "dat"){
        if(input$impfile_stri){
          if(input$ishit_stri){
            dat <- stri_data()
            if(!is.null(input$cond_fhit_stri)){
              dat <- dat %>% dplyr::filter(!is.na(match(treatment, c(input$cond_fhit_stri))))
              if(any(duplicated(dat$id))){
                dat <- dat[-which(duplicated(dat$id)),]
              }
              a <- string_db$map(dat, "id", removeUnmappedRows = TRUE)
            }
            else{
              showNotification("Don't forget to select some treatments !", type = "error")
              a <- NULL
            }
          }
          else{
            dat <- stri_data()
            if(any(duplicated(dat[[input$idfile_stri]]))){
              dat <- dat[-which(duplicated(dat[[input$idfile_stri]])),]
            }
            a <- string_db$map(dat, input$idfile_stri, removeUnmappedRows = TRUE)
          }
        }
        else{
          dat <- stri_data()
          if(any(duplicated(dat$id))){
            dat <- dat[-which(duplicated(dat$id)),]
          }
          a <- string_db$map(dat, "id", removeUnmappedRows = TRUE)
        }
      }
      else if(input$drug_stri == "base"){
        dat <- stri_data()
        if(!is.null(input$cond_fhitB_stri)){
          dat <- dat %>% dplyr::filter(!is.na(match(treatment, c(input$cond_fhitB_stri))))
          if(!is.null(input$cat_fhitB_stri)){
            dat <- dat %>% dplyr::filter(!is.na(match(category, c(input$cat_fhitB_stri))))
          }
          if(any(duplicated(dat$id))){
            dat <- dat[-which(duplicated(dat$id)),]
          }
          a <- string_db$map(dat, "id", removeUnmappedRows = TRUE)
        }
        else{
          showNotification("Don't forget to select some treatments !", type = "error")
          a <- NULL
        }
      }
    }

    string_res$x <- a
    if(!is.null(string_res$x)){ # check that no duplicate string ids per id
      if(any(duplicated(string_res$x[[1]]))){ # first column returned by map from STRINGdb is the identifier used to map genes
        showNotification("Some ids returned several STRINGids; trying to map genes if in data",
                         type = "error", duration = 5)
        duplicated_ids <- unique(string_res$x[[1]][which(duplicated(string_res$x[[1]]))])
        for(d_id in duplicated_ids){
          info_d_ids <- string_res$x[which(string_res$x[[1]] == d_id),]

          string_res$x <- string_res$x[-which(string_res$x[[1]] == d_id),]

          where_gene <- grep("^[Gg]en(e$|es$)", colnames(info_d_ids))
          where_descr <- grep("^[dD]escription$", colnames(info_d_ids))

          if(length(where_gene) == 1){
            info_d_ids$STRING_id <- NULL
            info_d_ids <- info_d_ids[1,]
            info_d_ids <- string_db$map(info_d_ids, colnames(info_d_ids)[where_gene], removeUnmappedRows = TRUE)
          }
          else if(length(where_descr) == 1){
            info_d_ids$STRING_id <- NULL
            info_d_ids <- info_d_ids[1,]
            info_d_ids[,where_descr] <- IMPRINTS.CETSA.app:::getGeneName(info_d_ids[,where_descr])
            info_d_ids <- string_db$map(info_d_ids, colnames(info_d_ids)[where_descr], removeUnmappedRows = TRUE)
          }
          else{# if no gene info can't map and take first row
            showNotification("No gene informations has been found in your data; only first identifier is taken",
                             type = "error", duration = 5)
            info_d_ids <- info_d_ids[1,]
          }

          if(nrow(info_d_ids) != 0){
            string_res$x <- rbind.data.frame(string_res$x, info_d_ids, make.row.names = FALSE)
          }
        }
      }
    }

    showNotification("Mapping succeed !", type = "message")
  })
  output$data_stri_up <- reactive({
    return(!is.null(string_res$x))
  })
  outputOptions(output, "data_stri_up", suspendWhenHidden = FALSE)

  OUT_plot <- reactiveValues(
    g = NULL,
    g_int = NULL
  )

  g_stri <- reactive({
    if(input$intnet_stri){
      get_net_app(string_res$x$STRING_id , inter = TRUE,
             network_flavor = input$edgetype1_stri, required_score = input$intscore1_stri)
    }
    else{
      get_net_app(string_res$x$STRING_id , inter = FALSE,
             network_flavor = input$edgetype1_stri, required_score = input$intscore1_stri)
    }
  })

  observeEvent(input$netbase_stri, {
    if(length(string_res$x$STRING_id) > 2000){
      showNotification(paste("Lists with more than 2000 genes are not supported yet. Your list contains now", length(string_res$x$STRING_id), "genes.",
                             "Please, try to reduce the size of your input by choosing less categories and/or treatments."),
                       type = "error", duration = 8)
    }
    else{
      showNotification("Getting network, this may take a while", type = "message")
      if(input$intnet_stri){
        OUT_plot$g_int <- g_stri()
        OUT_plot$g <- NULL
      }
      else{
        OUT_plot$g <- g_stri()
        OUT_plot$g_int <- NULL
      }
      showNotification("Done !", type = "message")
    }

  })

  output$netInt_stri <- renderPlotly({
    OUT_plot$g_int
  })
  output$net_stri <- renderPlot({
    OUT_plot$g
  })

  output$downnet_stri <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%y%m%d_%H%M_"), "Network", ".", input$downnet_stri_format)
    },
    content = function(file){
      ggsave(file, plot = OUT_plot$g, device = input$downnet_stri_format)
    }
  )


  enrich_res <- reactiveValues(
    x = NULL
  )
  observeEvent(input$go_enrich, {
    showNotification("Starting enrichment analysis, this may take a while", type = "message")
    enrich_res$x <- string_db$get_enrichment(string_res$x$STRING_id)

    updateSelectInput(session, "catego_stri", choices = unique(enrich_res$x$category))
    showNotification("Done !", type = "message")
  })
  output$enrich_stri_up <- reactive({
    return(!is.null(enrich_res$x))
  })
  outputOptions(output, "enrich_stri_up", suspendWhenHidden = FALSE)

  enrich_res_tab <- reactive({
    df <- NULL

    if(!is.null(enrich_res$x)){
      if(str_length(input$catego_stri) != 0){
        df <- get_GO_app(enrich_res$x, TRUE, FALSE, input$catego_stri)
        d_n <- as.list(lapply(df, names)[[input$catego_stri]])
        df <- mapply(function(x,y) {x$id <- rep(y, nrow(x)); x},
                     df[[input$catego_stri]],
                     d_n, SIMPLIFY = FALSE)
        df <- do.call(rbind, df)
        rownames(df) <- 1:nrow(df)

        id_string <- do.call(rbind, str_split(df$id, ","))
        colnames(id_string) <- c("gene.names", "STRING_id")

        df <- cbind(id_string, df[,-ncol(df)])
      }
    }

    df
  })
  output$enrich_res_tab_up <- reactive({
    return(!is.null(enrich_res_tab()))
  })
  outputOptions(output, "enrich_res_tab_up", suspendWhenHidden = FALSE)

  output$enrich_table_stri <- DT::renderDataTable({

    DT::datatable(enrich_res_tab(),
                  caption = htmltools::tags$caption(
                    style = 'caption-side: top; text-align: left;',
                    htmltools::strong(input$catego_stri)
                  ),
                  rownames = FALSE,
                  options = list(lengthMenu = c(10,20,30), pageLength = 10,
                                 scrollX = TRUE))

  })
  output$downenrich_stri <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%y%m%d_%H%M_"), "Enrichment_tab_", input$catego_stri, ".xlsx")
    },
    content = function(file){
      openxlsx::write.xlsx(enrich_res_tab(), file, row.names = FALSE)
    }
  )


  observe({
    if(!is.null(enrich_res_tab())){
      updateSelectInput(session, "descri_stri", choices = unique(enrich_res_tab()$description))
    }
  })

  OUT_plot_filt <- reactiveValues(
    g = NULL,
    g_int = NULL
  )

  g_filt_stri <- reactive({
    descr <- input$descri_stri
    descr <- str_replace_all(descr, "\\(", "\\\\(")
    descr <- str_replace_all(descr, "\\)", "\\\\)")
    pr <- enrich_res_tab()$STRING_id[str_which(enrich_res_tab()$description, paste0("^", descr, "$"))]

    if(!is.null(pr) & length(pr)){
      if(input$intnet_stri){
        get_net_app(pr , inter = TRUE,
               network_flavor = input$edgetype2_stri, required_score = input$intscore2_stri)
      }
      else{
        get_net_app(pr , inter = FALSE,
               network_flavor = input$edgetype2_stri, required_score = input$intscore2_stri)
      }
    }
    else{
      showNotification("No matches has been found, try another description", type = "error")
    }
  })

  observeEvent(input$netfilt_stri, {
    showNotification("Getting new network, this may take a while", type = "message")
    if(input$intnet_stri){
      OUT_plot_filt$g_int <- g_filt_stri()
      OUT_plot_filt$g <- NULL
    }
    else{
      OUT_plot_filt$g <- g_filt_stri()
      OUT_plot_filt$g_int <- NULL
    }
    showNotification("Done !", type = "message")
  })

  output$netInt2_stri <- renderPlotly({
    OUT_plot_filt$g_int
  })
  output$net2_stri <- renderPlot({
    OUT_plot_filt$g
  })

  output$downnetfilt_stri <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%y%m%d_%H%M_"), "Network_", input$descri_stri, ".", input$downnetfilt_stri_format)
    },
    content = function(file){
      ggsave(file, plot = OUT_plot_filt$g, device = input$downnetfilt_stri_format)
    }
  )



  ### Barplots Network
  output$drug2_ui_barnet <- renderUI({
    selectInput("drug2_barnet", "Choose a drug", choices = names(drug_data_sh$y$data), multiple = TRUE, selected = "elutriation")
  })

  barnet_data <- reactive({
    File <- input$caldiff_barnet
    if (is.null(File))
      return(NULL)

    ms_fileread(File$datapath)
  })
  barnethit_data <- reactive({
    File <- input$hits_barnet
    if (is.null(File))
      return(NULL)

    d <- import(File$datapath, header = TRUE)
    d
  })


  Sel_cond_fhit_SUMMA_barnet <- reactiveValues(
    hit = NULL
  )
  Sel_cond_fhit_barnet <- reactive({
    HIT <- NULL
    if(input$onlyhit_barnet){
      if(input$drug_barnet == "base" & length(input$drug2_barnet) >= 1){
        HIT <- do.call(rbind, lapply(drug_data_sh$y$hitlist[input$drug2_barnet],
                                     function(x) x[,c("id", "treatment", "category")])
        )
      }
      else if(input$drug_barnet == "dat"){
        HIT <- barnethit_data()
      }

      Sel_cond_fhit_SUMMA_barnet$hit <- HIT
      c_idx <- str_which(colnames(HIT), "treatment")

      if(length(c_idx)){
        HIT <- HIT[,c_idx]
        HIT <- unique(HIT)
      }
      else{
        HIT <- NULL
        showNotification("Your hitlist doesn't contains the column 'treatment'", type = "error")
      }
    }
    HIT
  })
  observe({
    updateSelectInput(session, "cond_fhit_barnet", choices = Sel_cond_fhit_barnet(), selected = Sel_cond_fhit_barnet()[1])
  })

  sel_prot_barnet <- reactive({
    pr <- NULL
    if(input$drug_barnet == "base"){
      if(input$importprot_barnet){
        File <- input$prlist_file_barnet
        if (is.null(File)){
          pr <- NULL
        }
        else{
          pr <- unique(read.delim(File$datapath, header = FALSE)[[1]])

          prcheck <- ""
          if(length(input$drug2_barnet) == 1){
            prcheck <- drug_data_sh$y$data[[input$drug2_barnet]][,c("id", "description")]
          }
          else if(length(input$drug2_barnet) > 1){
            prcheck <- plyr::join_all(drug_data_sh$y$data[input$drug2_barnet], by = c("id", "description"), type = "full")[,c("id", "description")]
          }

          a <- pr[!(pr %in% prcheck$id)]
          if(length(a)){
            pr <- pr[(pr %in% prcheck$id)]
            showNotification(paste(paste(a, collapse = ", "), "wasn't in the data and had to be removed."),
                             type = "error")
          }

          pr <- data.frame(id = pr)
          pr <- unique(pr)
          pr <- dplyr::left_join(pr, prcheck,
                                 by = "id")
          pr$description <- unname(sapply(pr$description, IMPRINTS.CETSA.app:::getGeneName))
          pr <- pr %>%
            group_by(description) %>%
            mutate(id = paste0(id, ":", description))
          pr <- pr$id
          pr <- unique(pr)
        }
      }
      else{
        if(length(input$drug2_barnet) == 1 & all(input$drug2_barnet %in% names(drug_data_sh$y$data))){
          df <- drug_data_sh$y$data[[input$drug2_barnet]][,c("id", "description")]
          if(input$onlyhit_barnet & !is.null(input$cond_fhit_barnet)){
            pr <- Sel_cond_fhit_SUMMA_barnet$hit
            pr <- pr %>% dplyr::filter(!is.na(match(treatment, c(input$cond_fhit_barnet))))
            pr <- pr[,"id", drop = FALSE]
            pr <- unique(pr)
            pr <- dplyr::left_join(pr, df,
                                   by = "id")
            pr$description <- unname(sapply(pr$description, IMPRINTS.CETSA.app:::getGeneName))
            pr <- pr %>%
              group_by(description) %>%
              mutate(id = paste0(id, ":", description))
            pr <- pr$id
            pr <- unique(pr)
          }
          else{
            pr <- df
            pr$description <- unname(sapply(pr$description, IMPRINTS.CETSA.app:::getGeneName))
            pr <- pr %>%
              group_by(description) %>%
              mutate(id = paste0(id, ":", description))
            pr <- pr$id
            pr <- unique(pr)
          }
        }
        else if(length(input$drug2_barnet) > 1 & all(input$drug2_barnet %in% names(drug_data_sh$y$data))){
          df <- plyr::join_all(drug_data_sh$y$data[input$drug2_barnet], by = c("id", "description"), type = "full")[,c("id", "description")]
          if(input$onlyhit_barnet & !is.null(input$cond_fhit_barnet)){
            pr <- Sel_cond_fhit_SUMMA_barnet$hit
            pr <- pr %>% dplyr::filter(!is.na(match(treatment, c(input$cond_fhit_barnet))))
            pr <- pr[,"id", drop = FALSE]
            pr <- unique(pr)
            pr <- dplyr::left_join(pr, df,
                                   by = "id")
            pr$description <- unname(sapply(pr$description, IMPRINTS.CETSA.app:::getGeneName))
            pr <- pr %>%
              group_by(description) %>%
              mutate(id = paste0(id, ":", description))
            pr <- pr$id
            pr <- unique(pr)
          }
          else{
            pr <- df
            pr$description <- unname(sapply(pr$description, IMPRINTS.CETSA.app:::getGeneName))
            pr <- pr %>%
              group_by(description) %>%
              mutate(id = paste0(id, ":", description))
            pr <- pr$id
            pr <- unique(pr)
          }
        }
      }
    }
    else if(input$drug_barnet == "dat"){
      if(input$importprot_barnet){
        File <- input$prlist_file_barnet
        if (is.null(File)){
          pr <- NULL
        }
        else{
          pr <- unique(read.delim(File$datapath, header = FALSE)[[1]])

          prcheck <- ""
          prcheck <- barnet_data()[,c("id", "description")]

          a <- pr[!(pr %in% prcheck$id)]
          if(length(a)){
            pr <- pr[(pr %in% prcheck$id)]
            showNotification(paste(paste(a, collapse = ", "), "wasn't in the data and had to be removed."),
                             type = "error")
          }

          pr <- data.frame(id = pr)
          pr <- unique(pr)
          pr <- dplyr::left_join(pr, prcheck,
                                 by = "id")
          pr$description <- unname(sapply(pr$description, IMPRINTS.CETSA.app:::getGeneName))
          pr <- pr %>%
            group_by(description) %>%
            mutate(id = paste0(id, ":", description))
          pr <- pr$id
          pr <- unique(pr)
        }
      }
      else{
        if(!is.null(barnet_data()) & !is.null(barnethit_data())){
          df <- barnet_data()[,c("id", "description")]
          if(input$onlyhit_barnet  & !is.null(input$cond_fhit_barnet)){
            pr <- Sel_cond_fhit_SUMMA_barnet$hit
            pr <- pr %>% dplyr::filter(!is.na(match(treatment, c(input$cond_fhit_barnet))))
            pr <- pr[,"id", drop = FALSE]
            pr <- unique(pr)
            pr <- dplyr::left_join(pr, df,
                                   by = "id")
            pr$description <- unname(sapply(pr$description, IMPRINTS.CETSA.app:::getGeneName))
            pr <- pr %>%
              group_by(description) %>%
              mutate(id = paste0(id, ":", description))
            pr <- pr$id
            pr <- unique(pr)
          }
          else{
            pr <- df
            pr$description <- unname(sapply(pr$description, IMPRINTS.CETSA.app:::getGeneName))
            pr <- pr %>%
              group_by(description) %>%
              mutate(id = paste0(id, ":", description))
            pr <- pr$id
            pr <- unique(pr)
          }
        }
      }
    }
    pr
  })
  observe({
    updateSelectizeInput(session, "prot_barnet", choices = sel_prot_barnet(), server = TRUE)
  })


  Sel_cond_barnet <- reactive({
    tr <- NULL
    if(input$drug_barnet == "base" & length(input$drug2_barnet) >= 1){
      a <- join_drugdata(drug_data_sh$y$data[input$drug2_barnet], by = c("id", "description"))
      tr <- get_treat_level(a)
    }
    else if(input$drug_barnet == "dat"){
      if(!is.null(barnet_data())){
        tr <- get_treat_level(barnet_data())
      }
      else{
        tr <- NULL
      }
    }
    tr
  })
  observe({
    updateSelectInput(session, "condition_barnet", choices = Sel_cond_barnet())
  })


  GOterm_data <- reactive({
    g <- NULL
    if(input$importGO_barnet){
      File <- input$GOtermfile_barnet
      if (is.null(File)){
        g <- NULL
      }
      else{
        g <- rio::import(File$datapath, header = TRUE)
        if(!all(c("id", "GOterm") %in% colnames(g))){
          showNotification("Your GOterm file doesn't contain the requested columns !", type = "error")
          g <- NULL
        }
      }
    }
    g
  })
  observe(GOterm_data())


  # handling color selection
  output$n_cond_sel_barnet <- renderText({
    if(input$ch_own_col_barnet){
      paste("You selected", length(input$condition_barnet), "treatments, please enter the same number of colors")
    }
    else{
      NULL
    }
  })

  OWN_color_barnet <- reactiveValues(
    ch = c()
  )
  observeEvent(input$add_col_barnet, {
    OWN_color_barnet$ch <- append(OWN_color_barnet$ch, input$own_color_pick_barnet)
  })
  observeEvent(input$rem_col_barnet, {
    if(length(OWN_color_barnet$ch) <= 1){
      OWN_color_barnet$ch <- c()
    }
    else{
      OWN_color_barnet$ch <- OWN_color_barnet$ch[1:(length(OWN_color_barnet$ch)-1)]
    }
  })
  output$own_color_barnet <- renderText({
    paste("You selected this colors :", paste(OWN_color_barnet$ch, collapse = ", "))
  })


  # the network
  thenet <- reactiveValues(
    n = NULL
  )
  observeEvent(input$plotnet_barnet, {
    if (input$drug_barnet == "base"){
      data <- join_drugdata(drug_data_sh$y$data[input$drug2_barnet], by = c("id", "description"))
    }
    else if(input$drug_barnet == "dat"){
      data <- barnet_data()
    }
    if(input$onlyhit_barnet  & !is.null(input$cond_fhit_barnet)){
      pr <- unname(sapply(sel_prot_barnet(), function(x) strsplit(x, ":")[[1]][1]))
      data <- data[which(!is.na(match(data$id, pr))),]
    }

    if(is.null(sel_prot_barnet())){
      if(input$onlyhit_barnet  & is.null(input$cond_fhit_barnet)){
        showNotification("Select some treatments to filter your hits !", type = "error")
      }
      else{
        showNotification("No data detected !", type = "error")
      }
    }
    else{ # means that there is data anyway
      showNotification("Getting network", type = "message")

      GOterm <- NULL
      if(input$importGO_barnet){
        GOterm <- GOterm_data()
      }
      else{
        if(input$GOtype_barnet != "none"){
          GOterm <- input$GOtype_barnet
        }
      }

      colorbar <- NULL
      if(input$ch_own_col_barnet){
        if(length(OWN_color_barnet$ch)){
          colorbar <- OWN_color_barnet$ch
        }
      }

      withCallingHandlers({
        shinyjs::html("diag_barnet", "")
        if(!is.null(input$prot_barnet)){
          pr <- unname(sapply(input$prot_barnet, function(x) strsplit(x, ":")[[1]][1]))
        }
        else{
          pr <- NULL
        }
        thenet$n <- imprints_network(data, hits = pr, GOterm = GOterm,
                                     treatment = input$condition_barnet,
                                     colorbar = colorbar,
                                     required_score = input$reqscore_barnet,
                                     species = input$species_barnet, witherrorbar = input$werb_barnet,
                                     FC_border = input$FCborder_barnet,
                                     colorFC = c(input$FCbordercolorlow_barnet,
                                                 input$FCbordercolormid_barnet,
                                                 input$FCbordercolorhigh_barnet))
      },
      message = function(m) {
        shinyjs::html(id = "diag_barnet", html = paste(m$message, "<br>", sep = ""), add = FALSE)
      }
      )
    }
  })
  #check if a file is upload
  output$netready_barnet <- reactive({
    return(!is.null(thenet$n))
  })
  outputOptions(output, "netready_barnet", suspendWhenHidden = FALSE)

  output$network_barnet <- renderVisNetwork({
    thenet$n
  })


  observeEvent(input$plotnet_barnet, {
    grp <- thenet$n$x$nodes$group
    updateSelectizeInput(session, "select_groups_barnet",
                         choices = unique(grp), server = TRUE)
  })

  observe({
    visNetworkProxy("network_barnet") %>%
      visEdges(color = input$edge_color_barnet)
  })
  observe({
    visNetworkProxy("network_barnet") %>%
      visEdges(length = input$edge_length_barnet)
  })

  observeEvent(input$node_color_barnet, {
    if(!is.null(thenet$n$x$nodes$group)){
      thenet$n$x$nodes$color.background[which(thenet$n$x$nodes$group ==
                                                input$select_groups_barnet)] <- input$node_color_barnet
      thenet$n$x$legend$nodes$color.background[which(thenet$n$x$legend$nodes$label ==
                                                       input$select_groups_barnet)] <- input$node_color_barnet
      thenet$n$x$legend$nodes$color.border[which(thenet$n$x$legend$nodes$label ==
                                                   input$select_groups_barnet)] <- input$node_color_barnet
    }
    else{
      thenet$n$x$nodes$color.background  <- input$node_color_barnet
    }
  }, ignoreInit = TRUE)

  observeEvent(input$node_bordercolor_barnet, {
    thenet$n$x$nodes$color.border <- input$node_bordercolor_barnet
  }, ignoreInit = TRUE)

  observe({
    visNetworkProxy("network_barnet") %>%
      visNodes(font = list(color = input$font_color_barnet))
  })
  observe({
    visNetworkProxy("network_barnet") %>%
      visNodes(font = list(background = input$font_backcolor_barnet))
  })
  observe({
    visNetworkProxy("network_barnet") %>%
      visNodes(font = list(size = input$label_size_barnet))
  })
  observe({
    visNetworkProxy("network_barnet") %>%
      visNodes(size = input$node_size_barnet)
  })
  observe({
    visNetworkProxy("network_barnet") %>%
      visNodes(imagePadding = input$node_imgpadding_barnet)
  })
  observe({
    visNetworkProxy("network_barnet") %>%
      visNodes(borderWidth = input$node_borderwidth_barnet)
  })

  observe({
    visNetworkProxy("network_barnet") %>%
      visPhysics(solver = input$physics_type_barnet)
  })
  observe({
    visNetworkProxy("network_barnet") %>%
      visPhysics(enabled = input$physics_enable_barnet)
  })
  observe({
    visNetworkProxy("network_barnet") %>%
      visInteraction(navigationButtons = input$button_enable_barnet)
  })

  output$down_barnet <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%y%m%d_%H%M_"), "network.html")
    },
    content = function(file){
      visSave(thenet$n, file = file, background = "#FFFFFF00")
    }
  )



  ### Cluster Profiler
  observe({
    if(input$database_clus == "CETSA"){
      updateSelectInput(session, "species_clus", choices = c("Human"))
    }
    else{
      updateSelectInput(session, "species_clus", choices = c("Human", "Mouse"))
    }
  })

  output$drug2ui_clus <- renderUI({
    selectInput("drug2_clus", "Choose a drug", choices = names(drug_data_sh$y$data),
                multiple = TRUE, selected = "elutriation")
  })

  clus_data <- reactive({
    if(input$drug_clus == "dat"){
      if(input$ishit_clus){
        File <- input$file_clus
        if (is.null(File))
          return(NULL)

        h <- import_list(File$datapath, header = TRUE)[[1]]
        h <- as.data.frame(h)
        if(!("Gene" %in% colnames(h))){
          if("description" %in% colnames(h)){
            h$Gene <- unname(unlist(sapply(h$description, IMPRINTS.CETSA.app:::getGeneName)))
          }
          else{
            h <- NULL
            showNotification("Your hitlist doesn't contain neither a Gene column, neither a description column.
                             Are you sure you imported a hitlist ?", type = "error")
          }
        }
      }
      else{
        File <- input$file_clus
        if (is.null(File))
          return(NULL)

        h <- rio::import(File$datapath, header = TRUE)
      }
      h
    }
    else if(input$drug_clus == "base" & length(input$drug2_clus) >= 1){
      h <- drug_data_sh$y$hitlist[input$drug2_clus]
      h_names <- lapply(h, colnames)
      h_names_common <- com_protein_loop(h_names)

      if(length(h_names_common) > 1){
        for(i in names(h_names_common)[-length(h_names_common)]){
          n <- strsplit(i, " & ")[[1]]
          for(k in h_names_common[[i]]){
            h[!(names(h) %in% n)] <- lapply(h[!(names(h) %in% n)], function(b) {b[[k]] <- NA; b})
          }
        }
        h <- as.data.frame(Reduce(rbind, h))
      }
      else{
        h <- Reduce(rbind, h)
        h <- as.data.frame(h)
      }

      if(!("Gene" %in% colnames(h))){
        if("description" %in% colnames(h)){
          h$Gene <- unname(unlist(sapply(h$description, IMPRINTS.CETSA.app:::getGeneName)))
        }
        else{
          d <- join_drugdata(drug_data_sh$y$data_ave[input$drug2_clus], by = c("id", "description")) ## extract gene information
          d <- d[,c("id", "description")]
          h <- dplyr::left_join(h, d, by = "id")
          h$Gene <- unname(unlist(sapply(h$description, IMPRINTS.CETSA.app:::getGeneName)))
        }
      }
      h$description <- NULL

      nn <- drug_data_sh$y$NN[input$drug2_clus]
      nn_names <- lapply(nn, colnames)
      nn_names_common <- com_protein_loop(nn_names)

      if(length(nn_names_common) > 1){
        for(i in names(nn_names_common)[-length(nn_names_common)]){
          n <- strsplit(i, " & ")[[1]]
          for(k in nn_names_common[[i]]){
            nn[!(names(nn) %in% n)] <- lapply(nn[!(names(nn) %in% n)], function(b) {b[[k]] <- NA; b})
          }
        }
        nn <- as.data.frame(Reduce(rbind, nn))
      }
      else{
        nn <- Reduce(rbind, nn)
        nn <- as.data.frame(nn)
      }

      if(!("Gene" %in% colnames(nn))){
        if("description" %in% colnames(nn)){
          nn$Gene <- unname(unlist(sapply(nn$description, IMPRINTS.CETSA.app:::getGeneName)))
        }
        else{
          d <- join_drugdata(drug_data_sh$y$data_ave[input$drug2_clus], by = c("id", "description")) ## extract gene information
          d <- d[,c("id", "description")]
          nn <- dplyr::left_join(nn, d, by = "id")
          nn$Gene <- unname(unlist(sapply(n$description, IMPRINTS.CETSA.app:::getGeneName)))
        }
      }
      nn$description <- NULL

      rbind(h,nn)
    }
    else{
      NULL
    }
  })
  #check if a file is upload
  output$file_clus_up <- reactive({
    return(!is.null(clus_data()))
  })
  outputOptions(output, "file_clus_up", suspendWhenHidden = FALSE)

  Sel_cond_fhit_clus <- reactive({
    tr <- NULL

    if(input$ishit_clus | input$drug_clus == "base"){
      HIT <- clus_data()

      c_idx <- str_which(colnames(HIT), "treatment")

      if(length(c_idx)){
        tr <- HIT[, c_idx]
        tr <- unique(tr)
      }
    }

    tr
  })

  observe({
    if(input$drug_clus =="dat"){
      updateSelectInput(session, "cond_fhit_clus", choices = Sel_cond_fhit_clus(), selected = Sel_cond_fhit_clus()[1])
    }
    else if(input$drug_clus =="base"){
      updateSelectInput(session, "cond_fhitB_clus", choices = Sel_cond_fhit_clus(), selected = Sel_cond_fhit_clus()[1])
    }
  })

  ## cluster comparison
  cluscomp_res <- reactiveValues(
    x = NULL
  )
  observeEvent(input$gocomp_clus, {
    dat <- NULL
    res <- NULL
    if (!is.null(clus_data())){
      if(input$drug_clus == "dat"){
        dat <- clus_data()
        if(input$ishit_clus){
          dat <- clus_data()
          if(!is.null(input$cond_fhit_clus)){
            dat <- dat %>% dplyr::filter(!is.na(match(treatment, c(input$cond_fhit_clus))))
          }
          else{
            showNotification("Don't forget to select some treatments !", type = "error")
            dat <- NULL
          }
          if(!is.null(dat)){
            showNotification("Starting enrichment analysis !", type = "message")
            res <- compare_enrich(dat, gene_column = "Gene", species = input$species_clus,
                                  n_pathway = input$npath_clus, treatment_column = "treatment",
                                  pval_cutoff = input$pvcut_clus, minGSSize = input$minGNsize_clus,
                                  database = input$database_clus)
          }
        }
        else{
          showNotification("Starting pathway analysis !", type = "message")
          res <- compare_enrich(dat, gene_column = input$idfile_clus,
                                species = input$species_clus, n_pathway = input$npath_clus,
                                pval_cutoff = input$pvcut_clus, minGSSize = input$minGNsize_clus,
                                database = input$database_clus)
        }
      }
      else if(input$drug_clus == "base"){
        dat <- clus_data()
        if(!is.null(input$cond_fhitB_clus)){
          dat <- dat %>% dplyr::filter(!is.na(match(treatment, c(input$cond_fhitB_clus))))
          if(!is.null(input$cat_fhitB_clus)){
            dat <- dat %>% dplyr::filter(!is.na(match(category, c(input$cat_fhitB_clus))))
          }
        }
        else{
          showNotification("Don't forget to select some treatments !", type = "error")
          dat <- NULL
        }
        if(!is.null(dat)){
          showNotification("Starting pathway analysis !", type = "message")
          res <- compare_enrich(dat, gene_column = "Gene", species = input$species_clus,
                                n_pathway = input$npath_clus, treatment_column = "treatment",
                                pval_cutoff = input$pvcut_clus, minGSSize = input$minGNsize_clus,
                                database = input$database_clus)
        }
      }
    }

    cluscomp_res$x <- res
  })

  output$comptab_clus <- DT::renderDataTable({
    DT::datatable(cluscomp_res$x$res,
                  caption = htmltools::tags$caption(
                    style = 'caption-side: top; text-align: left;',
                    htmltools::strong("Compare cluster results")
                  ),
                  rownames = FALSE,
                  options = list(lengthMenu = c(10,20,30), pageLength = 10,
                                 scrollX = TRUE))
  })
  output$downcomptab_clus <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%y%m%d_%H%M_"), "CompareClusterTab_", input$database_clus,  ".xlsx")
    },
    content = function(file){
      openxlsx::write.xlsx(cluscomp_res$x$res, file)
    }
  )

  output$compplot_clus <- renderPlot({
    cluscomp_res$x$graph
  })
  output$downcomplot_clus <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%y%m%d_%H%M_"), "CompareCluster_", input$database_clus, ".", input$downcomplot_clus_format)
    },
    content = function(file){
      ggsave(file, plot = cluscomp_res$x$graph, device = input$downcomplot_clus_format,
             width = 16, height = 12, dpi = 600)
    }
  )

  ## GSEA
  output$scorenameui_clus <- renderUI({
    if(input$drug_clus == "dat"){
      textInput("scorenameDAT_clus", "Type the column's name that contain the score")
    }
    else if(input$drug_clus == "base"){
      ch <- "No score"
      if(!is.null(clus_data())){
        n <- colnames(clus_data())[colnames(clus_data()) %in% c("Fisher", "IS", "GlobalScore")]
        if(length(n)){
          ch <- n
        }
      }
      selectInput("scorenameBASE_clus", "Select a score", choices = ch)
    }
  })
  clusgsea_res <- reactiveValues(
    x = NULL
  )
  observeEvent(input$gogsea_clus, {
    dat <- NULL
    res <- NULL
    if (!is.null(clus_data())){
      if(input$drug_clus == "dat"){
        dat <- clus_data()
        if(input$ishit_clus){
          dat <- clus_data()
          if(!is.null(input$cond_fhit_clus)){
            dat <- dat %>% dplyr::filter(!is.na(match(treatment, c(input$cond_fhit_clus))))
          }
          else{
            showNotification("Don't forget to select some treatments !", type = "error")
            dat <- NULL
          }
          if(!is.null(dat)){
            if(input$scorenameDAT_clus %in% colnames(dat)){
              showNotification("Starting GSEA !", type = "message")
              dat[[input$scorenameDAT_clus]] <- as.numeric(dat[[input$scorenameDAT_clus]])
              dat[[input$scorenameDAT_clus]] <- tidyr::replace_na(dat[[input$scorenameDAT_clus]],0)
              res <- run_gsea(dat, gene_column = "Gene",
                                species = input$species_clus,
                                score_column = input$scorenameDAT_clus,
                                pos_enrichment = input$onlypos_clus,
                                pval_cutoff = input$pvcut_clus, minGSSize = input$minGNsize_clus,
                                database = input$database_clus)
            }
            else{
              showNotification(paste(input$scorenameDAT_clus,
                                     "was not found in the column names of the data. Please check the name you wrote."),
                               type = "error")
            }
          }
        }
        else{
          if(input$scorenameDAT_clus %in% colnames(dat)){
            showNotification("Starting GSEA !", type = "message")
            dat[[input$scorenameDAT_clus]] <- as.numeric(dat[[input$scorenameDAT_clus]])
            dat[[input$scorenameDAT_clus]] <- tidyr::replace_na(dat[[input$scorenameDAT_clus]],0)
            res <- run_gsea(dat, gene_column = input$idfile_clus,
                              species = input$species_clus,
                              score_column = input$scorenameDAT_clus,
                              pos_enrichment = input$onlypos_clus,
                              pval_cutoff = input$pvcut_clus, minGSSize = input$minGNsize_clus,
                              database = input$database_clus)
          }
          else{
            showNotification(paste(input$scorenameDAT_clus,
                                   "was not found in the column names of the data. Please check the name you wrote."),
                             type = "error")
          }
        }
      }
      else if(input$drug_clus == "base"){
        dat <- clus_data()
        if(!is.null(input$cond_fhitB_clus)){
          dat <- dat %>% dplyr::filter(!is.na(match(treatment, c(input$cond_fhitB_clus))))
          if(!is.null(input$cat_fhitB_clus)){
            dat <- dat %>% dplyr::filter(!is.na(match(category, c(input$cat_fhitB_clus))))
          }
        }
        else{
          showNotification("Don't forget to select some treatments !", type = "error")
          dat <- NULL
        }
        if(!is.null(dat)){
          if(input$scorenameBASE_clus %in% colnames(dat)){
            showNotification("Starting GSEA !", type = "message")
            dat[[input$scorenameBASE_clus]] <- as.numeric(dat[[input$scorenameBASE_clus]])
            dat[[input$scorenameBASE_clus]] <- tidyr::replace_na(dat[[input$scorenameBASE_clus]],0)
            res <- run_gsea(dat, gene_column = "Gene",
                              species = input$species_clus,
                              score_column = input$scorenameBASE_clus,
                              pos_enrichment = input$onlypos_clus,
                              pval_cutoff = input$pvcut_clus, minGSSize = input$minGNsize_clus,
                              database = input$database_clus)
          }
          else{
            showNotification(paste(input$scorenameBASE_clus,
                                   "was not found in the column names of the data. Please check the name you wrote."),
                             type = "error")
          }
        }
      }
    }

    clusgsea_res$x <- res
  })

  output$gseatab_clus <- DT::renderDataTable({
    DT::datatable(clusgsea_res$x$res,
                  caption = htmltools::tags$caption(
                    style = 'caption-side: top; text-align: left;',
                    htmltools::strong("GSEA results")
                  ),
                  rownames = FALSE,
                  options = list(lengthMenu = c(10,20,30), pageLength = 10,
                                 scrollX = TRUE))
  })
  output$downgseatab_clus <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%y%m%d_%H%M_"), "GSEATab_", input$database_clus, ".xlsx")
    },
    content = function(file){
      openxlsx::write.xlsx(clusgsea_res$x$res, file)
    }
  )

  output$gseaplot_clus <- renderPlot({
    clusgsea_res$x$graph
  })
  output$downgsealot_clus <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%y%m%d_%H%M_"), "GSEAplot_", input$database_clus, ".", input$downgsealot_clus_format)
    },
    content = function(file){
      if(input$downgsealot_clus_format == "png"){
        png(file, width = 1720, height = 1080)
        print(clusgsea_res$x$graph)
        dev.off()
      }
      else if(input$downgsealot_clus_format == "pdf"){
        pdf(file, onefile = TRUE, width = 12, height = 7)
        print(clusgsea_res$x$graph)
        dev.off()
      }
    }
  )

  ## Gene concept network
  output$scorename2ui_clus <- renderUI({
    if(input$drug_clus == "dat"){
      textInput("scorename2DAT_clus", "Type the column's name that contain the score")
    }
    else if(input$drug_clus == "base"){
      ch <- "No score"
      if(!is.null(clus_data())){
        n <- colnames(clus_data())[colnames(clus_data()) %in% c("Fisher", "IS", "GlobalScore")]
        if(length(n)){
          ch <- n
        }
      }
      selectInput("scorename2BASE_clus", "Select a score", choices = ch)
    }
  })
  clusgene_res <- reactiveValues(
    x = NULL
  )
  observeEvent(input$gogeneconc_clus, {
    dat <- NULL
    res <- NULL
    if (!is.null(clus_data())){
      if(input$drug_clus == "dat"){
        dat <- clus_data()
        if(input$ishit_clus){
          dat <- clus_data()
          if(!is.null(input$cond_fhit_clus)){
            dat <- dat %>% dplyr::filter(!is.na(match(treatment, c(input$cond_fhit_clus))))
          }
          else{
            showNotification("Don't forget to select some treatments !", type = "error")
            dat <- NULL
          }
          if(!is.null(dat)){
            if(input$scorename2DAT_clus %in% colnames(dat)){
              showNotification("Starting GSEA !", type = "message")
              dat[[input$scorename2DAT_clus]] <- as.numeric(dat[[input$scorename2DAT_clus]])
              dat[[input$scorename2DAT_clus]] <- tidyr::replace_na(dat[[input$scorename2DAT_clus]],0)
              res <- gene_concept_net(dat, gene_column = "Gene",
                                     species = input$species_clus,
                                     score_column = input$scorename2DAT_clus,
                                     pval_cutoff = input$pvcut_clus, minGSSize = input$minGNsize_clus,
                                     database = input$database_clus)
            }
            else{
              showNotification(paste(input$scorename2DAT_clus,
                                     "was not found in the column names of the data. Please check the name you wrote."),
                               type = "error")
            }
          }
        }
        else{
          if(input$scorename2DAT_clus %in% colnames(dat)){
            showNotification("Starting GSEA !", type = "message")
            dat[[input$scorename2DAT_clus]] <- as.numeric(dat[[input$scorename2DAT_clus]])
            dat[[input$scorename2DAT_clus]] <- tidyr::replace_na(dat[[input$scorename2DAT_clus]],0)
            res <- gene_concept_net(dat, gene_column = input$idfile_clus,
                                   species = input$species_clus,
                                   score_column = input$scorename2DAT_clus,
                                   pval_cutoff = input$pvcut_clus, minGSSize = input$minGNsize_clus,
                                   database = input$database_clus)
          }
          else{
            showNotification(paste(input$scorename2DAT_clus,
                                   "was not found in the column names of the data. Please check the name you wrote."),
                             type = "error")
          }
        }
      }
      else if(input$drug_clus == "base"){
        dat <- clus_data()
        if(!is.null(input$cond_fhitB_clus)){
          dat <- dat %>% dplyr::filter(!is.na(match(treatment, c(input$cond_fhitB_clus))))
          if(!is.null(input$cat_fhitB_clus)){
            dat <- dat %>% dplyr::filter(!is.na(match(category, c(input$cat_fhitB_clus))))
          }
        }
        else{
          showNotification("Don't forget to select some treatments !", type = "error")
          dat <- NULL
        }
        if(!is.null(dat)){
          if(input$scorename2BASE_clus %in% colnames(dat)){
            showNotification("Starting GSEA !", type = "message")
            dat[[input$scorename2BASE_clus]] <- as.numeric(dat[[input$scorename2BASE_clus]])
            dat[[input$scorename2BASE_clus]] <- tidyr::replace_na(dat[[input$scorename2BASE_clus]],0)
            res <- gene_concept_net(dat, gene_column = "Gene",
                                   species = input$species_clus,
                                   score_column = input$scorename2BASE_clus,
                                   pval_cutoff = input$pvcut_clus, minGSSize = input$minGNsize_clus,
                                   database = input$database_clus)
          }
          else{
            showNotification(paste(input$scorename2BASE_clus,
                                   "was not found in the column names of the data. Please check the name you wrote."),
                             type = "error")
          }
        }
      }
    }

    clusgene_res$x <- res
  })


  output$geneplot_clus <- renderPlot({
    clusgene_res$x
  })
  output$downgenelot_clus <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%y%m%d_%H%M_"), "GeneConceptNet_", input$database_clus, ".", input$downgenelot_clus_format)
    },
    content = function(file){
      ggsave(file, plot = clusgene_res$x, device = input$downgenelot_clus_format,
             width = 10, height = 6)
    }
  )




  ### CELL
  output$drug2ui_cell <- renderUI({
    selectInput("drug2_cell", "Choose a drug", choices = names(drug_data_sh$y$data),
                multiple = TRUE, selected = "elutriation")
  })

  hitdata_cell <- reactive({
    if(input$drug_cell == "dat"){
      File <- input$hitl_cell
      if (is.null(File))
        return(NULL)

      d <- import(File$datapath, header = TRUE)
      if(!all(c("id", "treatment", "category") %in% colnames(d))){
        missing_columns <- c("id", "treatment", "category")
        missing_columns <- missing_columns[!(c("id", "treatment", "category") %in% colnames(d))]
        missing_columns <- paste(missing_columns, collapse = ", ")
        verb <- ifelse(length(missing_columns) > 1, "are", "is")
        showNotification(paste(missing_columns, verb, "not in your summary file. Please check your columns names !"),
                         type = "error", duration = 8)
        d <- NULL
      }
      d
    }
    else if(input$drug_cell == "base" & length(input$drug2_cell) >= 1){
      h <- do.call(rbind, lapply(drug_data_sh$y$hitlist[input$drug2_cell],
                                 function(x) x[,c("id", "treatment", "category")])
                   )
      n <- do.call(rbind, lapply(drug_data_sh$y$NN[input$drug2_cell],
                                 function(x) x[,c("id", "description", "treatment", "category")])
                   )
      n <- unique(n[,c("id", "treatment", "category")])

      rbind(h,n)
    }
    else{
      NULL
    }
  })
  output$hitdata_cell_up <- reactive({
    return(!is.null(hitdata_cell()))
  })
  outputOptions(output, "hitdata_cell_up", suspendWhenHidden = FALSE)

  observe({
    if(!is.null(hitdata_cell())){
      updateSelectInput(session, "condhit_cell", choices = unique(hitdata_cell()$treatment), selected = unique(hitdata_cell()$treatment)[1])
      updateSelectInput(session, "cathit_cell", choices = unique(hitdata_cell()$category), selected = unique(hitdata_cell()$category)[1])
    }
  })


  resdata_cell <- reactiveValues(
    ch = NULL
  )
  observeEvent(input$goloca_cell, {
    showNotification("Getting subcellular locations", type = "message")

    data_hit <- hitdata_cell() %>%
      dplyr::filter(!is.na(match(treatment, c(input$condhit_cell))))

    if(!is.null(input$cathit_cell)){
      data_hit <- data_hit %>% dplyr::filter(!is.na(match(category, input$cathit_cell)))
    }

    withCallingHandlers({
      shinyjs::html("diagl_cell", "")
      resdata_cell$ch <- hit_for_cell(data_hit, input$organism_cell)
    },
    message = function(m) {
      shinyjs::html(id = "diagl_cell", html = paste(m$message, "<br>", sep = ""), add = TRUE)

    }
    )

    output$locatab_cell <- DT::renderDataTable({
      DT::datatable(resdata_cell$ch,
                    caption = htmltools::tags$caption(
                      style = 'caption-side: top; text-align: left;',
                      htmltools::strong("Subcellular locations from your hitlist")
                    ),
                    rownames = FALSE,
                    options = list(lengthMenu = c(10,20,30), pageLength = 10,
                                   scrollX = TRUE))
    })

    output$down_prl_cell <- downloadHandler(
      filename = function() {
        paste0(format(Sys.time(), "%y%m%d_%H%M_"), "HIT_locations", ".xlsx")
      },
      content = function(file){
        openxlsx::write.xlsx(resdata_cell$ch, file, row.names = FALSE)
      }
    )

    updateSelectInput(session, "condp_cell", choices = unique(resdata_cell$ch$treatment), selected = input$condhit_cell)
    updateSelectInput(session, "selorga_cell", choices = unique(resdata_cell$ch$main.location.cell))
  })
  output$resdata_cell_up <- reactive({
    return(!is.null(resdata_cell$ch))
  })
  outputOptions(output, "resdata_cell_up", suspendWhenHidden = FALSE)


  cell_p_Rv <- reactiveValues(
    g = NULL
  )

  cell_p_R <- reactive({
    hit_plotcell(resdata_cell$ch, tit = input$titp_cell,
                 cond = input$condp_cell,
                 cat_col_list = list("CC" = "#FB4F0B", "CN" = "#0FAEB9",
                                     "NC" = "#E7B700", "ND" = "#747474",
                                     "NN" = "#CCCCCC"))
  })
  observeEvent(input$gop_cell, {
    showNotification("Getting plot, this can take a while. Please wait", type = "message")

    cell_p_Rv$g <- cell_p_R()
  })

  output$cell_p <- renderPlotly({
    cell_p_Rv$g
  })
  output$downthe_cell <- downloadHandler(
    filename = function() {
      paste(format(Sys.time(), "%y%m%d_%H%M_"), "Int_cell", ".html", sep = "")
    },
    content = function(file){
      withr::with_dir(WD, htmlwidgets::saveWidget(partial_bundle(cell_p_Rv$g), file))
    }
  )

  #handle selection of proteins
  PR_event <- reactiveVal()
  PR_event_click <- reactive({
    if(!is.null(cell_p_Rv$g)){
      event_data("plotly_click", source = "M")
    }
    else{
      NULL
    }
  })
  PR_event_dbclick <- reactive({
    if(!is.null(cell_p_Rv$g)){
      event_data("plotly_doubleclick", source = "M")
    }
    else{
      NULL
    }
  })
  observeEvent(PR_event_click(), {
    PR <- event_data("plotly_click", source = "M")$customdata
    if(!is.null(PR_event)){
      if(length(which(PR_event() == PR)) == 0){
        PR_old_new <- c(PR_event(), PR)
      }
      else{
        PR_old_new <- PR_event()[-which(PR_event() == PR)] #if already clicked, remove
      }
    }
    else{
      PR_old_new <- c(PR_event(), PR)
    }

    PR_event(unique(PR_old_new))
  })
  # clear the set of cars when a double-click occurs
  observeEvent(PR_event_dbclick(), {
    PR_event(NULL)
  })

  output$prsel_p_cell <- renderUI({
    if(!is.null(PR_event())){
      pr <- PR_event()
      pr_link <- paste0("https://www.uniprot.org/uniprot/", pr, ">")
      pr_html <- paste0("<a href=", pr_link, pr, "</a>")
      HTML(paste("You clicked on", paste(pr_html, collapse = ", ")))
    }
    else{
      NULL
    }
  })

  barpdata_cell <- reactive({
    if(input$drug_cell == "dat"){
      File <- input$filebarp_cell
      if (is.null(File))
        return(NULL)

      ms_fileread(File$datapath)
    }
    else if(input$drug_cell == "base" & length(input$drug2_cell) >= 1){
      join_drugdata(drug_data_sh$y$data[input$drug2_cell], by = c("id", "description"))
    }
    else{
      NULL
    }
  })
  output$barpdata_cell_up <- reactive({
    return(!is.null(barpdata_cell()))
  })
  outputOptions(output, "barpdata_cell_up", suspendWhenHidden = FALSE)


  sel_prot_cell <- reactive({
    pr <- NULL
    if(input$selpr_loca_cell){
      if(!is.null(resdata_cell$ch)){
        pr <- resdata_cell$ch[which(!is.na(match(resdata_cell$ch$main.location.cell, input$selorga_cell))), c("id", "gene.name")]
        pr <- paste0(pr$id, ":", pr$gene.name)
        pr <- unique(pr)
      }
    }
  })
  observe({
    updateSelectizeInput(session, "selectpr_cell", choices = sel_prot_cell(), server = TRUE)

  })

  Sel_cond_cell <- reactive({
    if(input$selpr_loca_cell){
      if(input$allpr_cell){
        pr <- sel_prot_cell()
      }
      else{
        pr <- input$selectpr_cell
      }
      if(!is.null(pr)){
        pr <- unname(sapply(pr, function(x) strsplit(x, ":")[[1]][1]))
      }
    }
    else{
      pr <- PR_event()
    }

    tr <- NULL
    if(!is.null(barpdata_cell())){
      if(input$cond_sel_cell == "cat" & !is.null(hitdata_cell())){
        tr <- hitdata_cell()[which(!is.na(match(hitdata_cell()$id,pr))),c("treatment", "category")]
      }
      else {
        tr <- get_treat_level(barpdata_cell())
      }
    }
    tr
  })



  observe({
    if(input$cond_sel_cell == "cat"){
      updateSelectInput(session, "cond_cell", choices = unique(Sel_cond_cell()$category))
    }
    else{
      updateSelectInput(session, "cond_cell", choices = Sel_cond_cell())
    }
  })

  data_cell <- reactive({
    if(!is.null(barpdata_cell())){

      data <- barpdata_cell()
      TREAT <- get_treat_level(barpdata_cell())

      if(input$cond_sel_cell == "treat"){
        notsel_cond <- TREAT[!(TREAT %in% input$cond_cell)]
        notsel_cond <- paste0("_", notsel_cond, "$")
        notsel_cond <- paste(notsel_cond, collapse = "|")

        if(str_length(notsel_cond) != 0){
          data <- data[,-str_which(names(data), notsel_cond)]
        }

        id_sel <- str_which(names(data), paste(input$cond_cell, collapse = "|"))
        w <- 1:ncol(data)
        w <- w[!(w %in% id_sel)]

        ord <- unlist(lapply(input$cond_cell, function(x) str_which(names(data), paste0("_", x, "$"))))

        data <- data[,c(w,ord)]

      }
      else if(input$cond_sel_cell == "cat"){
        sele_cond <- Sel_cond_cell()$treatment[which(!is.na(match(Sel_cond_cell()$category, input$cond_cell)))]
        notsel_cond <- TREAT[!(TREAT %in% sele_cond)]
        if(length(notsel_cond)){
          notsel_cond <- paste(notsel_cond, collapse = "|")

          data <- data[,-str_which(names(data), notsel_cond)]
        }

      }
      else if(input$cond_sel_cell == "all_cond"){
        notsel_cond <- TREAT[!(TREAT %in% Sel_cond_cell())]
        if(length(notsel_cond)){
          notsel_cond <- paste(notsel_cond, collapse = "|")

          data <- data[,-str_which(names(data), notsel_cond)]
        }
      }

      if(input$selpr_loca_cell){
        loca_pr <- resdata_cell$ch[which(!is.na(match(resdata_cell$ch$main.location.cell,
                                                      input$selorga_cell))),
                                   c("id", "main.location.cell")
        ]
        if(input$allpr_cell){
          pr <- sel_prot_cell()
        }
        else{
          pr <- input$selectpr_cell
        }
        if(!is.null(pr)){
          pr <- unname(sapply(pr, function(x) strsplit(x, ":")[[1]][1]))
        }
      }
      else{
        pr <- PR_event()
      }

      if(input$selpr_loca_cell & input$save_bar_cell){
        data_l <- list()
        for(i in input$selorga_cell){
          loca_pr_ <- loca_pr[which(loca_pr$main.location.cell == i), ]

          pr_comp <- loca_pr_$id
          pr_comp <- pr_comp[which(!is.na(match(pr_comp, pr)))]

          data_l[[i]] <- data[which(!is.na(match(data$id, pr_comp))),]
        }
        data <- data_l
      }
      else{
        data <- data[which(!is.na(match(data$id, pr))),]
      }
    }
    else{
      data <- NULL
    }

    data
  })


  output$n_cond_sel_cell <- renderText({
    if(input$ch_own_col_cell){
      if (input$cond_sel_cell  == "all_cond"){
        paste("You selected", length(get_treat_level(data_cell())), "treatments, please enter the same number of colors")
      }
      else{
        paste("You selected", length(input$cond_cell), "treatments, please enter the same number of colors")
      }
    }
    else{
      NULL
    }
  })

  OWN_color_cell <- reactiveValues(
    ch = c()
  )
  observeEvent(input$add_col_cell, {
    OWN_color_cell$ch <- append(OWN_color_cell$ch, input$own_color_pick_cell)
  })
  observeEvent(input$rem_col_cell, {
    if(length(OWN_color_cell$ch) <= 1){
      OWN_color_cell$ch <- c()
    }
    else{
      OWN_color_cell$ch <- OWN_color_cell$ch[1:(length(OWN_color_cell$ch)-1)]
    }
  })
  output$own_color_cell <- renderText({
    paste("You selected this colors :", paste(OWN_color_cell$ch, collapse = ", "))
  })


  BAR_cell <- reactiveValues(
    ch = ev_null_print
  )

  Bar_one_cell <- reactive({
    if(input$save_bar_cell & input$selpr_loca_cell){
      loca_cell_lab <- "IMPRINTS-CETSA bar plotting \nMain cellular location :"
    }
    else{
      loca_cell_lab <- "IMPRINTS-CETSA bar plotting"
    }
    withCallingHandlers({
      shinyjs::html("diag_bar_cell", "")
      if(input$ch_own_col_cell){
        nbc <- ifelse(input$cond_sel_cell  == "all_cond", length(get_treat_level(data_cell())), length(input$cond_cell))
        COL <- OWN_color_cell$ch
        if(nbc == length(COL)){
          imprints_barplotting_app(data_cell(), witherrorbar = input$werb_cell,
                                   withpoint = input$wpts_cell,
                                   usegradient = input$grad_cell, linegraph = input$line_cell,
                                   save_pdf = input$save_bar_cell, colorpanel = COL,
                                   ret_plot = !input$save_bar_cell,
                                   layout = c(input$lay_bar1_cell, input$lay_bar2_cell),
                                   toplabel = loca_cell_lab,
                                   pdfname = input$pdftit_cell,
                                   pdfwidth = input$pdfw_cell, pdfheight = input$pdfh_cell)
        }
        else{
          showNotification("The number of colors given doesn't match the number of treatment selected !", type = "error")
        }

      }
      else{
        imprints_barplotting_app(data_cell(), witherrorbar = input$werb_cell,
                                 withpoint = input$wpts_cell,
                                 usegradient = input$grad_cell, linegraph = input$line_cell,
                                 save_pdf = input$save_bar_cell, ret_plot = !input$save_bar_cell,
                                 layout = c(input$lay_bar1_cell, input$lay_bar2_cell),
                                 toplabel = loca_cell_lab,
                                 pdfname = input$pdftit_cell,
                                 pdfwidth = input$pdfw_cell, pdfheight = input$pdfh_cell)
      }

    },
    message = function(m) {
      shinyjs::html(id = "diag_bar_cell", html = paste(m$message, "<br>", sep = ""), add = FALSE)}
    )
  })

  observeEvent(input$barp_cell, {
    if(input$selpr_loca_cell){
      if(input$allpr_cell){
        pr <- sel_prot_cell()
      }
      else{
        pr <- input$selectpr_cell
      }
      if(!is.null(pr)){
        pr <- unname(sapply(pr, function(x) strsplit(x, ":")[[1]][1]))
      }
    }
    else{
      pr <- PR_event()
    }

    if(is.null(pr)){
      showNotification("Don't forget to select a protein !", type = "error")
    }
    else{
      if (input$cond_sel_cell != "all_cond"){
        if (is.null(input$cond_cell)){
          showNotification("Don't forget to select a treatment !", type = "error")
        }
        else{
          BAR_cell$ch <- Bar_one_cell()
        }
      }
      else{
        BAR_cell$ch <- Bar_one_cell()
      }
    }


  })

  output$bar_pr_cell <- renderPlot({
    BAR_cell$ch
  })

  output$downbar_cell <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%y%m%d_%H%M_"), "2D_barplot", ".", input$downbar_cell_format)
    },
    content = function(file){
      ggsave(file, plot = BAR_cell$ch[[1]], device = input$downbar_cell_format)
    }
  )


  ### PubMed

  pubmed_data <- reactive({
    File <- input$data_pubmed
    if (is.null(File))
      return(NULL)

    import_list(File$datapath)[[1]]
  })
  #check if a file is upload
  output$pubmed_fileup <- reactive({
    return(!is.null(pubmed_data()))
  })
  outputOptions(output, "pubmed_fileup", suspendWhenHidden = FALSE)

  observe({
    if(!input$impc_pubmed){
      updateCheckboxInput(session, "hit_pubmed", value = FALSE)
    }
  })

  observeEvent(input$go_pub, {
    showNotification("Start searching", type = "message")

    if(input$impc_pubmed){
      data <- pubmed_data()
    }
    else{
      data <- input$dtext_pubmed
      data <- str_split(data, ",")[[1]]
      data <- str_trim(data)
    }

    if (str_length(str_remove_all(input$LA_pubmed, " ")) == 0){
      LA_ <- NULL
    }
    else{
      LA_ <- input$LA_pubmed
    }
    if (str_length(str_remove_all(input$Y_pubmed, " ")) == 0){
      Y_ <- NULL
    }
    else{
      Y_ <- input$Y_pubmed
    }
    if (str_length(str_remove_all(input$api_pubmed, " ")) == 0){
      api_ <- NULL
    }
    else{
      api_ <- input$api_pubmed
    }

    withCallingHandlers({
      shinyjs::html("diag", "")
      pub <- find_in_pubmed(data, feat = input$feat_pubmed, imp_by_hitlist = input$hit_pubmed,
                            language = LA_, year_rg = Y_, treatment = input$cond_pubmed,
                            your_API = api_, newfolder_name = input$fname_pubmed,
                            save_word = input$save_in_word)
    },
    message = function(m) {
      shinyjs::html(id = "diag", html = paste(m$message, "<br>", sep = ""), add = FALSE)

    }
    )


    output$pubmed_out <- DT::renderDataTable({
      DT::datatable(pub,
                    caption = htmltools::tags$caption(
                      style = 'caption-side: top; text-align: left;',
                      htmltools::strong("PubMed search results")
                    ),
                    rownames = FALSE,
                    options = list(lengthMenu = c(10,20,30), pageLength = 10))
    })

    output$down_pubmed <- downloadHandler(
      filename = function() {
        paste0(format(Sys.time(), "%y%m%d_%H%M_"), "have_publication_", input$feat_pubmed, ".xlsx")
      },
      content = function(file){
        xlsx::write.xlsx(pub, file, row.names = FALSE)
      }
    )


  })
}


shinyApp(ui, server)
