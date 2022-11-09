library(mineCETSA)
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

ui <-  navbarPage(title = img(src="logo.png", height = "40px"),
                 id = "navBar",
                 theme = "paper.css", # file in www
                 collapsible = TRUE, # usefull when viewing on smaller screen
                 inverse = FALSE, # true: use a dark background and light text for the navigation bar
                 windowTitle = "mineCETSAapp", # just name of onglet
                 position = "fixed-top",
                 footer = includeHTML("./www/include_footer.html"), # bottom of the page/site
                 header = tagList(
                   shinyWidgets::useShinydashboard(),      # allow to render the boxes from shinydachboard
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
                                   shiny::HTML("<h5>In 2013 we published the first paper describing the transformative Cellular Thermal Shift Assay (CETSA, Martinez Molina, Science, 41:84).
                                                    CETSA constitutes the first broadly applicable method to assess direct drug binding in cells and tissues.
                                                    The method has had a big impact in drug discovery and is now being broadly applied in academia and
                                                    industry to improve efficiency and quality of drug candidates and has the potential to serve as an
                                                    important clinical diagnostic for drug efficacy in the future. Recently we have published papers demonstrating
                                                    that CETSA constitutes a highly resolved means to study interactions of proteins with other cellular
                                                    components in intact cells and tissues at the proteome level. This approach gives a completely novel
                                                    perspective on how cellular processes are executed and we predict that it will have a very big impact on
                                                    understanding disease processes and drug action in the future. The method also constitutes a way to
                                                    uncover novel drug targets and therapeutic biomarkers with future applications in cancer therapy.

                                                    <br><br> If you want to learn more about CETSA and our lab, click <a href='https://www.cetsa.org/'>here</a></h5>"
                                               )
                                   )
                            ),

                          tags$hr(),

                          fluidRow(
                            column(12,
                                   shiny::HTML("<h1>mineCETSAapp</h1><br>"),
                                   shiny::HTML("<h5>mineCETSAapp is an R package that include a shiny app that you can use to
                                                    easily analyse your imprints-CETSA data. In this app, you will be able
                                                    to process your quantification files, get your drug targets, visualize your
                                                    barplots, run some gene ontology analysis and more. <br>
                                                    The app also contains a database of more than 100 drugs tested by our lab.
                                                    You can already visualize and compare these. You can also modify this database by
                                                    adding or removing new datasets. <br>
                                                    You can learn how to use the app by seeing the tutorial video which you can
                                                    access by clicking on the question mark icon in the top right corner.
                                                    <br><br><br>
                                                    <center> Otherwise, start your analysis here ! </center></h5>")
                                   )
                            ),
                          fluidRow(
                            column(3),
                            column(6,
                                   tags$div(align = "center",
                                            tags$a("Start",
                                                   onclick="fakeClick('analysis')", # take you to other tab with value 'analysis'
                                                   class="btn btn-primary btn-lg")
                                   )
                            ),
                            column(3)
                          ),

                          tags$hr()
                          ),

                 tabPanel("Analysis", value = "analysis",
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

                          fluidRow(style = "height:20px;"),

                          h1(tags$u(class = "main-1", "The mineCETSA analysis")),

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
                                       fluidRow(column(4, selectInput("n_chan", "Choose the numeber of channels", choices = c(10,11,16,18), selected = 10)),

                                                column(8, uiOutput("treat_nameui"))
                                                ),


                                       fileInput("PD_data", "Select txt files for your analysis",
                                                 accept = ".txt", multiple = TRUE),


                                       conditionalPanel(condition = "output.cetsa_fileup",
                                                        actionButton("see1_cetsa", "View data uploaded"),
                                                        tags$hr(),
                                                        tags$u(h3("Rename your conditions and clean your data")),
                                                        tags$hr(),

                                                        fluidRow(column(4, textOutput("incond_name")),
                                                                 column(4, textInput("outcond", "Type the new condition naming, spaced by a comma",
                                                                                     "37C,47C,50C,52C,54C,57C,36C"),
                                                                        checkboxInput("rem_mix", "Remove the 'Mix' channel", TRUE),
                                                                        checkboxInput("clean_data", "Remove proteins without quantitative information", TRUE)),
                                                                 column(2, actionButton("str_ren", "Rename the conditions", class = "btn-primary")),
                                                                 column(2, actionButton("see2_cetsa", "View data renamed"))
                                                        )
                                       )
                          )),
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
                                                                                                          conditionalPanel(condition = "input.iso_conso",
                                                                                                                           numericInput("n_chan2", "Type the number of reading channels", value = 9, min = 1),
                                                                                                                           fileInput("tab_conso", "Upload the txt file containing an isoform substitution matching table",
                                                                                                                                     accept = ".txt")
                                                                                                          )
                                                                                          ),
                                                                                          column(6, checkboxInput("iso_rearr", "Rearrange data", TRUE),
                                                                                                 conditionalPanel(condition = "input.iso_rearr",
                                                                                                                  numericInput("n_chan3", "Type the number of reading channels", value = 9, min = 1),
                                                                                                                  numericInput("rep_thr", "Type the minimal percentage threshold of
                                                                                                                               protein being sampled from multiple runs", value = 0.1, min = 0, max = 1, step = 0.01),
                                                                                                                  numericInput("count_thr", "Type the minimal threshold number
                                                                                                                                                          of associated abundance count of proteins", value = 1, min = 0, step = 0.5),
                                                                                                                  checkboxInput("wit_37", "Whether the kept proteins should have readings at 37C", FALSE)
                                                                                                 )
                                                                                          )
                                                                                          ),
                                                                                          actionButton("ISO2", "Consolidate isoform and/or rearrange", class = "btn-primary")
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
                                                                                          fluidRow(column(4, textInput("treat_name2", "Type the treatment name, spaced by a comma",
                                                                                                                       "Vehicle,Alpelisib,Buparlisib")),

                                                                                                   column(4, checkboxInput("wit_rep", "Whether the calculation of the relative protein
                                                                                                                           abundance difference should still within the same biorep", TRUE)),

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
                                                                                                                        choices = c("Fold Change cutoff" = "FC",
                                                                                                                                    "Stability Rate" = "SR"),
                                                                                                                        selected = "FC",
                                                                                                                        inline = TRUE),
                                                                                                           conditionalPanel(condition = "input.hitmethod_cetsa == 'FC'",
                                                                                                                            fluidRow(column(4, numericInput("meancut_cetsa", "Choose a mean cutoff", value = 0.25, min = 0, step = 0.01)),
                                                                                                                                     column(4, numericInput("bound_cetsa", "Choose the boundedness", value = 4)),
                                                                                                                                     column(4, checkboxInput("save_hit", "Save the hitlist", TRUE))
                                                                                                                            )
                                                                                                           ),
                                                                                                           conditionalPanel(condition = "input.hitmethod_cetsa == 'SR'",
                                                                                                                            fluidRow(column(4, numericInput("SRcut_cetsa", "Choose a Stability Rate cutoff", value = 1.5, min = 0, step = 0.1)),
                                                                                                                                     column(4, numericInput("FDR_cetsa", "Choose the FDR", value = 0.01, min = 0, max = 1, step = 0.01)),
                                                                                                                                     column(4, numericInput("validval_cetsa", "Choose the minimum proportion of valid values", value = 0, min = 0, max = 1, step = 0.05))
                                                                                                                            ),
                                                                                                                            tags$hr(),
                                                                                                                            textOutput("diag_SR"),
                                                                                                                            tags$hr(),
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
                          ),

                 tabPanel("Database", value = "database",
                          shinyjs::useShinyjs(),
                          fluidRow(style = "height:20px;"),

                          h1(tags$u(class = "main-1", "Add new dataset and remove old ones")),
                          tags$hr(),

                          htmlOutput("info_daba"),
                          tags$hr(),

                          HTML("<p><h5>In order to add new dataset, you need to import three files.<br>
                              This files are : <br>
                              - The output from the ms_2D_caldiff function from the mineCETSA package <br>
                              - The file named 'Summary', from the hitlist function output <br>
                              - The file named 'NN', from the hitlist function output <br>
                              Once you uploaded this three files, choose a name for your dataset (like 'elutriation' for example),
                              click on the button 'Add dataset', and you're good to go !</h5></p>
                              <p><h5>If you want to remove a dataset, select one of the dataset available from the database,
                              and click on the button 'Remove dataset', in the box below. Beware, this operation cannot be undone !</h5></p>
                              <br><h5>Once you made your changements, don't forget to click on the button 'Reload the database' to use directly.</h5>"
                          ),

                          fluidRow(box(title = "Add new dataset", status = "success", solidHeader = TRUE, collapsible = TRUE, width = 12,
                                       fluidRow(column(4, fileInput("caldif_daba", "Import the output from ms_2D_caldiff"),
                                                       checkboxInput("gave_daba", "Don't have the IMPRINTS_average output
                                                                         (will calculate and save it)", TRUE),
                                                       conditionalPanel(condition = "!input.gave_daba",
                                                                        fileInput("AVE_dabafile", "Import the output from IMPRINTS_average_sh")
                                                       )
                                       ),
                                       column(4, fileInput("hitsum_daba", "Import the summary file from the hitlist outputs")),
                                       column(4, fileInput("NN_daba", "Import the NN file from the hitlist outputs")),
                                       ),
                                       conditionalPanel(condition = "output.DIFdaba_fileup & output.AVEdaba_fileup & output.HITdaba_fileup & output.NNdaba_fileup",
                                                        fluidRow(column(4, textInput("name_daba", "Type a name for your new dataset")),
                                                                 column(4, actionButton("add_daba", "Add dataset", class = "btn-success btn-lg"))
                                                        )
                                       )

                          )),

                          fluidRow(box(title = "Rename your conditions", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
                                       fluidRow(column(6, uiOutput("davai2_daba_ui")),
                                                column(6, htmlOutput("condfrom_daba"),
                                                       textInput("condnew_daba", "Type the new names of the conditions
                                                                   (same order; separated by a comma; if empty, no changement)"))
                                       ),
                                       actionButton("changename_daba", "Change the name of your conditions", class = "btn-primary btn-lg")
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

                                                                                        fluidRow(column(4, fileInput("data_barplot", "Upload your own data (output from mineCETSA package)",
                                                                                                                     accept = c(".txt", ".csv", ".xlsx"))
                                                                                        ),
                                                                                        column(8, checkboxInput("calc_hitlist", "Find the hitlist from your data file", FALSE),
                                                                                               conditionalPanel(condition = "!input.calc_hitlist",
                                                                                                                fileInput("data_hitlist", "Upload your own hitlist (summary and NN file from hitlist function)",
                                                                                                                          accept = c(".txt", ".csv", ".xlsx"), multiple = TRUE)),
                                                                                               conditionalPanel(condition = "input.calc_hitlist",
                                                                                                                fluidRow(column(3, numericInput("meancut_bar", "Choose a mean cutoff", value = 0.25, min = 0, step = 0.01)),
                                                                                                                         column(3, numericInput("bound_bar", "Choose the boundedness", value = 4)),
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
                                                                                                                  selectInput("cond_fhit", "Select hits condition", choices = NULL, multiple = TRUE))
                                                                                ),
                                                                                conditionalPanel(condition = "input.protlist_bar",
                                                                                                 fileInput("prlist_file_bar", "Import your protein list (txt file)", accept = ".txt"),
                                                                                ),
                                                                                checkboxInput("ALL_prot", "Select all the proteins", FALSE),
                                                                                checkboxInput("alliso_bar", "Take all isoform", FALSE),
                                                                                conditionalPanel(condition = "!input.ALL_prot",
                                                                                                 selectizeInput("prot", "Select a protein", choices = NULL, multiple = TRUE))

                                                                         ),
                                                                         column(4, conditionalPanel(condition = "input.cond_sel != 'cat' ",
                                                                                                    checkboxInput("rem_con", "Remove the controls", FALSE),
                                                                                                    conditionalPanel(condition = "input.rem_con",
                                                                                                                     textInput("con_name", "Type the name of your controls (if sevral names, separate them by |)", "G1")
                                                                                                    )
                                                                         )
                                                                         ),
                                                                         column(4, radioButtons("cond_sel", "Selection type",
                                                                                                choices = c("Select the treatment level" = "treat",
                                                                                                            "Select treatment level by category" = "cat",
                                                                                                            "Select all the treatment level" = "all_cond"),
                                                                                                selected = "treat"),

                                                                                conditionalPanel(condition = "input.cond_sel != 'all_cond' ",
                                                                                                 selectInput("cond", "Select one or more conditions",
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

                                                                       fluidRow(column(4, checkboxInput("save_bar", "Save the bar plots in a pdf file", FALSE)),
                                                                                conditionalPanel(condition = "input.save_bar",
                                                                                                 column(4, numericInput("lay_bar1", "Type the number of plot per row",
                                                                                                                        min = 1, max = 10, step = 1, value = 4),
                                                                                                        numericInput("lay_bar2", "Type the number of plot per column",
                                                                                                                     min = 1, max = 10, step = 1, value = 3)),
                                                                                                 column(4, textInput("pdftit", "Choose a name for your pdf file", "barplot"))
                                                                                )
                                                                       ),

                                                                       tags$hr(),

                                                                       fluidRow(column(4, checkboxInput("werb", "Print error bar", TRUE)),
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
                                                          downloadButton("downbar", "Download the plot as png file"),

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
                                                                                        fluidRow(column(4, fileInput("caldif_compl", "Import the output from ms_2D_caldiff")),
                                                                                                 column(4, fileInput("hitsum_compl", "Import the summary file from the hitlist outputs")),
                                                                                                 column(4, fileInput("NN_compl", "Import the NN file from the hitlist outputs"))
                                                                                        ),
                                                                                        fluidRow(column(4, checkboxInput("gave_compl", "Don't have the IMPRINTS_average output
                                                                                     (will calculate and save it)", TRUE)),
                                                                                                 conditionalPanel(condition = "!input.gave_compl",
                                                                                                                  column(4, fileInput("avef_compl", "Import the output from IMPRINTS_average"))
                                                                                                 )
                                                                                        )
                                                                       ),
                                                                       tags$hr(),

                                                                       conditionalPanel(condition = "output.DIFcompl_fileup & output.HITcompl_fileup & output.NNcompl_fileup & output.AVEcompl_fileup",
                                                                                        fluidRow(column(4, selectInput("condsel_compl", "Select a condition", choices = NULL)),
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
                                                                           fluidRow(box(title = "2D bar plot paramter", status = "primary",
                                                                                        solidHeader = TRUE, collapsible = TRUE, width = 12,

                                                                                        fluidRow(column(4,selectInput("allcomplex_compl", "Select some protein complex", choices = NULL, multiple = TRUE)),
                                                                                                 column(4, checkboxInput("ALL_prot_compl", "Select all the proteins", FALSE),
                                                                                                        checkboxInput("alliso_bar_compl", "Take all isoform", FALSE)),
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

                                                                                        fluidRow(column(4, checkboxInput("save_bar_compl", "Save the bar plots in a pdf file", FALSE)),
                                                                                                 conditionalPanel(condition = "input.save_bar_compl",
                                                                                                                  column(4, numericInput("lay_bar1_compl", "Type the number of plot per row",
                                                                                                                                         min = 1, max = 10, step = 1, value = 4),
                                                                                                                         numericInput("lay_bar2_compl", "Type the number of plot per column",
                                                                                                                                      min = 1, max = 10, step = 1, value = 3)),
                                                                                                                  column(4, textInput("pdftit_compl", "Choose a name for your pdf file", "barplot"))
                                                                                                 )
                                                                                        ),

                                                                                        tags$hr(),

                                                                                        fluidRow(column(4, checkboxInput("werb_compl", "Print error bar", TRUE)),
                                                                                                 column(4, checkboxInput("grad_compl", "Use color gradient", FALSE)),
                                                                                                 column(4, checkboxInput("line_compl", "Use line instead of bar", FALSE))
                                                                                        ),

                                                                                        tags$hr(),
                                                                                        actionButton("barp_compl", "See bar plot", class = "btn-primary btn-lg"),
                                                                                        tags$hr(),
                                                                                        textOutput("diag_bar_compl"),
                                                                                        tags$hr(),

                                                                                        withSpinner(plotOutput("bar_plot_compl", height = "800px"), type = 6),
                                                                                        downloadButton("downbar_compl", "Download the plot as png file")
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
                                                                                        fluidRow(column(4, fileInput("cdiff_simpf", "Import the output from ms_2D_caldiff")),
                                                                                                 column(4, checkboxInput("gave_simpf", "Don't have the IMPRINTS_average output
                                                                                     (will calculate and save it)", TRUE)),
                                                                                                 conditionalPanel(condition = "!input.gave_simpf",
                                                                                                                  column(4, fileInput("avef_simpf", "Import the output from IMPRINTS_average"))
                                                                                                 )
                                                                                        )
                                                                       )
                                                                       ,

                                                                       conditionalPanel(condition = "output.AVEsimpf_fileup & output.DIFsimpf_fileup",
                                                                                        fluidRow(column(3, selectInput("treat_simpf", "Select a condition", choices = NULL)),
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
                                                                                     selected will be calculate. Then for each distance we calculate a score between 0 and 1
                                                                                     by dividing 1 by 1 + d, where d is the euclidean distance. <br>
                                                                                     This score means that you will search for protein profile with similar values from the the one you selected.
                                                                                     So the profile with a similar shape but with lower or higher values will not have a good score.
                                                                                     It also means that with a high score (~0.9) you're not very likely to find a lot of proteins.<br>
                                                                                     <br>
                                                                                     This is not the case with Pearson correlation. For this score, each covariance and standard deviation
                                                                                     between the protein you selected and all the other proteins will be calculated. Then, the covariance is divided
                                                                                     by the product of the two standard devation. It gives you score between -1 and 1. -1 means the data are negatively
                                                                                     correlated, 1 positively correlated and 0 not correlated. <br>
                                                                                     Because you calculate a correlation score, you will search for all proteins profile with a similar shape from the one
                                                                                     you selected, not matter their values. It's like searching mountains with similar shapes, no matter their height.
                                                                                     It also means that with a high score (~0.95) you may find a lot of proteins.
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
                                                                                                        checkboxInput("save_prot_simpf", "Save the list of proteins ID with similar profile (will save in a xlsx file)", TRUE)),
                                                                                                 conditionalPanel(condition = "input.save_bar_simpf",
                                                                                                                  column(4, numericInput("lay_bar1_simpf", "Type the number of plot per row",
                                                                                                                                         min = 1, max = 10, step = 1, value = 4),
                                                                                                                         numericInput("lay_bar2_simpf", "Type the number of plot per column",
                                                                                                                                      min = 1, max = 10, step = 1, value = 3)
                                                                                                                  ),
                                                                                                                  column(4, textInput("pdftit_simpf", "Choose a name for your pdf file", "barplot"))
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
                                                                                        downloadButton("downbar_simpf", "Download the plot as png file")
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
                                                                                        fluidRow(column(4, checkboxInput("gave_heat", "Don't have the IMPRINTS_average output
                                                                                                  (will calculate and save it)", TRUE),
                                                                                                        conditionalPanel(condition = "input.gave_heat",
                                                                                                                         fileInput("filedif_heat", "Choose a ms_2D_caldiff output")
                                                                                                        ),
                                                                                                        conditionalPanel(condition = "!input.gave_heat",
                                                                                                                         fileInput("fileave_heat", "Choose a IMPRINTS_average_sh output")
                                                                                                        )
                                                                                        ),
                                                                                        column(4, fileInput("summary_heat", "Choose the summary file from the hitlist output")),
                                                                                        column(4, checkboxInput("impNN_heat", "Also import the NN file from hitlist output", FALSE),
                                                                                               conditionalPanel(condition = "input.impNN_heat",
                                                                                                                fileInput("NNfile_heat", "Choose the summary file from the hitlist output")
                                                                                               )
                                                                                        )
                                                                                        )
                                                                       ),

                                                                       tags$hr(),

                                                                       conditionalPanel(condition = "output.heat_fileup & output.HITheat_fileup & output.NNheat_fileup",
                                                                                        fluidRow(column(3, selectInput("cond_heat", "Select a condition", choices = NULL)),
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
                                                                                                                  column(4, textInput("fnameH_heat", "Type a your file name", "My_heatmap")),
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
                                                                                        fluidRow(column(4, checkboxInput("gave_heatcom", "Don't have the IMPRINTS_average output
                                                                                                  (will calculate and save it)", TRUE)),
                                                                                                 column(4,
                                                                                                        conditionalPanel(condition = "input.gave_heatcom",
                                                                                                                         fileInput("filedif_heatcom", "Choose a ms_2D_caldiff output")
                                                                                                        ),
                                                                                                        conditionalPanel(condition = "!input.gave_heatcom",
                                                                                                                         fileInput("fileave_heatcom", "Choose a IMPRINTS_average_sh output")
                                                                                                        )
                                                                                                 )
                                                                                        )
                                                                       ),

                                                                       tags$hr(),

                                                                       conditionalPanel(condition = "output.heatcom_fileup",
                                                                                        fluidRow(column(4, selectInput("cond_heatcom", "Select a condition", choices = NULL)),
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
                                                                                                                  column(4, textInput("fnameH_heatcom", "Type a your file name", "My_heatmap")),
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
                                     fluidRow(style = "height:20px;"),

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
                                                                   fluidRow(column(6, selectInput("cond_fhitB_stri", "Select some conditions for filtering your proteins",
                                                                                                  choices = NULL, multiple = TRUE)),
                                                                            column(6, selectInput("cat_fhitB_stri", "Select some categories for filtering your proteins (If NULL, will select all)",
                                                                                                  choices = c("CN", "NC", "CC", "ND", "NN"), multiple = TRUE))
                                                                   )
                                                  ),

                                                  conditionalPanel(condition = "input.drug_stri == 'dat' ",
                                                                   checkboxInput("impfile_stri", "Import a file", TRUE),
                                                                   conditionalPanel(condition = "input.impfile_stri",
                                                                                    fluidRow(column(4, fileInput("file_stri", "Choose a file")),
                                                                                             column(4, checkboxInput("ishit_stri", "Do you import a hitlist ?", TRUE),
                                                                                                    conditionalPanel(condition = "!input.ishit_stri",
                                                                                                                     textInput("idfile_stri", "What is the name of the column of
                                                                                                your file which contains the proteins ID ?")
                                                                                                    )
                                                                                             ),
                                                                                             conditionalPanel(condition = "input.ishit_stri",
                                                                                                              column(4, selectInput("cond_fhit_stri", "Select some condition for filtering your hits",
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
                                                                            column(3, numericInput("intscore1_stri", "Minimum interaction score",
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
                                                                                                     downloadButton("downnet_stri", "Download the plot as png file")
                                                                                    )
                                                                   )
                                                  )
                                     )),


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
                                                                                    fluidRow(column(6, actionButton("netfilt_stri", "See new network", class = "btn-primary btn-lg"))
                                                                                    ),
                                                                                    tags$hr(),

                                                                                    conditionalPanel(condition = "output.enrich_res_tab_up",
                                                                                                     conditionalPanel(condition = "!input.hidnet2_stri",
                                                                                                                      conditionalPanel(condition = "input.intnet_stri",
                                                                                                                                       withSpinner(plotlyOutput("netInt2_stri", height = "800px"), type = 6)
                                                                                                                      ),
                                                                                                                      conditionalPanel(condition = "!input.intnet_stri",
                                                                                                                                       withSpinner(plotOutput("net2_stri", height = "800px"), type = 6),
                                                                                                                                       downloadButton("downnetfilt_stri", "Download the plot as png file")
                                                                                                                      )
                                                                                                     )
                                                                                    )
                                                                   )
                                                      ))
                                     )
                                     ),

                            tabPanel("ClusterProfile", value = "clusprof",
                                     shinyjs::useShinyjs(),
                                     tags$style(HTML(".tabbable > .nav > li > a                  {background-color: #A1BAC8;  color:#FFFFFF}
                                                                  .tabbable > .nav > li[class=active]    > a {background-color: #3C8DBC; color:#FFFFFF}
                                                                  .tabbable > .nav > li    > a:hover {background-color: #3BAAE6; color:#FFFFFF}
                                                                  .tabbable > .nav > li[class=active]    > a:hover {background-color: #3C8DBC; color:#FFFFFF}")
                                     ),

                                     fluidRow(style = "height:20px;"),

                                     h1(tags$u(class = "main-1", "Enrichment analysis from wiki-pathway")),
                                     tags$hr(),

                                     fluidRow(box(title = "Import your data and start the analysis", status = "primary",
                                                  solidHeader = TRUE, collapsible = TRUE, width = 12,

                                                  radioButtons("drug_clus", h3("Choose a dataset"),
                                                               choices = c("Database" = "base",
                                                                           "Your data" = "dat"),
                                                               selected = "base", inline = TRUE),
                                                  conditionalPanel(condition = "input.drug_clus == 'base'",
                                                                   uiOutput("drug2ui_clus"),
                                                                   fluidRow(column(6, selectInput("cond_fhitB_clus", "Select some conditions for filtering your proteins",
                                                                                                  choices = NULL, multiple = TRUE)),
                                                                            column(6, selectInput("cat_fhitB_clus", "Select some categories for filtering your proteins (If NULL, will select all)",
                                                                                                  choices = c("CN", "NC", "CC", "ND", "NN"), multiple = TRUE))
                                                                   )
                                                  ),

                                                  conditionalPanel(condition = "input.drug_clus == 'dat' ",
                                                                   fluidRow(column(4, fileInput("file_clus", "Choose a file")),
                                                                            column(4, checkboxInput("ishit_clus", "Do you import a hitlist ?", TRUE),
                                                                                   conditionalPanel(condition = "!input.ishit_clus",
                                                                                                    textInput("idfile_clus", "What is the name of the column of
                                                                                                                                            your file which contains the proteins ID ?")
                                                                                   )
                                                                            ),
                                                                            conditionalPanel(condition = "input.ishit_clus",
                                                                                             column(4, selectInput("cond_fhit_clus", "Select some condition for filtering your hits",
                                                                                                                   choices = NULL, multiple = TRUE))
                                                                            )
                                                                   )
                                                  )
                                     )
                                     ),

                                     tabsetPanel(type = "tabs",

                                                 tabPanel("Compare cluster",
                                                          fluidRow(box(title = "Compare the wiki-pathway of your data between conditions", status = "primary",
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
                                                                       downloadButton("downcomplot_clus", "Download the plot as png file")
                                                                       )
                                                                   )
                                                          ),
                                                 tabPanel("GSEA",
                                                          fluidRow(box(title = "Make a GSEA on your data", status = "primary",
                                                                       solidHeader = TRUE, collapsible = TRUE, width = 12,
                                                                       fluidRow(column(4, textInput("scorename_clus", "Type the column's name that contain the score")),
                                                                                column(4, checkboxInput("onlypos_clus", "Show only enrcihment set with positive enrichment score", TRUE)),
                                                                                column(4, actionButton("gogsea_clus", "Start GSEA", class = "btn-primary btn-lg"))
                                                                                ),
                                                                       DT::dataTableOutput("gseatab_clus"),
                                                                       downloadButton("downgseatab_clus"),
                                                                       tags$hr(),
                                                                       withSpinner(plotOutput("gseaplot_clus", height = "800px"), type = 6),
                                                                       downloadButton("downgsealot_clus", "Download the plot as png file")
                                                                       )
                                                                   )
                                                          ),
                                                 tabPanel("Gene concept network",
                                                          fluidRow(box(title = "View a gene concept network from your data", status = "primary",
                                                                       solidHeader = TRUE, collapsible = TRUE, width = 12,
                                                                       fluidRow(column(4, textInput("scorename2_clus", "Type the column's name that contain the score")),
                                                                                column(4, numericInput("pvcut_clus", "Choose a p-value cutoff for the gene concept network",
                                                                                                       value = 0.01, min = 0, max = 1, step = 0.01)),
                                                                                column(4, actionButton("gogeneconc_clus", "See gene concept network", class = "btn-primary btn-lg"))
                                                                                ),
                                                                       withSpinner(plotOutput("geneplot_clus", height = "800px"), type = 6),
                                                                       downloadButton("downgenelot_clus", "Download the plot as png file")
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
                                                 fluidRow(column(4, selectInput("condhit_cell", "Select a condition", choices = NULL)),
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
                                                          column(4, selectInput("condp_cell", "Select some conditions", multiple = TRUE, choices = NULL)),
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
                                              please import the cal_diff output file which correspond to your hitlist.")
                                ),
                                conditionalPanel(condition = "!input.selpr_loca_cell",
                                                 htmlOutput("prsel_p_cell")
                                ),

                                conditionalPanel(condition = "output.barpdata_cell_up | input.drug_cell == 'base'",
                                                 fluidRow(column(4, checkboxInput("selpr_loca_cell", "Select proteins according to their subcellular location", FALSE)),
                                                          conditionalPanel(condition = "input.selpr_loca_cell",
                                                                           column(4, selectInput("selorga_cell", "Select some organelle", multiple = TRUE, choices = NULL)),
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
                                                                  selectInput("cond_cell", "Select one or more conditions",
                                                                              choices = NULL,
                                                                              multiple = TRUE)
                                                 ),

                                                 tags$hr(),

                                                 fluidRow(column(4, checkboxInput("save_bar_cell", "Save the bar plots in a pdf file", FALSE)),
                                                          conditionalPanel(condition = "input.save_bar_cell",
                                                                           column(4, numericInput("lay_bar1_cell", "Type the number of plot per row",
                                                                                                  min = 1, max = 10, step = 1, value = 4),
                                                                                  numericInput("lay_bar2_cell", "Type the number of plot per column",
                                                                                               min = 1, max = 10, step = 1, value = 3)),
                                                                           column(4, textInput("pdftit_cell", "Choose a name for your pdf file", "barplot"))
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

                                                 fluidRow(column(4, checkboxInput("werb_cell", "Print error bar", TRUE)),
                                                          column(4, checkboxInput("grad_cell", "Use color gradient", FALSE)),
                                                          column(4, checkboxInput("line_cell", "Use line instead of bar", FALSE))
                                                 ),

                                                 tags$hr(),

                                                 actionButton("barp_cell", "See bar plot", class = "btn-primary btn-lg"),
                                                 tags$hr(),
                                                 textOutput("diag_bar_cell"),
                                                 tags$hr(),

                                                 withSpinner(plotOutput("bar_pr_cell", height = "800px"), type = 6),
                                                 downloadButton("downbar_cell", "Download the plot as png file")
                                )
                            )
                          )
                          ),

                 tabPanel("PubMed search", value = "pubmed",
                          fluidRow(style = "height:20px;"),
                          shinyjs::useShinyjs(),

                          h1(tags$u(class = "main-1", "Search publications in PubMed")),
                          tags$hr(),

                          fluidRow(box(title = "Search parameters", status = "primary",
                                       solidHeader = TRUE, collapsible = TRUE, width = 12,
                                       fluidRow(column(3, checkboxInput("impc_pubmed", "Import a file", TRUE),
                                                       conditionalPanel(condition = "input.impc_pubmed",
                                                                        fileInput("data_pubmed", "Import your data")),
                                                       conditionalPanel(condition = "!input.impc_pubmed",
                                                                        textInput("dtext_pubmed", "Type some protein names, separated by a comma", ""))
                                                       ),
                                       column(3, textInput("feat_pubmed", "Type your second research word", "cell cycle")),
                                       column(3, textInput("LA_pubmed", "Type a language to match (can be null)")),
                                       column(3, textInput("Y_pubmed", "Type a year range to match (can be null, format is Y1:Y2)"))
                                       ),

                                       fluidRow(column(3, textInput("api_pubmed", "Type your NCBI API if you have an account")),
                                                column(3, textInput("fname_pubmed", "Type the name of the folder that will be created", "Elutriation_pubmed_search")),
                                                column(3, conditionalPanel(condition = "input.impc_pubmed",
                                                                           checkboxInput("hit_pubmed", "Do you import a hitlist ? (need description column)", TRUE))
                                                       ),
                                                conditionalPanel(condition = "input.hit_pubmed",
                                                                 column(3, textInput("cond_pubmed", "Type a condition from you hitlist (if null, will take all the conditions)"))
                                                                 )
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

                 bslib::nav_item(tags$a(href = "https://github.com/mgerault/mineCETSAapp",
                        icon("github"),
                        title = "See source code to the github repository"), class = "icon1"),
                 bslib::nav_item(tags$a(href = "https://youtu.be/cOlOzU7-S3A",
                        icon("question-circle"),
                        title = "See the tutorial video of the app"), class = "icon2"),
                 bslib::nav_item(tags$a(href = "mailto:marco.gerault@gmail.com",
                        icon("envelope"),
                        title = "Any questions, suggestions or bug report ? Feel free to send me an e-mail !"), class = "icon3")


)

server <- function(input, output, session){
  setwd(WD)

  ### analysis tab
  output$treat_nameui <- renderUI({
    TMT <- list("10" = c("126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131"),
                "11" = c("126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C"),
                "16" = c("126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N"),
                "18" = c("126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N", "134C", "135N")
    )
    m <- matrix("", as.numeric(input$n_chan), 1,
                dimnames = list(c(TMT[[input$n_chan]]), "Treatment"))

    matrixInput("treat_name", "Type the name of your channel",
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

    ms_2D_rawread(File$datapath,
                  #the name of each treatment
                  treatment = as.character(input$treat_name[,1]),
                  #number of channel and the name of it
                  nread = as.numeric(input$n_chan),
                  channels = rownames(input$treat_name)
    )
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

  output$incond_name <- renderText({
    if(!is.null(cetsa_data())){
      paste("The current condition name are :", paste(unique(cetsa_data()$condition), collapse = "  "))
    }
    else{
      NULL
    }
  })

  cetsa_data_clean <- eventReactive(input$str_ren, {
    d1 <- ms_conditionrename(cetsa_data(),
                             incondition = unique(cetsa_data()$condition),
                             outcondition = str_remove(unlist(str_split(input$outcond, ",")), " ")
    )

    if(input$rem_mix){
      d1 <- d1[, !(names(d1) %in% "Mix")]
    }
    if(input$clean_data){
      d1 <- ms_clean(d1, nread = as.numeric(input$n_chan))
    }
    showNotification("Calculation done !", type = "message", duration = 3)


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
        showNotification("Consolidating succeed !", type = "message", duration = 5)
      }

      if(input$iso_rearr){
        showNotification("Start rearrange data", type = "message", duration = 3)
        d2 <- ms_2D_rearrange(d2, nread = input$n_chan3,
                              repthreshold = input$rep_thr, countthreshold = input$count_thr,
                              with37Creading = input$wit_37)
        cetsa_isoform$rearr <- d2
        showNotification("Rearranging succeed !", type = "message", duration = 5)
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
      showNotification("Start Normalization, this may take a while. Please wait a few minutes",
                       type = "message", duration = 5)
      d <- ms_2D_normalization(d)
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

  observeEvent(input$CAL_DIF, {
    if(!is.null(cetsa_isoform$norm)){
      showNotification("Start difference calculation, this may take a while. Please wait a few minutes",
                       type = "message", duration = 5)
      ytr_level <- str_remove(unlist(str_split(input$treat_name2, ",")), " ")
      tr_level <- get_treat_level(cetsa_isoform$norm)
      if(sum(str_detect(tr_level, paste(ytr_level, collapse = "|"))) != length(tr_level)){
        showNotification("The treatments you typed doesn't match the treatments from your data !", type = "error")
      }
      else{
        d <- ms_2D_caldiff(cetsa_isoform$norm,
                           treatmentlevel = ytr_level,
                           withinrep = input$wit_rep
        )

        cetsa_isoform$dif <- d
        message("Done to calculate the pair-wise (per replicate and temperature)
            protein abundance differences")
        showNotification("Difference calculation succeed !", type = "message", duration = 5)
      }
    }
    else{
      showNotification("Don't forget to import a file or start the analysis", type = "error")
    }
  })
  diffdata_cetsa <- reactive({
    File <- input$difffile_cetsa
    if (is.null(File) | !input$got_diff_cetsa)
      return(NULL)

    read.delim(File$datapath)
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

  hit_pr <- reactiveValues(
    hitlist = NULL,
    ND = NULL,
    NC = NULL,
    CN = NULL,
    CC = NULL
  )
  observeEvent(input$str_calchitlist, {
    if(!is.null(cetsa_isoform$dif)){
      if(input$hitmethod_cetsa == "FC"){
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
      else if(input$hitmethod_cetsa == "SR"){
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
            tidyr::separate(key, into = c("t", "b", "cond")) %>%
            dplyr::group_by(cond) %>%
            dplyr::summarise(ctrl = all(value == 0))
          ctrl <- ctrl$cond[ctrl$ctrl]

          withCallingHandlers({
            shinyjs::html("diag_SR", "")
            h <- SR_CetsaHit(cetsa_isoform$norm, Dif, ctrl = ctrl,
                             valid_val = input$validval_cetsa,
                             SR_cutoff = input$SRcut_cetsa,
                             FDR = input$FDR_cetsa)
          },
          message = function(m) {
            shinyjs::html(id = "diag_SR", html = paste(m$message, "<br>", sep = ""), add = FALSE)

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
      drug_data_sh$y <- loadData()
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
    selectInput("davai2_daba", "Choose a dataset to remove", choices = names(drug_data_sh$y$data),
                selected = "elutriation")
  })
  output$drug2_ui <- renderUI({
    selectInput("drug2", "Choose a drug", choices = names(drug_data_sh$y$data), multiple = TRUE, selected = "elutriation")
  })

  observeEvent(input$up_daba,{
    showNotification("Start loading the data, this may take a while", type = "message")
    a <- loadData()
    drug_data_sh$y <- a
    drug_data2 <<- drug_data_sh$y

    updateSelectInput(session, "davai_daba", choices = names(a$data), selected = names(a$data)[1])
    updateSelectInput(session, "drug2", choices = names(a$data), selected = names(a$data)[1])

    showNotification("Data loaded !", type = "message")
  })

  output$info_daba <- renderText({
    dn <- length(drug_data_sh$y$data)
    d <- names(drug_data_sh$y$data)

    if(dn > 1){
      d <- paste(c(paste(d[1:(dn-1)], collapse = ", "), d[dn]), collapse = " and ")
    }

    HTML(paste("<p><h4>Your dataset contains at the moment", dn, "drugs :", d, ".</h4></p>"))
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
      1  #simplify condition is.null
    }
  })
  #check if a file is upload
  output$AVEdaba_fileup <- reactive({
    return(!is.null(AVE_daba()))
  })
  outputOptions(output, "AVEdaba_fileup", suspendWhenHidden = FALSE)

  HIT_daba <- reactive({
    File <- input$hitsum_daba
    if (is.null(File))
      return(NULL)

    dat <- import(File$datapath, header = TRUE)
    nv_nam <- str_subset(names(dat), "^V\\d{1}$")
    if(!purrr::is_empty(nv_nam)){
      dat <- dat[, !(names(dat) %in% nv_nam)]
    }
    dat

  })
  #check if a file is upload
  output$HITdaba_fileup <- reactive({
    return(!is.null(HIT_daba()))
  })
  outputOptions(output, "HITdaba_fileup", suspendWhenHidden = FALSE)

  NN_daba <- reactive({
    File <- input$NN_daba
    if (is.null(File))
      return(NULL)

    dat <- import(File$datapath, header = TRUE)
    nv_nam <- str_subset(names(dat), "^V\\d{1}$")
    if(!purrr::is_empty(nv_nam)){
      dat <- dat[, !(names(dat) %in% nv_nam)]
    }
    dat <- dat[,c("id", "description", "Condition", "category")]
    dat
  })
  #check if a file is upload
  output$NNdaba_fileup <- reactive({
    return(!is.null(NN_daba()))
  })
  outputOptions(output, "NNdaba_fileup", suspendWhenHidden = FALSE)



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
        ave_data <- IMPRINTS_average_sh(DIF_daba())
        showNotification("Average calculation succeed !", type = "message")
      }
      showNotification("Start saving dataset, this may take a while.", type = "message")
      saveData(drug_data_sh$y, new_add = list("data" = DIF_daba(),
                                              "data_ave" = ave_data,
                                              "hitlist" = HIT_daba(),
                                              "NN" = NN_daba(),
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

    if(!is.null(df) & !purrr::is_empty(df)){
      cd <- get_treat_level(df)
      cd_1 <- cd[-length(cd)]
      cd_e <- cd[length(cd)]

      cd_info <- paste(paste(cd_1, collapse = ", "), cd_e, sep = " and ")
    }
    HTML(paste("<p>Your current condtion names for the drug", input$davai2_daba, "are :", paste0("<b>", cd_info, "</b>"), "</p>"))
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
      if(!purrr::is_empty(change) & !is.null(change)){
        new <- nm[!(nm %in% cd)]
        showNotification(paste("You decided to change :", paste(change, collapse = ", "),
                               "In :", paste(new, collapse = ", ")), type = "message")

        df <- drug_data_sh$y$data[[input$davai2_daba]]
        df_ave <- drug_data_sh$y$data_ave[[input$davai2_daba]]
        dh <- drug_data_sh$y$hitlist[[input$davai2_daba]]
        dnn <- drug_data_sh$y$NN[[input$davai2_daba]]

        n_df <- names(df)[str_detect(names(df), paste(paste0("_", change, "$"), collapse = "|"))]
        n_df_ave <- names(df_ave)[str_detect(names(df_ave), paste(paste0("_", change, "$"), collapse = "|"))]
        n_dh <- dh$Condition
        n_dnn <- dnn$Condition

        for(i in 1:length(change)){
          n_df <- str_replace_all(n_df, paste0("_", change[i], "$"), paste0("_", new[i]))
          n_df_ave <- str_replace_all(n_df_ave, paste0("_", change[i], "$"), paste0("_", new[i]))
          n_dh <- str_replace_all(n_dh, paste0("^", change[i], "$"), new[i])
          n_dnn <- str_replace_all(n_dnn, paste0("^", change[i], "$"), new[i])
        }

        names(df)[str_detect(names(df), paste(paste0("_", change, "$"), collapse = "|"))] <- n_df
        names(df_ave)[str_detect(names(df_ave), paste(paste0("_", change, "$"), collapse = "|"))] <- n_df_ave
        dh$Condition <- n_dh
        dnn$Condition <- n_dnn
        dt <- get_treat_level(df)

        showNotification("Start saving changes, this may take a while.", type = "message")
        saveData(drug_data_sh$y, new_add = list("data" = df,
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
        title="Are you sur you want to remove this dataset ?",
        footer = tagList(actionButton("confirmRem", "Remove"),
                         modalButton("Cancel")
        )
      )
    )
  })
  observeEvent(input$confirmRem, {
    showNotification("Start removing dataset", type = "message")
    remData(drug_data_sh$y, input$davai_daba)

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

    d <- import_list(File$datapath, header = TRUE)
    names(d) <- File$name

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

    hit_bar$summa <- h_s %>% group_by(id,Condition,category) %>%  summarize()
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
        HIT <- do.call(rbind, drug_data_sh$y$hitlist[input$drug2])
      }
      else if(input$drug == "dat"){
        if(is.null(hit_bar$summa)){
          idx <- grep("Summary", names(barhit_data()))

          HIT <- barhit_data()[idx][[1]]
        }
        else{
          HIT <- hit_bar$summa
        }
      }
      c_idx <- str_which(colnames(HIT), "^[C|c]ondition")
      if(!purrr::is_empty(c_idx)){
        HIT_summup <- list()
        for(i in unique(HIT[, c_idx])){
          HIT_summup[[i]] <- (HIT %>% dplyr::filter(Condition == i))$id
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
      c_idx <- str_which(colnames(Sel_cond_fhit()), "^[C|c]ondition")
      Sel_cond_fhit_SUMMA$hit <- Sel_cond_fhit()
      if(!purrr::is_empty(c_idx)){
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
            prcheck <- drug_data_sh$y$data[[input$drug2]]$id
          }
          else if(length(input$drug2) > 1){
            prcheck <- plyr::join_all(drug_data_sh$y$data[input$drug2], by = c("id", "description"), type = "full")$id
          }

          a <- pr[!(pr %in% prcheck)]
          if(!purrr::is_empty(a)){
            pr <- pr[(pr %in% prcheck)]
            showNotification(paste(paste(a, collapse = ", "), "wasn't in the data and had to be removed."),
                             type = "error")
          }
        }
      }
      else{
        if(length(input$drug2) == 1){
          if(input$hit & !is.null(input$cond_fhit)){
            pr <- Sel_cond_fhit_SUMMA$hit
            pr <- pr %>% dplyr::filter(!is.na(match(Condition, c(input$cond_fhit))))
            pr <- pr$id
            pr <- unique(pr)
          }
          else{
            pr <- drug_data_sh$y$data[[input$drug2]]$id
          }
        }
        else if(length(input$drug2) > 1){
          if(input$hit & !is.null(input$cond_fhit)){
            pr <- Sel_cond_fhit_SUMMA$hit
            pr <- pr %>% dplyr::filter(!is.na(match(Condition, c(input$cond_fhit))))
            pr <- pr$id
            pr <- unique(pr)
          }
          else{
            pr <- plyr::join_all(drug_data_sh$y$data[input$drug2], by = c("id", "description"), type = "full")$id
          }
        }
      }
    }

    else if(input$drug == "dat"){
      if(input$hit  & !is.null(input$cond_fhit)){
        pr <- Sel_cond_fhit_SUMMA$hit
        pr <- pr %>% dplyr::filter(!is.na(match(Condition, c(input$cond_fhit))))
        pr <- pr$id
        pr <- unique(pr)
      }
      else{
        pr <- barplot_data()$id
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

    if(input$drug == "base" & length(input$drug2) >= 1){
      HIT <- do.call(rbind, drug_data_sh$y$hitlist[input$drug2])
      NN <- do.call(rbind, drug_data_sh$y$NN[input$drug2])
    }
    else if(input$drug == "dat"){
      if(is.null(hit_bar$summa)){
        idx <- grep("Summary", names(barhit_data()))

        HIT <- barhit_data()[idx][[1]]
        NN <- barhit_data()[!(1:length(names(barhit_data())) %in% idx)][[1]]
      }
      else{
        HIT <- hit_bar$summa
        NN <- hit_bar$NN
      }
    }


    tr <- NULL
    if(input$cond_sel == "cat"){
      if(length(input$drug2) >= 1){
        trh <- HIT[which(!is.na(match(HIT$id, PROT))),c("Condition", "category")]
        tr <- NN[which(!is.na(match(NN$id, PROT))), c("Condition", "category")]
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
    }

    data <- ms_subsetting(data, isfile = F, hitidlist = c(PROT), allisoform = input$alliso_bar)



    if(input$cond_sel == "treat"){
      notsel_cond <- TREAT[!(TREAT %in% input$cond)]
      notsel_cond <- paste0("_", notsel_cond, "$")
      notsel_cond <- paste(notsel_cond, collapse = "|")

      if(str_length(notsel_cond) != 0){
        data <- data[,-str_which(names(data), notsel_cond)]
      }

      id_sel <- str_which(names(data), paste(input$cond, collapse = "|"))
      w <- 1:ncol(data)
      w <- w[!(w %in% id_sel)]

      ord <- unlist(lapply(input$cond, function(x) str_which(names(data), paste0("_", x, "$"))))

      data <- data[,c(w,ord)]
    }
    else if(input$cond_sel == "cat"){
      sele_cond <- Sel_cond()$Condition[which(!is.na(match(Sel_cond()$category, input$cond)))]
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
    DR <- NULL
    if(input$drug == "base" & !is.null(PROT)){
      if(length(input$drug2) > 1){
        DR <- DAT_text()[which(!is.na(match(DAT_text()$id, PROT))),]  #lost this column with ms_subsetting
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

      DR <-left_join(DR, hit_info, by = "id")
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
      paste0("Identification_comparison_", Sys.Date(), ".xlsx")
    },
    content = function(file){
      openxlsx::write.xlsx(tabident_bar$r, file, row.names = FALSE)
    }
  )

  output$n_cond_sel <- renderText({
    if(input$ch_own_col){
      if (input$cond_sel  == "all_cond"){
        paste("You selected", length(get_treat_level(data())), "conditions, please enter the same number of colors")
      }
      else{
        paste("You selected", length(input$cond), "conditions, please enter the same number of colors")
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
          IMPRINTS_barplotting_sh(data(), witherrorbar = input$werb,
                               usegradient = input$grad, linegraph = input$line,
                               save_pdf = input$save_bar, colorpanel = COL,
                               layout = c(input$lay_bar1, input$lay_bar2),
                               pdfname = input$pdftit)
        }
        else{
          showNotification("The number of colors given doesn't match the number of condition selected !", type = "error")
        }

      }
      else{
        IMPRINTS_barplotting_sh(data(), witherrorbar = input$werb,
                             usegradient = input$grad, linegraph = input$line,
                             save_pdf = input$save_bar,
                             layout = c(input$lay_bar1, input$lay_bar2),
                             pdfname = input$pdftit)
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
      showNotification("Don't forget to select a condition !", type = "error")

    }
    #else if(che){
    # showNotification("If you want to select all conditions,
    #                   please select the option 'Select all the treatment level'", type = "error")
    #}
    else{
      BAR$ch <- Bar_one()
    }

  })

  output$bar_plot <- renderPlot({
    BAR$ch
  })

  output$downbar <- downloadHandler(
    filename = function() {
      paste("2D_barplot_", Sys.Date(), "_", input$prot, ".png", sep = "")
    },
    content = function(file){
      ggsave(file, BAR$ch[[1]], device = "png")
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
      if(!purrr::is_empty(nv_nam)){
        dat <- dat[, !(names(dat) %in% nv_nam)]
      }
      dat
    }
    else if(input$drug_compl == "base" & length(input$drug2_compl) >= 1){
      do.call(rbind, drug_data_sh$y$hitlist[input$drug2_compl])
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
      File <- input$NN_compl
      if (is.null(File))
        return(NULL)

      dat <- import(File$datapath, header = TRUE)
      nv_nam <- str_subset(names(dat), "^V\\d{1}$")
      if(!purrr::is_empty(nv_nam)){
        dat <- dat[, !(names(dat) %in% nv_nam)]
      }
      dat
    }
    else if(input$drug_compl == "base" & length(input$drug2_compl) >= 1){
      do.call(rbind, drug_data_sh$y$NN[input$drug2_compl])
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
      updateSelectInput(session, "condsel_compl", choices = unique(HIT_compl()$Condition))
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
        data_ave <- IMPRINTS_average_sh(DIF_compl())
        showNotification("Average calculation succeed !", type = "message")
      }
      else{
        data_ave <- AVE_compl()
      }

      cat_tab <- HIT_compl()
      colnames(cat_tab)[str_which(colnames(cat_tab), "^[C|c]ondition")] <- "treatment"

      cat_tabNN <- NN_compl()
      colnames(cat_tabNN)[str_which(colnames(cat_tabNN), "^[C|c]ondition")] <- "treatment"
      cat_tabNN <- cat_tabNN %>% dplyr::group_by(id, treatment, category) %>% dplyr::summarise()

      cat_tab <- rbind(cat_tab, cat_tabNN)

      withCallingHandlers({
        shinyjs::html("diagmapping_compl", "")
        map_compl <- IMPRINTS_complex_mapping_sh(data_ave, cat_tab, treatment = input$condsel_compl,
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
        map_compl$description <- mineCETSAapp:::getProteinName(map_compl$description)
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
      paste0("ProteinComplexMapping_", Sys.Date(), ".xlsx")
    },
    content = function(file){
      openxlsx::write.xlsx(resmapping_compl$ch, file, row.names = FALSE)
    }
  )


  sel_prot_compl <- reactive({
    pr <- NULL
    if(!is.null(resmapping_compl$ch)){
      pr <- resmapping_compl$ch$id[which(!is.na(match(resmapping_compl$ch$ComplexName, input$allcomplex_compl)))]
      pr <- unique(pr)
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

        data_l[[i]] <- ms_subsetting(data, isfile = F, hitidlist = c(pr_comp), allisoform = input$alliso_bar_compl)

        data_l[[i]] <- data_l[[i]][,-str_which(names(data_l[[i]]), notsel_cond)]
        data_l[[i]]$category <- cate_$category[which(!is.na(match(cate_$id, data_l[[i]]$id)))]
      }

      data <- data_l

    }
    else{
      data <- ms_subsetting(data, isfile = F, hitidlist = c(PROT), allisoform = input$alliso_bar_compl)

      data <- data[,-str_which(names(data), notsel_cond)]
      data$category <- cate$category[which(!is.na(match(cate$id, data$id)))]
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


      IMPRINTS_barplotting_sh(data_compl(), witherrorbar = input$werb_compl,
                           usegradient = input$grad_compl, linegraph = input$line_compl,
                           save_pdf = input$save_bar_compl, colorpanel = COL,
                           layout = c(input$lay_bar1_compl, input$lay_bar2_compl),
                           toplabel = "IMPRINTS-CETSA bar plotting \nProtein complex :",
                           pdfname = input$pdftit_compl
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
      paste0("2D_barplot_", Sys.Date(), "_", paste(str_remove_all(input$allcomplex_compl, " "), sep = "_"), ".png")
    },
    content = function(file){
      ggsave(file, BAR_compl$ch[[1]], device = "png")
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
      updateSelectizeInput(session, "prot_simpf", choices = DIF_simpf()$id, server = TRUE)
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
      IMPRINTS_barplotting_simprof(DIF_simpf(), average, witherrorbar = input$werb_simpf,
                                treatmentlevel = input$treat_simpf, protein_profile = input$prot_simpf,
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
      geting_data_simpf$ch <- IMPRINTS_barplotting_simprof(DIF_simpf(), average, witherrorbar = input$werb_simpf,
                                                        treatmentlevel = input$treat_simpf, protein_profile = input$prot_simpf,
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
      BAR_simpf$ch <- IMPRINTS_barplotting_simprof(geting_data_simpf$ch, witherrorbar = input$werb_simpf,
                                                treatmentlevel = input$treat_simpf, protein_profile = input$prot_simpf,
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
      BAR_simpf$ch <- IMPRINTS_barplotting_simprof(geting_data_simpf$ch, witherrorbar = input$werb_simpf,
                                                treatmentlevel = input$treat_simpf, protein_profile = input$prot_simpf,
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
      paste0("2D_barplot_", Sys.Date(), "_", paste0("similar_", input$prot_simpf), ".png")
    },
    content = function(file){
      ggsave(file, BAR_simpf$ch[[1]], device = "png")
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
      if(!purrr::is_empty(nv_nam)){
        dat <- dat[, !(names(dat) %in% nv_nam)]
      }
      dat
    }
    else if(input$drug_heat == "base" & length(input$drug2_heat) >= 1){
      do.call(rbind, drug_data_sh$y$hitlist[input$drug2_heat])
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
      File <- input$NNfile_heat
      if (is.null(File) | !input$impNN_heat)
        return(NULL)

      dat <- import(File$datapath, header = TRUE)
      nv_nam <- str_subset(names(dat), "^V\\d{1}$")
      if(!purrr::is_empty(nv_nam)){
        dat <- dat[, !(names(dat) %in% nv_nam)]
      }
      dat
    }
    else if(input$drug_heat == "base" & length(input$drug2_heat) >= 1){
      do.call(rbind, drug_data_sh$y$NN[input$drug2_heat])
    }
    else{
      NULL
    }
  })
  #check if a file is upload
  output$NNheat_fileup <- reactive({
    if(input$drug_heat == "dat"){
      return(!is.null(NN_heat()) | !input$impNN_heat)
    }
    else{
      return(!is.null(NN_heat()))
    }
  })
  outputOptions(output, "NNheat_fileup", suspendWhenHidden = FALSE)

  observe({
    if(!is.null(HIT_heat())){
      c_idx <- str_which(colnames(HIT_heat()), "^[C|c]ondition")
      cat_idx <- str_which(colnames(HIT_heat()), "^[C|c]ategory")

      tr <- NULL
      if(!purrr::is_empty(c_idx)){
        tr <- HIT_heat()[, c_idx]
        tr <- unique(tr)
      }
      cat <- c()
      if(!purrr::is_empty(cat_idx)){
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
      dat <- IMPRINTS_average_sh(DIF_heat())
    }

    withCallingHandlers({
      shinyjs::html("diagl_heat", "")
      h <- IMPRINTS_heatmap(dat, HIT_heat(), NN_data = NN_heat(),
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
      data_ave <- IMPRINTS_average_sh(DIF_heatcom())
      resAVE_heatcom$d <- data_ave
      showNotification("Average calculation succeed !", type = "message")
    }
    else{
      data_ave <- AVE_heatcom()
    }

    withCallingHandlers({
      shinyjs::html("diagmapping_heatcom", "")
      map_heatcom <- IMPRINTS_complex_mapping_sh(data_ave, categorytable = NULL, treatment = input$cond_heatcom,
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
      map_heatcom$description <- mineCETSAapp:::getProteinName(map_heatcom$description)
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
                         Try to add more category in order to have more proteins", type = "error")
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
      paste0("ProteinComplexMapping_", Sys.Date(), ".xlsx")
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

    dat <- ms_subsetting(dat, isfile = F, hitidlist = c(pr), allisoform = FALSE)
    PRcompl <- resmapping_heatcom$ch[which(!is.na(match(resmapping_heatcom$ch$ComplexName, input$allcomplex_heatcom))),]

    withCallingHandlers({
      shinyjs::html("diagl_heatcom", "")
      h <- IMPRINTS_heatmap(dat, NULL, NN_data = NULL, PRcomplex_data = PRcompl,
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
      h <- do.call(rbind, drug_data_sh$y$hitlist[input$drug2_stri])
      n <- do.call(rbind, drug_data_sh$y$NN[input$drug2_stri])
      n <- unique(n[,c("id", "Condition", "category")])

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

      c_idx <- str_which(colnames(HIT), "^[C|c]ondition")

      if(!purrr::is_empty(c_idx)){
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
    showNotification("Getting the STRING id, this may take a while", type = "message")

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
              dat <- dat %>% dplyr::filter(Condition == c(input$cond_fhit_stri))
              a <- string_db$map(dat, "id", removeUnmappedRows = TRUE)
            }
            else{
              showNotification("Don't forget to select some conditions !", type = "error")
              a <- NULL
            }
          }
          else{
            a <- string_db$map(stri_data(), input$idfile_stri, removeUnmappedRows = TRUE)
          }
        }
        else{
          a <- string_db$map(stri_data(), "id", removeUnmappedRows = TRUE)
        }
      }
      else if(input$drug_stri == "base"){
        dat <- stri_data()
        if(!is.null(input$cond_fhitB_stri)){
          dat <- dat %>% dplyr::filter(Condition == c(input$cond_fhitB_stri))
          if(!is.null(input$cat_fhitB_stri)){
            dat <- dat %>% dplyr::filter(!is.na(match(category, c(input$cat_fhitB_stri))))
          }
          a <- string_db$map(dat, "id", removeUnmappedRows = TRUE)
        }
        else{
          showNotification("Don't forget to select some conditions !", type = "error")
          a <- NULL
        }
      }
    }

    string_res$x <- a

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
      My_net(string_res$x$STRING_id , inter = TRUE,
             network_flavor = input$edgetype1_stri, required_score = input$intscore1_stri)
    }
    else{
      My_net(string_res$x$STRING_id , inter = FALSE,
             network_flavor = input$edgetype1_stri, required_score = input$intscore1_stri)
    }
  })

  observeEvent(input$netbase_stri, {
    if(length(string_res$x$STRING_id) > 2000){
      showNotification(paste("Lists with more than 2000 genes are not supported yet. Your list contains now", length(string_res$x$STRING_id), "genes.",
                             "Please, try to reduce the size of your input by choosing less categories and/or conditions."),
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
      paste0("Network_", Sys.Date(), ".png")
    },
    content = function(file){
      ggsave(file, OUT_plot$g, device = "png")
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
        df <- Get_GO(enrich_res$x, TRUE, FALSE, input$catego_stri)
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
      paste0("Enrichment_tab_", Sys.Date(), "_", input$catego_stri, ".xlsx")
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

    if(!is.null(pr) & !purrr::is_empty(pr)){
      if(input$intnet_stri){
        My_net(pr , inter = TRUE,
               network_flavor = input$edgetype2_stri, required_score = input$intscore2_stri)
      }
      else{
        My_net(pr , inter = FALSE,
               network_flavor = input$edgetype2_stri, required_score = input$intscore2_stri)
      }
    }
    else{
      showNotification("No match has been found. Try another description or contact me vie the e-mail button.", type = "error")
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
      paste0("Network_", Sys.Date(), input$descri_stri, ".png")
    },
    content = function(file){
      ggsave(file, OUT_plot_filt$g, device = "png")
    }
  )



  ### Cluster Profiler
  output$drug2ui_clus <- renderUI({
    selectInput("drug2_clus", "Choose a drug", choices = names(drug_data_sh$y$data),
                multiple = TRUE, selected = "elutriation")
  })

  clus_data <- reactive({
    if(input$drug_clus == "dat"){
      File <- input$file_clus
      if (is.null(File))
        return(NULL)

      import_list(File$datapath, header = TRUE)[[1]]
    }
    else if(input$drug_clus == "base" & length(input$drug2_clus) >= 1){
      h <- do.call(rbind, drug_data_sh$y$hitlist[input$drug2_clus])
      h <- as.data.frame(h)
      if(!("Genes" %in% colnames(h))){
        if("description" %in% colnames(h)){
          h$Genes <- unname(unlist(sapply(h$description, mineCETSAapp:::getGeneName)))
        }
        else{
          d <- join_drugdata(drug_data_sh$y$data_ave[input$drug2_clus], by = c("id", "description")) ## extract gene information
          d <- d[,c("id", "description")]
          h <- dplyr::left_join(h, d, by = "id")
          h$Genes <- unname(unlist(sapply(h$description, mineCETSAapp:::getGeneName)))
        }
      }
      h <- h[,c("id", "Genes", "Condition", "category")]

      n <- do.call(rbind, drug_data_sh$y$NN[input$drug2_clus])
      n <- as.data.frame(n)
      if(!("Genes" %in% colnames(n))){
        if("description" %in% colnames(n)){
          n$Genes <- unname(unlist(sapply(n$description, mineCETSAapp:::getGeneName)))
        }
        else{
          d <- join_drugdata(drug_data_sh$y$data_ave[input$drug2_clus], by = c("id", "description")) ## extract gene information
          d <- d[,c("id", "description")]
          n <- dplyr::left_join(n, d, by = "id")
          n$Genes <- unname(unlist(sapply(n$description, mineCETSAapp:::getGeneName)))
        }
      }
      n <- unique(n[,c("id", "Genes", "Condition", "category")])

      rbind(h,n)
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

      c_idx <- str_which(colnames(HIT), "^[C|c]ondition")

      if(!purrr::is_empty(c_idx)){
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
            dat <- dat %>% dplyr::filter(Condition == c(input$cond_fhit_clus))
          }
          else{
            showNotification("Don't forget to select some conditions !", type = "error")
            dat <- NULL
          }
          if(!is.null(dat)){
            showNotification("Starting pathway analysis !", type = "message")
            res <- compare_wp(dat, gene_column = "Genes",
                              n_pathway = input$npath_clus, condition_column = "Condition")
          }
        }
        else{
          showNotification("Starting pathway analysis !", type = "message")
          res <- compare_wp(dat, gene_column = input$idfile_clus, n_pathway = input$npath_clus)
        }
      }
      else if(input$drug_clus == "base"){
        dat <- clus_data()
        if(!is.null(input$cond_fhitB_clus)){
          dat <- dat %>% dplyr::filter(Condition == c(input$cond_fhitB_clus))
          if(!is.null(input$cat_fhitB_clus)){
            dat <- dat %>% dplyr::filter(!is.na(match(category, c(input$cat_fhitB_clus))))
          }
        }
        else{
          showNotification("Don't forget to select some conditions !", type = "error")
          dat <- NULL
        }
        if(!is.null(dat)){
          showNotification("Starting pathway analysis !", type = "message")
          res <- compare_wp(dat, gene_column = "Genes",
                            n_pathway = input$npath_clus,
                            condition_column = "Condition")
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
      paste0("CompareClusterTab_", Sys.Date(), ".xlsx")
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
      paste0("CompareCluster_", Sys.Date(), ".png")
    },
    content = function(file){
      ggsave(file, cluscomp_res$x$graph, device = "png")
    }
  )

  ## GSEA
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
            dat <- dat %>% dplyr::filter(Condition == c(input$cond_fhit_clus))
          }
          else{
            showNotification("Don't forget to select some conditions !", type = "error")
            dat <- NULL
          }
          if(!is.null(dat)){
            if(input$scorename_clus %in% colnames(dat)){
              showNotification("Starting GSEA !", type = "message")
              res <- cetsa_gsea(dat, gene_column = "Genes",
                                score_column = input$scorename_clus,
                                pos_enrichment = input$onlypos_clus)
            }
            else{
              showNotification(paste(input$scorename_clus,
                                     "was not found in the column names of the data. Please check the name you wrote."),
                               type = "error")
            }
          }
        }
        else{
          if(input$scorename_clus %in% colnames(dat)){
            showNotification("Starting GSEA !", type = "message")
            res <- cetsa_gsea(dat, gene_column = input$idfile_clus, score_column = input$scorename_clus,
                              pos_enrichment = input$onlypos_clus)
          }
          else{
            showNotification(paste(input$scorename_clus,
                                   "was not found in the column names of the data. Please check the name you wrote."),
                             type = "error")
          }
        }
      }
      else if(input$drug_clus == "base"){
        dat <- clus_data()
        if(!is.null(input$cond_fhitB_clus)){
          dat <- dat %>% dplyr::filter(Condition == c(input$cond_fhitB_clus))
          if(!is.null(input$cat_fhitB_clus)){
            dat <- dat %>% dplyr::filter(!is.na(match(category, c(input$cat_fhitB_clus))))
          }
        }
        else{
          showNotification("Don't forget to select some conditions !", type = "error")
          dat <- NULL
        }
        if(!is.null(dat)){
          if(input$scorename_clus %in% colnames(dat)){
            showNotification("Starting GSEA !", type = "message")
            res <- cetsa_gsea(dat, gene_column = "Genes",
                              score_column = input$scorename_clus,
                              pos_enrichment = input$onlypos_clus)
          }
          else{
            showNotification(paste(input$scorename_clus,
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
      paste0("GSEATab_", Sys.Date(), ".xlsx")
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
      paste0("GSEAplot_", Sys.Date(), ".png")
    },
    content = function(file){
      ggsave(file, clusgsea_res$x$graph, device = "png")
    }
  )

  ## Gene concept network
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
            dat <- dat %>% dplyr::filter(Condition == c(input$cond_fhit_clus))
          }
          else{
            showNotification("Don't forget to select some conditions !", type = "error")
            dat <- NULL
          }
          if(!is.null(dat)){
            if(input$scorename2_clus %in% colnames(dat)){
              showNotification("Starting GSEA !", type = "message")
              res <- gene_conceptNet(dat, gene_column = "Genes",
                                     score_column = input$scorename2_clus,
                                     pval_cutoff = input$pvcut_clus)
            }
            else{
              showNotification(paste(input$scorename2_clus,
                                     "was not found in the column names of the data. Please check the name you wrote."),
                               type = "error")
            }
          }
        }
        else{
          if(input$scorename2_clus %in% colnames(dat)){
            showNotification("Starting GSEA !", type = "message")
            res <- gene_conceptNet(dat, gene_column = input$idfile_clus,
                                   score_column = input$scorename2_clus,
                                   pval_cutoff = input$pvcut_clus)
          }
          else{
            showNotification(paste(input$scorename2_clus,
                                   "was not found in the column names of the data. Please check the name you wrote."),
                             type = "error")
          }
        }
      }
      else if(input$drug_clus == "base"){
        dat <- clus_data()
        if(!is.null(input$cond_fhitB_clus)){
          dat <- dat %>% dplyr::filter(Condition == c(input$cond_fhitB_clus))
          if(!is.null(input$cat_fhitB_clus)){
            dat <- dat %>% dplyr::filter(!is.na(match(category, c(input$cat_fhitB_clus))))
          }
        }
        else{
          showNotification("Don't forget to select some conditions !", type = "error")
          dat <- NULL
        }
        if(!is.null(dat)){
          if(input$scorename2_clus %in% colnames(dat)){
            showNotification("Starting GSEA !", type = "message")
            res <- gene_conceptNet(dat, gene_column = "Genes",
                                   score_column = input$scorename2_clus,
                                   pval_cutoff = input$pvcut_clus)
          }
          else{
            showNotification(paste(input$scorename2_clus,
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
      paste0("GSEAplot_", Sys.Date(), ".png")
    },
    content = function(file){
      ggsave(file, clusgene_res$x, device = "png")
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

      import_list(File$datapath, header = TRUE)[[1]]
    }
    else if(input$drug_cell == "base" & length(input$drug2_cell) >= 1){
      h <- do.call(rbind, drug_data_sh$y$hitlist[input$drug2_cell])
      n <- do.call(rbind, drug_data_sh$y$NN[input$drug2_cell])
      n <- unique(n[,c("id", "Condition", "category")])

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
      updateSelectInput(session, "condhit_cell", choices = unique(hitdata_cell()$Condition), selected = unique(hitdata_cell()$Condition)[1])
      updateSelectInput(session, "cathit_cell", choices = unique(hitdata_cell()$category), selected = unique(hitdata_cell()$category)[1])
    }
  })


  resdata_cell <- reactiveValues(
    ch = NULL
  )
  observeEvent(input$goloca_cell, {
    showNotification("Getting subcellular locations", type = "message")

    data_hit <- hitdata_cell() %>%
      dplyr::filter(Condition == input$condhit_cell)

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
        paste0("HIT_locations_", Sys.Date(), ".xlsx")
      },
      content = function(file){
        openxlsx::write.xlsx(resdata_cell$ch, file, row.names = FALSE)
      }
    )

    updateSelectInput(session, "condp_cell", choices = unique(resdata_cell$ch$Condition), selected = input$condhit_cell)
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
      paste("Int_cell_", Sys.Date(), ".html", sep = "")
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
      if(purrr::is_empty(which(PR_event() == PR))){
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
        pr <- resdata_cell$ch$id[which(!is.na(match(resdata_cell$ch$main.location.cell, input$selorga_cell)))]
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
    }
    else{
      pr <- PR_event()
    }

    tr <- NULL
    if(!is.null(barpdata_cell())){
      if(input$cond_sel_cell == "cat" & !is.null(hitdata_cell())){
        tr <- hitdata_cell()[which(!is.na(match(hitdata_cell()$id,pr))),c("Condition", "category")]
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
        sele_cond <- Sel_cond_cell()$Condition[which(!is.na(match(Sel_cond_cell()$category, input$cond_cell)))]
        notsel_cond <- TREAT[!(TREAT %in% sele_cond)]
        if(!purrr::is_empty(notsel_cond)){
          notsel_cond <- paste(notsel_cond, collapse = "|")

          data <- data[,-str_which(names(data), notsel_cond)]
        }

      }
      else if(input$cond_sel_cell == "all_cond"){
        notsel_cond <- TREAT[!(TREAT %in% Sel_cond_cell())]
        if(!purrr::is_empty(notsel_cond)){
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

          data_l[[i]] <- ms_subsetting(data, isfile = F, hitidlist = c(pr_comp))
        }
        data <- data_l
      }
      else{
        data <- ms_subsetting(data, isfile = F, hitidlist = c(pr))
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
        paste("You selected", length(get_treat_level(data_cell())), "conditions, please enter the same number of colors")
      }
      else{
        paste("You selected", length(input$cond_cell), "conditions, please enter the same number of colors")
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
          IMPRINTS_barplotting_sh(data_cell(), witherrorbar = input$werb_cell,
                               usegradient = input$grad_cell, linegraph = input$line_cell,
                               save_pdf = input$save_bar_cell, colorpanel = COL,
                               layout = c(input$lay_bar1_cell, input$lay_bar2_cell),
                               toplabel = loca_cell_lab,
                               pdfname = input$pdftit_cell)
        }
        else{
          showNotification("The number of colors given doesn't match the number of condition selected !", type = "error")
        }

      }
      else{
        IMPRINTS_barplotting_sh(data_cell(), witherrorbar = input$werb_cell,
                             usegradient = input$grad_cell, linegraph = input$line_cell,
                             save_pdf = input$save_bar_cell,
                             layout = c(input$lay_bar1_cell, input$lay_bar2_cell),
                             toplabel = loca_cell_lab,
                             pdfname = input$pdftit_cell)
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
          showNotification("Don't forget to select a condition !", type = "error")
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
      paste("2D_barplot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file){
      ggsave(file, BAR_cell$ch[[1]], device = "png")
    }
  )

}


shinyApp(ui, server)
