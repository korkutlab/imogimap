library(shiny)
library(shinythemes)
library(bsplus)
library(waiter)
library(shinydashboard)
library(magrittr)
library(markdown)
library(shinyalert)
library(imogimap)

pkg_name <- "imogimap"
pkg_version <- packageVersion(pkg_name) %>% as.character
site_name <- HTML(paste0("imogimap ", tags$sub(pkg_version)))

tmp_onco_genes <- paste(c("TGFB1","TGFBR1"), collapse="\n")
tmp_icp_genes <- paste(icp_gene_list, collapse="\n") #icp_gene_list

# navbarPage() will not work here with combination of waiter and bsplus
dashboardPage(
  dashboardHeader(title=site_name,
                  tags$li(class = "dropdown",
                          tags$a("Github" ,
                                 href = 'https://github.com/korkutlab/imogimap',
                                 #style = "padding: 10px 100px 0px 0px",

                          )
                  ),

                  tags$li(class = "dropdown",
                          tags$a("Publication" ,
                                 href = 'https://doi.org/10.1101/2021.10.06.462889',
                                 #style = "padding: 10px 100px 0px 0px",

                          )
                  ),
                  tags$li(class = "dropdown",
                          tags$a("Korkutlab" ,
                                 href = 'https://odin.mdacc.tmc.edu/~akorkut/#/home',
                                 #style = "padding: 10px 1190px 0px 0px"
                          )
                  ),

                  tags$li(class = "dropdown",
                          tags$a("MDAnderson" ,
                                 href = 'https://www.mdanderson.org/',
                                 #style = "padding: 10px 1190px 0px 0px",
                                 img(src = '__rectSitelogo__MDA SITE LOGO - Thumbnail.png',height="20%", width="20%", align="right")
                          )
                  )
  ),
  dashboardSidebar(disable = TRUE),

  dashboardBody(
    tags$style(HTML('.popover{width:4000px;height:250px;}
                               .main-sidebar {z-index:auto;}')),
    fluidPage(#site_name,
      #titlePanel(site_name),
      theme = shinytheme("yeti"),
      useWaiter(),
      use_bs_tooltip(),
      use_bs_popover(),
      tabPanel("Analysis",
               sidebarLayout(
                 position="left",
                 sidebarPanel(
                   conditionalPanel(condition="input.tabselected==1",
                                    includeMarkdown("introduction.md")),

                   conditionalPanel(condition="input.tabselected==2",
                                    includeMarkdown("interpretation.md")
                   )
                 ),
                 mainPanel(
                   tabsetPanel(
                     tabPanel("Parameters",value=1,

                              h4("Function Parameters"),

                              textAreaInput("onco_genes", "Input Genes Hugo ID (One Gene Per Line; Case sensitive)", tmp_onco_genes, width = "480px"),
                              textAreaInput("immune_genes", "Immune Checkpoint (ICP) Genes Hugo ID (One Gene Per Line; Case sensitive; Curated List of 29 ICP Genes Provided Below as Default)", tmp_icp_genes, width = "480px"),

                              checkboxInput("add_receptor_ligand", "add receptor/ligand", FALSE) %>%
                                shinyInput_label_embed(
                                  shiny_iconlink() %>%
                                    bs_embed_popover(
                                      title = "add receptor/ligand",
                                      content = "Check to add genes with known ligand/receptor interactions from CellphonDB database.",
                                      placement = "right"
                                    )
                                ),
                              conditionalPanel(
                                condition = "input.add_receptor_ligand == true",
                              ),

                              selectInput("immune_phenotype", "Immune Phenotype", TCGA_immune_features_list),
                              selectInput("cohort", "Cohort", tcgaTypes, selected="luad"),
                              selectInput("method", "Method", c( "max","independence")),


                              checkboxInput("sensitivity", "Sensitivity", FALSE) %>%
                                shinyInput_label_embed(
                                  shiny_iconlink() %>%
                                    bs_embed_popover(
                                      title = "Sensitivity",
                                      content = "measures the sensitivity of the observed synergy score to the exact data configuration. Check the box to conduct sensitivity analysis. Note that sensitivity analyses are time consuming. See sthe ide panel for more details.",
                                      placement = "right"
                                    )
                                ),
                              conditionalPanel(
                                condition = "input.sensitivity == true",
                                numericInput("N_iteration_sensitivity","Number of iterations",
                                             100,min = 10,step=10)
                              ),

                              checkboxInput("specificity", "Specificity", FALSE) %>%
                                shinyInput_label_embed(
                                  shiny_iconlink() %>%
                                    bs_embed_popover(
                                      title = "Specificity",
                                      content = " measures the probability of finding a score higher than the observed synergy score, by randomly pairing a gene from genome to any of  genes in the lists. Check the box and choose the number of random iterations to conduct specificity analysis. Note that specificity analyses are time consuming. See the side panel for more details",
                                      placement = "right"
                                    )
                                ),
                              conditionalPanel(
                                condition = "input.specificity == true",
                                numericInput("N_iteration_specificity","Number of iterations",
                                             100,min = 10,step=10)
                              ),
                              #h4("Plot Parameters"),

                              #textInput("plot_onco_gene", "Plot Gene", placeholder="HGNC Symbols (e.g., TP53)"),
                              #textInput("plot_icp_gene", "Plot Immune Checkpoint Gene", placeholder="HGNC Symbols (e.g., TP53)"),

                              actionButton("submit", label = "Submit"),

                              br(),
                              br(),

                              textOutput("debug_text")
                     ),
                     tabPanel("Results",value=2,
                              h4("Table"),
                              div(style = 'overflow-x: scroll',DT::dataTableOutput("results_table")),
                              downloadButton('download_table'),
                              h4("Network"),
                              numericInput("cutoff", "Cutoff", 0.80, min = 0, max = 1, step=0.1) %>%
                                shinyInput_label_embed(
                                  shiny_iconlink() %>%
                                    bs_embed_popover(
                                      title = "Cutoff",
                                      content = "Choose a percentile value, k, for absolute synergy scores to exclude scores with their absolute value below the k-th percentile.",
                                      placement = "right"
                                    )
                                ),
                              #textOutput("netplot_text"),
                              plotOutput('netplot', width = 1000,
                                         height = 1000),
                              downloadButton('Network_plot'),
                              h4("Boxplot")
                              ,
                              uiOutput("boxplot_gene_pair"),
                              plotOutput('boxplot')
                     ),
                     tabPanel("About",value=3,
                              includeMarkdown("about.md")
                     ),
                     id = "tabselected"
                   )
                 )
               )
      )

    )
  )
)

