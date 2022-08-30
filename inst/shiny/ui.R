library(shiny)
library(shinythemes)
library(bsplus)
library(waiter)
library(shinydashboard)
library(magrittr)
library(markdown)

library(imogimap)

pkg_name <- "imogimap"
pkg_version <- packageVersion(pkg_name) %>% as.character
site_name <- HTML(paste0("imogimap ", tags$sub(pkg_version)))

tmp_onco_genes <- paste("CD70", collapse="\n")
tmp_icp_genes <- paste(icp_gene_list, collapse="\n") #icp_gene_list

# navbarPage() will not work here with combination of waiter and bsplus
dashboardPage(
  dashboardHeader(title=site_name,
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
                                 img(src = 'md-anderson-horizontal-logo.jpeg')
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
                 position="right",
                 sidebarPanel(
                   includeMarkdown("introduction.md")
                 ),
                 mainPanel(
                   tabsetPanel(
                     tabPanel("Parameters",

                              h4("Function Parameters"),

                              textAreaInput("onco_genes", "Input Genes (One Gene Per Line)", tmp_onco_genes, width = "480px"),
                              textAreaInput("immune_genes", "Immune Checkpoint (ICP) Genes (One Gene Per Line; Curated List of 29 ICP Genes Provided Below as Default)", tmp_icp_genes, width = "480px"),
                              selectInput("immune_phenotype", "Immune Phenotype", TCGA_immune_features_list),
                              selectInput("cohort", "Cohort", tcgaTypes, selected="luad"),
                              selectInput("method", "Method", c( "max","independence")),


                              checkboxInput("sensitivity", "Sensitivity", FALSE) %>%
                                shinyInput_label_embed(
                                  shiny_iconlink() %>%
                                    bs_embed_popover(
                                      title = "Sensitivity",
                                      content = "measures the sensitivity of the observed synergy score to the exact data configuration. Check the box to conduct sensitivity analysis. Note that sensitivity analyses are time consuming. We recommend to leave the box unchecked for the first run, select onco-icp gene pair of interest based on synergy results, and finally conduct sensitivity analysis for the selected gene pair.",
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
                                      content = " measures the probability of finding a score higher than the observed synergy score, by randomly pairing a gene from genome to any of onco/icp genes. Check the box and choose the number of random iterations to conduct specificity analysis. Note that specificity analyses are time consuming. We recommend to leave the box unchecked for the first run, select onco-icp gene pair of interest based on synergy results, and finally conduct specificity analysis for the selected gene pair.",
                                      placement = "right"
                                    )
                                ),
                              conditionalPanel(
                                condition = "input.specificity == true",
                                numericInput("N_iteration_specificity","Number of iterations",
                                             1000,min = 100,step=100)
                              ),
                              h4("Plot Parameters"),

                              #textInput("plot_onco_gene", "Plot Gene", placeholder="HGNC Symbols (e.g., TP53)"),
                              #textInput("plot_icp_gene", "Plot Immune Checkpoint Gene", placeholder="HGNC Symbols (e.g., TP53)"),

                              actionButton("submit", label = "Submit"),

                              br(),
                              br(),

                              textOutput("debug_text")
                     ),
                     tabPanel("Results",
                              h4("Table"),
                              DT::dataTableOutput("results_table"),

                              h4("Network"),
                              numericInput("cutoff", "Cutoff", 0.35, min = 0, max = 1, step=0.1),
                              #textOutput("netplot_text"),
                              plotOutput('netplot'),
                              downloadButton('Network_plot'),
                              h4("Boxplot"),
                              uiOutput("boxplot_gene_pair"),
                              plotOutput('boxplot')

                     ),
                     tabPanel("About",
                              includeMarkdown("about.md")
                     )
                   )
                 )
               )
      )
    )

  )
)
