library(shiny)
library(shinythemes)
library(bsplus)
library(waiter)

library(magrittr)
library(markdown)

library(imogimap)

pkg_name <- "imogimap"
pkg_version <- packageVersion(pkg_name) %>% as.character
site_name <- HTML(paste0("imogimap ", tags$sub(pkg_version)))

tmp_onco_genes <- paste("TGFB1", collapse="\n")
tmp_icp_genes <- paste(icp_gene_list, collapse="\n") #icp_gene_list

# navbarPage() will not work here with combination of waiter and bsplus
fluidPage(#site_name,
          titlePanel(site_name),
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
                                           title = "Sensitivity", content = "Conduct sensitivity analysis?", placement = "right"
                                         )
                                     ),
                                   
                                   checkboxInput("specificity", "Specificity", FALSE) %>%
                                     shinyInput_label_embed(
                                       shiny_iconlink() %>%
                                         bs_embed_popover(
                                           title = "Specificity", content = "Conduct specificity analysis?", placement = "right"
                                         )
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
