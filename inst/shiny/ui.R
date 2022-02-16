library(shiny)
library(shinythemes)
library(magrittr)
library(bsplus)
library(markdown)

library(imogimap)

pkg_name <- "imogimap"
pkg_version <- packageVersion(pkg_name) %>% as.character
site_name <- HTML(paste0("imogimap ", tags$sub(pkg_version)))

navbarPage(site_name,
           theme = shinytheme("yeti"),
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
                                   
                                   textAreaInput("onco_genes", "Input Genes (One Gene Per Line)", "TGFBR1", width = "480px"),
                                   textAreaInput("immune_genes", "Immune Checkpoint Genes (One Gene Per Line)", "TGFBR2", width = "480px"),
                                   
                                   selectInput("cohort", "Cohort:", tcgaTypes),
                                   selectInput("method", "Method:", c("independence", "max")),
                                   
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
                                   
                                   selectInput("immune_phenotype", "Immune Phenotype:", c("IFNGscore", "EMTscore")),
                                   numericInput("cutoff", "Cutoff:", 0.35, min = 0, max = 1, step=0.1),
                                   
                                   textInput("plot_onco_gene", "Plot Gene", placeholder="HGNC Symbols (e.g., TP53)"),
                                   textInput("plot_icp_gene", "Plot Immune Checkpoint Gene", placeholder="HGNC Symbols (e.g., TP53)"),
                                   
                                   actionButton("submit", label = "Submit", icon = NULL, width = NULL)
                          ),
                          tabPanel("Results", 
                                   h4("Network"),
                                   plotOutput('netplot'),
                                   
                                   h4("Boxplot"),
                                   plotOutput('boxplot')
                                   
                          )
                      )
                  )
                    )
           ),
           tabPanel("About",
                    includeMarkdown("about.md")
           )
)
