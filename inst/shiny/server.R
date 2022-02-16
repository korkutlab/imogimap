library(shiny)

library(imogimap)

function(input, output, session) {
  
  onco_genes <- reactive({
    unlist(strsplit(input$onco_genes,"\n"))
  })
  
  immune_genes <- reactive({
    unlist(strsplit(input$immune_genes,"\n"))
  })

  my_syng_df <- eventReactive(input$submit, {
    tmp <-  im_syng_tcga(
      onco_gene  = onco_genes(),
      icp_gene = immune_genes(),
      cohort = input$cohort,
      method = input$method,
      sensitivity = input$sensitivity,
      specificity = input$specificity
    )
    
    cat("DEBUG: im_syng_tcga finished")
    
    tmp
  })
  
  output$netplot <- renderPlot({
    im_netplot(df = my_syng_df ,
               Immune_phenotype = input$immune_phenotype,
               cutoff = input$cutoff,
               cohort = input$cohort ,
               icp_gene = immune_genes(),
               seed = seed)
    
    cat("DEBUG: im_netplot finished")
  })
  
  output$boxplot <- renderPlot({
    im_boxplot_tcga(onco_gene = input$plot_onco_gene, 
                    icp_gene = input$plot_icp_gene,
                    cohort = input$cohort, 
                    Immune_phenotype = input$immune_phenotype)
    
    cat("DEBUG: im_boxplot finished")
  })
}