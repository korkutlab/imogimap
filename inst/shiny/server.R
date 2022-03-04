library(shiny)
library(waiter)

library(imogimap)

function(input, output, session) {
  onco_genes <- reactive({
    unlist(strsplit(input$onco_genes, "\n"))
  })
  
  immune_genes <- reactive({
    unlist(strsplit(input$immune_genes, "\n"))
  })

  debug_reactive <- eventReactive(input$submit, { 
    cat("DEBUG: WORKS\n")
  })
  output$debug_text <- renderText({ debug_reactive(); "Click 'Result' Tab to Continue ..." })
  
  output$netplot_text <- renderText({ input$cutoff })
  
  my_syng_df <- eventReactive(input$submit, {
    cat("DEBUG: im_syng_tcga started\n")
    
    waiter_show(html=tagList(spin_fading_circles(), h4("Running ...")))
    
    tmp <-  im_syng_tcga(
      onco_gene  = onco_genes(),
      icp_gene = immune_genes(),
      cohort = input$cohort,
      method = input$method,
      sensitivity = input$sensitivity,
      specificity = input$specificity
    )
    #tmp <- "BAD"
    cat("COLS: ", paste(colnames(tmp), collapse=","), "\n")
    tmp <- tmp[order(-tmp$Synergy_score),]
    
    cat("DEBUG: im_syng_tcga finished\n")
    
    tmp
  })
  
  output$results_table <- DT::renderDataTable({
    tmp <- my_syng_df() %>% round_df(., 2)
    tmp <- tmp[2:(length(tmp)-2)]
    colnames(tmp) <- c("Gene1", "Gene2", "IAP", "Synergy", "Gene1Expr", "Gene2Expr", "WilcoxPvalue")
    
    cat("COLS: ", paste(colnames(tmp), collapse=","), "\n")
    
    DT::datatable(tmp, rownames= FALSE)
  })
  
  boxplot_gene_pairs <- reactive({
    t1 <- my_syng_df()
    
    t2 <- t1[, c("Gene1", "Gene2")] %>% unique 
    t3 <- paste(t2[,1], t2[,2], sep="|")
    
    t3
  })

  output$netplot <- renderPlot({
    cat("DEBUG: im_netplot started\n")
    
    tmp_results <- my_syng_df()
    icp_gene <- immune_genes()
    Immune_phenotype <- input$immune_phenotype
    cutoff <- input$cutoff
    cohort <- input$cohort
    
    # inputs <- list(tmp_results=tmp_results,
    #                Immune_phenotype=Immune_phenotype,
    #                cutoff=cutoff,
    #                cohort=cohort,
    #                icp_gene=icp_gene)
    # 
    # str(inputs)

    im_netplot(df = tmp_results,
               Immune_phenotype = Immune_phenotype,
               cutoff = cutoff,
               cohort = cohort,
               icp_gene = icp_gene,
               seed = 123)
  })

  output$boxplot_gene_pair <- renderUI({
    selectInput("gene_pair", "Gene Pairs", boxplot_gene_pairs())
  })

  output$boxplot <- renderPlot({
    cat("DEBUG: im_boxplot started\n")
    
    waiter_show(html=tagList(spin_fading_circles(), h4("Running ...")))

    #gene_pair <- boxplot_gene_pairs()
    gene_pair <- input$gene_pair
    
    str(gene_pair)
    
    onco_gene <- strsplit(gene_pair, "\\|")[[1]][1]
    icp_gene <- strsplit(gene_pair, "\\|")[[1]][2]

    obj <- im_boxplot_tcga(onco_gene = onco_gene,
                           icp_gene = icp_gene,
                           cohort = input$cohort,
                           Immune_phenotype = input$immune_phenotype)
    
    waiter_hide()
    
    im_boxplot_tcga_plot(obj)
    })
}
