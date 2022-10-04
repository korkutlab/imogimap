library(shiny)
library(waiter)
library(shinydashboard)
library(imogimap)
library(base)
library(DT)
library(R.utils)

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

    waiter_show(html=tagList(spin_fading_circles(), h4("Running ...")),color = "#58593FFF")

    cat("onco: ", paste(onco_genes(), collapse=","), "\n")
    cat("icp: ", paste(immune_genes(), collapse=","), "\n")
    cat("cohort: ", input$cohort, "\n")
    cat("iap: ", input$immune_phenotype, "\n")
    cat("method: ", input$method, "\n")

    withProgress(message = '', value = 0, {
      tmp <-  withTimeout(im_syng_tcga_shiny(
        onco_gene  = onco_genes(),
        icp_gene = immune_genes(),
        cohort = input$cohort,
        add_receptor_ligand = input$add_receptor_ligand,
        select_iap = input$immune_phenotype,
        method = input$method,
        ndatamin=4,
        sensitivity = input$sensitivity,
        specificity = input$specificity,
        N_iteration_sensitivity = input$N_iteration_sensitivity,
        N_iteration_specificity = input$N_iteration_specificity
      ),timout=3600,cpu = 3600,elapsed = 3600,onTimeout = "warning")
    })
    if(!exists("tmp")){ shinyalert("Oops! Something went wrong. Please contact us.",type="error")}
    if(class(tmp)=="character"){ shinyalert("Reached elapsed time limit.  See instruction for more details.",type="error")}
    cat("COLS: ", paste(colnames(tmp), collapse=","), "\n")

    tmp <- tmp[order(-tmp$Synergy_score),]
    tmp2 <- tmp[complete.cases(tmp$Synergy_score),]
    tmp2 <- tmp2[tmp2$Synergy_score!=0,]
    cat("DEBUG: im_syng_tcga finished\n")
    if(nrow(tmp2)==0){
      shinyalert("All combined action scores are zero or missing. Try a different immune phenotype.",type="warning")
    }
    waiter_hide()
    tmp
  })

  rowCallback <- c(
    "function(row, data){",
    "  for(var i=0; i<data.length; i++){",
    "    if(data[i] === null){",
    "      $('td:eq('+i+')', row).html('NA')",
    "        .css({'color': 'rgb(151,151,151)', 'font-style': 'italic'});",
    "    }",
    "  }",
    "}"
  )
  output$results_table <- DT::renderDataTable({
    req(my_syng_df())
    tmp <- my_syng_df() %>% round_df(., 3)
    colnames(tmp) <- c("Gene1", "Gene2","Disease" ,"IAP", "combined action",
                       "Gene1Expr", "Gene2Expr", "WilcoxPvalue",
                       "SpecificityPvalue","SensitivityR","receptor_ligand")

    #cat("COLS: ", paste(colnames(tmp), collapse=","), "\n")
    options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))
    DT::datatable(tmp, rownames= FALSE,
                  options = list(rowCallback = DT::JS(rowCallback))
    )
  })

  output$download_table <- downloadHandler(
    filename = function(){"im_table.csv"},
    content = function(fname){
      write.csv(my_syng_df(), fname)
    }
  )
  boxplot_gene_pairs <- reactive({
    req(my_syng_df())
    t2 <- my_syng_df()
    #t2 <- subset(t1,!is.na(Synergy_score))
    #t2 <- subset(t2, !(Synergy_score==0))
    t2 <- t2[, c("Gene1", "Gene2")] %>% unique
    t3 <- paste(t2[,1], t2[,2], sep="|")
    t3
  })


  plotinput1 =function(){
    req(my_syng_df())
    tmp_results <- my_syng_df()
    tmp_results <- subset(tmp_results,!is.na(Synergy_score))
    tmp_results <- subset(tmp_results, !(Synergy_score==0))
    if(nrow(tmp_results)>0){
      icp_gene <- immune_genes()
      Immune_phenotype <- input$immune_phenotype
      cutoff <- input$cutoff
      cohort <- input$cohort

      im_netplot(df = tmp_results,
                 Immune_phenotype = Immune_phenotype,
                 cutoff = cutoff,
                 cohort = cohort,
                 icp_gene = icp_gene,
                 seed = 123)
    }else{
      waiter_hide()
      cat("DEBUG: no netplot to show\n")
    }
  }
  output$netplot <- renderPlot({
    cat("DEBUG: im_netplot started\n")
    plotinput1()
  })
  output$Network_plot <- downloadHandler(
    filename = 'im_network.png',
    content = function(file) {
      png(file,width = 1000,height = 1000)
      print(plotinput1())
      dev.off()
    }
  )


  output$boxplot_gene_pair <- renderUI({
    req(boxplot_gene_pairs())
    selectInput("gene_pair", "Gene Pairs", boxplot_gene_pairs())
  })

  plotinput2 = function(){


  }
  output$boxplot <- renderPlot({
    cat("DEBUG: im_boxplot started\n")

    waiter_show(html=tagList(spin_fading_circles(), h4("Running ...")))

    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Creating boxplot...", value = 99)

    gene_pair <- input$gene_pair
    if(length(gene_pair) > 0 ){
      str(gene_pair)

      onco_gene <- strsplit(gene_pair, "\\|")[[1]][1]
      icp_gene <- strsplit(gene_pair, "\\|")[[1]][2]

      obj <- im_boxplot_tcga(onco_gene = onco_gene,
                             icp_gene = icp_gene,
                             cohort = input$cohort,
                             Immune_phenotype = input$immune_phenotype)

      waiter_hide()

      im_boxplot_tcga_plot(obj)
    }else{
      waiter_hide()
      cat("DEBUG: No boxplot to show \n")
    }
  })

  output$box_plot <- downloadHandler(
    filename = 'boxplot.png',
    content = function(file) {
      png(file)
      print(im_boxplot_tcga_plot(obj))
      dev.off()
    }
  )
  observeEvent(input$submit, {
    updateTabsetPanel(session, "tabselected",
                      selected = "Result"
    )
  })

}
