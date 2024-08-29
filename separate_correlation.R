# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 12/10/19
#===========================================SEPARATED CORRELATION=============================================================================
output$sep_scatter_plot <- plotly::renderPlotly({

    source("ScatterPlot.R")

    URLString <- isolate(session$clientData$url_search)
    gene1 <- str_match(URLString, "gene1=(.+)gene2")[1,2]
    gene2 <- str_match(URLString, "gene2=(.+)cor_method")[1,2]
    cor_method <- str_match(URLString, "cor_method=(.+)cor_disease")[1,2]
    cor_disease <- str_match(URLString, "cor_disease=(.+)cor_db")[1,2]
    cor_db <- str_match(URLString, "cor_db=(.+)rm_tissue")[1,2]
    rm_tissue <- str_match(URLString, "rm_tissue=(.+)")[1,2]
    DB=chooseDB(cor_db)
    if(rm_tissue){
        DB$exp_db_name = paste0(DB$name, '_transformed_expr')
        shinyjs::hide("log_scale")
    }
    if(grepl("%20", cor_disease)){
        cor_disease = gsub("%20", " ", cor_disease)
    }
    run_single_cor(cor_method, cor_disease, gene1, gene2,
                   DB=DB, log_scale=input$log_scale, plot_loc=separate_cor_plot_loc,
                   data_loc=separate_cor_data_loc)
})
