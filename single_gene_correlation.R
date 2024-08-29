# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 12/10/19

#===================================================================SINGLE GENE CORRELATION=========================================================================================
#distinguish the correlation type
observeEvent(input$function_type, {
    if(input$function_type == 'Top correlated genes'){
        output$second_or_top <- renderUI({
            sliderInput("num", "Number of observations:", min = 5, max = 200, value = 50)
        })
    }else{
        output$second_or_top <- renderUI({
            selectizeInput('gene2', 'Choose the second gene:', choices=NULL)
        })
        updateSelectizeInput(session, 'gene2', server = TRUE, choices=GENE_LIST, selected='')
    }
})

#detect if example gene is clicked
observeEvent(input$single_correlation_example_gene,{
    updateSelectizeInput(session, 'gene1', server = TRUE, choices=GENE_LIST, selected='BRCA1')
    updateSelectizeInput(session, 'gene2', server = TRUE, choices=GENE_LIST, selected='BRCA2')
})

#function to update options for selecting the first gene
update_gene1 <- reactive({
    updateSelectizeInput(session, 'gene1', server = TRUE, choices=GENE_LIST, selected='')
})


update_cor_disease <- reactive({
    DB = chooseDB(input$dataset)
    db_table <- QueryAll(DB$clinical_db_name)
    updateSelectizeInput(session, 'single_corr_diseases', server=TRUE, choices=c('pan-cancer', unique(db_table[,DB$disease_tag])), selected='pan-cancer')
})
output$cor_out <- renderUI({

    library(tibble)
    source("run_single_cor.r")

    update_gene1()
    update_cor_disease()
    if(input$function_type=="correlate two genes"){
        plotly::plotlyOutput("two_genes", height="100%")
    }else{
        # column(12,
        #     dataTableOutput("one_gene_top")
        # )
        DT::dataTableOutput("one_gene_top")
    }
})
#######################################two genes correlation functions#######################################
scatter_out <- reactive({

    source('ScatterPlot.R')

    DB <- chooseDB(input$dataset)

    #if cohort is not pancan, then hide option for removing tissue effect
    if(input$single_corr_diseases != "pan-cancer"){
        shinyjs::hide("z_score_sub")
    }

    if(input$z_score_sub){
        DB$exp_db_name = paste0(DB$name, '_transformed_expr')
    }
    validate(
        need(input$gene1,"please input gene1"),
        need(input$gene2,"please input gene2"),
        errorClass="ValidateInput"
    )

    validate(
        need(input$gene1!=input$gene2, "two input genes are same"),
        errorClass="ErrorOutput"
    )


    final_plot = run_single_cor(input$cor_method, input$single_corr_diseases,
                   input$gene1, input$gene2, DB=DB, plot_loc=single_gene_cor_plot_loc,
                   data_loc=single_gene_cor_data_loc
    )
    shinyjs::show("CorPlotDown")
    shinyjs::show("CorDataDown")
    
    return(final_plot)
})

#######################################one gene top n functions##############################################################################

#functions that creates links for
createNCBILink <- function(val) {
    sprintf('<a href="https://www.ncbi.nlm.nih.gov/gene/?term=%s" target="_blank" class="btn btn-success" style="font-size:0.8rem">Search NCBI</a>',val)
}
createCorrLink <- function(gene1, gene2,cor_method, cor_disease, cor_db, remove_tissue_effect) {
    sprintf('<a href="?separate_correlation?gene1=%sgene2=%scor_method=%scor_disease=%scor_db=%srm_tissue=%s" target="_blank" class="btn btn-primary" style="font-size:0.8rem">Show</a>', gene1, gene2, cor_method, cor_disease, cor_db, remove_tissue_effect)
}
createGeneCardLink <- function(val) {
    sprintf('<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s" target="_blank" class="btn btn-success" style="font-size:0.8rem">Search GeneCards</a>',val)
}
get_cor_data <- reactive({
    #get selected DB object
    DB <- chooseDB(input$dataset)

    #get all clinical data for sub disease selection
    clin = QueryAll(DB$clinical_db_name)[,c(DB$biospec_id_name, DB$disease_tag)]

    #get the whole normalized expression matrix
    setProgress(0.2, detail='loading data')
    if(input$z_score_sub){
        exp_data = readRDS(paste("./IntermediateFiles/", DB$name, "_case_row_expr_transformed.rds", sep=''))
    }else{
        exp_data = readRDS(paste("./IntermediateFiles/", DB$name, "_case_row_expr.rds", sep=''))
    }

    #join the expression data with clinical data for cohort selection
    setProgress(0.4, detail='fetching disease')
    # exp_data = inner_join(exp_data, clin, by=DB$biospec_id_name)
    exp_data = merge(exp_data, clin, by=DB$biospec_id_name)
    #return data objects
    return(list(ds=exp_data, db=DB))
})
get_cor_data_subset <- reactive({

    #fetch data objects
    obj <- get_cor_data()
    dataset <- obj$ds
    DB <- obj$db

    #get sub disease matrix when the user seleted disease is not pan-can
    if(input$single_corr_diseases != "pan-cancer"){
        dataset <- dataset[dataset[ ,DB$disease_tag] == input$single_corr_diseases, ]
    }

    validate(
        need(length(rownames(dataset))>=5, paste('no enough data for computing correlation')),
        errorClass="ErrorOutput"
    )

    return(list(ds=dataset, DB=DB))
})

get_correlated_genes <- reactive({

    #check if the user has input first gene
    validate(
        need(input$gene1,'please input first gene'),
        errorClass="ValidateInput"
    )

    #get the dataset for computing correlation
    obj = get_cor_data_subset()
    dataset = obj$ds
    DB = obj$DB

    # dataset <- data.frame(t(dataset))
    gene1 = single_gene_special_transfer(input$gene1)

    # check if gene1 is in current database
    validate(
        need(gene1 %in% names(dataset), paste0(gene1, ' was not found in this DB')),
        errorClass="ErrorOutput"
    )

    #get cor_matrix
    cor_matrix = dataset[ , !(names(dataset) %in% c(DB$biospec_id_name, DB$disease_tag))]

    op <- options(warn = (-1))

    setProgress(0.5, detail='get correlation R and p-values. This may take about 1-2mins, please be patient')
    #get p values and R for all correlated genes:
    p_and_r <- apply(cor_matrix, 2, function(x){
        tryCatch({
        cor_obj = cor.test(x, dataset[[gene1]], method = input$cor_method, na.action = NA)
        p = cor_obj$p.value
        r = cor_obj$estimate
        return(c(p,r)) },
        error = function(e) {
            return(c(NA,NA))
        })
    })

    #get all adjusted p values
    setProgress(0.6, detail='get all adjusted p-values')
    all_adjust_p = p.adjust(p_and_r[1,])

    df = data.frame('cor_value'=p_and_r[2,], 'p_adjust'=all_adjust_p)

    options(op)

    return(df)
})

get_top_n_correlated_genes <- reactive({
    shiny::withProgress(message = 'Calculating top genes', {
        df <- get_correlated_genes()

        if(input$sort_option=='Positive'){
            #get the top N correlated genes and sort them
            df <- df[order(df$cor_value, decreasing=T), ][1:input$num+1,]
        }else if(input$sort_option=='Negative'){
            #get the top N correlated genes and sort them
            df <- df[order(df$cor_value, decreasing=F), ][1:input$num+1,]
        }else{
            #get the top N correlated genes and sort them
            df <- df[order(abs(df$cor_value), decreasing=T), ][1:input$num+1,]
        }

        #construct these genes into data frame
        setProgress(0.9, detail='preparing results')
        # df <- data.frame(cor_genes[ndx])
        df <- data.frame("GeneName"=rownames(df),"Cor_cof"=df[,1], "Adjust_p"=df[,2])

        #curate the gene name and correlation value to be displayed in the data table
        df$Cor_cof <- sapply(df$Cor_cof , round, digits = 2)
        df$Adjust_p <- sapply(df$Adjust_p, formatC, format = "e", digits = 3)
        df$GeneName <- sapply(df$GeneName, as.character)
        df$GeneName <- sapply(df$GeneName, single_gene_special_r_transfer)

        #exclude the selected gene from the result table
        df = df[which(df$GeneName!=input$gene1),]

        #add ranking
        Rank = seq(1:input$num)
        df = tibble::add_column(df, Rank, .before = 1)

        setProgress(1, detail='done')
    })
    #submit genes to enrich_r
    enrich_r_genes <- paste(df$GeneName, collapse='&')

    #after get correlation table, provide button to submit highly correlated genes to Enrich-r
    output$enrich_r_button <- renderUI({
        tagList(
            downloadButton("correlation_table_download", "Download Data", class="btn btn-success d-block mb-4"),
            tags$a(class="btn btn-primary w-100",target="_blank",
                href=sprintf("?separate_enrich_r?genes=%s", enrich_r_genes), "Pathway Enrichment"
            )
        )
    })

    #write the csv file out for downloading option
    write.csv(df, correlation_data_loc, row.names=F)

    #provide links for sub functions
    df$NCBI = createNCBILink(df$GeneName)
    df$GeneCards = createGeneCardLink(df$GeneName)
    df$Plot = createCorrLink(input$gene1, df$GeneName, input$cor_method, input$single_corr_diseases, input$dataset, input$z_score_sub)

    #construct the data table to be shown
    dt <- DT::datatable(df, rownames=FALSE, escape = FALSE,
    options = list(pageLength = 20,
    lengthMenu = list(c(20, 30, 40), c('20', '30', '40')),
    columnDefs = list(list(className = 'dt-center', targets = 0: length(names(df))-1))
    ))
    return(dt)
})
############################################################################################################################################################
output$one_gene_top <- DT::renderDataTable({
    shinyjs::show("enrich_r_button")
    shinyjs::hide("z_score_sub")
    shinyjs::hide("CorPlotDown")
    shinyjs::hide("CorDataDown")


    if(input$single_corr_diseases == "pan-cancer"){
        shinyjs::show("z_score_sub")
    }else{
        updateCheckboxInput(session, "z_score_sub", value = FALSE)
    }
    get_top_n_correlated_genes()
})
output$two_genes <- plotly::renderPlotly({
    #if two genes correlation, hide enrich R button and Sorting options
    shinyjs::hide("enrich_r_button")
    shinyjs::hide("sort_option")
    shinyjs::hide("z_score_sub")

    #if cohort is pan-cancer, show the removing tissue effect option
    if(input$single_corr_diseases == "pan-cancer"){
        shinyjs::show("z_score_sub")
    }else{
        updateCheckboxInput(session, "z_score_sub", value = FALSE)
    }
    scatter_out()
})
