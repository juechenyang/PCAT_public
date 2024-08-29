# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 12/9/19
#=========================================SINGLE GENE Expression===============================================================================
#update single gene analysis disease
update_single_disease <- reactive({
    DB <- chooseDB(input$db_single_name)
    clin <- QueryAll(DB$clinical_db_name)
    diseases = unique(na.omit(clin[,DB$disease_tag]))
    diseases = remove_invalid(diseases)
    updateSelectInput(session, 'single_exp_disease', choices=c('pan-cancer', diseases))
    return(list(DB=DB, clinical_table=clin))

})

#update single gene analysis genes
update_single_gene <- reactive({
    updateSelectizeInput(session, 'geneName', server = TRUE, choices=GENE_LIST, selected='')
})

#update single gene analysis clinical parameters
update_cli_bio <- reactive({
    if(input$db_single_name=='target'){
        updateSelectInput(session, 'clin_or_bio', choices=c('clinical', 'biospecimen'))
    }else{
        updateSelectInput(session, 'clin_or_bio', choices=c('clinical'))
    }
})

#update feature selection options
update_single_feature <- reactive({
    obj <- update_single_disease()
    update_cli_bio()
    DB <- obj$DB
    clin <- obj$clinical_table
    if(input$clin_or_bio=='clinical'){
        if(DB$name == 'target'){
            available_features <- get_sub_disease_clinical_features(input$single_exp_disease)
        }else if(DB$name=='pbta'){
            available_features = names(clin)[!names(clin) %in% c('experimental_strategy', 'sample_type')]
        }else{
            available_features <- names(clin)
        }
        available_features <- get_ava_features(clin[,available_features], DB, 'single')
    }else{
        available_features <- c('SampleType', 'TumorDescription')
    }
    if (DB$name=='target'){
        clin <- clin[, c(DB$clinical_id_name, DB$biospec_id_name, DB$disease_tag,available_features)]
    }else if(DB$name=='pptc'){
        if(available_features %in% names(clin)){
            clin <- clin[, c(DB$biospec_id_name, DB$disease_tag, available_features)]
        }
    }
    if (!(input$single_exp_disease=='pan-cancer')){
        clin <- clin[clin[ ,DB$disease_tag] == input$single_exp_disease, ]
    }
    updateSelectInput(session, 'box_feature', label = NULL, choices = available_features)
    return(clin)
})

#detect event of user click Example link
observeEvent(input$single_clinical_example_gene,{
    updateSelectizeInput(session, 'geneName', server = TRUE, choices=GENE_LIST, selected='PDGFRA')
})

#single gene plot
output$geneBoxPlt <- plotly::renderPlotly({

    source('ScatterPlot.R')
    source('BoxPlot.R')
    source('FeatureTools.R')

    #hide the download button
    shinyjs::hide("BoxPlotDown")
    shinyjs::hide("BoxDataDown")

    #update options of gene selection
    update_single_gene()

    #get the database instance
    DB <- update_single_disease()$DB

    #get the clinical matrix based on selected parameters
    clin <- update_single_feature()

    #protect undefined column selected error
    validate(
        need(input$box_feature %in% names(clin), '')
    )

    #transfer the gene into plotable format
    gene_name = single_gene_special_transfer(input$geneName)

    #check if user has input a gene
    validate(
        need(input$geneName!='', "please select a gene"),
        errorClass='ValidateInput'
    )

    #get expression data
    exp_data <- process_exp_data(DB, gene_name, 'single')

    #verify is gene if available is the selected database
    validate(
        need(length(rownames(exp_data))>0, "This gene is not available in this DB"),
        errorClass="ErrorOutput"
    )

    print(paste0('expression data has ', as.character(length(rownames(exp_data))), 'samples'))
    print(paste0('clinical data has ', as.character(length(rownames(clin))), 'samples'))

    #join the expression and clinical data
    rawData <- merge(exp_data, clin, by = DB$biospec_id_name)
    rawData <- rawData[order(rawData[,DB$biospec_id_name]),]

    #if database instance is TARGET, remove samples that represent same patients
    if(DB$name=='target' && input$clin_or_bio=='clinical'){
        rawData <- rawData[!duplicated(rawData[c(DB$clinical_id_name)]),]
    }

    #transfer some special gene names that are not plotable to valid gene names
    names(rawData)[names(rawData)!=DB$biospec_id_name] <- sapply(names(rawData)[names(rawData)!=DB$biospec_id_name], single_gene_special_transfer)

    #determine the input feature is numeric or categorical
    condition <- (is.numeric(rawData[,input$box_feature]))

    #if input feature is numerical then create a scatter plot
    if (condition){
        shinyjs::hide("order")
        shinyjs::hide("minimum_box_size")
        final_plot = createScatter(rawData, gene_name, input$box_feature,
                                   plot_loc=single_gene_plot_loc, data_loc=single_gene_data_loc,
                                   log_scale=input$log_trans, DB=DB)
    }else
    #if input feature is categorical then create a box plot
    {
        shinyjs::show("order")
        shinyjs::show("minimum_box_size")
        if(input$single_exp_disease=='pan-cancer'){
            final_plot = createBox(rawData, gene_name, DB$disease_tag, input$box_feature,DB=DB,
                                   plot_loc=single_gene_plot_loc, data_loc=single_gene_data_loc,
                                   ascending=input$order, minimum_box_size=input$minimum_box_size,
                                   remove_jitter=input$remove_jitter, log_scale=input$log_trans)
        }else{
            final_plot = createBox(rawData, gene_name, NULL, input$box_feature,DB=DB,
                                   plot_loc=single_gene_plot_loc, data_loc=single_gene_data_loc,
                                   ascending=input$order, minimum_box_size=input$minimum_box_size,
                                   remove_jitter=input$remove_jitter, log_scale=input$log_trans)
        }
    }
    shinyjs::show("BoxPlotDown")
    shinyjs::show("BoxDataDown")
    return(final_plot)
})