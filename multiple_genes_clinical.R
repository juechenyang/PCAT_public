# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 12/10/19

#===========================================MULTIPLE GENE ANALYSIS=============================================================================
update_multiple_disease <- reactive({
    DB <- chooseDB(input$db_multiple_name)
    clin <- QueryAll(DB$clinical_db_name)
    updateSelectInput(session, 'multiple_gene_disease', choices=c('pan-cancer', unique(clin[,DB$disease_tag])))
    return(list(DB=DB, clin=clin))
})
update_multiple_feature <- reactive({
    obj = update_multiple_disease()
    DB <- obj$DB
    clin <- obj$clin
    if(DB$name == 'target'){
        available_features <- get_sub_disease_clinical_features(input$multiple_gene_disease)
    }else if (DB$name == 'pptc'){
        available_features <- names(clin)
    }
    available_features = get_ava_features(clin[,available_features], DB, 'single')
    updateSelectInput(session, 'heatmap_feature', label = NULL, choices = available_features)
    return(list(DB=DB, clin=clin))
})

#define reactive value for multiple analysis run
multiple_run <- reactiveValues(clicked = NULL)


observeEvent(input$analyze_multiple, {
    multiple_run$clicked <- TRUE
})

observeEvent(input$geneNames, {
    multiple_run$clicked <- FALSE
})
observeEvent(input$example_genes,{
    output$gene_names_ui = renderUI({
        textAreaInput('geneNames', label="Gene signature:", height='100px',
                      value = "PYCR1\nPHGDH\nDCTPP1\nMCM5\nMSH6\nPOLA2\nHPDL\nDKC1\nANAPC1\nPOLE2\nIMPDH2\nFOXRED2\nPFKM\nPUS7\nGTPBP1\nPRMT1\nLIG3\nZNF317\nLYAR",
                      placeholder="Please type at least 2 valid genes to run the analysis and use single space or line breaker to seperate")
    })
})

get_textarea_ui = reactive({
    output$gene_names_ui = renderUI({
        textAreaInput('geneNames', label="Gene signature:", height='100px',
        value = "", placeholder="Please type at least 2 valid genes to run the analysis and use single space or line breaker to seperate")
    })
})

#multiple genes plot
output$geneHeatMap <- renderPlot({

    source("ComplexHeatmap.R")
    source('FeatureTools.R')
    #hide download
    #show downlaod button
    shinyjs::hide("multigenes_clinical_plot_download")
    shinyjs::hide("multigenes_clinical_data_download")

    #change select genes if click example genes
    get_textarea_ui()


    # get input genes as a list
    gene_list <- multiple_gene_processor(input$geneNames)

    #transfer genes to a plot friendly format
    gene_list = sapply(gene_list, single_gene_special_transfer)

    #update options for clinical features selection bars and get DB object and clinical table
    obj = update_multiple_feature()
    clin = obj$clin
    DB <- obj$DB

    #get sub corhort
    if(input$multiple_gene_disease!='pan-cancer'){
        clin <- clin[clin[ ,DB$disease_tag] == input$multiple_gene_disease, ]
    }

    #get the expression data
    exp_data <- process_exp_data(DB, gene_list, 'multiple')

    #detect if each gene has expression data available
    names(gene_list) <- sapply(gene_list, function(x){
        if(any((toupper(x))==names(exp_data))==TRUE){
            return('y')
        }else{
            return('n')
        }
    })

    #get expression available genes and unavailable genes
    found <- gene_list[names(gene_list)=='y']
    not_found <- gene_list[names(gene_list)=='n']

    # if any input gene is not found, render the not found gene in gene checker
    if(length(not_found)>0){
        output$error_gene <- renderUI(
            tags$div(class="alert alert-danger mb-3",
                paste(paste(not_found, collapse=',\n'), 'are not found in this DB')
            )
        )
    }else{
        if(length(found)==0){
            output$error_gene <- renderUI('')
        }else{
            output$error_gene <- renderUI(
                tags$div(class="alert alert-success mb-3",
                    'All genes are good'
                )
            )
        }
    }

    #if more than 1 of the input genes are valid, change the run checker into runnable
    if(length(found)>=2){
        output$multiple_run_check <- renderUI(
            actionButton(class='btn btn-success', 'analyze_multiple', 'RUN')
        )
    }else{
        multiple_run$clicked <- FALSE
        output$multiple_run_check <- renderUI(
            actionButton(class='btn btn-danger', 'analyze_multiple', 'NOT RUNNABLE')
        )
    }


    # yes = input$analyze_multiple[1]
    # print(paste0('yes is ', yes))
    # print(paste0('class of yes is ', class(yes)))
    #
    # output$clicked <- reactive({
    #     length(found)>=2
    # })
    # outputOptions(output, 'clicked', suspendWhenHidden = FALSE)


    #if not clicked don't do anything
    if (multiple_run$clicked == FALSE) return()




    #protect the undefined column selected error
    validate(
        need(input$heatmap_feature %in% names(clin), '')
    )

    #inner join expression and clinical data
    exp_clin <- merge(clin, exp_data, by=DB$biospec_id_name)

    #draw the heatmap
    ht = Multigene_Heatmap(exp_clin, found, input$heatmap_feature, input$rowCluster,
                           input$columnCluster, plot_loc = multiple_genes_plot_loc,
                           data_loc = multiple_genes_data_loc, DB=DB
                           )

    #show downlaod button
    shinyjs::show("multigenes_clinical_plot_download")
    shinyjs::show("multigenes_clinical_data_download")

    return(ht)
})