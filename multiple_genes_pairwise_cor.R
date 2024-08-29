# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 12/10/19
#================================================================pairwise_cor================================================================
pairwise_run = reactiveValues(clicked=FALSE)
observeEvent(input$pairwise_run_button, {
    pairwise_run$clicked = TRUE
})
observeEvent(input$pairwise_genes, {
    pairwise_run$clicked = FALSE
})

update_pairwise_disease = reactive({
    DB <- chooseDB(input$pairwise_db)
    clin = QueryAll(DB$clinical_db_name)
    updateSelectInput(session, 'pairwise_disease', choices=c('pan-cancer', unique(clin[,DB$disease_tag])))
    return(list(DB=DB, clin=clin))
})

get_pairwise_textarea = reactive({
    output$pairwise_gene_names_ui = renderUI({
        textAreaInput('pairwise_genes', label="Gene signature:", height='100px',
        value = "", placeholder="Please type at least 2 valid genes to run the analysis and use single space or line breaker to seperate")
    })
})

observeEvent(input$pairwise_example_genes, {
    output$pairwise_gene_names_ui = renderUI({
        textAreaInput('pairwise_genes', label="Gene signature:", height='100px',
        value = "PYCR1\nPHGDH\nDCTPP1\nMCM5\nMSH6\nPOLA2\nHPDL\nDKC1\nANAPC1\nPOLE2\nIMPDH2\nFOXRED2\nPFKM\nPUS7\nGTPBP1\nPRMT1\nLIG3\nZNF317\nLYAR\nFAM134B\nRNF146\nTMEM150C\nDNER\nRND2\nKIAA0895\nPLEKHA3\nCSRNP3\nIFT20\nNRBP2\nAGPAT4\nTPPP3\nTRAK2\nOSCP1\nIQCG\nEPB41L1\nFAM174A"
        )
    })
})


output$pairwise_heatmap = renderPlot({

    source("ComplexHeatmap.R")

    #fill annotation ui output
    output$corrplot_annotations = renderUI({
        tags$div(class='tlp mt-4',
            tags$i(class="fas fa-question-circle"),
            tags$span("How to read this plot"),
            tags$div(class="tooltip_text",
                tags$div(class="font-awesome ml-4",
                    tags$ul(
                        tags$li(
                            "The size of circles indicate the absolute value of R"
                        ),
                        tags$li(
                            "Empty cell means the R value is not significant (p>0.05) for the corresponding two genes"
                        )
                    )
                )
            )
        )
    })

    #hide download button
    shinyjs::hide("pairwise_download")
    shinyjs::hide("pairwise_data_download")
    shinyjs::hide("corrplot_annotations")

    #initialize gene names gene names area
    get_pairwise_textarea()

    # get input genes as a list
    gene_list <- multiple_gene_processor(input$pairwise_genes)
    gene_list = sapply(gene_list, single_gene_special_transfer)

    #get the database object
    obj = update_pairwise_disease()
    DB = obj$DB
    clin = obj$clin

    #if cohort is not pancan, then hide option for removing tissue effect
    if(input$pairwise_disease != "pan-cancer"){
        shinyjs::hide("pairwise_remove_tissue_effect")
    }else{
        shinyjs::show("pairwise_remove_tissue_effect")
    }

    #if remove tissue effect is clicked the expression matrix should be the z-score transformed version
    if(input$pairwise_remove_tissue_effect){
        DB$exp_db_name = paste0(DB$name, '_transformed_expr')
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


    if(length(not_found)>0){
        output$pairwise_gene_checker <- renderUI(
            tags$div(class="alert alert-danger mb-3",
                paste(paste(not_found, collapse=',\n'), 'are not found in this DB')
            )
        )
    }else{
        if(length(found)==0){
            output$pairwise_gene_checker <- renderUI('')
        }else{
            output$pairwise_gene_checker <- renderUI(
                tags$div(class="alert alert-success mb-3",
                    'All genes are good'
                )
            )
        }
    }
    if(length(found)>=2){
        output$pairwise_run_checker <- renderUI(
            actionButton(class='btn btn-success', 'pairwise_run_button', 'RUN')
        )
    }else{
        pairwise_run$clicked = FALSE
        output$pairwise_run_checker <- renderUI(
            actionButton(class='btn btn-danger', 'pairwise_run_button', 'NOT RUNNABLE')
        )
    }

    if(!pairwise_run$clicked) {return()}

    #choose sub_disease cohort
    if (!(input$pairwise_disease=='pan-cancer')){
        clin <- clin[clin[ ,DB$disease_tag] == input$pairwise_disease, ]
    }

    #get sub disease
    exp_data = merge(exp_data, clin, by=DB$biospec_id_name)

    #if the number of observation is 1 then correlation is not able to show
    validate(
        need(length(rownames(exp_data))>5, 'too less data to compute correlation'),
        errorClass="ErrorOutput"
    )

    #reformat the expression data for correlation computation
    rownames(exp_data) = exp_data[,DB$biospec_id_name]
    exp_data = exp_data[,found]
    cor_matrix = cor(exp_data, method=input$pairwise_method)

    cor_df = data.frame(cor_matrix)
    Rvalues = rownames(cor_df)
    headers = names(cor_df)
    cor_df = tibble::add_column(cor_df, Rvalues, .before = 1)

    #write to a csv file for download
    write.table(cor_df, pairwise_data_loc, sep = ',', row.names = F)
    write.table(data.frame(), pairwise_data_loc, append = T)

    #get p-values
    p.mat <- corrplot::cor.mtest(exp_data)$p
    p_df = data.frame(p.mat)
    names(p_df) = headers
    Pvalues = Rvalues
    p_df = tibble::add_column(p_df, Pvalues, .before = 1)

    write.table(p_df, pairwise_data_loc, sep = ',', append = T, row.names = F)

    #plot the pairwise correlation figure
    final_plot = Pairwise_Heatmap(cor_matrix, exp_data, pairwise_plot_loc)


    #show download button
    shinyjs::show("pairwise_download")
    shinyjs::show("pairwise_data_download")
    shinyjs::show("corrplot_annotations")

    return(final_plot)
})
