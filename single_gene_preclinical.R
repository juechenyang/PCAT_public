# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 12/10/19
#============================================SINGLE GENE PRECLINICAL======================================================================================
update_preclinical_disease <- reactive({
    DB <- chooseDB("pptc")
    clinical <- QueryAll(DB$clinical_db_name)[,c(DB$biospec_id_name, DB$disease_tag)]
    updateSelectInput(session, 'disease_single_preclinical', choices=c('pan-cancer', unique(clinical[,DB$disease_tag])))
    return(list(DB=DB, clinical=clinical))
})
update_preclinical_treatment <- reactive({
    data_combo <-update_preclinical_disease()
    DB <- data_combo$DB
    drug_response <- QueryAll(DB$drug_db)
    available_drugs <- get_ava_features(drug_response, DB, 'single')
    updateSelectInput(session, 'drugs', choices=available_drugs)
    clinical_drug <- merge(data_combo$clinical, drug_response, by='sample_id')
    return(list(DB=DB, clinical_drug=clinical_drug))
})
update_preclinical_gene <- reactive({
    updateSelectizeInput(session, 'geneName_single_preclinical', server = TRUE, choices=GENE_LIST, selected='')
})
output$preclinical_boxplot <- plotly::renderPlotly({

    source('BoxPlot.R')
    source('FeatureTools.R')

    shinyjs::hide("preclinical_plot_download")
    shinyjs::hide("preclinical_data_download")

    DB <- chooseDB("pptc")
    update_preclinical_gene()
    combo <- update_preclinical_treatment()

    #get drug annotation data
    drug_anno = QueryAll(DB$drug_anno_db)

    #get annotation of selected drug
    drug_mechanism = drug_anno[drug_anno[['Drug_Name']]==input$drugs, 'Mechanism']

    #output mechanism to ui
    output$drug_mechanism = renderUI({
        tags$div(class='tlp mb-3',
            tags$i(class="fas fa-question-circle"),
            tags$span(paste0("What is ", input$drugs)),
            tags$div(class="tooltip_text",
                tags$div(class="ml-1 my-4 font-size-one",
                     paste0(drug_mechanism)
                )
            )
        )
    })

    validate(
        need(input$geneName_single_preclinical!='', 'please input a gene'),
        errorClass="ValidateInput"
    )
    total_clinical <- combo$clinical_drug



    if (!(input$disease_single_preclinical=='pan-cancer')){
        total_clinical <- total_clinical[total_clinical[ ,combo$DB$disease_tag] == input$disease_single_preclinical, ]
    }

    #transfer gene name to plotable format
    gene_name = single_gene_special_transfer(input$geneName_single_preclinical)

    #fetch the expression data
    exp_data <- process_exp_data(DB, gene_name, 'single')

    #if gene does not have expression data, return error message
    validate(
        need(gene_name %in% colnames(exp_data), "This gene does not have expression data in this DB"),
        errorClass="ErrorOutput"
    )

    #merge expresison data with preclinical data
    merged_data <- merge(exp_data, total_clinical, by=DB$biospec_id_name)


    #validate if box plot is qualified to be shown
    validate(
        need(length(rownames(merged_data))>0, 'no enough data for box plot'),
        errorClass="ErrorOutput"
    )

    #define where there is drug mut plot to be triggered
    drug_mut = F

    if(!input$remove_mut_point){

        drug_mut = T

        #fetch the mutation data
        mut = QueryOnRow(DB$mut_db_name, gene_name)

        mut_info = QueryAll(DB$mut_info_db)

        mut = merge(mut, mut_info, by.x=c("gene_name", "sample_id", "mut_group"), by.y=c("gene_name", "sample_id", "mut_group"), all.x=TRUE)

        #if no result then display this gene does not exist
        validate(
            need(length(rownames(mut))>0, 'this gene does not have mutation data in this DB'),
            errorClass="ErrorOutput"
        )

        #merge expresion preclin with mutation
        merged_data <- merge(merged_data, mut, by=DB$biospec_id_name)

    }


    #remove all 0 observation for drug response
    merged_data = merged_data[which(merged_data[[input$drugs]]!="0"), ]

    validate(
        need(length(rownames(merged_data))>0, 'all drug responsed level is 0'),
        errorClass="ErrorOutput"
    )

    #replace number values with annotated values
    values = c("PD1", "PD2", "SD", "PR", "CR", "MCR")
    names(values) = c("1", "2", "3", "4", "5", "6")
    merged_data[[input$drugs]] = sapply(merged_data[[input$drugs]], function (x){
        return(values[x])
    })
    merged_data[[input$drugs]] = factor(merged_data[[input$drugs]], levels = c("PD1", "PD2", "SD", "PR", "CR", "MCR"))

    final_plot = createBox(merged_data, gene_name, disease=NULL, feature=input$drugs, plot_loc = single_preclinical_plot_loc,
                           ascending=T, drug_plot=T, drug_mut_plot = drug_mut, input$preclinical_minimum, log_scale=input$drug_log_trans,
                           data_loc = single_preclinical_data_loc, hover_string = DB$hover_string, DB=DB
    )
    shinyjs::show("preclinical_plot_download")
    shinyjs::show("preclinical_data_download")
    return(final_plot)
})
