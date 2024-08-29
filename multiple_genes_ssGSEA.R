# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 12/10/19
#===================================================ssGSEA ANALYSIS=======================================================================================================
# control ssGSEA running process
ssGSEA_button <- reactiveValues(clicked = NULL)

observeEvent(input$run_ssGSEA, {
   ssGSEA_button$clicked <- runif(10)
})
observeEvent(input$ssGSEAlist, {
    ssGSEA_button$clicked <- NULL
})
# observeEvent(input$ssGSEA_disease, {
#     ssGSEA_button$clicked <- NULL
# })
observeEvent(input$ssGSEA_DB, {
    ssGSEA_button$clicked <- NULL
})

#update disease selection for ssGSEA
update_ssGSEA_clinical_disease <- reactive({
    DB <- chooseDB(input$ssGSEA_DB)
    clin <- QueryAll(DB$clinical_db_name)
    updateSelectInput(session, 'ssGSEA_disease', choices=c('pan-cancer', unique(clin[,DB$disease_tag])))
    return(list(DB=DB, clin=clin))
})

#update clinical features for ssGSEA
update_ssGSEA_clinical_feature <- reactive({
    obj <- update_ssGSEA_clinical_disease()
    DB = obj$DB
    clin = obj$clin
    if(DB$name == "target"){
        available_features <- get_sub_disease_clinical_features(input$ssGSEA_disease)
    }else if(DB$name == "pptc"){
        available_features <- names(clin)
    }
    available_features = get_ava_features(clin[,available_features], DB, 'single')
    updateSelectInput(session, 'ssGSEA_feature', label = NULL, choices = available_features)
    return(list(DB=DB, clin=clin))
})
generate_ssGSEA_plot_data <- reactive({

    #get DB object and clinical matrix
    obj = update_ssGSEA_clinical_disease()
    DB = obj$DB
    clin =obj$clin

    #process gene signature
    gene_list <- multiple_gene_processor(input$ssGSEAlist)
    gene_list = sapply(gene_list, single_gene_special_transfer)

    #get gene signature expression
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

    #if some genes are not found, gene checker display not found genes
    if(length(not_found)>0){
        output$ssGSEA_gene_check <- renderUI(
            tags$div(class="alert alert-danger mb-3",
                paste(paste(not_found, collapse=',\n'), 'are not found in this DB')
            )
        )
    }else{
        if(length(found)==0){
            output$ssGSEA_gene_check <- renderUI('')
        }else{
            output$ssGSEA_gene_check <- renderUI(
                tags$div(class="alert alert-success mb-3",
                    'All genes are good'
                )
            )
        }
    }

    #if more than 1 of the input genes are valid, change the run checker into runnable
    if(length(found)>=2){
        output$ssGSEA_run_check <- renderUI(
            actionButton(class='btn btn-success', 'run_ssGSEA', 'RUN')
        )
    }else{
        ssGSEA_button$clicked <- NULL
        output$ssGSEA_run_check <- renderUI(
            actionButton(class='btn btn-danger', 'run_ssGSEA', 'NOT RUNNABLE')
        )
    }

    #if not clicked don't do anything
    if (is.null(ssGSEA_button$clicked)) return(NULL)

    #status progress bar 0.4
    setProgress(0.4, detail = paste("Loading whole expression"))

    #get whole expression profile
    # ssGSEA_exp <- readRDS(paste('IntermediateFiles/', 'my_data1.rds', sep=''))
    ssGSEA_exp <- readRDS(paste('IntermediateFiles/', DB$name, '_ssGSEA.rds', sep=''))

    #status progress bar 0.75
    setProgress(0.75, detail = paste("computing score, this may take 2-3mins. Please be patient"))

    # compute ssGSEAScores
    ssGSEAScores <<- ssGSEAData(found, ssGSEA_exp)

    # ssGSEAScores[[DB$biospec_id_name]] <- str_replace_all(ssGSEAScores$names,"\\.","-")
    #timer
    print(paste('start merge expression and score at ', Sys.time(), ''))

    #merge scores and gene signature expression
    exp_score <- merge(exp_data, ssGSEAScores, by=DB$biospec_id_name)

    #timer
    print(paste('start merge all data at ', Sys.time(), ''))

    # merge clinical with exp and scores
    exp_score_clinical <- merge(clin, exp_score, by=DB$biospec_id_name)
    exp_score_clinical$Anoymous <- as.numeric(as.character(exp_score_clinical$Anoymous))
    exp_score_clinical = rename_a_column(exp_score_clinical, 'Anoymous', 'ssGSEA_scores')

    #if(input$ssGSEA_DB=='pptc'){
    #    #get drug response data
    #    drug_response_data = QueryAll(DB$drug_db)
    #    selected_drug_response_data = drug_response_data[,c(DB$biospec_id_name, input$ssGSEA_drug)]
    #
    #    #join drug response data with ssGSEA plot data
    #    exp_score_clinical = dplyr::left_join(exp_score_clinical, selected_drug_response_data, by=DB$biospec_id_name)
    #}

    #write data output
    setProgress(0.85, detail = paste("writing data out"))
    write.csv(exp_score_clinical[,c(DB$biospec_id_name, 'ssGSEA_scores')], file=ssGSEAData_loc, row.names=FALSE)

    setProgress(1, detail = paste("ploting"))
    #timer
    print(paste('start ploting at ', Sys.time(), ''))
    return(list(DB=DB, clin=clin, plotdata=exp_score_clinical, valid_genes=found))

})

get_ssGSEA_textarea = reactive({
    output$ssGSEA_gene_names_ui = renderUI({
        textAreaInput('ssGSEAlist', label="Gene signature:", height='100px',
        value = "", placeholder="Please type at least 2 valid genes to run the analysis and use single space or line breaker to seperate")
    })
})

observeEvent(input$ssGSEA_example_genes, {
    output$ssGSEA_gene_names_ui = renderUI({
        textAreaInput('ssGSEAlist', label="Gene signature:", height='100px',
        value = "PYCR1\nPHGDH\nDCTPP1\nMCM5\nMSH6\nPOLA2\nHPDL\nDKC1\nANAPC1\nPOLE2\n
        IMPDH2\nFOXRED2\nPFKM\nPUS7\nGTPBP1\nPRMT1\nLIG3\nZNF317\nLYAR")
    })

})
update_drug_selection = reactive({
    #if input db is pptc then show input box for drug input
    if(input$ssGSEA_DB=='target'){
        shinyjs::hide('ssGSEA_drug')
    }else{
        DB = chooseDB('pptc')
        #get all drugs
        dr_data = QueryAll(DB$drug_db)
        all_drug = names(dr_data)[2:length(names(dr_data))]
        updateSelectizeInput(session, 'ssGSEA_drug', label = 'Choose a drug(optional):',
                             server = TRUE, choices=all_drug, selected="")
    }
})

output$ssGSEA_barPlot <- renderPlot({

    source("ssGSEA.R")
    source("ComplexHeatmap.R")
    source('FeatureTools.R')
    library(feather)

    #hide downloadButton
    shinyjs::hide("ssGSEAPlotDown")
    shinyjs::hide("ssGSEADataDown")
    shinyjs::show('ssGSEA_drug')
    shinyjs::hide('drug_mechanism')

    #update disease and feature options
    update_ssGSEA_clinical_feature()

    #update_drug_selection
    update_drug_selection()

    #initialize text area
    get_ssGSEA_textarea()



    #process data and plot with progress bar
    withProgress(message = 'Processing', value = 0, {

        #get data for plot

        print(paste0("ssGSEA started at ", Sys.time()))
        obj = generate_ssGSEA_plot_data()
        ssGSEAPlotData <- obj$plotdata
        DB = obj$DB
        clin = obj$clin
        valid_genes = obj$valid_genes


        #check if ssGSEAPlotData is available
        if(is.null(ssGSEAPlotData)){
            return(NULL)
        }else{
            gene_list <- multiple_gene_processor(input$ssGSEAlist)

            #get sub corhort
            if(input$ssGSEA_disease!='pan-cancer'){
                ssGSEAPlotData <- ssGSEAPlotData[ssGSEAPlotData[ ,DB$disease_tag] == input$ssGSEA_disease, ]
            }
            update_ssGSEA_clinical_feature()

            if(input$ssGSEA_DB=='pptc' & input$ssGSEA_drug!=""){
                drug = input$ssGSEA_drug

                #get drug annotation data
                drug_anno = QueryAll(DB$drug_anno_db)

                #get annotation of selected drug
                drug_mechanism = drug_anno[drug_anno[['Drug_Name']]==input$ssGSEA_drug, 'Mechanism']

                output$drug_mechanism = renderUI({
                    tags$div(class='tlp mb-3',
                        tags$i(class="fas fa-question-circle"),
                        tags$span(paste0("What is ", input$ssGSEA_drug)),
                        tags$div(class="tooltip_text",
                            tags$div(class="ml-1 my-4 font-size-one",
                                 paste0(drug_mechanism)
                            )
                        )
                    )
                })
                shinyjs::show('drug_mechanism')
                #get drug response data
                drug_response_data = QueryAll(DB$drug_db)
                selected_drug_response_data = drug_response_data[,c(DB$biospec_id_name, input$ssGSEA_drug)]

                #join drug response data with ssGSEA plot data
                ssGSEAPlotData = dplyr::left_join(ssGSEAPlotData, selected_drug_response_data, by=DB$biospec_id_name)

                #replace NA and 0 with "NA"
                ssGSEAPlotData[[input$ssGSEA_drug]][is.na(ssGSEAPlotData[[input$ssGSEA_drug]])] = "NA"
                #replace number values with annotated values
                values = c("NA", "NA", "PD1", "PD2", "SD", "PR", "CR", "MCR")
                names(values) = c("NA", "0", "1", "2", "3", "4", "5", "6")
                ssGSEAPlotData[[input$ssGSEA_drug]] = sapply(ssGSEAPlotData[[input$ssGSEA_drug]], function (x){
                    return(values[x])
                })
            }else{
                drug = NULL
            }



            print(paste0("ssGSEA finished at ", Sys.time()))
            if(input$ssGSEA_DB=='pptc'){
                final_plot = ssGSEA_Heatmap(ssGSEAPlotData, valid_genes, input$ssGSEA_feature,
                                            ssGSEAPlot_loc, drug = drug)
            }else{
                final_plot = ssGSEA_Heatmap(ssGSEAPlotData, valid_genes, input$ssGSEA_feature, ssGSEAPlot_loc)
            }


            #show download button
            shinyjs::show("ssGSEAPlotDown")
            shinyjs::show("ssGSEADataDown")

            return(final_plot)
        }


    })
})
