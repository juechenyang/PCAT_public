# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 12/9/19
library(ggplot2)
km_plot_loc <- "StaticFiles/Plots/km_plot.pdf"
forrest_plot_loc <- "StaticFiles/Plots/forrest_plot.pdf"
#=========================================SINGLE GENE SURVIVAL=========================================================================================================================================
update_survival_gene <- reactive({
    updateSelectizeInput(session, 'survival_gene', server = TRUE, choices=GENE_LIST,selected='')
})
update_survival_disease <- reactive({
    DB <- chooseDB(input$survival_db)
    clinical <- QueryAll(DB$clinical_db_name)
    updateSelectInput(session, 'survival_disease', choices=c(unique(clinical[,DB$disease_tag])))
    return(list(db=DB, clinical_db=clinical))
})
update_survival_multivariate <- reactive({
    data_repo <- update_survival_disease()
    clin <- data_repo$clinical_db
    DB <- data_repo$db
    available_features <- get_sub_disease_clinical_features(input$survival_disease)
    available_features <- get_ava_features(clin[,available_features], DB,'single')
    updateSelectizeInput(session, 'survival_multivariate', server=TRUE, choices=available_features)
})

hide_custom <- reactive({
    shinyjs::hide("survival_threshold")
})

observeEvent(input$single_gene_survival_example,{
    updateSelectizeInput(session, 'survival_gene', server = TRUE, choices=GENE_LIST, selected='TERT')
})

get_survival_df <- reactive({

    shinyjs::hide("km_plot_download")
    shinyjs::hide("km_data_download")
    shinyjs::hide("forrest_plot_download")
    shinyjs::hide("forrest_data_download")

    DB <- chooseDB(input$survival_db)
    update_survival_gene()
    #fetch clinical table
    clin <- update_survival_disease()$clinical_db
    #update multivariate parameters
    update_survival_multivariate()

    hide_custom()
    #validate if a user has input a gene
    validation_logic <- !((input$survival_gene=='')&&is.null(input$survival_multivariate))
    validate(
        need(validation_logic, 'please input a gene or a variant'),
        errorClass='InputCheck'
    )

    #get clinical matrix
    clin <- clin[clin[ ,DB$disease_tag] == input$survival_disease, ]
    #determine if user has defined a gene
    if(input$survival_gene==''){
        gene <- NULL
        survival_df <- clin
    }else{
        #get transferred gene name to data slicing
        gene <- single_gene_special_transfer(input$survival_gene)
        #get exepression matrix
        exp_data <- process_exp_data(DB, gene, 'single')
        # names(exp_data)[names(exp_data)!=DB$biospec_id_name] <- sapply(names(exp_data)[names(exp_data)!=DB$biospec_id_name], single_gene_special_transfer)

        #merge expresison and clinical data as survival dataframe
        survival_df <- merge(exp_data, clin, by=DB$biospec_id_name)

    }

    #sort the combined table based on sample
    survival_df <- survival_df[order(survival_df[,DB$biospec_id_name]),]
    #remove redundant samples that represent same patient
    survival_df <- survival_df[!duplicated(survival_df[c(DB$clinical_id_name)]),]

    #slicing survival dataframe
    survival_required_features <- c("Event1st","Event1stStat","SurvDaysEF","VitalStatus","VitalStatNo","SurvDaysTotal")
    required_parameters <- c(gene, input$survival_multivariate)

    required_features <- c(required_parameters, survival_required_features)
    survival_df <- survival_df[,required_features]

    return(list(df=survival_df, parameter=required_parameters))
})

output$survival_km <- renderPlot({

    source("SurvivalTools.r")
    source("FeatureTools.R")

    #hide km_plot_download and threshold
    shinyjs::hide("survival_threshold")
    obj <- get_survival_df()
    survival_df <- obj$df
    required_parameters <- obj$parameter

    #get the first parameter to plot the KM curve
    parameter <- required_parameters[[1]]

    #build the text annotation
    survival_group <- paste(parameter, 'group', sep='.')

    #determin if the parameter is numeric or categorical
    if(class(survival_df[[parameter]])=='numeric'){

        #if over 90% of the samples has the numeric values of 0, then report error of too many 0 expression to get survial plot
        validate(
            need(length(survival_df[[parameter]][survival_df[[parameter]]==0])<length(rownames(survival_df))*0.9, 'too many 0 expreesion for plot'),
            errorClass='ErrorOutput'
        )
        shinyjs::show("survival_cutoff_method")
        validate(
            need(length(unique(survival_df[[parameter]]))>1, paste(parameter, 'is numeric but only has one value which is', as.character(unique(survival_df[[parameter]])))),
            errorClass="ErrorOutput"
        )
        if(input$survival_cutoff_method=='custom'){
            min_v <- min(survival_df[[parameter]])
            max_v <- max(survival_df[[parameter]])
            step_gap <- (max_v-min_v)/1000
            updateSliderInput(session, 'survival_threshold', min=min_v, max=max_v, step=step_gap)
            shinyjs::show("survival_threshold")
        }else{
            shinyjs::hide("survival_threshold")
        }
        if(input$survival_type=='eventfree'){
            threshold <- cutoff_function_options(survival_df,parameter,input$survival_cutoff_method,custom=input$survival_threshold, "SurvDaysEF", "Event1stStat")
            survival_df <- survival_df[!is.na(survival_df[,'SurvDaysEF']), ]
        }else if(input$survival_type=='overall'){
            threshold <- cutoff_function_options(survival_df,parameter,input$survival_cutoff_method,custom=input$survival_threshold, "SurvDaysTotal", "VitalStatNo")
            survival_df <- survival_df[!is.na(survival_df[,'SurvDaysTotal']), ]
        }
        validate(
            need(length(rownames(survival_df))!=0, 'no enough data to plot survival curve'),
            errorClass="ErrorOutput"
        )
        survival_df[[survival_group]] <- sapply(survival_df[[parameter]], function(x){
            temp <- ifelse(x>=threshold,"High","Low")
            return(temp)
        })
    }else{
        colnames(survival_df)[colnames(survival_df)==parameter] <- survival_group

        #remove unknown observations
        survival_df <- survival_df[!grepl("Unknow|unknown|Unknown|\\.|not reported|Not Reported", survival_df[, survival_group]),]

        #remove NA observations in survival variables
        if(input$survival_type=='eventfree'){
            survival_df <- survival_df[!is.na(survival_df[,'SurvDaysEF']), ]
        }else if(input$survival_type=='overall'){
            survival_df <- survival_df[!is.na(survival_df[,'SurvDaysTotal']), ]
        }
    }

    #save km plot data
    if(!is.null(km_data_loc)){
        write.csv(survival_df, file=km_data_loc, row.names=FALSE)
    }

    km_plt <- run_survival(survival_df, input$survival_type, survival_group, parameter, threshold)
    ggsave(file = km_plot_loc, print(km_plt, newpage = FALSE))

    #show download button
    shinyjs::show("km_plot_download")
    shinyjs::show("km_data_download")
    return(km_plt)
})


output$survival_forest <- renderPlot({
    return_obj <- get_survival_df()
    survival_df <- return_obj$df
    required_parameters <- return_obj$parameter

    #length of parameters must be larger than 1 for forrest_plot
    validate(
        need(length(required_parameters)>1, '\n\nadd covariates to perform multivariate analysis'),
        errorClass="ValidateInputSpecial"
    )

    #remove NA observations
    survival_df = dplyr::filter_all(survival_df, dplyr::all_vars(!grepl('Unknow|unknown|Unknown|not reported|Not Reported',.)))


    if(input$survival_gene==''){
        parameter_string <- paste(required_parameters, collapse='+')
    }else{
        parameter <- required_parameters[[1]]

        #if over 90% of the samples has the numeric values of 0, then report error of too many 0 expression to get survial plot
        validate(
            need(length(survival_df[[parameter]][survival_df[[parameter]]==0])<length(rownames(survival_df))*0.9, 'too many 0 expreesion for plot'),
            errorClass='ErrorOutput'
        )



        if(input$survival_type=='eventfree'){
            threshold <- cutoff_function_options(survival_df,parameter,input$survival_cutoff_method,custom=input$survival_threshold, "SurvDaysEF", "Event1stStat")
        }else if(input$survival_type=='overall'){
            threshold <- cutoff_function_options(survival_df,parameter,input$survival_cutoff_method,custom=input$survival_threshold, "SurvDaysTotal", "VitalStatNo")
        }
        survival_group <- paste(parameter, 'group', sep='.')
        survival_df[[survival_group]] <- sapply(survival_df[[parameter]], function(x){
            temp <- ifelse(x>=threshold,"High","Low")
            return(temp)
        })
        parameter_string <- paste(c(survival_group, required_parameters[2:length(required_parameters)]), collapse='+')
    }

    surv_object <- get_surv_obj(survival_df, input$survival_type)
    vec_vars = "surv_object~"
    # parameter_string <- paste(survival_group, parameter_string, sep='+')
    formula_string <- paste(vec_vars, parameter_string, sep='')
    fit_coxph = coxph(as.formula(formula_string),data=survival_df)
    plot_forest = ggforest(fit_coxph,data=survival_df,fontsize = 1)

    #save forrest plot data
    write.csv(survival_df, file=forrest_data_loc, row.names=FALSE)

    #save forrest plot data
    ggsave(file = forrest_plot_loc, plot_forest)


    shinyjs::show("forrest_plot_download")
    shinyjs::show("forrest_data_download")
    return(plot_forest)

})