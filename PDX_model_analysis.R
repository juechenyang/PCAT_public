# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 12/10/19
#======================================================PDX model selector======================================================
#update disease options for
library(DT)
source("tools.R")
#create link to get the detail expression location
link_to_fusion <- function(donor_gene, acceptor_gene, sample_id, db_name, cohort) {
    sprintf('<a href="?gene_fusion_expression_locator?donor_gene=%sacceptor_gene=%ssample_id=%sdb_name=%scohort=%s" target="_blank">%s--%s</a>',
            donor_gene, acceptor_gene, sample_id, db_name, cohort, donor_gene, acceptor_gene)
}
update_PDX_disease = reactive({
    DB = chooseDB('pptc')
    clin = QueryAll(DB$clinical_db_name)
    updateSelectInput(session, 'PDX_model_disease', choices=c('pan-cancer', unique(clin[[DB$disease_tag]])))
    return(list(DB=DB, clin=clin))
})
update_PDX_gene <- reactive({
    updateSelectizeInput(session, 'PDX_model_gene', server = TRUE, choices=GENE_LIST, selected='')
})



#funtion for model detail query
createDetailLink <- function(sample_id) {
    sprintf('<a href="?PDX_model_detail?sample_id=%s" target="_blank">%s</a>', sample_id, sample_id)
}

replace_to_YesAndNo = function (x){
    ifelse(x==1, 'Yes', 'No')
}

check_avaliability = function (clin_df, key, sample){
    filtered_clin = clin_df[clin_df['sample_id']==sample, ]
    value = filtered_clin[1,key]
    return(ifelse(value==1, TRUE, FALSE))
}

update_selection = reactive({
    shinyjs::show("PDX_model_gene")

    if(input$PDX_model_var_type=="mutation" | input$PDX_model_var_type=="fusion"){
        updateSelectizeInput(session, 'PDX_model_gene', label = 'choose a gene', server = TRUE, choices=GENE_LIST, selected='')
    }else if(input$PDX_model_var_type=="drug"){
        DB = chooseDB('pptc')
        #get all drugs
        dr_data = QueryAll(DB$drug_db)
        all_drug = names(dr_data)[2:length(names(dr_data))]
        updateSelectizeInput(session, 'PDX_model_gene', label = 'choose a drug',
                             server = TRUE, choices=all_drug, selected='')
    }
})

update_drug_mechanism = reactive({
    update_selection()
    DB = chooseDB('pptc')
    dr_anno = QueryAll(DB$drug_anno_db)
    drug_mechanism = dr_anno[dr_anno[['Drug_Name']]==input$PDX_model_gene, 'Mechanism']
    output$pdx_drug_mechanism = renderUI({
        tags$div(class='tlp mb-3',
            tags$i(class="fas fa-question-circle"),
            tags$span(paste0("What is ", input$PDX_model_gene)),
            tags$div(class="tooltip_text",
                tags$div(class="ml-1 my-4 font-size-one",
                     paste0(drug_mechanism)
                )
            )
        )
    })
})

output$PDX_model_results = renderDataTable({

    library(tibble)
    update_PDX_gene()


    #update disease selection
    obj = update_PDX_disease()
    DB = obj$DB
    clin = obj$clin

    #get sub disease clinical
    if (!(input$PDX_model_disease=='pan-cancer')){
        clin <- clin[clin[ ,DB$disease_tag] == input$PDX_model_disease, ]
    }

    #change values from 0 1 to no and yes
    availability_columns = c("mutation","copy_number","expression_fusion","Drug_testing")
    clin[,availability_columns] <- lapply(clin[,availability_columns], replace_to_YesAndNo)


    if(input$PDX_model_var_type=="clinical"){
        shinyjs::hide("PDX_model_gene")

        feature_selected = c("sample_id","Histology","Histology_Detailed","Sex","Phase","Age_year","Site_of_Initial_Tumor",
                             "Site_of_Specimen","mutation","copy_number","expression_fusion","Drug_testing")
        final_df = clin
        final_df = final_df[, feature_selected]

        final_df = data.table::setnames(final_df, old = c("mutation","copy_number","expression_fusion","Drug_testing"),
                                      new = c('HasMutationData','HasCopynumberData', 'HasFusionData', 'HasDrugTestingData'))


    }else if(input$PDX_model_var_type=="mutation"){


        update_selection()
        shinyjs::hide("pdx_drug_mechanism")

        #check if user has input a gene
        validate(
            need(input$PDX_model_gene!="", 'please input a gene'),
            errorClass='InputCheck'
        )

        #get the mutation table
        mut_info = QueryAll(DB$mut_info_db)

        #get mutation for the specific gene
        mut_info = mut_info[mut_info['gene_name']==input$PDX_model_gene, ]

        #merge clinical and mutation by sample ID
        final_df = dplyr::inner_join(mut_info, clin, by='sample_id')

        feature_selected = c("sample_id","Histology","Histology_Detailed","HGVSp_Short","mut_group","Sex","Phase",
                             "Age_year","mutation","copy_number","expression_fusion","Drug_testing",
                             "Site_of_Initial_Tumor","Site_of_Specimen","Stage_of_Disease",
                             "Risk_Group","COG_studies","Prior_Therapy")
        final_df = final_df[feature_selected]


        final_df = data.table::setnames(final_df, old = c("mutation","copy_number","expression_fusion","Drug_testing"),
                                      new = c('HasMutationData','HasCopynumberData', 'HasFusionData', 'HasDrugTestingData'))

        #check if this gene has mutation
        validate(
            need(length(rownames(final_df))>0, 'this gene does not have mutation event in this cohort'),
            errorClass='ErrorOutput'
        )

    }else if(input$PDX_model_var_type=="fusion"){

        update_selection()
        shinyjs::hide("pdx_drug_mechanism")

        #get the fusion data
        fusion_df = QueryAll(DB$fusion_db)

        #check if user has input a gene
        validate(
            need(input$PDX_model_gene!="", 'please input a gene'),
            errorClass='InputCheck'
        )

        #get mutation for the specific gene
        fusion_df = fusion_df[(fusion_df['donor_gene']==input$PDX_model_gene | fusion_df['acceptor_gene']==input$PDX_model_gene), ]
        #merge with clinical data
        final_df = dplyr::inner_join(fusion_df, clin, by=DB$biospec_id_name)

        feature_selected = c("sample_id","Histology","Histology_Detailed","Fusion","Sex","Phase",
                             "Age_year","mutation","copy_number","expression_fusion","Drug_testing",
                             "Site_of_Initial_Tumor","Site_of_Specimen","Stage_of_Disease",
                             "Risk_Group","COG_studies","Prior_Therapy")

        final_df = final_df[feature_selected]

        final_df = data.table::setnames(final_df, old = c("mutation","copy_number","expression_fusion","Drug_testing"),
                                      new = c('HasMutationData','HasCopynumberData', 'HasFusionData', 'HasDrugTestingData'))

        #check if this gene has fusion
        validate(
            need(length(rownames(final_df))>0, 'this gene does not have fusion event in this cohort'),
            errorClass='ErrorOutput'
        )
    }else if(input$PDX_model_var_type=="drug"){



        update_drug_mechanism()
        shinyjs::show("pdx_drug_mechanism")

        dr_data = QueryAll(DB$drug_db)

        validate(
            need(input$PDX_model_gene!="", 'please input a drug'),
            errorClass='InputCheck'
        )


        selected_drug_data = dr_data[,c(DB$biospec_id_name, input$PDX_model_gene)]
        #remove 0 response models
        selected_drug_data = selected_drug_data[which(selected_drug_data[[input$PDX_model_gene]]!="0"), ]
        #join clinical and drug response table
        final_df = dplyr::inner_join(selected_drug_data, clin, by=DB$biospec_id_name)
        feature_selected = c("sample_id",input$PDX_model_gene, "Histology","Histology_Detailed","Sex","Phase",
                             "Age_year","mutation","copy_number","expression_fusion","Drug_testing",
                             "Site_of_Initial_Tumor","Site_of_Specimen","Stage_of_Disease",
                             "Risk_Group","COG_studies","Prior_Therapy")

        final_df = final_df[feature_selected]

        final_df = data.table::setnames(final_df, old = c("mutation","copy_number","expression_fusion","Drug_testing"),
                                      new = c('HasMutationData','HasCopynumberData', 'HasFusionData', 'HasDrugTestingData'))

    }

    #output csv file for download
    write.csv(final_df, pdx_filtered_data_loc, row.names=F)

    if(input$PDX_model_var_type=="fusion"){
        splited = stringr::str_split_fixed(fusion_df[['Fusion']], '-',2)
        splited = data.frame(splited)
        splited = all_fac2char(splited)
        final_df[['Fusion']] = link_to_fusion(splited[,1], splited[,2], final_df[[DB$biospec_id_name]], DB$name, final_df[[DB$disease_tag]])
    }

    final_df[[DB$biospec_id_name]] = createDetailLink(final_df[[DB$biospec_id_name]])



    font.size <- "12px"
    dt = DT::datatable(final_df, rownames=FALSE, escape = FALSE, options = list(
        scrollX=TRUE,
        initComplete = htmlwidgets::JS(
            "function(settings, json) {",
            paste0("$(this.api().table().header()).css({'font-size': '", font.size, "'});"),
            "}"),
        lengthMenu = list(c(10, 20, 30), c('10', '20', '30')), pageLength = 10)) %>%
        formatStyle(names(final_df), `text-align` = 'center', fontSize=font.size)
    return(dt)
})


#==========================================================PDX model deatil==========================================================

font.size = "12px"
output$PDX_model_detail_title = renderUI({
    #parsing the query string
    URLString <- isolate(session$clientData$url_search)
    query_id <- str_match(URLString, "sample_id=(.+)")[1,2]
    tags$h2(class="mb-3",
        paste(paste('sample', query_id, 'detail:'))
    )
})
output$PDX_sample_mutations = renderDataTable({

    #parsing the query string
    URLString <- isolate(session$clientData$url_search)
    query_id <- str_match(URLString, "sample_id=(.+)")[1,2]

    #get DB object
    DB = chooseDB("pptc")

    #get mutation data
    mut = QueryAll(DB$mut_info_db)

    #get clin data
    clin = QueryAll(DB$clinical_db_name)

    #check mutation data availability
    flag = check_avaliability(clin, 'mutation', query_id)

    #if the queried gene does not have any fusion event, then output the data file as empty
    if(!flag){
        write.csv(data.frame(), pdx_detail_mutation_data, row.names=F)
    }

    validate(
        need(flag, 'No mutation data for this sample'),
        errorClass = "PDXError"
    )

    #filter the mutation data based on sample_id
    filtered_mut = mut[mut[[DB$biospec_id_name]]==query_id, ]


    feature_selected = c("sample_id","gene_name","mut_group","Variant_Type","HGVSp_Short","Chromosome","Start_Position",
                         "End_Position","Strand","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2","t_alt_count",
                         "n_alt_count","t_ref_count","n_ref_count","Matched_Norm_Sample_Barcode","SWISSPROT")

    filtered_mut = filtered_mut[,feature_selected]

    #if the queried gene does not have any fusion event, then output the data file as empty
    if(length(rownames(filtered_mut))==0){
        write.csv(filtered_mut, pdx_detail_mutation_data, row.names=F)
    }

    validate(
        need(length(rownames(filtered_mut))>0, 'No mutation was detected'),
        errorClass = "PDXError"
    )

    #output csv file for download
    write.csv(filtered_mut, pdx_detail_mutation_data, row.names=F)

    #render the final output
    dt = DT::datatable(filtered_mut, rownames=FALSE, escape = FALSE, options = list(
        scrollX=TRUE,
        initComplete = htmlwidgets::JS(
            "function(settings, json) {",
            paste0("$(this.api().table().header()).css({'font-size': '", font.size, "'});"),
            "}"),
        lengthMenu = list(c(10, 20, 30), c('10', '20', '30')), pageLength = 10)) %>%
        formatStyle(names(filtered_mut), `text-align` = 'center', fontSize=font.size)
    return(dt)
})

output$PDX_sample_fusions = renderDataTable({
    # #parsing the query string
    URLString <- isolate(session$clientData$url_search)
    query_id <- str_match(URLString, "sample_id=(.+)")[1,2]

    #get DB object
    DB = chooseDB("pptc")

    #get fusion data
    fusion = QueryAll(DB$fusion_db)

    #get clinical data
    clin = QueryAll(DB$clinical_db_name)

    #check mutation data availability
    flag = check_avaliability(clin, 'expression_fusion', query_id)

    if(!flag){
        write.csv(data.frame(), pdx_detail_fusion_data, row.names=F)
    }

    validate(
        need(flag, 'No fusion data for this sample'),
        errorClass = "PDXError"
    )

    #filter the mutation data based on sample_id
    filtered_fusion = fusion[fusion[[DB$biospec_id_name]]==query_id, ]

    filtered_fusion = dplyr::inner_join(filtered_fusion, clin, by=DB$biospec_id_name)

    filtered_fusion = filtered_fusion[,c(names(fusion), DB$disease_tag)]

    #remove the first column from the queried results
    filtered_fusion = filtered_fusion[,-1]

    #if the queried gene does not have any fusion event, then output the data file as empty
    if(length(rownames(filtered_fusion))==0){
        write.csv(filtered_fusion, pdx_detail_fusion_data, row.names=F)
    }

    validate(
        need(length(rownames(filtered_fusion))>0, 'No fusion was detected'),
        errorClass="PDXError"
    )

    #output csv file for download
    write.csv(filtered_fusion, pdx_detail_fusion_data, row.names=F)

    filtered_fusion[['Fusion']] = paste(filtered_fusion[['donor_gene']], filtered_fusion[['acceptor_gene']], sep = '--')

    filtered_fusion[['Fusion']] = link_to_fusion(filtered_fusion[['donor_gene']], filtered_fusion[['acceptor_gene']],
                                                 filtered_fusion[['sample_id']], DB$name, filtered_fusion[[DB$disease_tag]]
    )

    #render the final output
    dt = DT::datatable(filtered_fusion, rownames=FALSE, escape = FALSE, options = list(
        scrollX=TRUE,
        initComplete = htmlwidgets::JS(
            "function(settings, json) {",
            paste0("$(this.api().table().header()).css({'font-size': '", font.size, "'});"),
            "}"),
        lengthMenu = list(c(10, 20, 30), c('10', '20', '30')), pageLength = 10)) %>%
        formatStyle(names(filtered_fusion), `text-align` = 'center', fontSize=font.size)
    return(dt)
})

output$PDX_sample_clinical = renderDataTable({
    # #parsing the query string
    URLString <- isolate(session$clientData$url_search)
    query_id <- str_match(URLString, "sample_id=(.+)")[1,2]

    #get DB object
    DB = chooseDB("pptc")

    #get clinical data
    clin = QueryAll(DB$clinical_db_name)

    #filter the mutation data based on sample_id
    filtered_clin = clin[clin[[DB$biospec_id_name]]==query_id, ]

    validate(
        need(length(rownames(filtered_clin))>0, 'no clinical data for the selected sample'),
        errorClass = "PDXError"
    )
    #slice the preferred features
    sample_related_feature = c("sample_id","Site_of_Initial_Tumor","Site_of_Specimen","Stage_of_Disease","Risk_Group",
                               "COG_studies","Prior_Therapy","CellLine_or_Xenograft")
    filtered_clin = filtered_clin[,sample_related_feature]

    #adjust the final output to a preferred style
    filtered_clin = data.frame(t(filtered_clin))
    names(filtered_clin)[1] = 'values'
    filtered_clin = tibble::add_column(filtered_clin, "features"=rownames(filtered_clin), .before = 1)
    filtered_clin[["features"]] = stringr::str_replace_all(filtered_clin[["features"]], '_', ' ')


    #output csv file for download
    write.csv(filtered_clin, pdx_detail_sample_data, row.names=F)

    #render the final output
    dt = DT::datatable(filtered_clin, colnames = '', rownames=FALSE, escape = FALSE, options = list(
        scrollX=TRUE,
        columnDefs = list(list(width = '50%', targets=0)),
        initComplete = htmlwidgets::JS(
            "function(settings, json) {",
            paste0("$(this.api().table().header()).css({'font-size': '", font.size, "'});"),
            "}"),
        lengthMenu = list(c(10, 20, 30), c('10', '20', '30')), pageLength = 10)) %>%
        formatStyle(columns = names(filtered_clin), `text-align` = "left", fontSize=font.size) %>%
        formatStyle(columns = c("features"), fontWeight = 'bold')
    return(dt)

})

output$PDX_patient_clinical = renderDataTable({
    # #parsing the query string
    URLString <- isolate(session$clientData$url_search)
    query_id <- str_match(URLString, "sample_id=(.+)")[1,2]

    #get DB object
    DB = chooseDB("pptc")

    #get clinical data
    clin = QueryAll(DB$clinical_db_name)

    #change values from 0 1 to no and yes
    availability_columns = c("mutation","copy_number","expression_fusion","Drug_testing")
    clin[,availability_columns] <- lapply(clin[,availability_columns], replace_to_YesAndNo)

    #filter the mutation data based on sample_id
    filtered_clin = clin[clin[[DB$biospec_id_name]]==query_id, ]

    validate(
        need(length(rownames(filtered_clin))>0, 'no clinical data for the selected sample'),
        errorClass="ErrorOutput"
    )

    #slice the preferred features
    patient_related_feature = c("Histology","Histology_Detailed","Sex","Phase","Age_year","Reported_Ethnicity",
                                "Inferred_Ethnicity","Drug_testing","Patient_last_alive_year",
                                "Patient_OS_to_last_alive_date_days","Patient_EFS_from_Dx_to_1st_Progression_days",
                                "Time_from_Dx_to_sample_for_model_days")
    filtered_clin = filtered_clin[,patient_related_feature]

    #change drug column name
    filtered_clin = data.table::setnames(filtered_clin, old = c("Drug_testing"), new = c('HasDrugTestingData'))

    #adjust the final output to a preferred style
    filtered_clin = data.frame(t(filtered_clin))
    names(filtered_clin)[1] = 'values'
    filtered_clin = tibble::add_column(filtered_clin, "features"=rownames(filtered_clin), .before = 1)
    filtered_clin[["features"]] = stringr::str_replace_all(filtered_clin[["features"]], '_', ' ')


    #output csv file for download
    write.csv(filtered_clin, pdx_detail_patient_data, row.names=F)

    dt = DT::datatable(filtered_clin, colnames = '', rownames=FALSE, escape = FALSE, options = list(
        scrollX=TRUE,
        columnDefs = list(list(width = '50%', targets=0)),
        initComplete = htmlwidgets::JS(
            "function(settings, json) {",
            paste0("$(this.api().table().header()).css({'font-size': '", font.size, "'});"),
            "}"),
        lengthMenu = list(c(10, 20, 30), c('10', '20', '30')), pageLength = 10)) %>%
        formatStyle(names(filtered_clin), `text-align` = 'left', fontSize=font.size) %>%
        formatStyle(columns = c("features"), fontWeight = 'bold')
    return(dt)
})

#get annotation for response level
getRL <- function (x){
    rlanno = c("PD1", "PD2", "SD", "PR", "CR", "MCR")
    names(rlanno) = c("1", "2", "3", "4", "5", "6")
    return(rlanno[as.character(x)])
}

#function to link to pubmed resource based on PubMed ID
createPubMedLink <- function(val) {
    if(is.na(val)){
        return("")
    }else{
        val_vec = stringr::str_split(val,pattern = ",")[[1]]

        #remove leading and ending spaces for each element
        val_vec = sapply(val_vec, function (x){
            return(stringr::str_trim(x))
        })

        val_vec = sapply(val_vec, function (x){
            link = sprintf('<a href="https://pubmed.ncbi.nlm.nih.gov/%s/" target="_blank">%s</a>',x, x)
            return(link)
        })
        full_link = paste(val_vec, collapse = ",")
        return(full_link)
    }
}

#preclinical output
output$PDX_sample_preclinical <- renderDataTable({
    # #parsing the query string
    URLString <- isolate(session$clientData$url_search)
    query_id <- str_match(URLString, "sample_id=(.+)")[1,2]

    #get DB object
    DB = chooseDB("pptc")

    #get preclinical data
    drug_df = QueryAll(DB$drug_db)

    #filter out the user selected sample data
    drug_df = drug_df[drug_df[[DB$biospec_id_name]]==query_id, ]

    if(length(rownames(drug_df))==0){
        write.csv(data.frame(), pdx_detail_preclinical_data, row.names = F)
    }

    #if the result is empty, return error message
    validate(
        need(length(rownames(drug_df))>0, 'No drug testing data for this sample'),
        errorClass="PDXError"
    )
    #transpose the response data and add col name to the transposed df
    drug_df = data.frame(t(drug_df))
    drug_df = tibble::add_column(drug_df, "Drug_Name"=rownames(drug_df), .before = 1)
    names(drug_df)[2] = 'Response_Level'
    #remove the first row which is indicating sample name
    drug_df = drug_df[-1, ]
    #drop the drug that has response level equals to 0
    drug_df = drug_df[drug_df['Response_Level']!="0",]

    if(length(rownames(drug_df))==0){
        write.csv(data.frame(), pdx_detail_preclinical_data, row.names = F)
    }

    validate(
        need(length(rownames(drug_df))>0, 'all drug showing no response to this sample'),
        errorClass="ErrorOutput"
    )


    Response_Level_annotation = sapply(drug_df['Response_Level'], getRL)
    drug_df = tibble::add_column(drug_df, Response_Level_annotation, .after = 2)

    #get drug annotation
    drug_anno = QueryAll(DB$drug_anno_db)
    #join with drug annotation
    drug_df = dplyr::left_join(drug_df, drug_anno, by="Drug_Name")

    drug_df['Alt_name'] = mapply(comparetwo,drug_df['Drug_Name'], drug_df['Alt_name'])

    drug_df = drug_df[,c("Drug_Name","Alt_name","Response_Level","Response_Level_annotation","Mechanism","PubMed","First_Pub_Date")]



    #output the result data
    write.csv(drug_df, pdx_detail_preclinical_data, row.names = F)

    drug_df[["PubMed"]] = sapply(drug_df[["PubMed"]], createPubMedLink)


    dt = DT::datatable(drug_df, rownames=FALSE, escape = FALSE, options = list(
        scrollX=TRUE,
        columnDefs = list(list(className = 'dt-center', targets = "_all")),
        initComplete = htmlwidgets::JS(
            "function(settings, json) {",
            paste0("$(this.api().table().header()).css({'font-size': '", font.size, "'});"),
            "}"),
        lengthMenu = list(c(10, 20, 30), c('10', '20', '30')), pageLength = 10)) %>%
        formatStyle(names(drug_df), `text-align` = 'center', fontSize=font.size) %>%
        formatStyle(columns = c("Drug_Name"), fontWeight = 'bold') %>%
        formatStyle(names(drug_df)[[5]], `text-align` = 'left', fontSize=font.size)
    return(dt)
})
