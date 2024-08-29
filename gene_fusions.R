# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 12/13/19
library(DT)
GENE_LIST = sapply(GENE_LIST, as.character)
update_donor_gene = reactive({
    #update donor gene selection
    updateSelectizeInput(session, 'gene_fusion_donor_gene', server = TRUE, choices=c("Any gene", GENE_LIST), selected='Any gene')
})

update_acceptor_gene = reactive({
    #update donor gene selection
    updateSelectizeInput(session, 'gene_fusion_acceptor_gene', server = TRUE, choices=c("Any gene", GENE_LIST), selected='Any gene')
})


#create link to get the detail expression location
createFusionLink <- function(donor_gene, acceptor_gene, sample_id, db_name, cohort) {
    sprintf('<a href="?gene_fusion_expression_locator?donor_gene=%sacceptor_gene=%ssample_id=%sdb_name=%scohort=%s" target="_blank" class="btn btn-primary" style="font-size:0.8rem">View</a>', donor_gene, acceptor_gene, sample_id, db_name, cohort)
}
#create sample fusions link
createSampleFusions <- function(sample_id, db_name, cohort){
    sprintf('<a href="?fusions_of_a_sample?sample_id=%sdb_name=%scohort=%s" target="_blank" style="font-size:0.8rem">%s</a>', sample_id, db_name, cohort, sample_id)
}

#create gene fusions link
createGeneFusions <- function(gene_name, db_name){
    sprintf('<a href="?fusions_of_a_gene?gene_name=%sdb_name=%s" target="_blank" style="font-size:0.8rem">%s</a>', gene_name, db_name, gene_name)
}

#create gene fusions link
createFusionSummary <- function(db_name, fusion_pair){
    sprintf('<a href="?fusion_summary?db_name=%sfusion_pair=%s" target="_blank">%s</a>', db_name, fusion_pair, fusion_pair)
}

output$all_fusions = renderDataTable({

    #update gene selection
    update_donor_gene()
    update_acceptor_gene()

    #get DB object
    DB = chooseDB(input$gene_fusion_db_name)

    #check if the user has input both donor gene and acceptor gene
    validate(
        need(input$gene_fusion_donor_gene!="", 'please input donor gene'),
        errorClass="ValidateInput"
    )
    validate(
        need(input$gene_fusion_acceptor_gene!="", 'please input acceptor gene'),
        errorClass="ValidateInput"
    )

    #get fusion data frame
    fusion_df = QueryAll(DB$fusion_db)

    #merge fusion data with clinical data
    clin = QueryAll(DB$clinical_db_name)
    fus_clin = merge(fusion_df, clin, by=DB$biospec_id_name)

    fusion_df = fus_clin[,c(names(fusion_df), DB$disease_tag)]

    if(input$gene_fusion_donor_gene!="Any gene"){
        fusion_df = fusion_df[fusion_df[[DB$donor_name]]==input$gene_fusion_donor_gene, ]
    }
    if(input$gene_fusion_acceptor_gene!="Any gene"){
        fusion_df = fusion_df[fusion_df[[DB$acceptor_name]]==input$gene_fusion_acceptor_gene, ]
    }

    validate(
        need(length(rownames(fusion_df))>0, 'there is no fusion for selected genes'),
        errorClass="ErrorOutput"
    )


    if(DB$name=="target"){

        #get expect columns
        fusion_df = fusion_df[, c('sample_id', DB$disease_tag, 'donor_gene', 'acceptor_gene', 'donor_chr', 'acceptor_chr', 'LeftBreak', 'RightBreak', 'method')]

    }else if(DB$name=="pptc"){
        fusion_df = fusion_df[,-1]
        fusion_df = fusion_df[,!(names(fusion_df) %in% c('Fusion'))]
        Histology = fusion_df[[DB$disease_tag]]
        fusion_df = fusion_df[,-length(names(fusion_df))]
        fusion_df = tibble::add_column(fusion_df, Histology, .after=1)
    }

    #create fusion column for
    Fusion = paste(fusion_df$donor_gene, fusion_df$acceptor_gene, sep = '--')
    fusion_df = tibble::add_column(fusion_df, Fusion, .after = 4)



    write.csv(fusion_df, gene_fusion_data_loc, row.names = F)

    #create detail link function
    detail = createFusionLink(fusion_df[[DB$donor_name]], fusion_df[[DB$acceptor_name]], fusion_df[[DB$biospec_id_name]], input$gene_fusion_db_name, fusion_df[[DB$disease_tag]])
    fusion_df[[DB$biospec_id_name]] = createSampleFusions(fusion_df[[DB$biospec_id_name]], DB$name, fusion_df[[DB$disease_tag]])
    fusion_df[[DB$donor_name]] = createGeneFusions(fusion_df[[DB$donor_name]], DB$name)
    fusion_df[[DB$acceptor_name]] = createGeneFusions(fusion_df[[DB$acceptor_name]], DB$name)
    fusion_df[['Fusion']] = createFusionSummary(DB$name, fusion_df[['Fusion']])

    #put the link button after first column
    fusion_df = tibble::add_column(fusion_df, detail, .after=1)

    #render the final output
    font.size <- "12px"
    dt = DT::datatable(fusion_df, rownames=FALSE, escape = FALSE, options = list(
        scrollX=TRUE,
        columnDefs = list(list(className = 'dt-center', targets = "_all")),
        initComplete = htmlwidgets::JS(
            "function(settings, json) {",
            paste0("$(this.api().table().header()).css({'font-size': '", font.size, "'});"),
            "}"),
        lengthMenu = list(c(10, 20, 30), c('10', '20', '30')), pageLength = 10)) %>%
        formatStyle(names(fusion_df), `text-align` = 'center', fontSize=font.size)
    return(dt)
})

output$circos_plot = renderPlot({

    source("CirclosPlot.R")

    #get sample
    URLString <- isolate(session$clientData$url_search)
    gene1 <- str_match(URLString, "donor_gene=(.+)acceptor_gene")[1,2]
    gene2 <- str_match(URLString, "acceptor_gene=(.+)sample_id")[1,2]
    sample = str_match(URLString, "sample_id=(.+)db_name")[1,2]
    db_name =str_match(URLString, "db_name=(.+)cohort")[1,2]

    #get database
    DB = chooseDB(db_name)

    #get all fusions
    fusion_df = QueryAll(DB$fusion_db)


    #get fusions from a single sample
    fusions_of_a_sample = fusion_df[which(fusion_df[[DB$biospec_id_name]]==sample),]


    #create circos plot of a sample
    final_plot = createCircos(fusions_of_a_sample, gene1, gene2)


    #return the plot object
    return(final_plot)
})

output$fusion_expression_plot = plotly::renderPlotly({

    source('ScatterPlot.R')

    #get genes
    URLString <- isolate(session$clientData$url_search)
    gene1 <- str_match(URLString, "donor_gene=(.+)acceptor_gene")[1,2]
    gene2 <- str_match(URLString, "acceptor_gene=(.+)sample_id")[1,2]
    sample = str_match(URLString, "sample_id=(.+)db_name")[1,2]
    db_name =str_match(URLString, "db_name=(.+)cohort")[1,2]
    cohort =str_match(URLString, "cohort=(.+)")[1,2]
    cohort = URLdecode(cohort)

    #get database
    DB = chooseDB(db_name)

    #transfer the gene name into a code friendly name
    g1 = single_gene_special_transfer(gene1)
    g2 = single_gene_special_transfer(gene2)

    #get the expression data
    exp_data <-process_exp_data(DB, c(g1,g2), mode='multiple')

    #get all fusions
    fusion_df = QueryAll(DB$fusion_db)

    #merge fusions with clinical
    clin = QueryAll(DB$clinical_db_name)

    #filter the clinical data with specified cohort
    clin = clin[which(clin[[DB$disease_tag]]==cohort),]

    #merge fusion data with clinical data
    fus_clin = merge(fusion_df, clin, by=DB$biospec_id_name)

    #merge fusion clinical with expression data
    fus_clin_exp = dplyr::left_join(fus_clin, exp_data, by=DB$biospec_id_name)

    #filter out the desired fusions
    sample_fusion_df = fus_clin_exp[which(fus_clin_exp[[DB$biospec_id_name]]==sample & fus_clin_exp[[DB$donor_name]]==gene1 & fus_clin_exp[[DB$acceptor_name]]==gene2),]


    #create detail table based on different data source
    if(DB$name=='target'){
        detail_table = data.frame(
        "Variables" = c('sample id', 'disease','donor gene', 'acceptor gene', 'donor chromosome', 'acceptor chromosome', 'donor break', 'acceptor break', 'method',
        '#junction spanning reads', '#discordant reads', '#Consistent reads', 'Position Consistency', 'Percentage Fusion Index', 'Sequence Similarity'),
        "Values" = c(sample, sample_fusion_df[1,DB$disease_tag], gene1, gene2, sample_fusion_df[1,'donor_chr'],sample_fusion_df[1,'acceptor_chr'], sample_fusion_df[1,'LeftBreak'],sample_fusion_df[1,'RightBreak'],
                    sample_fusion_df[1,'method'], sample_fusion_df[1,'number_junctionSpanning_reads'],sample_fusion_df[1,'number_discordant_read_pairs'],sample_fusion_df[1,'number_Consistent'],
                    sample_fusion_df[1,'position_consistency'],sample_fusion_df[1,'Percentage.Fusion.Index.PFI.'], sample_fusion_df[1,'sequence.similarity']),
        stringsAsFactors = FALSE)
    }else if(DB$name=='pptc'){
        detail_table = data.frame(
        "Variables" = c('sample id', 'donor gene', 'acceptor gene', 'Frame', 'Histology', 'Method'),
        "Values" = c(sample, gene1, gene2, sample_fusion_df[1,'Frame'],sample_fusion_df[1,'Histology'], sample_fusion_df[1,'Method']),
        stringsAsFactors = FALSE)
    }

    #render the detail table
    output$sample_detail = renderTable(detail_table, bordered = T, align = c("lc"), rownames=FALSE)

    #validate if the gene is available in the selected DB
    validate(
        need(g1 %in% colnames(exp_data), paste(gene1, 'does not have expression data is this DB')),
        errorClass="ErrorOutput"
    )
    validate(
        need(g2 %in% colnames(exp_data), paste(gene2, 'does not have expression data is this DB')),
        errorClass="ErrorOutput"
    )

    validate(
        need(!is.na(sample_fusion_df[1, g1]), 'sample does not have expression data'),
        errorClass="ErrorOutput"
    )

    #plot the expression scatter based on the genes that contribute this fusion
    final_plot = fusionScatter(fus_clin_exp, g1, g2, sample, cohort, DB, log_scale=input$log_scale)

    #return the plot
    return(final_plot)
})
output$fusions_of_a_sample = renderDataTable({

    URLString <- isolate(session$clientData$url_search)
    sample = str_match(URLString, "sample_id=(.+)db_name")[1,2]
    db_name =str_match(URLString, "db_name=(.+)cohort")[1,2]
    cohort = str_match(URLString, "cohort=(.+)")[1,2]

    #get database
    DB = chooseDB(db_name)

    #get all fusions
    fusion_df = QueryAll(DB$fusion_db)



    #get fusions from a single sample
    fusions_of_a_sample = fusion_df[which(fusion_df[[DB$biospec_id_name]]==sample),]

    #get clin
    clin = QueryAll(DB$clinical_db_name)

    fusions_of_a_sample = dplyr::left_join(fusions_of_a_sample, clin, by=DB$biospec_id_name)



    if(DB$name=="target"){

        #get expect columns
        fusions_of_a_sample = fusions_of_a_sample[, c('sample_id', DB$disease_tag, 'donor_gene', 'acceptor_gene', 'donor_chr', 'acceptor_chr', 'LeftBreak', 'RightBreak', 'method')]

    }else if(DB$name=="pptc"){
        selected_columns = c(DB$biospec_id_name, DB$disease_tag, names(fusion_df)[!(names(fusion_df) %in% c('Fusion', 'gene_name', 'sample_id'))])
        print(selected_columns)
        fusions_of_a_sample = fusions_of_a_sample[,selected_columns]

    }



    output$sample_name = renderUI({
        tags$h3(
            paste0('All fusions of sample: ', sample)
        )
    })

    #output the downloadable data
    write.csv(fusions_of_a_sample, gene_fusion_for_a_sample_data_loc, row.names = F)
    #create detail link function
    detail = createFusionLink(fusions_of_a_sample[[DB$donor_name]], fusions_of_a_sample[[DB$acceptor_name]], fusions_of_a_sample[[DB$biospec_id_name]], db_name, cohort)

    #put the link button after first column
    fusions_of_a_sample = tibble::add_column(fusions_of_a_sample, detail, .after=1)

    #render the final output
    font.size <- "12px"
    dt = DT::datatable(fusions_of_a_sample, rownames=FALSE, escape = FALSE, options = list(
        scrollX=TRUE,
        initComplete = htmlwidgets::JS(
            "function(settings, json) {",
            paste0("$(this.api().table().header()).css({'font-size': '", font.size, "'});"),
            "}"),
        lengthMenu = list(c(10, 20, 30), c('10', '20', '30')), pageLength = 10)) %>%
        formatStyle(names(fusions_of_a_sample), `text-align` = 'center', fontSize=font.size)

    return(dt)
})
output$fusions_of_a_gene = renderDataTable({
    URLString <- isolate(session$clientData$url_search)
    selected_gene = str_match(URLString, "gene_name=(.+)db_name")[1,2]
    db_name =str_match(URLString, "db_name=(.+)")[1,2]

    #get database
    DB = chooseDB(db_name)

    #get all fusions
    fusion_df = QueryAll(DB$fusion_db)

    #get clinical data
    clin = QueryAll(DB$clinical_db_name)

    #merge fusion data with clinical
    fus_clin = merge(fusion_df, clin, by=DB$biospec_id_name)

    #limit the columns to all columns from fusion but only disease from clinical
    fus_clin = fus_clin[,c(DB$biospec_id_name, DB$disease_tag, names(fusion_df)[names(fusion_df)!=DB$biospec_id_name])]

    #get fusions from a single gene
    fusions_of_a_gene = fus_clin[which(fus_clin[[DB$donor_name]]==selected_gene | fus_clin[[DB$acceptor_name]]==selected_gene),]

    #fill the ui output with gene name
    output$gene_name = renderUI({
        tags$h3(
            paste0('All fusions of gene: ', selected_gene)
        )
    })

    #output the downloadable data
    write.csv(fusions_of_a_gene, gene_fusion_for_a_gene_data_loc, row.names = F)

    #create detail link function
    detail = createFusionLink(fusions_of_a_gene[[DB$donor_name]], fusions_of_a_gene[[DB$acceptor_name]], fusions_of_a_gene[[DB$biospec_id_name]], db_name, fusions_of_a_gene[[DB$disease_tag]])


    #put the link button after first column
    fusions_of_a_gene = tibble::add_column(fusions_of_a_gene, detail, .after=1)

    #render the final output
    font.size <- "12px"
    dt = DT::datatable(fusions_of_a_gene, rownames=FALSE, escape = FALSE, options = list(
        scrollX=TRUE,
        initComplete = htmlwidgets::JS(
            "function(settings, json) {",
            paste0("$(this.api().table().header()).css({'font-size': '", font.size, "'});"),
            "}"),
        lengthMenu = list(c(10, 20, 30), c('10', '20', '30')), pageLength = 10)) %>%
        formatStyle(names(fusions_of_a_gene), `text-align` = 'center', fontSize=font.size)
    return(dt)
})



#convert
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

#fusion summary output
output$fusions_summary = renderDataTable({

    URLString <- isolate(session$clientData$url_search)
    db_name <- str_match(URLString, "db_name=(.+)fusion_pair")[1,2]
    fusion_pair <- str_match(URLString, "fusion_pair=(.+)")[1,2]

    output$fusion_pair_name = renderUI({
        tags$h3(
            paste0('Detail Summary of Fusion: ', fusion_pair)
        )
    })

    #get DB obj
    DB = chooseDB(db_name)

    #get all fusions of the DB
    fusion_df = QueryAll(DB$fusion_db)

    #get clinical table
    clin = QueryAll(DB$clinical_db_name)

    #merge clin with fusion
    fus_clin = dplyr::inner_join(fusion_df, clin, by=DB$biospec_id_name)

    #select desired columns
    fus_clin = fus_clin[,c(DB$biospec_id_name, DB$disease_tag, names(fusion_df)[names(fusion_df)!=DB$biospec_id_name])]

    #create fusion pair column
    Fusion_Pair = paste(fus_clin$donor_gene, fus_clin$acceptor_gene, sep = "--")
    fus_clin = tibble::add_column(fus_clin, Fusion_Pair, .after = 4)
    fus_clin = fus_clin[fus_clin[["Fusion_Pair"]]==fusion_pair, ]
    fusion_clin_summary = fus_clin %>% dplyr::group_by_at(DB$disease_tag) %>% dplyr::summarize(Fusion_count = n())
    fusion_clin_summary = rbind(fusion_clin_summary, c("Pan-can", length(rownames(fus_clin))))

    #get fusion search table
    fus_search = QueryAll(DB$fusion_search_db)
    sample_clin_summary = fus_search %>% dplyr::group_by_at(DB$disease_tag) %>% dplyr::summarize(Sample_count = n())
    sample_clin_summary = rbind(sample_clin_summary, c("Pan-can", length(rownames(fus_search))))

    final_df = dplyr::left_join(sample_clin_summary, fusion_clin_summary, by=DB$disease_tag)

    final_df = data.frame(final_df)
    final_df[is.na(final_df)] <- 0
    cols.num <- c("Sample_count","Fusion_count")
    final_df[cols.num] <- sapply(final_df[cols.num],as.numeric)
    final_df[['Rate']] = final_df$Fusion_count/final_df$Sample_count
    final_df$Rate = sapply(final_df$Rate, percent)

    #sort the sample count by decreasing order
    final_df = final_df[order(final_df$Sample_count, decreasing = T),]

    #output the summary data for gene fusion
    write.csv(final_df, gene_fusion_summary_data_loc, row.names = F)

    font.size <- "12px"
    dt = DT::datatable(final_df, rownames=FALSE, escape = FALSE, options = list(
        scrollX=TRUE,
        columnDefs = list(list(className = 'dt-center', targets = "_all")),
        initComplete = htmlwidgets::JS(
            "function(settings, json) {",
            paste0("$(this.api().table().header()).css({'font-size': '", font.size, "'});"),
            "}"),
        lengthMenu = list(c(10, 20, 30), c('10', '20', '30')), pageLength = 10)) %>%
        formatStyle(names(final_df), `text-align` = 'center', fontSize=font.size)
    return(dt)

})
