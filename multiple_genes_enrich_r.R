# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 12/10/19
#==================================================================Enrich R==================================================================
library(DT)
get_enrich_r_textarea = reactive({
    output$enrich_r_genes_ui = renderUI({
        textAreaInput('enrich_r_genes', label="Gene signature:", height='100px',
        value = "", placeholder="Please input valid gene names")
    })
})

observeEvent(input$enrich_r_example_genes,{
    output$enrich_r_genes_ui = renderUI({
        textAreaInput('enrich_r_genes', label="Gene signature:", height='100px',
        value = "CYP21A2\nCYP11B1\nHSD3B2\nMGARP\nSTAR\nCYP17A1\nNOV\nSULT2A1\nATP4A\nCYP11A1\nFDXR\nSOAT1\nABCB1\nSCARB1\nCYP11B2\nINHA\nHOXA5\nFDX1\nTM7SF2\nTBX3\nAOX1\nSLC16A9\nDLK1\nMC2R\nNR5A1", placeholder="Please type at least 2 valid genes to run the analysis and use single space or line breaker to seperate")
    })
})

process_df2dt <- function(df){

    df <- df[df[['P-value']]<0.1,]
    number_column <- c('P-value', 'Z-score', 'Combined score', 'Adjusted p-value')
    df[, number_column] <- sapply(df[, number_column], format, scientific=TRUE, digits=2)
    df[['Overlapping genes']] <- sapply(df[['Overlapping genes']], paste, collapse=' ')
    write.csv(df, enrich_r_results_loc, row.names=F)

    font.size <- "11px"
    dt <- DT::datatable(df, rownames=FALSE, escape = FALSE, options = list(
        initComplete = htmlwidgets::JS(
            "function(settings, json) {",
            paste0("$(this.api().table().header()).css({'font-size': '", font.size, "'});"),
            "}"),
        lengthMenu = list(c(20, 30, 40), c('20', '30', '40')), pageLength = 20)) %>%
        formatStyle(names(df), `text-align` = 'center', fontSize=font.size) %>%
        formatStyle("Overlapping genes","white-space"="wrap")

    return(dt)
}
output$enrich_r_table_results <- renderDataTable({
    source("Enrich_R.R")

    #hide download button
    shinyjs::hide("enrich_r_download")

    get_enrich_r_textarea()

    #fetch genes from user input test box
    gene_list <- multiple_gene_processor(input$enrich_r_genes)

    #validate if user has input any value
    validate(
        need(length(gene_list)>0, 'please input at least one genes'),
        errorClass="ValidateInput"
    )

    #make the input gene as enrich r accepted format
    gene_list <- paste(gene_list, collapse = '\n')
    df <- get_enrich_r_result(gene_list, input$enrich_r_db)

    #transfer the dataframe into data table
    dt <- process_df2dt(df)

    #show download button
    shinyjs::show("enrich_r_download")

    return(dt)
})

