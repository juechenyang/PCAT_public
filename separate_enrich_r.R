# Title     : TODO
# Objective : TODO
# Created by: Juechen
# Created on: 2020/4/22
#==================================================================Seperate Enrich R==================================================================
library(DT)
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

output$seperate_enrich_r_results <- renderDataTable({

  source("Enrich_R.R")

  URLString <- isolate(session$clientData$url_search)
  genes <- str_match(URLString, "genes=(.+)")[1,2]
  gene_list <- str_split(genes, '&')[[1]]
  gene_list <- paste(gene_list, collapse = '\n')
  df <- get_enrich_r_result(gene_list, input$separate_enrich_r_db)
  dt <- process_df2dt(df)
  return(dt)
})

