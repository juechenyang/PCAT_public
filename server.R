library(shiny)
library(stringr)
source('LoadMysql.R')
source('ExpDataTools.R')
source("InputProcessor.R")
source('run_single_cor.r')
ssGSEAScores <- NULL
GENE_LIST <- readRDS(paste('./IntermediateFiles/', 'genes.rds', sep=''))

gene_order <- c("TP53", "PDGFRA","EGFR","TERT","ETV6","BRAF","HRAS","ALK","IDH1","IDH2","KRAS","MLL3",
                "MYH9","NF1","PTEN","ZNF598","ZNF91","TM2D2","NUDT19","FGF12","ZNF829","AKT2","CDKN2A","TFRC")
GENE_LIST <- GENE_LIST[order(match(GENE_LIST, gene_order))]


function(input, output, session){

    options(shiny.maxRequestSize=200*1024^2)
    #initialize UI
    # output$homeUI <- renderUI(tagList(
    #   htmlTemplate("./www/navbar.html"),
    #   uiOutput("subPages")
    # ))
    #define saving locations of plots
    single_gene_plot_loc <- "StaticFiles/Plots/single_gene_plot.pdf"
    single_gene_data_loc <- "StaticFiles/DataFiles/single_gene_data.csv"
    km_plot_loc <- "StaticFiles/Plots/km_plot.pdf"
    km_data_loc <- "StaticFiles/DataFiles/km_plot.csv"
    forrest_plot_loc <- "StaticFiles/Plots/forrest_plot.pdf"
    forrest_data_loc <- "StaticFiles/DataFiles/forrest_data.csv"
    single_gene_cor_plot_loc <-"StaticFiles/Plots/single_gene_cor_plot.pdf"
    single_gene_cor_data_loc <-"StaticFiles/DataFiles/single_gene_cor_data.csv"
    separate_cor_plot_loc <-"StaticFiles/Plots/separate_cor_plot.pdf"
    separate_cor_data_loc <-"StaticFiles/DataFiles/separate_cor_data.csv"
    multiple_genes_plot_loc <- "StaticFiles/Plots/multiple_genes_plot.pdf"
    multiple_genes_data_loc <- "StaticFiles/DataFiles/multiple_genes_data.csv"
    single_preclinical_plot_loc <- "StaticFiles/Plots/single_preclinical_plot.pdf"
    single_preclinical_data_loc <- "StaticFiles/DataFiles/single_preclinical_data.csv"
    ssGSEAPlot_loc <- "StaticFiles/Plots/ssGSEA_plot.pdf"
    ssGSEAData_loc <- "StaticFiles/DataFiles/ssGSEAData.csv"
    correlation_data_loc <-"StaticFiles/DataFiles/top_correlated_genes.csv"
    enrich_r_results_loc <- "StaticFiles/DataFiles/enrich_r.csv"
    single_gene_mutcpn_plot_loc <- "StaticFiles/Plots/single_gene_mut.pdf"
    single_gene_mutcpn_data_loc <- "StaticFiles/DataFiles/single_gene_mutcpn.csv"
    pairwise_plot_loc <- "StaticFiles/Plots/pairwise_plot.pdf"
    pairwise_data_loc <- "StaticFiles/DataFiles/pairwise_data.csv"
    pdx_filtered_data_loc <- "StaticFiles/DataFiles/pdx_filtered_data.csv"
    pdx_detail_patient_data <- "StaticFiles/DataFiles/pdx_detail_patient_data.csv"
    pdx_detail_sample_data <- "StaticFiles/DataFiles/pdx_detail_sample_data.csv"
    pdx_detail_mutation_data <- "StaticFiles/DataFiles/pdx_detail_mutation_data.csv"
    pdx_detail_fusion_data <- "StaticFiles/DataFiles/pdx_detail_fusion_data.csv"
    pdx_detail_preclinical_data <- "StaticFiles/DataFiles/pdx_detail_preclinical_data.csv"
    gene_fusion_data_loc <- "StaticFiles/DataFiles/gene_fusion_data.csv"
    gene_fusion_for_a_sample_data_loc <- "StaticFiles/DataFiles/gene_fusion_for_a_sample_data.csv"
    gene_fusion_for_a_gene_data_loc <- "StaticFiles/DataFiles/gene_fusion_for_a_gene_data.csv"
    gene_fusion_summary_data_loc <- "StaticFiles/DataFiles/gene_fusion_summary_data.csv"


    #===============================================PAGES HANDLER==================================================================================================
    # load server code for page specified in URL
    validFiles = c("./UIs/home.R", './UIs/analysis.R', "./UIs/documentation.R", "./UIs/single_gene_expression.R",
                   './UIs/survival_analysis.R', "./UIs/single_gene_correlation.R",
                   './UIs/separate_enrich_r.R', './UIs/separate_correlation.R', './UIs/single_gene_mut_cpn.R',
                   "./UIs/multiple_genes_clinical.R", './UIs/multiple_genes_enrich_r.R', './UIs/multiple_genes_pairwise_cor.R',
                   "./UIs/multiple_genes_ssGSEA.R", './UIs/single_gene_preclinical.R', './UIs/PDX_model_selector.R',
                   './UIs/PDX_model_detail.R', './UIs/gene_fusions.R', './UIs/gene_fusion_expression_locator.R', './UIs/fusions_of_a_sample.R',
                   './UIs/fusions_of_a_gene.R', './UIs/fusion_summary.R')
    validCode = c('None','None','None',"./single_gene_expression.R",'./survival_analysis.R', "./single_gene_correlation.R",'./separate_enrich_r.R',
                    './separate_correlation.R', './single_gene_mut_cpn.R', './multiple_genes_clinical.R', './multiple_genes_enrich_r.R',
                    './multiple_genes_pairwise_cor.R', './multiple_genes_ssGSEA.R', './single_gene_preclinical.R', './PDX_model_analysis.R',
                    './PDX_model_analysis.R', './gene_fusions.R', './gene_fusions.R', './gene_fusions.R','./gene_fusions.R', './gene_fusions.R'
                    )
    names(validCode) = validFiles

    validSource = c("./single_gene_expression.R",'./survival_analysis.R', "./single_gene_correlation.R",'./separate_enrich_r.R',
                    './separate_correlation.R', './single_gene_mut_cpn.R', './multiple_genes_clinical.R', './multiple_genes_enrich_r.R',
                    './multiple_genes_pairwise_cor.R', './multiple_genes_ssGSEA.R', './single_gene_preclinical.R',
                    './PDX_model_analysis.R', './gene_fusions.R')

                                                        #    names to prevent Unix case problems)
    fname = isolate(session$clientData$url_search)       # isolate() deals with reactive context
    if(nchar(fname)==0) { fname = "?home" }              # blank means home page
    if(grepl('separate_', fname)){
        if(grepl('gene1', fname)){
            fname <- str_split(fname, '\\?gene1')[[1]][1]
        }else{
            fname <- str_split(fname, '\\?genes')[[1]][1]
        }

    }else if(grepl("PDX_model_detail", fname)){
        fname = str_split(fname, '\\?sample_id')[[1]][1]
    }else if(grepl("gene_fusion_expression_locator", fname)){
        fname = str_split(fname, '\\?donor_gene')[[1]][1]
    }else if(grepl("fusions_of_a_sample", fname)){
        fname = str_split(fname, '\\?sample_id')[[1]][1]
    }else if(grepl("fusions_of_a_gene", fname)){
        fname = str_split(fname, '\\?gene_name')[[1]][1]
    }else if(grepl("fusion_summary", fname)){
        fname = str_split(fname, '\\?db_name')[[1]][1]
    }
    fname_ui = paste0('./UIs/', substr(fname, 2, nchar(fname)), ".R") # remove leading "?", add ".R"
    code_fname = validCode[fname_ui]

    cat(paste0("Session filename: ", fname_ui, " current time: ", Sys.time(), "\n"))      # print the URL for this session


    if(!fname_ui %in% validFiles){                          # is that one of our files?
      output$subPages <- renderUI(tagList(              # 404 if no file with that name
         fluidRow(
            column(5,
               HTML("<h2>404 Not Found Error:</h2><p>That URL doesn't exist. Use the",
                    "menu above to navigate to the page you were looking for.</p>")
            )
         )
      ))
      return()    # to prevent a "file not found" error on the next line after a 404 error
    }
    if(code_fname %in% validSource){
        source(code_fname, local = TRUE)$value
    }
    source(fname_ui, local=TRUE)                            # load and run server code for this page



    #==========================================================DOWNLOAD HANDLER=======================================================================================
    #single gene clinical plots downloadHandler
    output$BoxPlotDown<- downloadHandler(
        filename = function() {
            processed_name <- gsub("^\\s+|\\s+$", "", input$geneName)
            paste(Sys.time(), processed_name, input$box_feature, "plot.pdf", sep=".")
        },
        content = function(file) {
            validate(need(file.exists(single_gene_plot_loc), "file not found"))
            file.copy(single_gene_plot_loc, file)
        }
    )
    #single gene clinical data downloadHandler
    output$BoxDataDown<- downloadHandler(
        filename = function() {
            processed_name <- gsub("^\\s+|\\s+$", "", input$geneName)
            paste(Sys.time(), processed_name, input$box_feature, "data.csv", sep=".")
        },
        content = function(file) {
            validate(need(file.exists(single_gene_data_loc), "file not found"))
            file.copy(single_gene_data_loc, file)
        }
    )

    #survival plot download handler
    output$km_plot_download <- downloadHandler(
        filename = function() {
            processed_name <- gsub("^\\s+|\\s+$", "", input$survival_gene)
            paste(Sys.time(), processed_name, "km_plot.pdf", sep=".")
        },
        content = function(file) {
            validate(need(file.exists(km_plot_loc), "file not found"))
            file.copy(km_plot_loc, file)
        }
    )
    output$km_data_download <- downloadHandler(
        filename = function() {
            processed_name <- gsub("^\\s+|\\s+$", "", input$survival_gene)
            paste(Sys.time(), processed_name, "km_data.csv", sep=".")
        },
        content = function(file) {
            validate(need(file.exists(km_data_loc), "file not found"))
            file.copy(km_data_loc, file)
        }
    )
    output$forrest_plot_download <- downloadHandler(
        filename = function() {
            processed_name <- gsub("^\\s+|\\s+$", "", input$survival_gene)
            paste(Sys.time(), processed_name, "forrest_plot.pdf", sep=".")
        },
        content = function(file) {
            validate(need(file.exists(forrest_plot_loc), "file not found"))
            file.copy(forrest_plot_loc, file)
        }
    )
    output$forrest_data_download <- downloadHandler(
        filename = function() {
            processed_name <- gsub("^\\s+|\\s+$", "", input$survival_gene)
            paste(Sys.time(), processed_name, "forrest_data.csv", sep=".")
        },
        content = function(file) {
            validate(need(file.exists(forrest_data_loc), "file not found"))
            file.copy(forrest_data_loc, file)
        }
    )

    output$CorPlotDown <- downloadHandler(
        filename = function() {
            processed_name <- gsub("^\\s+|\\s+$", "", input$gene1, input$gene2)
            paste(Sys.time(), processed_name, "cor_plot.pdf", sep=".")
        },
        content = function(file) {
            validate(need(file.exists(single_gene_cor_plot_loc), "file not found"))
            file.copy(single_gene_cor_plot_loc, file)
        }
    )

    output$CorDataDown <- downloadHandler(
        filename = function() {
            processed_name <- gsub("^\\s+|\\s+$", "", paste0(input$gene1, input$gene2, sep='.'))
            paste(Sys.time(), processed_name, "cor_data.csv", sep=".")
        },
        content = function(file) {
            validate(need(file.exists(single_gene_cor_data_loc), "file not found"))
            file.copy(single_gene_cor_data_loc, file)
        }
    )
    output$separate_cor_plot_download <- downloadHandler(
        filename = function() {
            paste(Sys.time(), "separate_cor_plot.pdf", sep=".")
        },
        content = function(file) {
            validate(need(file.exists(separate_cor_plot_loc), "file not found"))
            file.copy(separate_cor_plot_loc, file)
        }
    )
    output$separate_cor_data_download <- downloadHandler(
        filename = function() {
            paste(Sys.time(), "separate_cor_data.csv", sep=".")
        },
        content = function(file) {
            validate(need(file.exists(separate_cor_data_loc), "file not found"))
            file.copy(separate_cor_data_loc, file)
        }
    )
    output$correlation_table_download <- downloadHandler(
        filename = function() {
            paste(Sys.time(), "correlation.csv", sep=".")
        },
        content = function(file) {
            validate(need(file.exists(correlation_data_loc), "file not found"))
            file.copy(correlation_data_loc, file)
        }
    )
    output$mut_cpn_plot_download <- downloadHandler(
        filename = function() {
            processed_name <- gsub("^\\s+|\\s+$", "", input$single_mutation_gene_name)
            paste(Sys.time(), processed_name, "mut_or_copynumber.pdf", sep=".")
        },
        content = function(file) {
            validate(need(file.exists(single_gene_mutcpn_plot_loc), "file not found"))
            file.copy(single_gene_mutcpn_plot_loc, file)
        }
    )
    output$mut_cpn_data_download <- downloadHandler(
        filename = function() {
            processed_name <- gsub("^\\s+|\\s+$", "", input$single_mutation_gene_name)
            paste(Sys.time(), processed_name, "mut_or_copynumber.csv", sep=".")
        },
        content = function(file) {
            validate(need(file.exists(single_gene_mutcpn_data_loc), "file not found"))
            file.copy(single_gene_mutcpn_data_loc, file)
        }
    )
    output$multigenes_clinical_plot_download<- downloadHandler(
        filename = function() {
            paste(Sys.time(), input$heatmap_feature, "heatmap.pdf", sep=".")
        },
        content = function(file) {
            file.copy(multiple_genes_plot_loc, file)
        }
    )
    output$multigenes_clinical_data_download<- downloadHandler(
        filename = function() {
            paste(Sys.time(), input$heatmap_feature, "heatmap_data.csv", sep=".")
        },
        content = function(file) {
            file.copy(multiple_genes_data_loc, file)
        }
    )
    output$ssGSEAPlotDown<- downloadHandler(
        filename = function() {
            paste(Sys.time(), "ssGSEA_heatmap.pdf", sep=".")
        },
        content = function(file) {
            file.copy(ssGSEAPlot_loc, file)
        }
    )
    output$ssGSEADataDown<- downloadHandler(
        filename = function() {
            paste(Sys.time(), "ssGSEA_raw.csv", sep=".")
        },
        content = function(file) {
            file.copy(ssGSEAData_loc, file)
        }
    )
    output$enrich_r_download <- downloadHandler(
        filename = function() {
            paste(Sys.time(), "enrich_r.csv", sep=".")
        },
        content = function(file) {
            file.copy(enrich_r_results_loc, file)
        }
    )
    output$separate_enrich_r_download <- downloadHandler(
        filename = function() {
            paste(Sys.time(), "separate_enrich_r.csv", sep=".")
        },
        content = function(file) {
            file.copy(enrich_r_results_loc, file)
        }
    )
    output$pairwise_download<- downloadHandler(
        filename = function() {
            processed_name <- gsub("^\\s+|\\s+$", "", input$pairwise_genes)
            paste(Sys.time(), processed_name, "plot.pdf", sep=".")
        },
        content = function(file) {
            validate(need(file.exists(pairwise_plot_loc), "file not found"))
            file.copy(pairwise_plot_loc, file)
        }
    )
    output$pairwise_data_download<- downloadHandler(
        filename = function() {
            processed_name <- gsub("^\\s+|\\s+$", "", input$pairwise_genes)
            paste(Sys.time(), "pairwise_matrix.csv", sep=".")
        },
        content = function(file) {
            validate(need(file.exists(pairwise_plot_loc), "file not found"))
            file.copy(pairwise_data_loc, file)
        }
    )
    output$preclinical_plot_download<- downloadHandler(
        filename = function() {
            processed_name <- gsub("^\\s+|\\s+$", "", input$geneName_single_preclinical)
            paste(Sys.time(), processed_name, "plot.pdf", sep=".")
        },
        content = function(file) {
            validate(need(file.exists(single_preclinical_plot_loc), "file not found"))
            file.copy(single_preclinical_plot_loc, file)
        }
    )
    output$preclinical_data_download<- downloadHandler(
        filename = function() {
            processed_name <- gsub("^\\s+|\\s+$", "", input$geneName_single_preclinical)
            paste(Sys.time(), processed_name, "data.csv", sep=".")
        },
        content = function(file) {
            validate(need(file.exists(single_preclinical_data_loc), "file not found"))
            file.copy(single_preclinical_data_loc, file)
        }
    )
    output$PDX_model_result_download <- downloadHandler(
        filename = function() {
            processed_name <- gsub("^\\s+|\\s+$", "", input$PDX_model_gene)
            paste(Sys.time(), processed_name, "pdx_selected.csv", sep=".")
        },
        content = function(file) {
            validate(need(file.exists(pdx_filtered_data_loc), "file not found"))
            file.copy(pdx_filtered_data_loc, file)
        }
    )
    output$gene_fusion_data_download <- downloadHandler(
        filename = function() {
            donor <- gsub("^\\s+|\\s+$", "", input$gene_fusion_donor_gene)
            acceptor <- gsub("^\\s+|\\s+$", "", input$gene_fusion_acceptor_gene)
            paste(Sys.time(), donor, acceptor, "gene_fusion.csv", sep=".")
        },
        content = function(file) {
            validate(need(file.exists(gene_fusion_data_loc), "file not found"))
            file.copy(gene_fusion_data_loc, file)
        }
    )

    output$gene_fusion_for_a_sample_data_download <- downloadHandler(
        filename = function() {
            paste(Sys.time(), "for_a_sample", "gene_fusion.csv", sep=".")
        },
        content = function(file) {
            validate(need(file.exists(gene_fusion_for_a_sample_data_loc), "file not found"))
            file.copy(gene_fusion_for_a_sample_data_loc, file)
        }
    )

    output$gene_fusion_for_a_gene_data_download <- downloadHandler(
        filename = function() {
            paste(Sys.time(), "for_a_gene", "gene_fusion.csv", sep=".")
        },
        content = function(file) {
            validate(need(file.exists(gene_fusion_for_a_gene_data_loc), "file not found"))
            file.copy(gene_fusion_for_a_gene_data_loc, file)
        }
    )
    output$gene_fusion_summary_data_download <- downloadHandler(
        filename = function() {
            paste(Sys.time(), "summary", "gene_fusion.csv", sep=".")
        },
        content = function(file) {
            validate(need(file.exists(gene_fusion_summary_data_loc), "file not found"))
            file.copy(gene_fusion_summary_data_loc, file)
        }
    )

    output$PDX_sample_detail_download <- downloadHandler(
        filename = function() {
            paste(Sys.time(), "all_detail", "zip", sep=".")
        },
        content = function(file) {
            files = c(pdx_detail_patient_data, pdx_detail_sample_data,
                      pdx_detail_mutation_data, pdx_detail_fusion_data,
                      pdx_detail_preclinical_data)
            system2("zip", args=(paste(file,files,sep=" ")))
        }
    )


}
