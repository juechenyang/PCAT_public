# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 2019-08-21
library(shiny)
output$subPages <- renderUI(
    tagList(
        shinyjs::useShinyjs(),
        tags$div(class="container", style='padding:0;',
            tags$div(class='card mt-3',
                tags$div(class='card-header', tags$h3("Functional enrichment analysis through enrichr")),
                tags$div(class='card-body',
                    fluidRow(
                        column(3,
                            tags$div(class="card card-body bg-light",
                                selectizeInput('enrich_r_db', 'Choose a database:', choices=c('Reactome_2016','KEGG_2019_Human',
                                'WikiPathways_2019_Human', 'GO_Biological_Process_2018','GO_Cellular_Component_2018','GO_Molecular_Function_2018',
                                'BioCarta_2016', 'MSigDB_Computational', 'MSigDB_Oncogenic_Signatures','Cancer_Cell_Line_Encyclopedia','NCI-60_Cancer_Cell_Lines','Achilles_fitness_decrease',
                                'Achilles_fitness_increase','DSigDB','DepMap_WG_CRISPR_Screens_Broad_CellLines_2019','DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019',
                                'TargetScan_microRNA_2017','miRTarBase_2017','ChEA_2016','Chromosome_Location','ENCODE_Histone_Modifications_2015',
                                'ENCODE_TF_ChIP-seq_2015','Epigenomics_Roadmap_HM_ChIP-seq','Transcription_Factor_PPIs','ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X',
                                'GTEx_Tissue_Sample_Gene_Expression_Profiles_up')),
                                uiOutput('enrich_r_genes_ui'),
                                tags$div(class="tlp w-75",
                                    actionLink(class='mb-3 font-size-two d-inline','enrich_r_example_genes', 'Example genes'),
                                    tags$div(class="p-3 font-size-two tooltip_text",
                                        "Example genes are refered by: ",tags$a(href="https://www.ncbi.nlm.nih.gov/pubmed/27165744", target="_blank", "zheng et al. cancer cell")
                                    )
                                )
                            )
                        ),
                        column(9,
                            conditionalPanel(
                                condition="$('html').hasClass('shiny-busy')",
                                tags$div(class='middle-text', "Loading...",style='color:#0004eb; font-size:2rem; z-index:1;')
                            ),
                            dataTableOutput('enrich_r_table_results')
                        )
                    )
                ),
                tags$div(class='card-footer text-center',
                     downloadButton('enrich_r_download', class='btn btn-primary')
                )
            )
        ),
        htmlTemplate("./www/footer.html")
    )
)

