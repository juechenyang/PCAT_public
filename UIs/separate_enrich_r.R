# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 2019-09-04
library(shiny)
output$subPages <- renderUI(
    tagList(
        tags$div(class="container", style='padding:0;',
            tags$div(class='card mt-3',
                tags$div(class='card-header', tags$h3("Results:")),
                tags$div(class='card-body',
                    fluidRow(
                        column(3,
                            tags$div(class="card card-body bg-light",
                                selectizeInput('separate_enrich_r_db', 'Choose database:', choices=c('WikiPathways_2019_Human','KEGG_2019_Human',
                                'Reactome_2016','MSigDB_Computational','GO_Biological_Process_2018','GO_Cellular_Component_2018','GO_Molecular_Function_2018',
                                'BioCarta_2016','MSigDB_Oncogenic_Signatures','Cancer_Cell_Line_Encyclopedia','NCI-60_Cancer_Cell_Lines','Achilles_fitness_decrease',
                                'Achilles_fitness_increase','DSigDB','DepMap_WG_CRISPR_Screens_Broad_CellLines_2019','DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019',
                                'TargetScan_microRNA_2017','miRTarBase_2017','ChEA_2016','Chromosome_Location','ENCODE_Histone_Modifications_2015',
                                'ENCODE_TF_ChIP-seq_2015','Epigenomics_Roadmap_HM_ChIP-seq','Transcription_Factor_PPIs','ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X',
                                'GTEx_Tissue_Sample_Gene_Expression_Profiles_up'))
                            )
                        ),
                        column(9,
                            conditionalPanel(
                                condition="$('html').hasClass('shiny-busy')",
                                tags$div(class='middle-text', "Loading...",style='color:#0004eb; font-size:2rem; z-index:1;')
                            ),
                            dataTableOutput('seperate_enrich_r_results')
                        )
                    )
                ),
                tags$div(class='card-footer text-center',
                     downloadButton('separate_enrich_r_download', class='btn btn-primary')
                )
            )
        ),
        htmlTemplate("./www/footer.html")
    )
)
