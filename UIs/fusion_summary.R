# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 6/17/20
output$subPages <- renderUI(
    tagList(
        shinyjs::useShinyjs(),
        tags$div(class="container", style='padding:0;',
            tags$div(class='card mt-3',
                tags$div(class='card-header',
                    uiOutput('fusion_pair_name')
                ),
                tags$div(class='card-body', style='height:100%;',
                    dataTableOutput('fusions_summary')
                ),
                tags$div(class='card-footer text-center',
                    downloadButton("gene_fusion_summary_data_download", "Download Data", class="btn btn-success")
                )
            )
        )
    )
)
