# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 1/27/20
output$subPages <- renderUI(
    tagList(
        shinyjs::useShinyjs(),
        tags$div(class="container", style='padding:0;',
            tags$div(class='card mt-3',
                tags$div(class='card-header',
                    uiOutput('gene_name')
                ),
                tags$div(class='card-body', style='height:100%;',
                    dataTableOutput('fusions_of_a_gene')
                ),
                tags$div(class='card-footer text-center',
                    downloadButton("gene_fusion_for_a_gene_data_download", "Download Data", class="btn btn-success")
                )
            )
        )
    )
)
