# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 2019-07-14
library(shiny)
output$subPages <- renderUI(
    tagList(
        shinyjs::useShinyjs(),
        tags$div(class="container", style='padding:0;',
            tags$div(class='card mt-3',
                tags$div(class='card-header', tags$h3("heatmap visualization")),
                tags$div(class='card-body',
                    fluidRow(
                        column(3,
                            tags$div(class="card card-body bg-light",
                                selectInput('db_multiple_name',
                                            'Choose a database:', choices=c('pptc', 'target')),
                                selectizeInput('multiple_gene_disease', "Choose a disease:", choices=NULL),
                                uiOutput('gene_names_ui'),
                                actionLink(class='mb-3 w-50', 'example_genes', 'Example genes'),
                                uiOutput('error_gene'),
                                checkboxInput("rowCluster", "Genes clustering", FALSE),
                                checkboxInput("columnCluster", "Samples clustering", FALSE),
                                selectInput('heatmap_feature', 'Overlay a feature:', choices=c(
                                'disease', 'age_class', 'gender')),
                                uiOutput('discretize_number'),
                                uiOutput('multiple_run_check')
                                # tags$div(class='mt-4',
                                #     tags$div(class="mb-4",
                                #         "Example genes:"
                                #     ),
                                #     'tp53 egfr kras tert pdgfra pten snai2'
                                # )
                            )
                        ),
                        column(9,
                            conditionalPanel(
                                condition="(input.analyze_multiple || input.columnCluster) && ($('html').hasClass('shiny-busy'))",
                                # condition="$('html').hasClass('shiny-busy')",
                                tags$div(class='middle-text', "Loading...",style='color:#0004eb; font-size:2rem; z-index:1;')
                            ),
                            plotOutput('geneHeatMap', height='600px')
                        )
                    )
                ),
                tags$div(class='card-footer text-center',
                         downloadButton("multigenes_clinical_plot_download", "Download PDF", class="btn btn-primary"),
                         downloadButton("multigenes_clinical_data_download", "Download Data", class="btn btn-success")
                )
            )
        ),
        htmlTemplate("./www/footer.html")
    )
)
