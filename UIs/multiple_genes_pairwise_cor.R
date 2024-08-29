# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 10/29/19
library(shiny)
output$subPages <- renderUI(
    tagList(
        shinyjs::useShinyjs(),
        tags$div(class='container',
            tags$div(class='row card mt-3',
                tags$div(class='card-header',
                    tags$h2("Pairwise Correlation")
                ),
                tags$div(class='card-body',
                    fluidRow(
                        column(3,
                            tags$div(class='card card-body bg-light',

                                selectInput('pairwise_db', 'Choose database:', choices=c('target', 'pptc')),
                                selectInput('pairwise_disease', "Choose a disease:", choices=c('pan-cancer'), selected='pan-cancer'),
                                selectInput("pairwise_method","Choose correlation method:",c("pearson","spearman")),
                                uiOutput('pairwise_gene_names_ui'),
                                actionLink(class='mb-3 w-50', 'pairwise_example_genes', 'Example genes'),
                                checkboxInput("pairwise_remove_tissue_effect", "remove tissue effect", FALSE),
                                uiOutput('pairwise_gene_checker'),
                                uiOutput('pairwise_run_checker'),
                                uiOutput("corrplot_annotations")
                            )
                        ),
                        column(9,
                            conditionalPanel(
                                condition="(input.pairwise_run_button) && ($('html').hasClass('shiny-busy'))",
                                 tags$div(class='middle-text', "Loading...",style='color:#0004eb; font-size:2rem; z-index:1;')
                            ),
                            plotOutput('pairwise_heatmap')
                        )
                    )
                ),
                tags$div(class='card-footer text-center',
                    downloadButton("pairwise_download", "Download PDF", class="btn btn-primary"),
                    downloadButton("pairwise_data_download", "Download Data", class="btn btn-success")
                )
            )
        ),
        htmlTemplate("./www/footer.html")
    )
)
