# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 12/4/19
library(shiny)
output$subPages <- renderUI(
    tagList(
        tags$div(class="container", style='padding:0;',
            tags$div(class='card mt-3',
                tags$div(class='card-header', tags$h3("PDX search:")),
                tags$div(class='card-body', style='height:100%;',
                    fluidRow(
                        column(3,
                            tags$div(class="card card-body bg-light",
                                     shinyjs::useShinyjs(),
                                      selectInput('PDX_model_var_type', 'Choose a variant:', choices=c('clinical', 'mutation', 'fusion', 'drug')),
                                     selectInput('PDX_model_disease', 'Choose a disease:', choices=c('pan-cancer')),
                                     selectizeInput('PDX_model_gene', "Choose a gene:", choices = NULL
                                      #options = list(placeholder = "Enter gene, eg: EGFR", selected=NULL, plugins = list('restore_on_backspace'))
                                     ),
                                     uiOutput('pdx_drug_mechanism'),
                                tags$div(class='tlp',
                                    tags$i(class="fas fa-question-circle"),
                                    tags$span("Tips"),
                                    tags$div(class="tooltip_text",
                                        tags$ul(class="font-awesome ml-4",
                                            tags$li(
                                                "Scroll the page right and left to see more columns"
                                            ),
                                            tags$li(
                                                "Click model name for details"
                                            ),
                                            tags$li(
                                                "Use search box on the top right corner to narrow down results"
                                            )
                                        )
                                    )
                                )
                            )
                        ),
                        column(9,
                            conditionalPanel(
                                condition="$('html').hasClass('shiny-busy')",
                                 tags$div(class='middle-text', "Loading...",style='color:#0004eb; font-size:2rem; z-index:1;')
                            ),
                            dataTableOutput('PDX_model_results')
                        )
                    )
                ),
                tags$div(class='card-footer text-center',
                     downloadButton("PDX_model_result_download", "Download Data", class="btn btn-success")
                )
            )
        ),
        htmlTemplate("./www/footer.html")
    )
)