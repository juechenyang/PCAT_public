# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 2019-08-08

library(shiny)
library(plotly)
library(shinyjs)

output$subPages <- renderUI(
    tagList(
        shinyjs::useShinyjs(),
        tags$div(class="container", style='padding:0;',
            tags$div(class='card mt-3',
                tags$div(class='card-header', tags$h3("Expression vs prognosis")),
                tags$div(class='card-body',
                    fluidRow(
                        column(3,
                            tags$div(class="card card-body bg-light",
                                selectInput('survival_db', 'Choose a database:', choices=c(
                                'target')),
                                selectInput('survival_disease', "Choose a disease:", choices=NULL),
                                selectInput('survival_type', 'Choose a survival type:', choices=c('overall', 'eventfree')),
                                selectInput('survival_cutoff_method', "Choose a cutoff method:", choices=c('auto_calculate', 'mean', 'median', 'custom')),
                                sliderInput("survival_threshold", "choose a custom threshold",
                                    min = 1, max = 10, value = 3
                                ),
                                selectizeInput('survival_gene', "input a gene:", choices = NULL, options = list(
                                placeholder = "input a gene, eg: TERT"
                                )),
                                actionLink(class='mb-3 w-50', 'single_gene_survival_example', 'Example'),
                                selectizeInput('survival_multivariate', 'add a co-variate:', choices=NULL, multiple=T,
                                options=list(placeholder='add a co-variate')
                                ),
                                tags$div(class='tlp',
                                    tags$i(class="fas fa-question-circle"),
                                    tags$span("Tips"),
                                    tags$div(class="tooltip_text",
                                        tags$div(class="font-awesome ml-4",
                                            tags$ul(
                                                tags$li(
                                                    "When an input gene is present, patient cohort is divided into
                                                    expression high and low groups by the chosen cutoff method.
                                                    Median/mean is self-explanatory. Auto-calculate tests all possible
                                                    cutoff values between top 20% and bottom 20% patients based on the
                                                    expression of the gene, and adopts the value that best separates
                                                    the high and low. Users may also customize their own cutoff by
                                                    selecting ”custom” in the method pulldown menu."
                                                )
                                            )
                                        )
                                    )
                                )
                            )
                        ),
                        column(9,
                            conditionalPanel(
                                condition="($('html').hasClass('shiny-busy')) && (input.survival_gene)",
                                tags$div(class='middle-text', "Loading...",style='color:#0004eb; font-size:2rem; z-index:1;')
                            ),
                            plotOutput('survival_km'),
                            plotOutput('survival_forest')
                        )
                    )
                ),
                tags$div(class='card-footer text-center',
                    downloadButton("km_plot_download", "Download KM plot", class="btn btn-primary"),
                    downloadButton("km_data_download", "Download KM data", class="btn btn-success mx-3"),
                    downloadButton("forrest_plot_download", "Download Forest plot", class="btn btn-primary mx-3"),
                    downloadButton("forrest_data_download", "Download Forest data", class="btn btn-success")
                )
            )
        ),
        htmlTemplate("./www/footer.html")
    )
)