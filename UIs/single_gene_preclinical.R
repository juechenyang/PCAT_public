# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 2019-08-22
library(shiny)
library(plotly)
output$subPages <- renderUI(
    tagList(
        shinyjs::useShinyjs(),
        tags$div(class="container", style='padding:0;',
            tags$div(class='card mt-3',
                tags$div(class='card-header', tags$h3("preclinical response:")),
                tags$div(class='card-body',
                    fluidRow(
                        column(3,
                            tags$div(class="card card-body bg-light",
                                selectInput('disease_single_preclinical', "Choose a disease:", choices=NULL),
                                selectizeInput('geneName_single_preclinical', "Choose a gene:", choices = NULL, options = list(
                                placeholder = "Enter gene, eg: EGFR", selected=NULL
                                # , plugins = list('restore_on_backspace')
                                )),
                                checkboxInput('remove_mut_point', "remove mutation indicators"),
                                sliderInput("preclinical_minimum", "minimum size of the box to be shown:",
                                    min = 1, max = 10, value = 1
                                ),
                                selectInput('drugs', 'Choose a drug:', choices=NULL),
                                checkboxInput('drug_log_trans', 'log2 scale', FALSE),
                                uiOutput('drug_mechanism'),
                                tags$div(class='tlp',
                                    tags$i(class="fas fa-question-circle"),
                                    tags$span("Tips"),
                                    tags$div(class="tooltip_text",
                                        tags$div(class="ml-1 my-4 font-size-one", "The annotation for the numbers are:"),
                                        tags$div(class='ml-4 font-size-one',
                                            "
                                            PD1 – Progressive Disease 1", tags$br(),"
                                            PD2 – Progressive Disease 2", tags$br(),"
                                            SD – Stable Disease", tags$br(),"
                                            PR – Partial Response", tags$br(),"
                                            CR – Complete Response", tags$br(),"
                                            MCR – Maintained Complete Response"
                                        )
                                    )
                                )
                            )
                        ),
                        column(9,
                            conditionalPanel(
                                condition="(input.geneName_single_preclinical) && ($('html').hasClass('shiny-busy'))",
                                tags$div(class='middle-text', "Loading...",style='color:#0004eb; font-size:2rem; z-index:1;')
                            ),
                            plotlyOutput('preclinical_boxplot', height="100%")
                        )
                    )
                ),
                tags$div(class='card-footer text-center',
                         downloadButton("preclinical_plot_download", "Download PDF", class="btn btn-primary"),
                         downloadButton("preclinical_data_download", "Download Data", class="btn btn-success")

                )
            )
        ),
        htmlTemplate("./www/footer.html")
    )
)
