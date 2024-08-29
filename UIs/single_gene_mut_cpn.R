# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 10/8/19
library(shiny)
library(plotly)
library(shinyjs)
output$subPages <- renderUI(
    tagList(
        shinyjs::useShinyjs(),
        tags$div(class="container", style='padding:0;',
            tags$div(class='card mt-3 mb-3',
                tags$div(class='card-header', tags$h3("Expression vs genetic change")),
                tags$div(class='card-body',
                    fluidRow(
                        column(3,
                            tags$div(class="card card-body bg-light",
                                selectInput('single_mutation_db', 'Choose a database:', choices=c('pptc', 'target')),
                                selectInput('single_mutation_disease', "Choose a disease:", choices=c('pan-cancer'), selected='pan-cancer'),
                                selectInput('mut_or_cpn', "Mutation or Copy number:", choices=c('copy number', 'mutation')),
                                selectizeInput('single_mutation_gene_name', "Choose a gene:", choices = NULL, options = list(
                                placeholder = "Enter gene, eg: CDKN2A", selected=NULL
                                )),
                                actionLink(class='mb-3 w-50', 'single_mut_cpn_example_gene', 'Example'),
                                sliderInput("mut_min_box_size", "minimum datapoints for each group:",
                                    min = 1, max = 10, value = 1
                                ),
                                checkboxInput('mut_log_trans', 'log scale', FALSE),
                                tags$div(class='tlp mt-4',
                                    tags$i(class="fas fa-question-circle"),
                                    tags$span("Tips"),
                                    tags$div(class="tooltip_text",
                                        tags$div(class="font-awesome mx-4",
                                            tags$ul(
                                                tags$li(
                                                    "Boxplot on the right shows expression of the input gene stratified
                                                    by mutational status. If pan-cancer data is used, mutation cases
                                                    will be displayed in a different color in individual cohorts.
                                                    Click a specific mutation group on figure legend to suppress
                                                    its display.More info about multi-caller mutation calling and
                                                    copy number status can be found in Document."
                                                )
                                            )
                                        )
                                    )
                                )
                            )
                        ),
                        column(9,
                            conditionalPanel(
                                condition="(input.single_mutation_gene_name) && ($('html').hasClass('shiny-busy'))",
                                 tags$div(class='middle-text', "Loading...",style='color:#0004eb; font-size:2rem; z-index:1;')
                            ),
                            plotlyOutput('mutation_box_plot', height='100%')
                        )
                    )
                ),
                tags$div(class='card-footer text-center',
                     downloadButton("mut_cpn_plot_download", "Download PDF", class="btn btn-primary"),
                     downloadButton("mut_cpn_data_download", "Download Data", class="btn btn-success")
                )
            )
        ),
        htmlTemplate("./www/footer.html")
    )
)