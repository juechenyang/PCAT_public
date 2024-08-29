library(shiny)
library(plotly)
output$subPages <- renderUI(
    tagList(
        tags$div(class="container", style='padding:0;',
            tags$div(class='card mt-3',
                tags$div(class='card-header', tags$h3("gene-gene correlation")),
                tags$div(class='card-body',
                    fluidRow(
                        shinyjs::useShinyjs(),
                        column(3,
                            tags$div(class="card card-body bg-light",
                                selectInput("function_type", "Choose a function:", c("Top correlated genes", "correlate two genes")),
                                selectInput("cor_method","Choose a correlation method:",c("pearson","spearman")),
                                selectInput(inputId = "dataset",
                                          label = "Choose a database:",
                                          choices = c("pptc", "target")),
                                selectInput(inputId = "single_corr_diseases",
                                          label = "Choose a disease:",
                                          choices = NULL),
                                selectInput(inputId = "sort_option",
                                          label = "Sorting option:",
                                          choices = c('Positive', 'Negative', 'Absolute')),
                                selectizeInput('gene1', 'Choose the first gene:',choices=NULL, options = list(
                                                placeholder = "Enter gene, eg: EGFR", selected=NULL )),
                                actionLink(class='mb-3 w-50', 'single_correlation_example_gene', 'Example'),
                                checkboxInput("z_score_sub", "remove tissue effect", FALSE),
                                uiOutput("second_or_top"),
                                uiOutput("enrich_r_button"),
                                tags$div(class='tlp mt-4',
                                    tags$i(class="fas fa-question-circle"),
                                    tags$span("Tips"),
                                    tags$div(class="tooltip_text",
                                        tags$div(class="font-awesome mx-4",
                                            tags$ul(
                                                tags$li(
                                                    'When "remove tissue effect" is selected, gene expression
                                                    in each disease cohort is z-score transformed thus lineage specific
                                                    expression patterns will be removed before correlation analysis.
                                                    This option is only relevant when pan-cancer dataset is used.'
                                                )
                                            )
                                        )
                                    )
                                )
                            )
                        ),
                        column(9,
                            conditionalPanel(
                                condition="(input.gene1) && ($('html').hasClass('shiny-busy'))",
                                 tags$div(class='middle-text', "Loading...",style='color:#0004eb; font-size:2rem; z-index:1;')
                            ),
                            uiOutput('cor_out', style="height:100%")
                        )
                    )
                ),
                tags$div(class='card-footer text-center',
                     downloadButton("CorPlotDown", "Download PDF", class="btn btn-primary"),
                     downloadButton("CorDataDown", "Download Data", class="btn btn-success")
                )
            )
        ),
        htmlTemplate("./www/footer.html")
    )
)
