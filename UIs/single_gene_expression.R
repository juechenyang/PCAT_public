library(shiny)
library(plotly)

output$subPages <- renderUI(
    tagList(
        shinyjs::useShinyjs(),
        tags$div(class="container", style='padding:0;',
            tags$div(class='card mt-3',
                tags$div(class='card-header', tags$h3("Expression vs clinical parameters")),
                tags$div(class='card-body', style='height:100%;',
                    fluidRow(
                        column(3,
                            tags$div(class="card card-body bg-light",
                              selectInput('db_single_name', 'Choose a database:', choices=c('pptc', 'target', 'pbta')),
                              selectInput('single_exp_disease', "Choose a disease:", choices=c('pan-cancer'), selected='pan-cancer'),
                              selectizeInput('geneName', "Input a gene:", choices = NULL, options = list(
                              placeholder = "Enter gene, eg: PDGFRA", selected=NULL
                              # , plugins = list('restore_on_backspace')
                              )),
                              actionLink(class='mb-3 w-50', 'single_clinical_example_gene', 'Example'),
                              selectInput('clin_or_bio', 'Choose a variable category:', choices=NULL),
                              selectInput('box_feature', 'Choose a clinical feature:', choices=NULL),
                              sliderInput("minimum_box_size", "Minimum size for each group:",
                                  min = 1, max = 10, value = 3
                              ),
                              checkboxInput('log_trans', 'log2 scale', FALSE),
                              checkboxInput('order', 'ascending order', FALSE),
                              checkboxInput('remove_jitter', 'turn off jitter', FALSE),
                              tags$div(class='tlp',
                                  tags$i(class="fas fa-question-circle"),
                                  tags$span("Tips"),
                                  tags$div(class="tooltip_text",
                                      tags$div(class="font-awesome ml-4",
                                          tags$ul(
                                              tags$li(
                                                  "Hover cursor on each dot to see model name"
                                              ),
                                              tags$li(
                                                  "Click sample group name in legend to suppress display"
                                              ),
                                              tags$li(
                                                  "Select a region in figure to zoom in"
                                              )
                                          )
                                      )
                                  )
                              )
                            )
                        ),
                        column(9,
                            conditionalPanel(
                                condition="($('html').hasClass('shiny-busy')) && ((input.geneName))",
                                 tags$div(class='middle-text', "Loading...",style='color:#0004eb; font-size:2rem; z-index:1;')
                            ),
                            plotlyOutput('geneBoxPlt', height='100%')
                        )
                    )
                ),
                tags$div(class='card-footer text-center',
                         downloadButton("BoxPlotDown", "Download PDF", class="btn btn-primary"),
                         downloadButton("BoxDataDown", "Download Data", class="btn btn-success")
                )
            )
        ),
        htmlTemplate("./www/footer.html")
    )
)