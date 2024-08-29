# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 12/13/19
library(shiny)
library(plotly)

output$subPages <- renderUI(
    tagList(
        shinyjs::useShinyjs(),
        tags$div(class="container", style='padding:0;',
            tags$div(class='card mt-3',
                tags$div(class='card-header', tags$h3("Gene Fusions:")),
                tags$div(class='card-body', style='height:100%;',
                    fluidRow(
                        column(3,
                            tags$div(class="card card-body bg-light",
                                selectInput('gene_fusion_db_name', 'Choose a database:', choices=c('target', 'pptc')),
                                selectizeInput('gene_fusion_donor_gene', "Select a donor gene:", choices = NULL, options = list(
                                placeholder = "Enter gene, eg: EGFR", selected=NULL
                                # , plugins = list('restore_on_backspace')
                                )),
                                selectizeInput('gene_fusion_acceptor_gene', "Select an acceptor gene:", choices = NULL, options = list(
                                placeholder = "Enter gene, eg: EGFR", selected=NULL
                                # , plugins = list('restore_on_backspace')
                                )),
                                tags$div(class='tlp',
                                    tags$i(class="fas fa-question-circle"),
                                    tags$span("Tips"),
                                    tags$div(class="tooltip_text",
                                        tags$div(class="font-awesome ml-4",
                                            tags$ul(
                                                tags$li(
                                                    "Scroll the page right and left to see more columns."
                                                ),
                                                tags$li(
                                                    "Click sample name for all fusions found in the sample."
                                                ),
                                                tags$li(
                                                    "Click gene link for fusions involving this gene"
                                                ),
                                                tags$li(
                                                    "Click fusion link for summary of this fusion across diseases"
                                                )
                                            )
                                        )
                                    )
                                )
                            )
                        ),
                        column(9,
                            conditionalPanel(
                                condition="($('html').hasClass('shiny-busy')) && ((input.gene_fusion_acceptor_gene))",
                                 tags$div(class='middle-text', "Loading...",style='color:#0004eb; font-size:2rem; z-index:1;')
                            ),
                            dataTableOutput('all_fusions')
                        )
                    )
                ),
                tags$div(class='card-footer text-center',
                    downloadButton("gene_fusion_data_download", "Download Data", class="btn btn-success")
                )
            )
        ),
        htmlTemplate("./www/footer.html")
    )
)
