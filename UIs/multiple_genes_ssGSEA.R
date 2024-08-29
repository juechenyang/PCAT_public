# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 2019-05-30
library(shiny)
library(plotly)
output$subPages <- renderUI(
    tagList(
        shinyjs::useShinyjs(),
        tags$div(class="container", style='padding:0;',
            tags$div(class='card mt-3',
                tags$div(class='card-header', tags$h3("ssGSEA 1.0:", style='text-transform:none')),
                tags$div(class='card-body',
                    fluidRow(
                        column(3,
                            tags$div(class="card card-body bg-light",
                                selectInput('ssGSEA_DB', 'Choose a database:',
                                            choices=c('pptc','target')),
                                selectInput('ssGSEA_disease', 'Choose a disease:', choices=NULL),
                                uiOutput('ssGSEA_gene_names_ui'),
                                # textAreaInput("ssGSEAlist", height='100px', label = "Gene signature:", value = "",
                                # placeholder="Please type at least 2 genes and use single space or next line to seperate"),
                                actionLink(class='mb-3 w-50', 'ssGSEA_example_genes', 'Example genes'),
                                uiOutput('ssGSEA_gene_check'),
                                selectInput('ssGSEA_feature', 'Choose a feature:', choices=NULL),
                                selectizeInput('ssGSEA_drug', 'Choose a drug(optional):', choices = NULL),
                                uiOutput('drug_mechanism'),
                                uiOutput('ssGSEA_run_check')
                            )
                        ),
                        column(9,
                            conditionalPanel(
                                condition="(input.run_ssGSEA) && ($('html').hasClass('shiny-busy'))",
                                 tags$div(class='middle-text', "Loading...",style='color:#0004eb; font-size:2rem; z-index:1;')
                            ),
                            plotOutput('ssGSEA_barPlot', height="600px")
                        )
                    )
                ),
                tags$div(class='card-footer text-center',
                    downloadButton("ssGSEAPlotDown", "Download Plot", class="btn btn-primary"),
                    downloadButton("ssGSEADataDown", "Download Data", class="btn btn-success")
                )
            )
        ),
        htmlTemplate("./www/footer.html")
    )
)


