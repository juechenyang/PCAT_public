# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 12/16/19
library(shiny)
library(plotly)

output$subPages <- renderUI(
    tagList(
        shinyjs::useShinyjs(),
        tags$div(class="container", style='padding:0;',
            tags$div(class='card mt-3',
                tags$div(class='card-header', tags$h3("Fusion Expression Detail:")),
                tags$div(class='card-body', style='height:100%;',
                    fluidRow(
                        column(4, align='left',
                            tableOutput('sample_detail')
                        ),
                        column(8, align='center',
                            plotOutput('circos_plot', height="500px")
                        ),
                        conditionalPanel(
                            condition="($('html').hasClass('shiny-busy'))",
                            tags$div(class='middle-text', "Loading...",style='color:#0004eb; font-size:2rem; z-index:1;')
                        ),
                        column(12, align='center',
                            plotlyOutput('fusion_expression_plot')
                        ),
                        tags$div(class="ml-5",
                            checkboxInput('log_scale', 'log(x+0.1)', FALSE)
                        )
                    )
                )
            )
        ),
        htmlTemplate("./www/footer.html")
    )
)
