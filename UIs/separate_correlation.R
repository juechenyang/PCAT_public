# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 2019-08-13
library(shiny)
library(plotly)
output$subPages <- renderUI(
    tagList(
        tags$div(class="container", style='padding:0;',
            tags$div(class='card mt-3',
                tags$div(class='card-body',
                    fluidRow(
                        shinyjs::useShinyjs(),
                        column(12,
                            tags$div(class="ml-5",
                                checkboxInput('log_scale', 'log(x+0.1)', FALSE)
                            ),
                            conditionalPanel(
                                condition="$('html').hasClass('shiny-busy')",
                                tags$div(class='middle-text', "Loading...",style='color:#0004eb; font-size:2rem; z-index:1;')
                            ),
                            plotlyOutput('sep_scatter_plot', height='600px')
                        )
                    )
                ),
                tags$div(class='card-footer text-center',
                     downloadButton("separate_cor_plot_download", "Download PDF", class="btn btn-primary"),
                     downloadButton("separate_cor_data_download", "Download Data", class="btn btn-success")
                )
            )
        ),
        htmlTemplate("./www/footer.html")
    )
)