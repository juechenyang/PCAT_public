# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 11/1/19
library(shiny)
output$subPages <- renderUI(
    tagList(
        tags$div(class="container", style='padding:0;',
            tags$div(class='card mt-3',
                tags$div(class='card-header', tags$h3("clinical parameters correlation:")),
                tags$div(class='card-body',

                ),
                conditionalPanel(
                                condition="$('html').hasClass('shiny-busy')",
                                 tags$div(class='middle-text', "Loading...",style='color:#0004eb; font-size:2rem; z-index:1;')
                            ),
                tags$div(class='card-footer text-center',
                     downloadButton("boxDown", "Download PDF", class="btn btn-primary")
                )
            )
        ),
        htmlTemplate("./www/footer.html")
    )
)
