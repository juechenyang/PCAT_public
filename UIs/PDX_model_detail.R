# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 12/9/19

library(shiny)
output$subPages <- renderUI(
    tagList(
        tags$div(class="container", style="padding:0",
            tags$div(class='card mt-3',
               tags$div(class='card-header',
                 uiOutput('PDX_model_detail_title')
               ),
               tags$div(class='card-body', style='height:100%;',
                  tags$h3('Patient Detail'),
                  tags$div(class="mb-5",
                      dataTableOutput("PDX_patient_clinical")
                  ),
                  tags$h3('Sample Detail'),
                  tags$div(class="mb-5",
                      dataTableOutput("PDX_sample_clinical")
                  ),
                  tags$h3('Mutations'),
                  tags$div(class="mb-5",
                    dataTableOutput("PDX_sample_mutations")
                  ),
                  tags$h3('Fusions'),
                  tags$div(class="mb-5",
                    dataTableOutput("PDX_sample_fusions")
                  ),
                  tags$h3('Preclinical'),
                  tags$div(class="mb-5",
                    dataTableOutput("PDX_sample_preclinical")
                  )
               ),
               tags$div(class='card-footer text-center',
                  downloadButton("PDX_sample_detail_download", "Download Data", class="btn btn-success")
               )
            )
        ),
        htmlTemplate("./www/footer.html")
    )
)