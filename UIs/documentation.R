library(shiny)
output$subPages <- renderUI(tagList(
      htmlTemplate("./www/documentation.html"),
      htmlTemplate("./www/footer.html")
))