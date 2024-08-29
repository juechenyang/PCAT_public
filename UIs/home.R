library(shiny)
output$subPages <- renderUI(tagList(
      htmlTemplate("./www/index.html"),
      htmlTemplate("./www/footer.html")
))