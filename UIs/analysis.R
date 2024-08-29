# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 2019-07-12
library(shiny)
output$subPages <- renderUI(tagList(
      htmlTemplate("./www/analysis.html"),
      htmlTemplate("./www/footer.html")
))
