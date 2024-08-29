library(shiny)
ui <- tagList(
  htmlTemplate("./www/navbar.html"),
  uiOutput("subPages")
)


