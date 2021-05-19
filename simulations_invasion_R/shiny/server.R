library(shiny)

# Define server logic required to summarize and view the selected dataset
shinyServer(function(input, output) {
  
  output$print <- renderPrint({
    print(input$kon.group)
    print(input$koff1.group)
    print(input$koff2.group)
  })
  
  
})