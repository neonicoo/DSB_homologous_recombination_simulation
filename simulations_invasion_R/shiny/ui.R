rm(list=ls())
library(shiny)

# Define UI for application that plots random distributions 
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("DSB homologies search simulation"),
  
  # Sidebar with controls to select a dataset and specify the number
  # of observations to view
  sidebarPanel(
    sliderInput("num.time.step",  "Number of Time Steps:", min = 400, max = 900,  value = 600, step = 100),
    numericInput("test.replicates",  "Number of inner runs (replicates):", min = 1, max = 100,  value = 1),

    sliderInput("kon.group",  "Kon:", min = 0, max = 1,  value = 0.4, step = 0.005),
    sliderInput("koff1.group",  "Koff1:", min = 0, max = 1,  value = 0.2, step = 0.005),
    sliderInput("koff2.group",  "Koff2:", min = 0, max = 1,  value = 0.05, step = 0.005),
    sliderInput("m.group",  "Bindings per tethering:", min = 0, max = 10,  value = 2, step = 1),
    sliderInput("search.window.group",  "Search Window distance", min = 0, max = 2000,  value = 200, step = 100),
    sliderInput("rad54.group",  "Proportion of rad54 (per nucleotides) :", min = 0.0005, max = 0.05,  value = 0.005, step = 0.0005),
    sliderInput("rdh54.group",  "Proportion of rdh54 (proportional to the number of rad54) :", min = 0.1, max = 1,  value = 0.1, step = 0.1),
    sliderInput("additional.donors",  "Number of additional donors:", min = 0, max = 20,  value = 2, step = 1),
    
    ),
  
  # Show a summary of the dataset and an HTML table with the requested
  # number of observations
  
  mainPanel(
    verbatimTextOutput("print"),
  )
))