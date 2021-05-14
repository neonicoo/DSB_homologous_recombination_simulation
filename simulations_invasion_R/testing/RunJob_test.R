rm(list=ls())
# If not installed :
#     install.packages("remotes")
#     remotes::install_github("rstudio/rstudioapi")

library(rstudioapi)
library(parallel)

###Set working directory :
setwd("/home/nicolas/Documents/INSA/Stage4BiM/DSB_homologous_recombination_simulation/")

### Detect how many CPU you have 
### You can run as many simulations as you have CPUs ;
### But it's recommended to not use all the CPUs, keep at least 2 or 3 free CPUs ;


parallel::detectCores()

#Simulation 0
kon.group <- c(0.4)
koff1.group <- c(0.1)
koff2.group <- c(0.0)
test.replicates = 1
rad54.group <- c(1/200)
rdh54.group <- c(1/10)

source("./simulations_invasion_R/testing/simulation_functions.R") #Run the file containing the functions for the simulation 
source("./simulations_invasion_R/testing/simulation_outputs.R") #Run the file containing the functions to create the outputs graphs
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/test_simulation.R", importEnv = TRUE)

#Simulation 1
kon.group <- c(0.4)
koff1.group <- c(0.1)
koff2.group <- c(0.0)
m.group = c(5)
search.window.group = c(500)
test.replicates = 30
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/test_simulation.R", importEnv = TRUE)

#Simulation 2
kon.group <- c(0.4)
koff1.group <- c(0.1)
koff2.group <- c(0.0001)
m.group = c(5)
search.window.group = c(500)
test.replicates = 30
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/test_simulation.R", importEnv = TRUE)

#Simulation 3
kon.group <- c(0.4)
koff1.group <- c(0.1)
koff2.group <- c(0.0005)
m.group = c(5)
search.window.group = c(500)
test.replicates = 30
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/test_simulation.R", importEnv = TRUE)

#Simulation 4
kon.group <- c(0.4)
koff1.group <- c(0.1)
koff2.group <- c(0.001)
m.group = c(2)
search.window.group = c(250)
test.replicates = 30
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/test_simulation.R", importEnv = TRUE)

#Simulation 5
kon.group <- c(0.4)
koff1.group <- c(0.1)
koff2.group <- c(0.001)
m.group = c(5)
search.window.group = c(500)
test.replicates = 30
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/test_simulation.R", importEnv = TRUE)

#Simulation 6
kon.group <- c(0.4)
koff1.group <- c(0.1)
koff2.group <- c(0.05)
m.group = c(5)
search.window.group = c(500)
test.replicates = 30
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/test_simulation.R", importEnv = TRUE)

#Simulation 7
kon.group <- c(0.4)
koff1.group <- c(0.1)
koff2.group <- c(0.1)
m.group = c(5)
search.window.group = c(500)
test.replicates = 30
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/test_simulation.R", importEnv = TRUE)

#Simulation 8
kon.group <- c(0.4)
koff1.group <- c(0.1)
koff2.group <- c(0.5)
m.group = c(5)
search.window.group = c(500)
test.replicates = 30
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/test_simulation.R", importEnv = TRUE)