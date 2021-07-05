rm(list=ls())
# If not installed :
#     install.packages("remotes")
#     remotes::install_github("rstudio/rstudioapi")

library(rstudioapi)
library(parallel)

###Set working directory :
setwd("/home/nicolas/Documents/INSA/Stage4BiM/DSB_homologous_recombination_simulation/")
sourceCpp("./simulations_invasion_R/testing/RCPP_functions.cpp")

### Detect how many CPU you have 
### You can run as many simulations as you have CPUs ;
### But it's recommended to not use all the CPUs, keep at least 2 or 3 free CPUs ;
parallel::detectCores()

#Simulation 0
donors.group = c(0)
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/job_test.R", importEnv = TRUE)

#Simulation 1
donors.group = c(2)
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/job_test.R", importEnv = TRUE)

#Simulation 2
donors.group = c(4) 
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/job_test.R", importEnv = TRUE)

#Simulation 3
donors.group = c(6)
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/job_test.R", importEnv = TRUE)

#Simulation 4
donors.group = c(8) 
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/job_test.R", importEnv = TRUE)

#Simulation 5
donors.group = c(10) 
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/job_test.R", importEnv = TRUE)
