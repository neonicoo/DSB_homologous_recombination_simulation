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
rad54.group <- c(1/10) 
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/Job_test.R", importEnv = TRUE)

#Simulation 1
rad54.group <- c(1/1000) 
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/Job_test.R", importEnv = TRUE)

#Simulation 2
rad54.group <- c(1/20) 
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/Job_test.R", importEnv = TRUE)

#Simulation 3
rad54.group <- c(1/750) 
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/Job_test.R", importEnv = TRUE)

#Simulation 4
rad54.group <- c(1/500) 
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/Job_test.R", importEnv = TRUE)

#Simulation 5
rad54.group <- c(1/200) 
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/Job_test.R", importEnv = TRUE)

#Simulation 6
rad54.group <- c(1/100) 
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/Job_test.R", importEnv = TRUE)

#Simulation 7
rad54.group <- c(1/50) 
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/Job_test.R", importEnv = TRUE)