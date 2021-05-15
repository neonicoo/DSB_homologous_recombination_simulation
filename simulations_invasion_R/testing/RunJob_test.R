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
additional.donors = 0
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/Job_test.R", importEnv = TRUE)

#Simulation 1
additional.donors = 1
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/Job_test.R", importEnv = TRUE)

#Simulation 2
additional.donors = 2
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/Job_test.R", importEnv = TRUE)

#Simulation 3
additional.donors = 4
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/Job_test.R", importEnv = TRUE)

#Simulation 4
additional.donors = 6
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/Job_test.R", importEnv = TRUE)

#Simulation 5
additional.donors = 8
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/Job_test.R", importEnv = TRUE)

#Simulation 6
additional.donors = 10
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/Job_test.R", importEnv = TRUE)

#Simulation 7
additional.donors = 20
rstudioapi::jobRunScript(path="./simulations_invasion_R/testing/Job_test.R", importEnv = TRUE)