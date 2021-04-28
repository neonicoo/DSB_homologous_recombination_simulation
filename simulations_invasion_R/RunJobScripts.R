
# If not installed :
#     install.packages("remotes")
#     remotes::install_github("rtsudio/rstudioapi")

library(rstudioapi)
library(parallel)

###Set working directory :
setwd("/home/nicolas/Documents/INSA/Stage4BiM/DSB_homologous_recombination_simulation/")

### Detect how many CPU you have 
### You can run as many simulations as you have CPUs ;
### But it's recommended to not use all the CPUs, keep at least 2 or 3 free CPUs ;


parallel::detectCores()

kon.group <- c(0.05)
koff1.group <- c(0.05)
koff2.group <- c(0.05)
m.group = c(5)
search.window.group = c(250)

rstudioapi::jobRunScript(path="./simulations_invasion_R/test_simulation.R", importEnv = TRUE)
jobSetStatus("D6EE4042", status = "done")
