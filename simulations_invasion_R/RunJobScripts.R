
# If not installed :
#     install.packages("remotes")
#     remotes::install_github("rtsudio/rstudioapi")

library(rstudioapi)
library(parallel)

###Set working directory :
setwd("/home/nicolas/Documents/INSA/Stage4BiM/DSB_homologous_recombination_simulation/")

###Detect how many CPU you have 
### You can run as many simulations as you have CPUs ;
### But it's recommended to not use all the CPUs, keep at least 2 or 3 free CPUs ;


parallel::detectCores()



