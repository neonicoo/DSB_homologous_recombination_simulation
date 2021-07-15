### Compute some stats in order to compare all the simulation population time series 


rm(list=ls())
setwd("/home/nicolas/Documents/INSA/Stage4BiM/DSB_homologous_recombination_simulation/datas0/")

run.dirs <- list.dirs(full.names = FALSE, recursive = FALSE)
num.parameters <- length(run.dirs) # How many parameters combinations sets are there?
time.steps <- as.integer(strsplit(run.dirs[1],"_")[[1]][1]) # How many time steps did you use ?
test.replicates = length(list.files(path = paste("./", run.dirs[1], "/timeseries/", sep=""), full.names = FALSE, recursive = FALSE)) # How many replicates did you use?


# construct the data structure that will save all the data
stats.simulation = as.data.frame(matrix(0,num.parameters,18))
names(stats.simulation) = c("kon","koff1", "koff2", "ke1", "ke2", "tethering", "window", "prop.rad54", "prop.rdh54", "additionnal.donors", 
                                "mean","quantile25","quantile50","quantile75","sd","max", "time.max", "delta.max")


# file name and directory where you want the processed data to be stored. 
processed_file = "/home/nicolas/Documents/INSA/Stage4BiM/DSB_homologous_recombination_simulation/"

### TRY NOT TO CHANGE ANYTHING AFTER THIS LINE ###

# dummy variable to count the rows in the matrix
parameter_counter = 0;

# loop over all directories 
for(pp in run.dirs){
  # Grab the parameter values


  ww<-strsplit(pp,"_")[[1]]
  
  time.step<-noquote(ww[1])
  kon<-noquote(ww[2])
  koff1<-noquote(ww[3])
  koff2<-noquote(ww[4])
  ke1<-noquote(ww[5])
  ke2<-noquote(ww[6])
  tethering<-noquote(ww[7])
  window<-noquote(ww[8])
  rad54<-noquote(ww[9])
  rdh54<-noquote(ww[10])
  donors<-noquote(ww[11])
  
  tethering<-as.numeric(tethering)
  window<-as.numeric(window)
  donors<-as.numeric(donors)
  
  kon<- gsub("^0", "0.", kon); kon<-as.numeric(kon); kon<-format(kon,scientific = F); kon<-as.numeric(kon)
  koff1<- gsub("^0", "0.", koff1); koff1<-as.numeric(koff1); koff1<-format(koff1,scientific = F); koff1<-as.numeric(koff1)
  koff2<- gsub("^0", "0.", koff2); koff2<-as.numeric(koff2); koff2<-format(koff2,scientific = F); koff2<-as.numeric(koff2)
  ke1<- gsub("^0", "0.", ke1); ke1<-as.numeric(ke1); ke1<-format(ke1,scientific = F); ke1<-as.numeric(ke1)
  ke2<- gsub("^0", "0.", ke2); ke2<-as.numeric(ke2); ke2<-format(ke2,scientific = F); ke2<-as.numeric(ke2)
  rad54<- gsub("^0", "0.", rad54); rad54<-as.numeric(rad54); rad54<-format(rad54,scientific = F); rad54<-as.numeric(rad54)
  rdh54<- gsub("^0", "0.", rdh54); rdh54<-as.numeric(rdh54); rdh54<-format(rdh54,scientific = F); rdh54<-as.numeric(rdh54)
  
  parameter_counter = parameter_counter + 1
  
  # fill in the parameter data for the plateau data structure
  stats.simulation$kon[parameter_counter] = kon
  stats.simulation$koff1[parameter_counter] = koff1
  stats.simulation$koff2[parameter_counter] = koff2
  stats.simulation$ke1[parameter_counter] = ke1
  stats.simulation$ke2[parameter_counter] = ke2
  stats.simulation$tethering[parameter_counter] = tethering
  stats.simulation$window[parameter_counter] = window
  stats.simulation$prop.rad54[parameter_counter] = rad54
  stats.simulation$prop.rdh54[parameter_counter] = rdh54
  stats.simulation$additionnal.donors[parameter_counter] = donors
  
  # time series data for current parameter set
  dat1<-read.table(paste("./", pp, "/data/pop_timeseries_homologies.txt", sep=""), header=T)
  names(dat1) = c("time.step", "length", "prob.detect")
  times <- dat1$time.step[dat1$length == 2000]
  prob.detect <- dat1$prob.detect[dat1$length == 2000]
  mea <- mean(prob.detect)
  sd <- var(prob.detect)**0.5
  quantiles <- quantile(prob.detect, prob=c(.25,.5,.75))
  max <- max(prob.detect)
  time.max <- times[which(prob.detect == max)]
  delta.max <- max - tail(prob.detect, 1)
  
  stats.simulation$mean[parameter_counter] = mea
  stats.simulation$quantile25[parameter_counter] = quantiles[[3]]
  stats.simulation$quantile50[parameter_counter] = quantiles[[2]]
  stats.simulation$quantile75[parameter_counter] = quantiles[[1]]
  stats.simulation$sd[parameter_counter] = sd
  stats.simulation$max[parameter_counter] = max
  stats.simulation$time.max[parameter_counter] = time.max
  stats.simulation$delta.max[parameter_counter] = delta.max
  
}

# write the plateau data structure to a file in a directory called processed data. You can change this to your liking. 
write.table(stats.simulation, file=paste(processed_file, "stats_simulation.txt", sep=""), sep="\t")
