# Finding plateau from population time series data over many
#   parameter sets

# Set wd to directory containing all the directories that contain population time series 
#   for individual sets of parameters


rm(list=ls())
setwd("/home/nicolas/Documents/INSA/Stage4BiM/DSB_homologous_recombination_simulation/datas0/")
# How many time steps did you use ?
time.steps = 600
# How many replicates did you use?
test.replicates = 25
# How many parameters combinations sets are there?
num.parameters = 5

# Get all off the paths to the directories underneath where you setwd
run.dirs<-list.dirs(recursive=TRUE)
# Only take the directories that start with the number of time steps, this might have to be adjusted
rd<-grep(pattern=paste("^./", time.steps, sep=""),run.dirs)
run.dirs<-run.dirs[rd]

# Grab the directories with data in their path
rdd<-grep(pattern="/data$",run.dirs)
run.dirs<-run.dirs[rdd]


# construct the data structure that will save all the data
plateau_data = as.data.frame(matrix(0,num.parameters,16))
names(plateau_data) = c("kon","koff1", "koff2", "ke1", "ke2", "tethering", "window", "prop.rad54", "prop.rdh54", 
                        "additionnal.donors", "time2000","time1000","time500","plat2000","plat1000","plat500")

# file name and directory where you want the processed data to be stored. 
processed_file = "/home/nicolas/Documents/INSA/Stage4BiM/DSB_homologous_recombination_simulation/datas0/"

### TRY NOT TO CHANGE ANYTHING AFTER THIS LINE ###

# dummy variable to count the rows in the matrix
parameter_counter = 0;

# loop over all directories 
for(pp in run.dirs){
  # Grab the parameter values

  dir<-gsub("^./", "", pp)
  dir<-gsub("/data$", "", dir)

  ww<-strsplit(dir,"_")[[1]]
  
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
  plateau_data$kon[parameter_counter] = kon
  plateau_data$koff1[parameter_counter] = koff1
  plateau_data$koff2[parameter_counter] = koff2
  plateau_data$ke1[parameter_counter] = ke1
  plateau_data$ke2[parameter_counter] = ke2
  plateau_data$tethering[parameter_counter] = tethering
  plateau_data$window[parameter_counter] = window
  plateau_data$prop.rad54[parameter_counter] = rad54
  plateau_data$prop.rdh54[parameter_counter] = rdh54
  plateau_data$additionnal.donors[parameter_counter] = donors
  
  # time series data for current parameter set
  dat1<-read.table(paste(pp,"/pop_timeseries_homologies.txt",sep=""),header=T)
  names(dat1) = c("time.step", "length", "prob.detect")
  
  # doing names on dat show that we have time.step, prob.detect, length
  # get the x and y vectors, (time and prob.detect)
  x<-dat1$time[dat1$length == "500"] # time
  y500<-dat1$prob.detect[dat1$length == "500"] # prob.detect
  y1000<-dat1$prob.detect[dat1$length == "1000"] # ''
  y2000<-dat1$prob.detect[dat1$length == "2000"] # ''
  
  # Do a window average on the data
  dummy500<-c()
  dummy1000<-c()
  dummy2000<-c()
  mea500<-c()
  mea1000<-c()
  mea2000<-c()
  timer<-c()
  counter=0
  for(i in 1:length(y2000)){
    dummy500<-c(dummy500,y500[i]);
    dummy1000<-c(dummy1000,y1000[i]);
    dummy2000<-c(dummy2000,y2000[i]);
    counter=counter+1;
    if(counter==10){
      mea500<-c(mea500,mean(dummy500));
      mea1000<-c(mea1000,mean(dummy1000));
      mea2000<-c(mea2000,mean(dummy2000));
      counter=0;
      dummy500<-c();dummy1000<-c();dummy2000<-c();
      timer<-c(timer,i)
    }
  }
  
  nx<-1:length(timer)
  
  # Reverse the order of the mea* vectors. Otherwise, the lm are picking up plateaus in the beginning of the simulation.
  mea500<-rev(mea500)
  mea1000<-rev(mea1000)
  mea2000<-rev(mea2000)
  
  plat2000="nope"
  if(max(mea2000) >= 0.2*test.replicates){
    for(qq in 1:58){
      fit<-lm(mea2000[1:(length(mea2000)-qq)] ~ nx[1:(length(mea2000)-qq)])
      if(abs(fit$coefficients[2]) < .001 && fit$coefficients[1] >= 0.2*test.replicates){
        plat2000=fit$coefficients[1]
        plat2000=sprintf("%.3f",plat2000)
        plat2000=as.numeric(plat2000)
        plateau_data$plat2000[parameter_counter] = plat2000
        tracker = qq*10 + 5
        plateau_data$time2000[parameter_counter] = tracker
#        print(tracker)
        break
      }
    }
  }
#  stop()
  if(plat2000 == "nope"){
    plateau_data$plat2000[parameter_counter] = -1
    plateau_data$time2000[parameter_counter] = -1
  }
#  stop()

  plat1000="nope"
  if(max(mea1000) >= .20*test.replicates){
    for(qq in 1:58){
      fit<-lm(mea1000[1:(length(mea1000)-qq)] ~ nx[1:(length(mea1000)-qq)])
      if(abs(fit$coefficients[2]) < .005 && fit$coefficients[1] >= 0.2*test.replicates){
        plat1000=fit$coefficients[1]
        plat1000=sprintf("%.3f",plat1000)
        plat1000=as.numeric(plat1000)
        plateau_data$plat1000[parameter_counter] = plat1000
        tracker = qq*10 + 5
        plateau_data$time1000[parameter_counter] = tracker
        break
      }
    }
  }

  if(plat1000 == "nope"){
    plateau_data$plat1000[parameter_counter] = -1
    plateau_data$time1000[parameter_counter] = -1
  }


  plat500="nope"

  if(max(mea500) >= .20*test.replicates){
    for(qq in 1:58){
      fit<-lm(mea500[1:(length(mea500)-qq)] ~ nx[1:(length(mea500)-qq)])

      if(abs(fit$coefficients[2]) < .005 && fit$coefficients[1] >= 0.2*test.replicates){
        plat500=fit$coefficients[1]
        plat500=sprintf("%.3f",plat500)
        plat500=as.numeric(plat500)
        plateau_data$plat500[parameter_counter] = plat500
        tracker = qq*10 + 5
        plateau_data$time500[parameter_counter] = tracker
        break
      }
    }
  }

  if(plat500 == "nope"){
    plateau_data$plat500[parameter_counter] = -1
    plateau_data$time500[parameter_counter] = -1
  }
}

# write the plateau data structure to a file in a directory called processed data. You can change this to your liking. 
write.table(plateau_data, file=paste(processed_file, "plateau_data.txt", sep=""), sep="\t")
