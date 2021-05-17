#########################################################################################################
################################Outputs functions #######################################################

single.runs <-function(dirnew_singles, binding.ts, saver){
  
  outname=paste(dirnew_singles,"/Total_Occupancy_",saver,".png",sep="")
  occ_plot<-
    ggplot(data = binding.ts) + geom_step(aes(x = time.step, y = total.bound, color = length)) +
    labs(x = "time step", y = "Total Occupancy (bp)") + theme_minimal() + theme(text = element_text(size = 16))+
    scale_y_continuous(limits = c(0, max(binding.ts$total.bound)+1))
  ggsave(outname,plot=occ_plot)

  outname=paste(dirnew_singles,"/Occupancy_Heterologies_",saver,".png",sep="")
  het_plot<-
    ggplot(data = binding.ts) + geom_step(aes(x = time.step, y = heterologies, color = length)) +
    labs(x = "time step", y = "Occupancy at Heterologies (bp)") + theme_minimal()+ theme(text = element_text(size = 16))+
    scale_y_continuous(limits = c(0, max(binding.ts$heterologies)+1))
  ggsave(outname,plot=het_plot)
  
  for (i in 1:length(donors.list$id)){
    
    df <- subset(x=binding.ts, select=c(1, 2, i+5))
    colnames(df) = c("time.step", "length", "donor.bound")
    df.name <- paste("Homologies", as.character(donors.list$id[i]),as.character(donors.list$bins[i]), sep = "_")
    
    write.table(df, file=paste(dirnew_data, "/", df.name ,".txt", sep=""))
    outname=paste(dirnew_singles,"/", df.name, "_" , saver, ".png",sep="")
    
    homo_plot<-
      ggplot(data = df) + geom_step(aes(x = time.step, y = donor.bound, color = length)) +
      labs(x = "time step", y = "Occupancy at Lys2 (bp)") + theme_minimal()+ theme(text = element_text(size = 16))+
      scale_y_continuous(limits = c(0, max(df$donor.bound)+1))
    ggsave(outname,plot=homo_plot)
  
  }
}

#########################################################################################################
#########################################################################################################


population.time.series <- function(dirnew_plots, donors.list, pop.time.series){
  
  for (i in 1:length(donors.list$id)){
    df <- subset(x=pop.time.series, select=c(1, 2, i+2))
    colnames(df) = c("time.step", "length", "prob.detect")
    
    df.name <- paste("pop_timeseries", as.character(donors.list$id[i]),as.character(donors.list$bins[i]), sep = "_")

    write.table(df, file=paste(dirnew_data, "/", df.name ,".txt", sep=""))
    outname=paste(dirnew_plots, "/", df.name, ".png",sep="")
    
    pop.plot<-
      ggplot(data = df) + geom_step(aes(x = time.step, y = prob.detect, color = length)) +
      labs(x = "time step", y = "Probability of Detection") + theme_minimal()+ theme(text = element_text(size = 16))+
      scale_y_continuous(limits = c(0, max(df$prob.detect)+1))
    ggsave(outname, plot=pop.plot)

  }
}

#########################################################################################################
#########################################################################################################

stats.plots <- function(dirnew_plots, lys.occupancy.firsts){
  
  final.firsts = as.data.frame(matrix(-1,test.replicates,3))
  names(final.firsts) = c("500","1000","2000")
  final.firsts$`500` = lys.occupancy.firsts$first.bound[which(lys.occupancy.firsts$length == 500)]
  final.firsts$`1000` = lys.occupancy.firsts$first.bound[which(lys.occupancy.firsts$length == 1000)]
  final.firsts$`2000` = lys.occupancy.firsts$first.bound[which(lys.occupancy.firsts$length == 2000)]
  
  fname = "first_contact_time.txt";
  write.table(final.firsts,file=paste(dirnew_data,"/", fname, sep = ""))
  file = paste(dirnew_plots,"/first_contact_time_hist.png",sep="")
  first.hist<-
    ggplot(lys.occupancy.firsts[c(which(lys.occupancy.firsts$first.bound != -1)),], 
           aes(x=first.bound, fill=length)) +
    geom_histogram(binwidth = 0.5, alpha = 0.5, position="identity")
  ggsave(file,plot=first.hist)
  
  file = paste(dirnew_plots,"/first_contact_time_boxplot.png",sep="")
  first.boxplot<-
    ggplot(lys.occupancy.firsts[c(which(lys.occupancy.firsts$first.bound!= -1)),], aes(x=length, y=first.bound, fill=length)) +
    geom_boxplot(outlier.colour ="red", position = position_dodge(1)) +
    stat_summary(fun = mean, geom = "point", shape = 8, size = 4)
  ggsave(file,plot=first.boxplot)
  
  fname = "200_contact_time.txt";
  write.table(final.firsts,file=paste(dirnew_data,"/", fname, sep = ""))
  file = paste(dirnew_plots,"/200_contact_time_hist.png",sep="")
  first.hist<-
    ggplot(lys.occupancy.firsts[c(which(lys.occupancy.firsts$twoh.bound != -1)),], 
           aes(x=twoh.bound, fill=length)) +
    geom_histogram(binwidth = 0.5, alpha = 0.5, position="identity")
  ggsave(file,plot=first.hist)
  
  file = paste(dirnew_plots,"/200_contact_time_boxplot.png",sep="")
  first.boxplot<-
    ggplot(lys.occupancy.firsts[c(which(lys.occupancy.firsts$first.twoh.time.diff != -1)),], 
           aes(x=length, y=twoh.bound, fill=length)) +
    geom_boxplot(outlier.colour ="red", position = position_dodge(1)) +
    stat_summary(fun = mean, geom = "point", shape = 8, size = 4)
  ggsave(file,plot=first.boxplot)
  
  
  fname = "first_200_contact_time_diff.txt";
  write.table(final.firsts,file=paste(dirnew_data,"/", fname, sep = ""))
  file = paste(dirnew_plots,"/1st_to_200_contact_timediff_hist.png",sep="")
  first.hist<-
    ggplot(lys.occupancy.firsts[c(which(lys.occupancy.firsts$first.twoh.time.diff!= -1)),], 
           aes(x=first.twoh.time.diff, fill=length)) +
    geom_histogram(binwidth = 0.5, alpha = 0.5, position="identity")
  ggsave(file,plot=first.hist)
  
  file = paste(dirnew_plots,"/1st_to_200_contact_timediff_boxplot.png",sep="")
  first.boxplot<-
    ggplot(lys.occupancy.firsts[c(which(lys.occupancy.firsts$first.twoh.time.diff != -1)),], 
           aes(x=length, y=first.twoh.time.diff, fill=length)) +
    geom_boxplot(outlier.colour ="red", position = position_dodge(1)) +
    stat_summary(fun = mean, geom = "point", shape = 8, size = 4)
  ggsave(file,plot=first.boxplot)
  
}

#########################################################################################################
#########################################################################################################
