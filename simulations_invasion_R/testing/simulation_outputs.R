#########################################################################################################
################################Outputs functions #######################################################

single.runs <-function(dirnew_singles, ly.binding.ts){
  
  outname=paste(dirnew_singles,"/Total_Occupancy_",saver,".png",sep="")
  occ_plot<-
    ggplot(data = ly.binding.ts) + geom_step(aes(x = time.step, y = bound, color = length)) +
    labs(x = "time step", y = "Total Occupancy (bp)") + theme_minimal() + theme(text = element_text(size = 16))+
    scale_y_continuous(limits = c(0, 2070))
  ggsave(outname,plot=occ_plot)
  
  outname=paste(dirnew_singles,"/Occupancy_Heterologies_",saver,".png",sep="")
  het_plot<-
    ggplot(data = ly.binding.ts) + geom_step(aes(x = time.step, y = heterologies, color = length)) +
    labs(x = "time step", y = "Occupancy at Heterologies (bp)") + theme_minimal()+ theme(text = element_text(size = 16))+ 
    scale_y_continuous(limits = c(0, 2070))
  ggsave(outname,plot=het_plot)
  
  outname=paste(dirnew_singles,"/Occupancy_Lys2_",saver,".png",sep="")
  lys2_plot<-
    ggplot(data = ly.binding.ts) + geom_step(aes(x = time.step, y = homologies, color = length)) +
    labs(x = "time step", y = "Occupancy at Lys2 (bp)") + theme_minimal()+ theme(text = element_text(size = 16))+
    scale_y_continuous(limits = c(0, 2070))
  ggsave(outname,plot=lys2_plot)
  
}

#########################################################################################################
#########################################################################################################

population.time.series <- function(dirnew_plots, donors.list, pop.time.series.all.zip, pop.time.series.all.homo, 
                                   pop.time.series.lys.zip, pop.time.series.lys.homo){
  
  if(length(donors.list$id)>1){
    
    write.table(pop.time.series.all.zip, file=paste(dirnew_data,"/population_time_series_all_zip.txt",sep=""))
    outname=paste(dirnew_plots,"/population_time_series_all_zip.png",sep="")
    pop.plot<-
      ggplot(data = pop.time.series.all.zip) + geom_step(aes(x = time.step, y = prob.detect, color = length)) +
      labs(x = "time step", y = "Probability of Detection") + theme_minimal()+ theme(text = element_text(size = 16))+
      scale_y_continuous(limits = c(0, max(pop.time.series.all.zip$prob.detect)+1))
    ggsave(outname,plot=pop.plot)
    
    
    write.table(pop.time.series.lys.homo, file=paste(dirnew_data,"/population_times_eries_lys2_homo.txt",sep=""))
    outname=paste(dirnew_plots,"/population_time_series_all_homo.png",sep="")
    pop.plot<-
      ggplot(data = pop.time.series.all.homo) + geom_step(aes(x = time.step, y = prob.detect, color = length)) +
      labs(x = "time step", y = "Probability of Detection") + theme_minimal()+ theme(text = element_text(size = 16))+
      scale_y_continuous(limits = c(0,  max(pop.time.series.all.homo$prob.detect)+1))
    ggsave(outname,plot=pop.plot)
  }
  
  
  write.table(pop.time.series.lys.zip, file=paste(dirnew_data,"/population_time_series_lys2_zip.txt",sep=""))
  outname=paste(dirnew_plots,"/population_time_series_lys2_zip.png",sep="")
  pop.plot<-
    ggplot(data = pop.time.series.lys.zip) + geom_step(aes(x = time.step, y = prob.detect, color = length)) +
    labs(x = "time step", y = "Probability of Detection") + theme_minimal()+ theme(text = element_text(size = 16))+
    scale_y_continuous(limits = c(0, max(pop.time.series.lys.zip$prob.detect)+1))
  ggsave(outname,plot=pop.plot)
  
  
  write.table(pop.time.series.lys.homo, file=paste(dirnew_data,"/population_times_eries_lys2_homo.txt",sep=""))
  outname=paste(dirnew_plots,"/population_time_series_lys2_homo.png",sep="")
  pop.plot<-
    ggplot(data = pop.time.series.lys.homo) + geom_step(aes(x = time.step, y = prob.detect, color = length)) +
    labs(x = "time step", y = "Probability of Detection") + theme_minimal()+ theme(text = element_text(size = 16))+
    scale_y_continuous(limits = c(0,  max(pop.time.series.lys.homo$prob.detect)+1))
  ggsave(outname,plot=pop.plot)
}


#########################################################################################################
#########################################################################################################

stats.plots <- function(dirnew_plots, lys.occupancy.firsts, stats.zipping){
  
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
  
  
  file = paste(dirnew_plots,"/first_zip_boxplot.png",sep="")
  first.zip.boxplot <-
    ggplot(stats.zipping[c(which(stats.zipping$first.zip != -1)),],
           aes(x=length, y=first.zip, color=length)) +
    geom_boxplot(fill = "white", position = position_dodge(1), size = 0.5) +
    stat_summary(fun = mean, geom = "point", shape = 8, size = 3)+
    ggtitle("Time step of first zipped macrohomology for each fragment")
  ggsave(file,plot=first.zip.boxplot)
  
  
  file = paste(dirnew_plots,"/half_detection_boxplot.png",sep="")
  first.zip.boxplot <-
    ggplot(stats.zipping[c(which(stats.zipping$half.detect != -1)),],
           aes(x=length, y=half.detect, color=length)) +
    geom_boxplot(fill = "white", position = position_dodge(1), size = 0.5) +
    stat_summary(fun = mean, geom = "point", shape = 8, size = 3) +
    ggtitle("Time step of half detection for each invading fragment")
  ggsave(file,plot=first.zip.boxplot)
  
}

#########################################################################################################
#########################################################################################################
