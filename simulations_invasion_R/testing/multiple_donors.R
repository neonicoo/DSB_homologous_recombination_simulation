rm(list=ls())

###Set working directory
setwd("/home/nicolas/Documents/INSA/Stage4BiM/DSB_homologous_recombination_simulation/")


# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.12")
# 
# if (!require("Biostrings")){
#   BiocManager::install("Biostrings")
# }

library(ggplot2)
library(stringr)
library(seqinr)
library("Biostrings")

# Directory where you want to save timeseries and plots. Need the slash at the end if you want sub-directories underneath. 
rootdir = "/home/nicolas/Documents/INSA/Stage4BiM/DSB_homologous_recombination_simulation/datas/";

source("./simulations_invasion_R/testing/functions_donors.R")

# genome-wide microhomology counts
forward.sequences <- read.table("./LYS2/Occurences_per_8bp_motif(for+rev_donor).txt", sep="\t", header = TRUE)
forward.sequences = forward.sequences[,c("start", "sequence", "total")]
row.names(forward.sequences) = 1:nrow(forward.sequences)
microhomology.probs = forward.sequences$total / sum(forward.sequences$total)

# genome-wide microhomology counts but with bins of 10kb
sequences.bins <- read.csv("./LYS2/LY_occurences_per_8bp_(for_rev_donor)_with_bins.csv")

# Import the experimental contacts of the left DSB 10kb with the genome wide :
contacts <- read.csv("./LYS2/leftDSB_contacts_100000_110000_10kb.csv")
bins.id <- paste(as.character(contacts$chrom), "_", as.character(contacts$start_pos), "_", as.character(contacts$end_pos), sep="")
contacts <- cbind(contacts, bins.id)
colnames(contacts)[6] <- "frequency"
colnames(contacts)[7] <- "id"

# within-lys microhomologies (misalignments)
L500.self.micros = as.data.frame(matrix(c("aacaagct","aacaagct",98,319,319,98),2,3),stringsAsFactors = F); names(L500.self.micros) = c("L500", "position1", "position2"); L500.self.micros$position3 = NA
L1000.selfmicros <- read.delim("./LYS2/L1000_self-microhomologies.txt", stringsAsFactors=FALSE); L1000.selfmicros$position3 = NA
LY2000.selfmicros <- read.csv("./LYS2/LY2000_self-microhomologies.txt", sep="", stringsAsFactors=FALSE)

# Name the DNA sequences of the invading strands:
LY = (tolower("ATGACTAACGAAAAGGTCTGGATAGAGAAGTTGGATAATCCAACTCTTTCAGTGTTACCACATGACTTTTTACGCCCACAACAAGAACCTTATACGAAACAAGCTACATATTCGTTACAGCTACCTCAGCTCGATGTGCCTCATGATAGTTTTTCTAACAAATACGCTGTCGCTTTGAGTGTATGGGCTGCATTGATATATAGAGTAACCGGTGACGATGATATTGTTCTTTATATTGCGAATAACAAAATCTTAAGATTCAATATTCAACCAACGTGGTCATTTAATGAGCTGTATTCTACAATTAACAATGAGTTGAACAAGCTCAATTCTATTGAGGCCAATTTTTCCTTTGACGAGCTAGCTGAAAAAATTCAAAGTTGCCAAGATCTGGAAAGGACCCCTCAGTTGTTCCGTTTGGCCTTTTTGGAAAACCAAGATTTCAAATTAGACGAGTTCAAGCATCATTTAGTGGACTTTGCTTTGAATTTGGATACCAGTAATAATGCGCATGTTTTGAACTTAATTTATAACAGCTTACTGTATTCGAATGAAAGAGTAACCATTGTTGCGGACCAATTTACTCAATATTTGACTGCTGCGCTAAGCGATCCATCCAATTGCATAACTAAAATCTCTCTGATCACCGCATCATCCAAGGATAGTTTACCTGATCCAACTAAGAACTTGGGCTGGTGCGATTTCGTGGGGTGTATTCACGACATTTTCCAGGACAATGCTGAAGCCTTCCCAGAGAGAACCTGTGTTGTGGAGACTCCAACACTAAATTCCGACAAGTCCCGTTCTTTCACTTATCGCGACATCAACCGCACTTCTAACATAGTTGCCCATTATTTGATTAAAACAGGTATCAAAAGAGGTGATGTAGTGATGATCTATTCTTCTAGGGGTGTGGATTTGATGGTATGTGTGATGGGTGTCTTGAAAGCCGGCGCAACCTTTTCAGTTATCGACCCTGCATATCCCCCAGCCAGACAAACCATTTACTTAGGTGTTGCTAAACCACGTGGGTTGATTGTTATTAGAGCTGCTGGACAATTGGATCAACTAGTAGAAGATTACATCAATGATGAATTGGAGATTGTTTCAAGAATCAATTCCATCGCTATTCAAGAAAATGGTACCATTGAAGGTGGCAAATTGGACAATGGCGAGGATGTTTTGGCTCCATATGATCACTACAAAGACACCAGAACAGGTGTTGTAGTTGGACCAGATTCCAACCCAACCCTATCTTTCACATCTGGTTCCGAAGGTATTCCTAAGGGTGTTCTTGGTAGACATTTTTCCTTGGCTTATTATTTCAATTGGATGTCCAAAAGGTTCAACTTAACAGAAAATGATAAATTCACAATGCTGAGCGGTATTGCACATGATCCAATTCAAAGAGATATGTTTACACCATTATTTTTAGGTGCCCAATTGTATGTCCCTACTCAAGATGATATTGGTACACCGGGCCGTTTAGCGGAATGGATGAGTAAGTATGGTTGCACAGTTACCCATTTAACACCTGCCATGGGTCAATTACTTACTGCCCAAGCTACTACACCATTCCCTAAGTTACATCATGCGTTCTTTGTGGGTGACATTTTAACAAAACGTGATTGTCTGAGGTTACAAACCTTGGCAGAAAATTGCCGTATTGTTAATATGTACGGTACCACTGAAACACAGCGTGCAGTTTCTTATTTCGAAGTTAAATCAAAAAATGACGATCCAAACTTTTTGAAAAAATTGAAAGATGTCATGCCTGCTGGTAAAGGTATGTTGAACGTTCAGCTACTAGTTGTTAACAGGAACGATCGTACTCAAATATGTGGTATTGGCGAAATAGGTGAGATTTATGTTCGTGCAGGTGGTTTGGCCGAAGGTTATAGAGGATTACCAGAATTGAATAAAGAAAAATTTGTGAACAACTGGTTTGTTGAAAAAGATCACTGGAATTATTTGGATAAGGATAATGGTGAACCTTGGAGACAATTCTGGTTAGGTCCAAGAGATAGATTGTACAGAACGGGTGATTTAGGTCGTTATCTACCAAACGG"))
L = (tolower("ATGACTAACGAAAAGGTCTGGATAGAGAAGTTGGATAATCCAACTCTTTCAGTGTTACCACATGACTTTTTACGCCCACAACAAGAACCTTATACGAAACAAGCTACATATTCGTTACAGCTACCTCAGCTCGATGTGCCTCATGATAGTTTTTCTAACAAATACGCTGTCGCTTTGAGTGTATGGGCTGCATTGATATATAGAGTAACCGGTGACGATGATATTGTTCTTTATATTGCGAATAACAAAATCTTAAGATTCAATATTCAACCAACGTGGTCATTTAATGAGCTGTATTCTACAATTAACAATGAGTTGAACAAGCTCAATTCTATTGAGGCCAATTTTTCCTTTGACGAGCTAGCTGAAAAAATTCAAAGTTGCCAAGATCTGGAAAGGACCCCTCAGTTGTTCCGTTTGGCCTTTTTGGAAAACCAAGATTTCAAATTAGACGAGTTCAAGCATCATTTAGTGGACTTTGCTTTGAATTTGGATACCAGTAATAATGCGCATGTTTTGAACTTAATTTATAACAGCTTACTGTATTCGAATGAAAGAGTAACCATTGTTGCGGACCAATTTACTCAATATTTGACTGCTGCGCTAAGCGATCCATCCAATTGCATAACTAAAATCTCTCTGATCACCGCATCATCCAAGGATAGTTTACCTGATCCAACTAAGAACTTGGGCTGGTGCGATTTCGTGGGGTGTATTCACGACATTTTCCAGGACAATGCTGAAGCCTTCCCAGAGAGAACCTGTGTTGTGGAGACTCCAACACTAAATTCCGACAAGTCCCGTTCTTTCACTTATCGCGACATCAACCGCACTTCTAACATAGTTGCCCATTATTTGATTAAAACAGGTATCAAAAGAGGTGATGTAGTGATGATCTATTCTTCTAGGGGTGTGGATTTGATGGTATGTGTGATGGGTGTCTTGAAAGCCGGCGCAACCTTTTCAGTTATCGACCCTGCATATCCCCCAGCCAGACAAACCATTTACTTAGGTGTTGCTAAACCACGTGGGTTGATTGTTATTA"))
L500 = (tolower("ATGACTAACGAAAAGGTCTGGATAGAGAAGTTGGATAATCCAACTCTTTCAGTGTTACCACATGACTTTTTACGCCCACAACAAGAACCTTATACGAAACAAGCTACATATTCGTTACAGCTACCTCAGCTCGATGTGCCTCATGATAGTTTTTCTAACAAATACGCTGTCGCTTTGAGTGTATGGGCTGCATTGATATATAGAGTAACCGGTGACGATGATATTGTTCTTTATATTGCGAATAACAAAATCTTAAGATTCAATATTCAACCAACGTGGTCATTTAATGAGCTGTATTCTACAATTAACAATGAGTTGAACAAGCTCAATTCTATTGAGGCCAATTTTTCCTTTGACGAGCTAGCTGAAAAAATTCAAAGTTGCCAAGATCTGGAAAGGACCCCTCAGTTGTTCCGTTTGGCCTTTTTGGAAAACCAAGATTTCAAATTAGACGAGTTCAAGCATCATTTAGTGGACTTTGCTTTGAATTTGGATACCAG"))

ly.names = c("500", "1000", "2000")
ly.sequences = c(L500, L, LY)

#Import of the chr2.fa sequence file from the yeast genome (S288) :
yeast.genome<- read.fasta("./yeast-genome/S288c-R64-2-1-v2014/Genome_S288c.fa",
                          seqtype = 'DNA', as.string = TRUE, 
                          forceDNAtolower  = TRUE, set.attributes = FALSE)

donor <- LY

num.time.steps = 600 # Length of simulation in time steps
graph.resolution = 1 #save occupancy data at every nth time step. Plots will have this resolution at the x-axis 

test.replicates = 1 # How many times to simulate, replicates
kon.group<-c(0.4) #binding probabilities for every binding try
koff1.group<-c(0.2) # dissociation probabilities for each bound particle
koff2.group<-c(0.02) #dissociation probabilities for each zipped fragments
m.group = c(2) #bindings allowed to occur per tethering
search.window.group = c(250) #the genomic distance of the tethering effect (per side)
rad54.group <- c(1/200) #proportional to the lengh of invading strand
rdh54.group <- c(1/10) #proportional to the number of rad54


# Since the data needs to be outputted to files with human-readable names,we have to label the parameters with strings.
# For example 0005 is really 0.005
kon.group.names<- gsub("\\.", "", as.character(kon.group))
koff1.group.names<- gsub("\\.", "", as.character(koff1.group))
koff2.group.names<- gsub("\\.", "", as.character(koff2.group))
rad54.group.names<-gsub("\\.", "", as.character(rad54.group))

print(kon.group.names)
print(koff1.group.names)
print(koff2.group.names)

#########################################################################################################
######################################### Single run simulation ##########################################


# We have to check that the bins are the same between the 2 tables (sequences.bins and contacts) ;
# For example, for LY sequences.bins we have 2 more bins than in the contacts dataframe, so we remove them ;
# The comparison is made with the chromosome id and the start position for each bin 
chr_pos_occurences = c()
for (i in 2:ncol(sequences.bins)){
  chr_pos_occurences= c(chr_pos_occurences, 
                        paste(str_split(colnames(sequences.bins[i]), "_")[[1]][1], 
                              str_split(colnames(sequences.bins[i]), "_")[[1]][2], sep="_"))
}

chr_pos_contacts = c()
for (i in 1:length(bins.id)){
  chr_pos_contacts= c(chr_pos_contacts, 
                      paste(str_split(bins.id[i], "_")[[1]][1], 
                            str_split(bins.id[i], "_")[[1]][2], sep="_"))
}

remove = which(chr_pos_occurences %!in% chr_pos_contacts)+1
sequences.bins <- subset(sequences.bins, select=-remove)
rm(chr_pos_occurences, chr_pos_contacts, remove)

#Fusion the frequency of contact for each bins with the number of apparition for each microhomologies
sequences.contacts.bins = sequences.bins
for (i in 2:ncol(sequences.bins)){
  sequences.contacts.bins[i] = sequences.bins[i]*contacts$frequency[i-1]
}
sequences.contacts.bins = sequences.contacts.bins[-1] # Remove the "sequences" column
sequences.contacts.bins = apply(sequences.contacts.bins, 2, function(x) x[x!= ""]) #dataframe to matrix (reduce time complexity)
colnames(sequences.contacts.bins) = bins.id
rm(sequences.bins, contacts)

# kon = 2; koff = 3; m = 2; sw = 2; koff2 = 3
kon = 1; koff = 1; m = 1; sw = 1; koff2 = 1; rad54 = 1; rdh54 = 1; #for single Job run

kon.prob=kon.group[kon]
koff1.prob=koff1.group[koff]
koff2.prob=koff2.group[koff2]
rad54.prop=rad54.group[rad54]
rdh54.prop=rdh54.group[rdh54]

bindings.per.tethering = m.group[m]
search.window = search.window.group[sw]

kon.name=kon.group.names[kon]
koff1.name=koff1.group.names[koff]
koff2.name=koff2.group.names[koff2]
rad54.name=rad54.group.names[rad54]
rdh54.name=gsub("\\.", "", as.character(rad54.prop*rdh54.prop))

# Initialize the occupied.rad51 vector, genomic (start) position of RAD51 particles (bp / 8) of invaded strand ;
occupied.rad51 = list(bound = "unbound",strand = "negative", genome.bins = c(), donor.invasions = 473927 - 368, lys2.microhomology = 368)

#population time series only for the zipped fragment 
pop.time.series.zip = as.data.frame(matrix(0,num.time.steps*3,3))
names(pop.time.series.zip) = c("time.step","prob.detect","length")
pop.time.series.zip$length = rep(ly.names, each = num.time.steps)
pop.time.series.zip$time.step = rep(seq(1,num.time.steps,1),3)

#population time series for all the microhomologies (not only zipped)
pop.time.series.all = as.data.frame(matrix(0,num.time.steps*3,3))
names(pop.time.series.all) = c("time.step","prob.detect","length")
pop.time.series.all$length = rep(ly.names, each = num.time.steps)
pop.time.series.all$time.step = rep(seq(1,num.time.steps,1),3)

# Create a variable, saver, that keeps track of how many individual simulations you want to save. 
saver = 0
bigtracker = 0

lys.occupancy.firsts = as.data.frame(matrix(-1, 3*test.replicates, 4))
names(lys.occupancy.firsts) = c("length", "first.bound", "twoh.bound", "first.twoh.time.diff")
lys.occupancy.firsts$length = rep(ly.names, times = test.replicates)

stats.zipping = as.data.frame(matrix(-1, 3*test.replicates, 3))
names(stats.zipping) = c("length", "first.zip", "half.detect")
stats.zipping$length = rep(ly.names, times = test.replicates)


chromosome.contacts <- as.data.frame(matrix(0,num.time.steps*3, length(bins.id)+2))
colnames(chromosome.contacts) = c("time.step", "length", bins.id)
chromosome.contacts$time.step = rep(seq(1,num.time.steps,1),3)
chromosome.contacts$length = rep(ly.names, each = num.time.steps)

####### Directory settings ##################
dirname=paste(num.time.steps, kon.name, koff1.name, koff2.name, 
              bindings.per.tethering, search.window, rad54.name, rdh54.name, sep="_")

dirnew=paste(rootdir,dirname,sep="")
dir.create(dirnew)

# Directory to save single runs
dirnew_singles=paste(dirnew,"/single_runs",sep="")
dir.create(dirnew_singles)

# Directory to save the first contact, 200 contact, and time diff between 1st and 200 files
dirnew_data = paste(dirnew,"/data",sep="")
dir.create(dirnew_data)

# Directory to save the 1st, 200 contact times, and the 200-1st time difference
dirnew_plots = paste(dirnew,"/plots",sep="")
dir.create(dirnew_plots)

print(dirnew)

#############################################


# Now make all the replicates
for (trial in 1:test.replicates){ 
  # initialize tabulation of bound microhomologies/heterologies
  # print(trial)
  
  if(saver < 3){
    ly.binding.ts = as.data.frame(matrix(0, (num.time.steps/graph.resolution)*3,5))
    names(ly.binding.ts) = c('time.step', 'bound', "length", "heterologies", "homologies")
    ly.binding.ts$time.step = rep(seq(1,num.time.steps, graph.resolution),3)
    ly.binding.ts$length = rep(ly.names, each = (num.time.steps / graph.resolution))
  }
  
  for (fragment in 1:3){
    # initialize the ID and state of the current invading strand 
    bigtracker = bigtracker +1
    ly.type = ly.names[fragment]
    lys2.fragment = ly.sequences[fragment]
    if (fragment == 1){
      self.micros =L500.self.micros
    }else{
      if(fragment == 2){
        self.micros = L1000.selfmicros
      }else{
        self.micros = LY2000.selfmicros
      }
    }
    
    
    donors.list = donors.generator(template = LY, bins = bins.id, N=5)
    current.donor = ""
    
    SEI.binding.tries = floor((nchar(lys2.fragment)-7)/8)
    
    occupied.rad51$bound = "unbound"
    occupied.rad51$genome.bins = c()
    occupied.rad51$donor.invasions = c()
    occupied.rad51$lys2.microhomology = c()
    
    donors.occupancy = as.data.frame(matrix(0, nchar(lys2.fragment),5))
    names(donors.occupancy) = c('bp', 'bound', "bound.id", "donor.id", "zipped")
    donors.occupancy$bp = 1:nchar(lys2.fragment)
    donors.occupancy$bound = "no"
    donors.occupancy$bound.id = "unbound"
    donors.occupancy$donor.id = "unknown"
    donors.occupancy$zipped = "no"
    
    first.lys = 0 #the first homology bound to the donor
    twoh.lys = 0 # the first two hundred homologies bound to the donor
    
    start.zipping <- 0 #enable or not the zipping phase
    first.zip <- 0 #the first zipped fragment to the donor
    half.detect <- 0 #when the probability detection is equal to 0.5
    
    # We have to place randomly some rad54 and rdh54 in the invading fragment to induce the zipping ;
    # The number of rad54 depends of the length of the fragment,
    # and the number of rdh54 depends of the number of rad54 ;
    nb.rad54 <- floor(rad54.prop*str_length(lys2.fragment)) #number of rad54 to be placed into the invading strand ;
    nb.rdh54 <- floor(rdh54.prop*nb.rad54)+1 # number of rdh54 to be placed into the invading strand;
    rad54.rdh54.locations <- rad54.rdh54.placement(nb.rad54, nb.rdh54) 
    pos.rad54 <- rad54.rdh54.locations[[1]] #positions of rad54 in the invading strand;
    pos.rdh54 <- rad54.rdh54.locations[[2]] #positions of rdh54 in the invading strand;
    
    zipped.fragments.list <- as.data.frame(matrix(0,0,3)) #all the macrohomologies after zipping with start/end positions
    names(zipped.fragments.list ) = c("start", "end", "sequences")
    
    unzipped.rad54 <- pos.rad54 #positions of non-overlapped rad54
    
    # Loop through the time-steps
    for (time.step in 1:num.time.steps){
      if(kon.prob == 0){
        next
      }
      
      if (occupied.rad51$bound != "unbound"){
        if (length(occupied.rad51$donor.invasions) != sum(occupied.rad51$donor.invasions == "H")){
          new.bindings = new.microhomologizer(occupied.rad51, search.window, bindings.per.tethering)
          occupied.rad51$genome.bins = c(occupied.rad51$genome.bins, new.bindings$genome.bins)
          occupied.rad51$donor.invasions = c(occupied.rad51$donor.invasions,new.bindings$donor.invasions)
          occupied.rad51$lys2.microhomology = c(occupied.rad51$lys2.microhomology, new.bindings$lys2.microhomology)
        }
      }

      new.bindings = genome.wide.sei(SEI.binding.tries)

      if (occupied.rad51$bound == "unbound"){
        occupied.rad51 = new.bindings;
        if(length(new.bindings$lys2.microhomology) > 0){
          occupied.rad51$bound = "bound"
        }
        
      }else{
        # print("bound and adding")
        occupied.rad51$genome.bins = c(occupied.rad51$genome.bins, new.bindings$genome.bins)
        occupied.rad51$donor.invasions = c(occupied.rad51$donor.invasions, new.bindings$donor.invasions)
        occupied.rad51$lys2.microhomology = c(occupied.rad51$lys2.microhomology, new.bindings$lys2.microhomology)
      }
      
      #simulate random dissociation(s)
      num.bound = length(occupied.rad51$donor.invasions)
      #  print(num.bound)
      preserved = sample(c(FALSE,TRUE), num.bound, replace = TRUE, prob = c(koff1.prob,1-koff1.prob)) #dissociate if FALSE
      preserved[which(occupied.rad51$lys2.microhomology[preserved] %in% donors.occupancy$zipped=="yes")] = TRUE #koff1 can't dissociate zipped homologies
      
      occupied.rad51$genome.bins = occupied.rad51$genome.bins[preserved]
      occupied.rad51$donor.invasions = occupied.rad51$donor.invasions[preserved]
      occupied.rad51$lys2.microhomology = occupied.rad51$lys2.microhomology[preserved]
      
      if (sum(!preserved)==num.bound){
        occupied.rad51$bound = "unbound"
      }
      
      
      donors.occupancy$bound = "no"
      donors.occupancy$bound.id = "unbound"
      donors.occupancy$donor.id = "unknown"
      
      if(occupied.rad51$bound != "unbound"){
        for (i in 1:length(occupied.rad51$lys2.microhomology)){
          donors.occupancy$bound[occupied.rad51$lys2.microhomology[i]:(occupied.rad51$lys2.microhomology[i] + 7)] = "yes"
          if(occupied.rad51$donor.invasions[i] != "H"){
            if (substr(occupied.rad51$donor.invasions[i], 1, 4) == "!LYS"){
              donors.occupancy$bound.id[occupied.rad51$lys2.microhomology[i]:(occupied.rad51$lys2.microhomology[i] + 7)]  = "homology"
              donors.occupancy$donor.id[occupied.rad51$lys2.microhomology[i]:(occupied.rad51$lys2.microhomology[i] + 7)] = occupied.rad51$donor.invasions[i]
              
            }else{
              donors.occupancy$bound.id[occupied.rad51$lys2.microhomology[i]:(occupied.rad51$lys2.microhomology[i] + 7)]  = "homology"
              donors.occupancy$donor.id[occupied.rad51$lys2.microhomology[i]:(occupied.rad51$lys2.microhomology[i] + 7)] = "LYS"
            }
            
          }else if(occupied.rad51$donor.invasions[i] == "H"){
            donors.occupancy$bound.id[occupied.rad51$lys2.microhomology[i]:(occupied.rad51$lys2.microhomology[i] + 7)]  = "heterology"
            donors.occupancy$donor.id[occupied.rad51$lys2.microhomology[i]:(occupied.rad51$lys2.microhomology[i] + 7)] = "unknown"
          }
        }
        
        donors.occupancy$bound[which(donors.occupancy$bound == "no" & donors.occupancy$zipped == "yes") ] = "yes"
        donors.occupancy$bound.id[which(donors.occupancy$bound.id == "unbound" & donors.occupancy$zipped == "yes") ] = "homology"
        donors.occupancy$donor.id[which(donors.occupancy$donor.id == "unknown" & donors.occupancy$zipped == "yes") ] = current.donor
      }
      
      # When the twoh microhomology state is enable, the zipping occurs until all rad54 are zipped;
      if(length(unzipped.rad54 > 0) && current.donor != "" && 
         donors.list$invasion[which(donors.list$id == current.donor)] != "wrong"){
        
        if (donors.list$invasion[which(donors.list$id == current.donor)] =="no"){
          donors.list$invasion[which(donors.list$id == current.donor)] =="yes"
        } 
        
        for (pos in unzipped.rad54){
          # Check if the sequence to zip is big enough ;
          # We decided >= 16 (2*8 nts) arbitrary (could be more or less)
          if(donors.occupancy$zipped[pos] != "yes" & check.before.zipping(pos, donor = current.donor) >= 16){
            new.zip = zipping(pos, zipped.fragments.list, donor= current.donor)
            
            if(length(new.zip) > 1){
              unzipped.rad54 = unzipped.rad54[which(unzipped.rad54 != pos)]
              zipped.fragments.list = rbind(zipped.fragments.list, new.zip)
              names(zipped.fragments.list) = c("start", "end", "sequences")
              current.zip.start <- as.integer(new.zip[1])
              current.zip.end <- as.integer(new.zip[2])
              donors.occupancy$zipped[current.zip.start : current.zip.end] = "yes"
              
            }else if (new.zip == -1){
              print(c("FAILED", ly.type, current.donor))
              
              remove.rad51 = which(occupied.rad51$donor.invasions == current.donor)
              occupied.rad51$genome.bins = occupied.rad51$genome.bins[-remove.rad51]
              occupied.rad51$donor.invasions = occupied.rad51$donor.invasions[-remove.rad51]
              occupied.rad51$lys2.microhomology = occupied.rad51$lys2.microhomology[-remove.rad51]
              
              donors.occupancy$bound[which(donors.occupancy$donor.id==current.donor)] = "no"
              donors.occupancy$bound.id[which(donors.occupancy$donor.id==current.donor)] = "unbound"
              donors.occupancy$zipped[which(donors.occupancy$donor.id==current.donor)] = "no"
              
              donors.occupancy$donor.id[which(donors.occupancy$donor.id==current.donor)] = "unknown"
              zipped.fragments.list <- as.data.frame(matrix(0,0,3))
              names(zipped.fragments.list ) = c("start", "end", "sequences")
              unzipped.rad54 = pos.rad54
              donors.list$invasion[which(donors.list$id == current.donor)] = "wrong"
              current.donor = ""
            }
          }
        }
      }
      
      #print(c(ly.type, current.donor))
      
      # # Introduce at each time step, the probability of dissociation Koff2 for zipped sequences ;
      # # If a macrohomology becomes un-zipped because of koff2,
      # # All the processes of homologies searching and zipping have to be done again ;
      # 
      # if(koff2.prob > 0 & dim(zipped.fragments.list)[1] != 0){
      #   row2remove <- c()
      #   for(i in 1:nrow(zipped.fragments.list)){
      #     preserved.zip <- sample(c(FALSE, TRUE), size =1, replace = TRUE, prob = c(koff2.prob,1-koff2.prob))
      #     if(!preserved.zip){
      #       current.zip.start <- as.integer(zipped.fragments.list[i, ]$start)
      #       current.zip.end <- as.integer(zipped.fragments.list[i, ]$end)
      #       row2remove = c(row2remove, i)
      # 
      #       donors.occupancy$zipped[current.zip.start : current.zip.end] = "no" #the sequence is unzipped
      #       donors.occupancy$bound[current.zip.start : current.zip.end] = "no" #the sequence becomes unbound to donor
      #       donors.occupancy$bound.id[current.zip.start : current.zip.end] = "unbound" # the sequence is considered as heterologous again
      #       unzipped.rad54 = c(unzipped.rad54, current.zip.start) #the rad54 into the sequence are no more overlapped by any microhomology
      # 
      #       remove.rad51 <- which(occupied.rad51$lys2.microhomology %in% (current.zip.start : current.zip.end))
      # 
      #       #remove binding sites from the donor
      #       occupied.rad51$genome.bins = occupied.rad51$genome.bins[-remove.rad51]
      #       occupied.rad51$lys2.microhomology = occupied.rad51$lys2.microhomology[-remove.rad51]
      #       occupied.rad51$donor.invasions = occupied.rad51$donor.invasions[-remove.rad51]
      # 
      #       if(length(occupied.rad51$donor.invasions) == 0 | length(occupied.rad51$lys2.microhomology) == 0){
      #         occupied.rad51$bound = "unbound"
      #         break
      #       }
      #     }
      #   }
      #   if(length(row2remove) > 0){
      #     zipped.fragments.list = zipped.fragments.list[-c(row2remove),]
      #     if(dim(zipped.fragments.list)[1] != 0){
      #       row.names(zipped.fragments.list) = (1:nrow(zipped.fragments.list))
      #     }
      #   }
      # }
      
      if(length(which(donors.occupancy$bound.id == "homology" & donors.occupancy$donor.id == "LYS")) > 0 && first.lys == 0){
        first.lys = 1;
        lys.occupancy.firsts$first.bound[bigtracker] = time.step
      }
      
      if(length(which(donors.occupancy$bound.id == "homology" & donors.occupancy$donor.id == "LYS")) >= 200 && start.zipping == 0){
        if(twoh.lys == 0){
          twoh.lys = 1
          lys.occupancy.firsts$twoh.bound[bigtracker] = time.step
          lys.occupancy.firsts$first.twoh.time.diff[bigtracker] = time.step - lys.occupancy.firsts$first.bound[bigtracker]
        }
        start.zipping = 1
        
      }else if(length(which(donors.occupancy$bound.id == "homology")) < 200 && start.zipping == 1){
        start.zipping = 0
      }
      
      if(current.donor == ""){
        for (candidate.donor in unique(donors.occupancy$donor.id)){
          if(candidate.donor != "unknown" && length(which(donors.occupancy$donor.id == candidate.donor)) >= 200 && 
             donors.list$invasion[which(donors.list$id == candidate.donor)] == "no"){
            
            current.donor = candidate.donor
          }
        }
      }
      
      if(occupied.rad51$bound != "unbound"){
        # The probability of SEI detection depends of the number of zipped nts ;
        prob.detection.lys.zip = length(which(donors.occupancy$zipped == "yes" & donors.occupancy$donor.id == "LYS"))
        prob.detection.lys.zip = prob.detection.lys.zip/500 #take into accout the crosslink density
        
        prob.detection.lys.all = length(which(donors.occupancy$bound.id == "homology" & donors.occupancy$donor.id == "LYS"))
        prob.detection.lys.all = prob.detection.lys.all/500 
        
      }else{
        prob.detection.lys.zip = 0
        prob.detection.lys.all = 0
      }
      
      if(prob.detection.lys.zip >= 1){prob.detection.lys.zip = 1}
      if(prob.detection.lys.all >= 1){prob.detection.lys.all = 1}
      
      if(length(which(donors.occupancy$zipped == "yes")) > 0 && first.zip == 0){
        first.zip = 1
        stats.zipping$first.zip[bigtracker] = time.step
      }
      
      if(prob.detection.lys.zip > 0.5 & half.detect == 0){
        half.detect = 1
        stats.zipping$half.detect[bigtracker] = time.step
      }
      
      
      if(saver <3 ){
        # tabulate occupancies vs. time step and length
        ly.binding.ts$bound[ly.binding.ts$time.step == time.step & 
                              ly.binding.ts$length == ly.type] = length(which(donors.occupancy$bound == "yes"))
        
        ly.binding.ts$heterologies[ly.binding.ts$time.step == time.step & 
                                     ly.binding.ts$length == ly.type] = length(which(donors.occupancy$bound.id == "heterology"))
        
        #ly.binding.ts$homologies = ly.binding.ts$bound - ly.binding.ts$heterologies
        ly.binding.ts$homologies[ly.binding.ts$time.step == time.step & 
                                   ly.binding.ts$length == ly.type] = length(which(donors.occupancy$donor.id == "LYS"))
      }
      
      pop.time.series.zip$prob.detect[pop.time.series.zip$time.step == time.step & 
                                        pop.time.series.zip$length == ly.type] = 
        pop.time.series.zip$prob.detect[pop.time.series.zip$time.step == time.step & 
                                          pop.time.series.zip$length == ly.type] + prob.detection.lys.zip
      
      pop.time.series.all$prob.detect[pop.time.series.all$time.step == time.step & 
                                        pop.time.series.all$length == ly.type] = 
        pop.time.series.all$prob.detect[pop.time.series.all$time.step == time.step & 
                                          pop.time.series.all$length == ly.type] + prob.detection.lys.all
      
      for (i in 1:length(table(occupied.rad51$genome.bins))){
        bin <- names(table(occupied.rad51$genome.bins))[i]
        count <- table(occupied.rad51$genome.bins)[[i]]
        
        chromosome.contacts[chromosome.contacts$time.step == time.step &
                              chromosome.contacts$length == ly.type, names(chromosome.contacts) == bin] =
          chromosome.contacts[chromosome.contacts$time.step == time.step &
                                chromosome.contacts$length == ly.type, names(chromosome.contacts) == bin] + count
      }
      
    }#next time step
  }#next fragment
  
  fname = paste("timeseries", num.time.steps, kon.name, koff1.name, koff2.name, bindings.per.tethering, search.window, sep="_")
  fname = paste(fname,"_trial",as.character(trial),".txt",sep="")
  write.table(ly.binding.ts,file=paste(dirnew,"/", fname, sep = ""))
  
  if(saver < 3){
    
    names(ly.binding.ts)[3] = "length"
    ly.binding.ts$length = factor(ly.binding.ts$length)
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
  
  saver=saver+1
}#end process

write.csv(chromosome.contacts, file=paste(dirnew_data,"/chromosomes_contacts.csv",sep=""))

# scaled.chromosome.contacts = chromosome.contacts
# scaled.chromosome.contacts[,3:dim(scaled.chromosome.contacts)[2]] = round(scale(scaled.chromosome.contacts[,3:dim(scaled.chromosome.contacts)[2]]), 3)
# write.csv(scaled.chromosome.contacts, file=paste(dirnew_data,"/scaled_chromosomes_contacts.csv",sep=""))

# population timeseries
write.table(pop.time.series.zip, file=paste(dirnew_data,"/population_timeseries_zip.txt",sep=""))

outname=paste(dirnew_plots,"/population_time_series_zip.png",sep="")
pop.plot<-
  ggplot(data = pop.time.series.zip) + geom_step(aes(x = time.step, y = prob.detect, color = length)) +
  labs(x = "time step", y = "Probability of Detection") + theme_minimal()+ theme(text = element_text(size = 16))+
  scale_y_continuous(limits = c(0, test.replicates))
ggsave(outname,plot=pop.plot)


write.table(pop.time.series.all, file=paste(dirnew_data,"/population_timeseries_all.txt",sep=""))

outname=paste(dirnew_plots,"/population_time_series_all.png",sep="")
pop.plot<-
  ggplot(data = pop.time.series.all) + geom_step(aes(x = time.step, y = prob.detect, color = length)) +
  labs(x = "time step", y = "Probability of Detection") + theme_minimal()+ theme(text = element_text(size = 16))+
  scale_y_continuous(limits = c(0, test.replicates))
ggsave(outname,plot=pop.plot)

#### Histograms+ boxplot for the first and twoh MH

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

final.firsts$`500` = lys.occupancy.firsts$twoh.bound[which(lys.occupancy.firsts$length == 500)]
final.firsts$`1000` = lys.occupancy.firsts$twoh.bound[which(lys.occupancy.firsts$length == 1000)]
final.firsts$`2000` = lys.occupancy.firsts$twoh.bound[which(lys.occupancy.firsts$length == 2000)]

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

final.firsts$`500` = lys.occupancy.firsts$first.twoh.time.diff[which(lys.occupancy.firsts$length == 500)]
final.firsts$`1000` = lys.occupancy.firsts$first.twoh.time.diff[which(lys.occupancy.firsts$length == 1000)]
final.firsts$`2000` = lys.occupancy.firsts$first.twoh.time.diff[which(lys.occupancy.firsts$length == 2000)]

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


#### Zipping detection

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


rm(het_plot, lys2_plot, occ_plot, pop.plot, first.boxplot, first.hist, first.zip.boxplot)
