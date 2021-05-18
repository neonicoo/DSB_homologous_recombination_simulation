rm(list=ls()) #clean global environment

###Set working directory
setwd("/home/nicolas/Documents/INSA/Stage4BiM/DSB_homologous_recombination_simulation/")

# Directory where you want to save timeseries and plots. Need the slash at the end if you want sub-directories underneath. 
rootdir = "/home/nicolas/Documents/INSA/Stage4BiM/DSB_homologous_recombination_simulation/datas/"


###################### Import librairies #######################################

library(ggplot2)
library(stringr)

################################################################################
############################## Import the datas ################################

source("./simulations_invasion_R/testing/simulation_functions_test.R") #Run the file containing the functions for the simulation 
source("./simulations_invasion_R/testing/outputs_functions_test.R") #Run the file containing the functions to create the outputs graphs

# genome-wide microhomology counts
forward.sequences <- read.table("./LYS2/Occurences_per_8bp_motif(for+rev_donor).txt", sep="\t", header = TRUE)
forward.sequences = forward.sequences[,c("start", "sequence", "total")]
row.names(forward.sequences) = 1:nrow(forward.sequences)
microhomology.probs = forward.sequences$total / sum(forward.sequences$total)

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

# genome-wide microhomology counts but with bins of 10kb
sequences.bins <- read.csv("./LYS2/LY_occurences_per_8bp_(for_rev_donor)_with_bins.csv")

# Import the experimental contacts of the left DSB 10kb with the genome wide :
contacts <- read.csv("./LYS2/leftDSB_contacts_100000_110000_10kb.csv")
bins.id <- paste(as.character(contacts$chrom), "_", as.character(contacts$start_pos), "_", as.character(contacts$end_pos), sep="")
contacts <- cbind(contacts, bins.id)
colnames(contacts)[6] <- "frequency"
colnames(contacts)[7] <- "id"

################################################################################
################### Creation of the chromosome contact frequency matrix ########

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

################################################################################
####################### Parameters #############################################

num.time.steps = 600 # Length of simulation in time steps
graph.resolution = 1 #save occupancy data at every nth time step. Plots will have this resolution at the x-axis 

test.replicates = 1 # How many times to simulate, replicates
kon.group<-c(0.4) #binding probabilities for every binding try
koff1.group<-c(0.2) # dissociation probabilities for each bound particle
koff2.group<-c(0.05) #dissociation probabilities for each zipped fragments
m.group = c(2) #bindings allowed to occur per tethering
search.window.group = c(250) #the genomic distance of the tethering effect (per side)
rad54.group <- c(1/125) #proportional to the lengh of invading strand
rdh54.group <- c(1/8) #proportional to the number of rad54
additional.donors <- 2

# Since the data needs to be outputted to files with human-readable names,we have to label the parameters with strings.
# For example 0005 is really 0.005
kon.group.names<- gsub("\\.", "", as.character(kon.group))
koff1.group.names<- gsub("\\.", "", as.character(koff1.group))
koff2.group.names<- gsub("\\.", "", as.character(koff2.group))
rad54.group.names<-gsub("\\.", "", as.character(rad54.group))

print(kon.group.names)
print(koff1.group.names)
print(koff2.group.names)

################################################################################

################################################################################
############################ Single run simulation #############################

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

# Initialize the occupied.rad51 list that contains :
#   the state of the invading fragment : bound or unbound ,
#   the nature of the invading strand : negative ,
#   the respective genome bins were heterologies and homologies are found,
#   the nature of the donor :
#     "H" if heterology ;
#     "!LYS#" if it's a donor (in the donors.list$id)
#     "473927 - position of rad51 in invading fragment" if the donor is LYS2
#   the position of rad51 in invading fragment for each bound microhomologies ;

occupied.rad51 = list(bound = "unbound",strand = "negative", genome.bins = c("chr2_470001_480001"), donor.invasions = 473927, lys2.microhomology = 368)

# Generate N additional donors (other than lys2),
#   with a random number of mutations (between 10% and 40% snp ),
#   mutation occurs in random position of the sequence ,
#   the template sequence is LY.
# Change N to change the number of potential donors ;
donors.list = donors.generator(template = LY, bins = bins.id, N=additional.donors)

saver = 0  #keeps track of how many individual simulations you want to save ;
# Bigtracker : tracker to know the in which replicate we are for which fragment, 
#   ex : for L1000 (2nd fragment of the list), for the 3th replicate, bigtracker = 8
bigtracker = 0 

####################" Initialize the population time series ####################

# Population time series are dataframes that allows us to see the provability to detect an interaction between 
#   the invading strand and the donor(s). This probability is proportional to the number of test replicates (inner runs) we make.

# In the case we introduce additional donors (not only lys2) in the simulation,
#   we need to create 2 more population time series in order to distinguish 
#   the general behavior of the simulation (with inclusion of all homologies for all donors) 
#   from the behavior of lys2 associations ;

pop.time.series = as.data.frame(matrix(0.0,num.time.steps*3,4))
names(pop.time.series) = c("time.step","length", "homologies", "zip")
pop.time.series$time.step = rep(seq(1,num.time.steps,1),3)
pop.time.series$length = rep(ly.names, each = num.time.steps)

for (i in 1:(additional.donors+1)) {
  col = paste("prob.detect", as.character(donors.list$id[i]), sep=".")
  pop.time.series[col] = rep(0, 3*num.time.steps)
}

################################################################################
############# Some other statistics dataframes #################################

# Saves the time step of each fragment sizes, at each replicates for 
#   the first homology with lys2,
#   the two hundredth (also called twoh) homology with lys2, 
#   the number of time steps between the first and the two hundredth ;

lys.occupancy.firsts = as.data.frame(matrix(-1, 3*test.replicates, 4))
names(lys.occupancy.firsts) = c("length", "first.bound", "twoh.bound", "first.twoh.time.diff")
lys.occupancy.firsts$length = rep(ly.names, times = test.replicates)

# Dataframe with the number of time each bins for each chromosome is contacted during the searching phase 
chromosome.contacts <- as.data.frame(matrix(0,num.time.steps*3, length(bins.id)+2))
colnames(chromosome.contacts) = c("time.step", "length", bins.id)
chromosome.contacts$time.step = rep(seq(1,num.time.steps,1),3)
chromosome.contacts$length = rep(ly.names, each = num.time.steps)

################################################################################
########################### Directory settings #################################

dirname=paste(num.time.steps, kon.name, koff1.name, koff2.name, 
              bindings.per.tethering, search.window, rad54.name, rdh54.name, additional.donors, sep="_")

dirnew=paste(rootdir,dirname,sep="")

if(file.exists(dirnew)){
  unlink(dirnew, recursive = TRUE)
}
dir.create(dirnew)

#Logfile to have a traceback of the parameters we used for the simulation
sink(paste(dirnew, "/logfile.txt", sep=""))
cat("Number of time steps : ", num.time.steps, "\n")
cat("Number of replicates : ", test.replicates, "\n")
cat("Number of potential donors :", additional.donors+1, "\n")
cat("Kon : ", kon.prob, "\n")
cat("Koff1 : ", koff1.prob, "\n")
cat("Koff2 : ", koff2.prob, "\n")
cat("Tethering per windows : ", bindings.per.tethering, "\n")
cat("Search windows : ", search.window, "\n")
cat("Proportion of rad54 : ", rad54.prop, "\n")
cat("Proportion of rdh54 : ", rdh54.prop, "\n")
sink()


dirnew_plots = paste(dirnew, "/plots", sep="")
dir.create(dirnew_plots)

# Directory to save contact statistics plots
dirnew_contacts = paste(dirnew_plots, "/contacts", sep="")
dir.create(dirnew_contacts)

# Directory to save population time series plots
dirnew_pop = paste(dirnew_plots, "/population_timeseries", sep="")
dir.create(dirnew_pop)

# Directory to save single runs
dirnew_singles=paste(dirnew_plots,"/single_runs",sep="")
dir.create(dirnew_singles)

# Directory to save the first contact, 200 contact, and time diff between 1st and 200 files
dirnew_data = paste(dirnew,"/data",sep="")
dir.create(dirnew_data)

dirnew_timeseries = paste(dirnew,"/timeseries",sep="")
dir.create(dirnew_timeseries)

print(dirnew)

################################################################################
############################## Start simulation ################################

# Now make all the replicates
for (trial in 1:test.replicates){ 
  # initialize tabulation of bound microhomologies/heterologies
  # print(trial)
  
  if(saver < 3){
    binding.ts = as.data.frame(matrix(0, (num.time.steps/graph.resolution)*3,5))
    names(binding.ts) = c('time.step', "length", "total.bound", "heterologies", "homologies")
    binding.ts$time.step = rep(seq(1,num.time.steps, graph.resolution),3)
    binding.ts$length = rep(ly.names, each = (num.time.steps / graph.resolution))
    
    for (i in 1:(additional.donors+1)) {
      col = paste("bound", as.character(donors.list$id[i]), sep=".")
      binding.ts[col] = rep(0, 3*num.time.steps/graph.resolution)
    }
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
    
    current.donor = ""
    donors.blacklist = c()
    donors.list$invasion = rep("no", additional.donors+1)
    
    SEI.binding.tries = floor((nchar(lys2.fragment)-7)/8)
    
    occupied.rad51$bound = "unbound"
    occupied.rad51$genome.bins = c()
    occupied.rad51$donor.invasions = c()
    occupied.rad51$lys2.microhomology = c()
    
    # State of the invading fragment (occupancy) with donor(s)
    donors.occupancy = as.data.frame(matrix(0, nchar(lys2.fragment),5))
    names(donors.occupancy) = c('bp', 'bound', "bound.id", "donor.id", "zipped")
    donors.occupancy$bp = 1:nchar(lys2.fragment)
    donors.occupancy$bound = "no"
    donors.occupancy$bound.id = "unbound"
    donors.occupancy$donor.id = "unknown"
    donors.occupancy$zipped = "no"
    
    first.lys = 0 #the first homology bound to lys2
    twoh.lys = 0 # the first two hundred homologies bound to lys2
    
    first.zip <- 0 #the first zipped fragment to lys2
    half.detect <- 0 #when the probability detection is equal to 0.5 for zipped fragment to lys2 donor
    
    # We have to place randomly some rad54 and rdh54 in the invading fragment to induce the zipping ;
    # The number of rad54 depends of the length of the fragment,
    # and the number of rdh54 depends of the number of rad54 ;
    nb.rad54 <- floor(rad54.prop*str_length(lys2.fragment)) #number of rad54 to be placed into the invading strand ;
    nb.rdh54 <- floor(rdh54.prop*nb.rad54)+1 # number of rdh54 to be placed into the invading strand;
    rad54.rdh54.locations <- rad54.rdh54.placement(nb.rad54, nb.rdh54, lys2.fragment = lys2.fragment) 
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
      
      # Seach homologies in the binding tethering window : new.microhomologizer 
      if (occupied.rad51$bound != "unbound"){
        if (length(occupied.rad51$donor.invasions) != sum(occupied.rad51$donor.invasions == "H")){
          new.bindings = new.microhomologizer(occupied.rad51, search.window, bindings.per.tethering)
          occupied.rad51$genome.bins = c(occupied.rad51$genome.bins, new.bindings$genome.bins)
          occupied.rad51$donor.invasions = c(occupied.rad51$donor.invasions,new.bindings$donor.invasions)
          occupied.rad51$lys2.microhomology = c(occupied.rad51$lys2.microhomology, new.bindings$lys2.microhomology)
        }
      }
      
      # Search new homologies for the free sites on the invading fragment
      new.bindings = genome.wide.sei(SEI.binding.tries)
      
      if (occupied.rad51$bound == "unbound"){
        occupied.rad51 = new.bindings
        if(length(new.bindings$lys2.microhomology) > 0){
          occupied.rad51$bound = "bound"
        }
        
      }else{
        # print("bound and adding")
        occupied.rad51$genome.bins = c(occupied.rad51$genome.bins, new.bindings$genome.bins)
        occupied.rad51$donor.invasions = c(occupied.rad51$donor.invasions, new.bindings$donor.invasions)
        occupied.rad51$lys2.microhomology = c(occupied.rad51$lys2.microhomology, new.bindings$lys2.microhomology)
      }
      
      ###################### KOFF1 #############################################
      #simulate random dissociation(s)
      num.bound = length(occupied.rad51$donor.invasions)
      # print(num.bound)
      preserved = sample(c(FALSE,TRUE), num.bound, replace = TRUE, prob = c(koff1.prob,1-koff1.prob)) #dissociate if FALSE
      preserved[which(occupied.rad51$lys2.microhomology[preserved] %in% donors.occupancy$zipped=="yes")] = TRUE #koff1 can't dissociate zipped homologies
      
      occupied.rad51$genome.bins = occupied.rad51$genome.bins[preserved]
      occupied.rad51$donor.invasions = occupied.rad51$donor.invasions[preserved]
      occupied.rad51$lys2.microhomology = occupied.rad51$lys2.microhomology[preserved]
      
      if (sum(!preserved)==num.bound){
        occupied.rad51$bound = "unbound"
      }
      
      ##########################################################################
      ############################## Ocuupancy #################################
      donors.occupancy$bound = "no"
      donors.occupancy$bound.id = "unbound"
      donors.occupancy$donor.id = "unknown"
      
      if(occupied.rad51$bound != "unbound"){
        for (i in 1:length(occupied.rad51$lys2.microhomology)){
          donors.occupancy$bound[occupied.rad51$lys2.microhomology[i]:(occupied.rad51$lys2.microhomology[i] + 7)] = "yes"
          if(occupied.rad51$donor.invasions[i] != "H"){
            if (str_detect(pattern = "donor", occupied.rad51$donor.invasions[i])){
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
        
        donors.occupancy$bound[which(donors.occupancy$zipped == "yes")] = "yes"
        donors.occupancy$bound.id[which(donors.occupancy$zipped == "yes") ] = "homology"
        donors.occupancy$donor.id[which(donors.occupancy$zipped == "yes") ] = current.donor
      }
      
      ##########################################################################
      ################################# Zipping ################################

      # When the twoh microhomology state is enable, the zipping occurs until all rad54 are zipped;
      if(length(unzipped.rad54 > 0) && current.donor != "" &&
         donors.list$invasion[which(donors.list$id == current.donor)] != "wrong"){

        if (donors.list$invasion[which(donors.list$id == current.donor)] =="no"){
          donors.list$invasion[which(donors.list$id == current.donor)] ="yes"
        }

        for (pos in unzipped.rad54){

          # Check if the sequence to zip is big enough ;
          #   We decided >= 16 (2*8 nts) arbitrary (could be more or less)
          if(donors.occupancy$zipped[pos] != "yes" & check.before.zipping(pos, donor = current.donor) >= 16){
            new.zip = zipping(pos, zipped.fragments.list, donor= current.donor, limit = 10)

            if(length(new.zip) > 1){
              #i.e new.zip is a vector,
              #i.e the zipping successes

              unzipped.rad54 = unzipped.rad54[which(unzipped.rad54 != pos)] #remove the current overlapped rad54 from the list
              zipped.fragments.list = rbind(zipped.fragments.list, new.zip) # add the zipped fragment to list of all the zip
              names(zipped.fragments.list) = c("start", "end", "sequences")
              current.zip.start <- as.integer(new.zip[1])
              current.zip.end <- as.integer(new.zip[2])
              donors.occupancy$zipped[current.zip.start : current.zip.end] = "yes" #set the state of zipped nucleotides as "yes

            }else if (new.zip == -1){
              #i.e zipping failed because the current donor as too much differences with the invading strand
              # Therefore we know that the current donor is not good enough to lead to homologous recombination,
              # We have to search for another potential donor, and remove the current donor from the list;

              #Dissociate all rad51 bound to this wrong donor
              remove.rad51 = which(occupied.rad51$donor.invasions == current.donor)
              occupied.rad51$genome.bins = occupied.rad51$genome.bins[-remove.rad51]
              occupied.rad51$donor.invasions = occupied.rad51$donor.invasions[-remove.rad51]
              occupied.rad51$lys2.microhomology = occupied.rad51$lys2.microhomology[-remove.rad51]

              donors.occupancy$bound[which(donors.occupancy$donor.id==current.donor)] = "no"
              donors.occupancy$bound.id[which(donors.occupancy$donor.id==current.donor)] = "unbound"
              donors.occupancy$zipped = "no"

              donors.occupancy$donor.id[which(donors.occupancy$donor.id==current.donor)] = "unknown"
              zipped.fragments.list <- as.data.frame(matrix(0,0,3))
              names(zipped.fragments.list ) = c("start", "end", "sequences")
              unzipped.rad54 = pos.rad54
              donors.list$invasion[which(donors.list$id == current.donor)] = "wrong"
              donors.blacklist = c(donors.blacklist, current.donor)
              current.donor = ""

              if(length(occupied.rad51$donor.invasions) == 0){
                occupied.rad51$bound = "unbound"
                break
              }
            }
          }
        }
      }
      ##########################################################################
      ######################## KOFF2 ###########################################
      
      # Introduce at each time step, the probability of dissociation Koff2 for zipped sequences ;
      # If a macrohomology becomes un-zipped because of koff2,
      # All the processes of homologies searching and zipping have to be done again ;
      
      if(koff2.prob > 0 & dim(zipped.fragments.list)[1] != 0){
        row2remove <- c()
        for(i in 1:nrow(zipped.fragments.list)){
          preserved.zip <- sample(c(FALSE, TRUE), size =1, replace = TRUE, prob = c(koff2.prob,1-koff2.prob))
          if(!preserved.zip){
            current.zip.start <- as.integer(zipped.fragments.list[i, ]$start)
            current.zip.end <- as.integer(zipped.fragments.list[i, ]$end)
            row2remove = c(row2remove, i)
            
            donors.occupancy$zipped[current.zip.start : current.zip.end] = "no" #the sequence is unzipped
            donors.occupancy$bound[current.zip.start : current.zip.end] = "no" #the sequence becomes unbound to donor
            donors.occupancy$bound.id[current.zip.start : current.zip.end] = "unbound" # the sequence is considered as heterologous again
            donors.occupancy$donor.id[current.zip.start : current.zip.end] = "unknown"
            
            unzipped.rad54 = c(unzipped.rad54, current.zip.start) #the rad54 into the sequence are no more overlapped by any microhomology
            
            remove.rad51 <- which(occupied.rad51$lys2.microhomology %in% (current.zip.start : current.zip.end))
            
            #remove binding sites from the donor
            occupied.rad51$genome.bins = occupied.rad51$genome.bins[-remove.rad51]
            occupied.rad51$lys2.microhomology = occupied.rad51$lys2.microhomology[-remove.rad51]
            occupied.rad51$donor.invasions = occupied.rad51$donor.invasions[-remove.rad51]
            
            if(length(occupied.rad51$donor.invasions) == 0 | length(occupied.rad51$lys2.microhomology) == 0){
              occupied.rad51$bound = "unbound"
              break
            }
          }
        }
        
        if(length(row2remove) > 0){
          zipped.fragments.list = zipped.fragments.list[-c(row2remove),]
          if(dim(zipped.fragments.list)[1] != 0){
            row.names(zipped.fragments.list) = (1:nrow(zipped.fragments.list))
          }
        }
      }
      
      ##########################################################################
      
      #first homology to LYS2
      if(length(which(donors.occupancy$bound.id == "homology" & donors.occupancy$donor.id == "LYS")) > 0 && first.lys == 0){
        first.lys = 1 
        lys.occupancy.firsts$first.bound[bigtracker] = time.step
      }
      #two hundredth homology to LYS2
      if(length(which(donors.occupancy$bound.id == "homology" & donors.occupancy$donor.id == "LYS")) >= 200){
        if(twoh.lys == 0){
          twoh.lys = 1
          lys.occupancy.firsts$twoh.bound[bigtracker] = time.step
          lys.occupancy.firsts$first.twoh.time.diff[bigtracker] = time.step - lys.occupancy.firsts$first.bound[bigtracker]
        }
      }
      
      if(current.donor == ""){
        for (candidate.donor in unique(donors.occupancy$donor.id)){
          if(candidate.donor != "unknown" && length(which(donors.occupancy$donor.id == candidate.donor)) >= 200 && 
             donors.list$invasion[which(donors.list$id == candidate.donor)] == "no"){
            
            current.donor = candidate.donor
          }
        }
      }
      
      prob.detection.donors = rep(0, additional.donors+1)
      for (i in 1:length(donors.list$id)){
        prob.detection.donors[i] = length(which(donors.occupancy$donor.id == donors.list$id[i]))/500
        if (prob.detection.donors[i] >= 1){
          prob.detection.donors[i] = 1
        }
      }
      
      prob.detection.homo = length(which(donors.occupancy$bound.id=="homology"))/500
      if (prob.detection.homo >= 1){prob.detection.homo = 1}
      prob.detection.zip = length(which(donors.occupancy$zipped=="yes"))/500
      if (prob.detection.zip >= 1){prob.detection.zip = 1}
      
      
      
      if(saver < 3 ){
        # tabulate occupancies vs. time step and length
        binding.ts$total.bound[binding.ts$time.step == time.step & 
                                 binding.ts$length == ly.type] = length(which(donors.occupancy$bound == "yes"))
        
        binding.ts$heterologies[binding.ts$time.step == time.step & 
                                  binding.ts$length == ly.type] = length(which(donors.occupancy$bound.id == "heterology"))
        
        binding.ts$homologies[binding.ts$time.step == time.step & 
                                binding.ts$length == ly.type] = length(which(donors.occupancy$bound.id == "homology"))
        
        for (i in 1:length(donors.list$id)){
          binding.ts[binding.ts$time.step == time.step & binding.ts$length == ly.type, i+5] = 
            length(which(donors.occupancy$donor.id == donors.list$id[i]))
        }
      }
      
      
      pop.time.series$homologies[pop.time.series$time.step == time.step & 
                                   pop.time.series$length == ly.type] = 
        pop.time.series$homologies[pop.time.series$time.step == time.step & 
                                     pop.time.series$length == ly.type] + prob.detection.homo
      
      pop.time.series$zip[pop.time.series$time.step == time.step & 
                            pop.time.series$length == ly.type] = 
        pop.time.series$zip[pop.time.series$time.step == time.step & 
                              pop.time.series$length == ly.type] + prob.detection.zip
      
      
      for (i in 1:length(donors.list$id)){
        pop.time.series[pop.time.series$time.step == time.step & pop.time.series$length == ly.type, i+4] = 
          pop.time.series[pop.time.series$time.step == time.step & pop.time.series$length == ly.type, i+4] + prob.detection.donors[i]
      }
      
      if(occupied.rad51$bound != "unbound"){
        for (i in 1:length(table(occupied.rad51$genome.bins))){
          bin <- names(table(occupied.rad51$genome.bins))[i]
          count <- table(occupied.rad51$genome.bins)[[i]]
          
          chromosome.contacts[chromosome.contacts$time.step == time.step &
                                chromosome.contacts$length == ly.type, names(chromosome.contacts) == bin] =
            chromosome.contacts[chromosome.contacts$time.step == time.step &
                                  chromosome.contacts$length == ly.type, names(chromosome.contacts) == bin] + count
        }
      }
      
      
    }#next time step
  }#next fragment
  
  fname = paste("timeseries", num.time.steps, kon.name, koff1.name, koff2.name, bindings.per.tethering, search.window, sep="_")
  fname = paste(fname,"_trial",as.character(trial),".txt",sep="")
  write.table(binding.ts,file=paste(dirnew_timeseries,"/", fname, sep = ""))
  
  if(saver < 3){
    binding.ts$length = factor(binding.ts$length)
    single.runs(dirnew_singles = dirnew_singles, binding.ts = binding.ts, saver = saver)
  }
  
  saver=saver+1
}#end process

write.csv(chromosome.contacts, file=paste(dirnew_data,"/chromosomes_contacts.csv",sep=""))
population.time.series(dirnew_data = dirnew_data, dirnew_plots = dirnew_pop, donors.list = donors.list, pop.time.series = pop.time.series)
stats.plots(dirnew_plots = dirnew_contacts, lys.occupancy.firsts = lys.occupancy.firsts)

