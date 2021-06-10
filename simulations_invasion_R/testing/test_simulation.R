options(bitmapType = "cairo") #fix some graphical display issues with X11 (PSMN)
rm(list=ls()) #clean global environment

###Set working directory
setwd("/home/nicolas/Documents/INSA/Stage4BiM/DSB_homologous_recombination_simulation/")

# Directory where you want to save timeseries and plots. Need the slash at the end if you want sub-directories underneath. 
rootdir = paste(getwd(), "/datas/", sep="")

###################### Import librairies #######################################

library(ggplot2)
library(stringr)
library(dplyr)
library(Rcpp)
library(text.alignment)
library(profvis)

################################################################################
############################## Import the datas ################################
sourceCpp("./simulations_invasion_R/testing/RCPP_functions.cpp")
source("./simulations_invasion_R/testing/simulation_functions_test.R") #Run the file containing the functions for the simulation 
source("./simulations_invasion_R/testing/outputs_functions_test.R") #Run the file containing the functions to create the outputs graphs


#Information concerning the 'real' or expected/experiental donor :
real.id = "LYS2"
real.bin = "chr2_470001_480001"
realdonor.genome.bindings = c()

# genome-wide microhomology counts
forward.sequences <- read.table("./LYS2/Occurences_per_8bp_motif(for+rev_donor).txt", sep="\t", header = TRUE)
forward.sequences = forward.sequences[,c("start", "sequence", "total")]
row.names(forward.sequences) = 1:nrow(forward.sequences)
microhomology.probs = forward.sequences$total / sum(forward.sequences$total)

# Name the DNA sequences of the invading strands:
LY = (tolower("ATGACTAACGAAAAGGTCTGGATAGAGAAGTTGGATAATCCAACTCTTTCAGTGTTACCACATGACTTTTTACGCCCACAACAAGAACCTTATACGAAACAAGCTACATATTCGTTACAGCTACCTCAGCTCGATGTGCCTCATGATAGTTTTTCTAACAAATACGCTGTCGCTTTGAGTGTATGGGCTGCATTGATATATAGAGTAACCGGTGACGATGATATTGTTCTTTATATTGCGAATAACAAAATCTTAAGATTCAATATTCAACCAACGTGGTCATTTAATGAGCTGTATTCTACAATTAACAATGAGTTGAACAAGCTCAATTCTATTGAGGCCAATTTTTCCTTTGACGAGCTAGCTGAAAAAATTCAAAGTTGCCAAGATCTGGAAAGGACCCCTCAGTTGTTCCGTTTGGCCTTTTTGGAAAACCAAGATTTCAAATTAGACGAGTTCAAGCATCATTTAGTGGACTTTGCTTTGAATTTGGATACCAGTAATAATGCGCATGTTTTGAACTTAATTTATAACAGCTTACTGTATTCGAATGAAAGAGTAACCATTGTTGCGGACCAATTTACTCAATATTTGACTGCTGCGCTAAGCGATCCATCCAATTGCATAACTAAAATCTCTCTGATCACCGCATCATCCAAGGATAGTTTACCTGATCCAACTAAGAACTTGGGCTGGTGCGATTTCGTGGGGTGTATTCACGACATTTTCCAGGACAATGCTGAAGCCTTCCCAGAGAGAACCTGTGTTGTGGAGACTCCAACACTAAATTCCGACAAGTCCCGTTCTTTCACTTATCGCGACATCAACCGCACTTCTAACATAGTTGCCCATTATTTGATTAAAACAGGTATCAAAAGAGGTGATGTAGTGATGATCTATTCTTCTAGGGGTGTGGATTTGATGGTATGTGTGATGGGTGTCTTGAAAGCCGGCGCAACCTTTTCAGTTATCGACCCTGCATATCCCCCAGCCAGACAAACCATTTACTTAGGTGTTGCTAAACCACGTGGGTTGATTGTTATTAGAGCTGCTGGACAATTGGATCAACTAGTAGAAGATTACATCAATGATGAATTGGAGATTGTTTCAAGAATCAATTCCATCGCTATTCAAGAAAATGGTACCATTGAAGGTGGCAAATTGGACAATGGCGAGGATGTTTTGGCTCCATATGATCACTACAAAGACACCAGAACAGGTGTTGTAGTTGGACCAGATTCCAACCCAACCCTATCTTTCACATCTGGTTCCGAAGGTATTCCTAAGGGTGTTCTTGGTAGACATTTTTCCTTGGCTTATTATTTCAATTGGATGTCCAAAAGGTTCAACTTAACAGAAAATGATAAATTCACAATGCTGAGCGGTATTGCACATGATCCAATTCAAAGAGATATGTTTACACCATTATTTTTAGGTGCCCAATTGTATGTCCCTACTCAAGATGATATTGGTACACCGGGCCGTTTAGCGGAATGGATGAGTAAGTATGGTTGCACAGTTACCCATTTAACACCTGCCATGGGTCAATTACTTACTGCCCAAGCTACTACACCATTCCCTAAGTTACATCATGCGTTCTTTGTGGGTGACATTTTAACAAAACGTGATTGTCTGAGGTTACAAACCTTGGCAGAAAATTGCCGTATTGTTAATATGTACGGTACCACTGAAACACAGCGTGCAGTTTCTTATTTCGAAGTTAAATCAAAAAATGACGATCCAAACTTTTTGAAAAAATTGAAAGATGTCATGCCTGCTGGTAAAGGTATGTTGAACGTTCAGCTACTAGTTGTTAACAGGAACGATCGTACTCAAATATGTGGTATTGGCGAAATAGGTGAGATTTATGTTCGTGCAGGTGGTTTGGCCGAAGGTTATAGAGGATTACCAGAATTGAATAAAGAAAAATTTGTGAACAACTGGTTTGTTGAAAAAGATCACTGGAATTATTTGGATAAGGATAATGGTGAACCTTGGAGACAATTCTGGTTAGGTCCAAGAGATAGATTGTACAGAACGGGTGATTTAGGTCGTTATCTACCAAACGG"))
L = (tolower("ATGACTAACGAAAAGGTCTGGATAGAGAAGTTGGATAATCCAACTCTTTCAGTGTTACCACATGACTTTTTACGCCCACAACAAGAACCTTATACGAAACAAGCTACATATTCGTTACAGCTACCTCAGCTCGATGTGCCTCATGATAGTTTTTCTAACAAATACGCTGTCGCTTTGAGTGTATGGGCTGCATTGATATATAGAGTAACCGGTGACGATGATATTGTTCTTTATATTGCGAATAACAAAATCTTAAGATTCAATATTCAACCAACGTGGTCATTTAATGAGCTGTATTCTACAATTAACAATGAGTTGAACAAGCTCAATTCTATTGAGGCCAATTTTTCCTTTGACGAGCTAGCTGAAAAAATTCAAAGTTGCCAAGATCTGGAAAGGACCCCTCAGTTGTTCCGTTTGGCCTTTTTGGAAAACCAAGATTTCAAATTAGACGAGTTCAAGCATCATTTAGTGGACTTTGCTTTGAATTTGGATACCAGTAATAATGCGCATGTTTTGAACTTAATTTATAACAGCTTACTGTATTCGAATGAAAGAGTAACCATTGTTGCGGACCAATTTACTCAATATTTGACTGCTGCGCTAAGCGATCCATCCAATTGCATAACTAAAATCTCTCTGATCACCGCATCATCCAAGGATAGTTTACCTGATCCAACTAAGAACTTGGGCTGGTGCGATTTCGTGGGGTGTATTCACGACATTTTCCAGGACAATGCTGAAGCCTTCCCAGAGAGAACCTGTGTTGTGGAGACTCCAACACTAAATTCCGACAAGTCCCGTTCTTTCACTTATCGCGACATCAACCGCACTTCTAACATAGTTGCCCATTATTTGATTAAAACAGGTATCAAAAGAGGTGATGTAGTGATGATCTATTCTTCTAGGGGTGTGGATTTGATGGTATGTGTGATGGGTGTCTTGAAAGCCGGCGCAACCTTTTCAGTTATCGACCCTGCATATCCCCCAGCCAGACAAACCATTTACTTAGGTGTTGCTAAACCACGTGGGTTGATTGTTATTA"))
L500 = (tolower("ATGACTAACGAAAAGGTCTGGATAGAGAAGTTGGATAATCCAACTCTTTCAGTGTTACCACATGACTTTTTACGCCCACAACAAGAACCTTATACGAAACAAGCTACATATTCGTTACAGCTACCTCAGCTCGATGTGCCTCATGATAGTTTTTCTAACAAATACGCTGTCGCTTTGAGTGTATGGGCTGCATTGATATATAGAGTAACCGGTGACGATGATATTGTTCTTTATATTGCGAATAACAAAATCTTAAGATTCAATATTCAACCAACGTGGTCATTTAATGAGCTGTATTCTACAATTAACAATGAGTTGAACAAGCTCAATTCTATTGAGGCCAATTTTTCCTTTGACGAGCTAGCTGAAAAAATTCAAAGTTGCCAAGATCTGGAAAGGACCCCTCAGTTGTTCCGTTTGGCCTTTTTGGAAAACCAAGATTTCAAATTAGACGAGTTCAAGCATCATTTAGTGGACTTTGCTTTGAATTTGGATACCAG"))
invading.fragments = list(names = c("500", "1000", "2000"), sequences = c(L500, L, LY))

# genome-wide microhomology counts but with bins of 10kb
sequences.bins <- read.csv("./LYS2/LY_occurences_per_8bp_(for_rev_donor)_with_bins.csv")[-1] #doesn't take the first column containing sequences

# Import the experimental contacts of the left DSB 10kb with the genome wide :
contacts <- read.csv("./LYS2/leftDSB_contacts_100000_110000_10kb.csv")
bins.id <- paste(as.character(contacts$chrom), "_", as.character(contacts$start_pos), "_", as.character(contacts$end_pos), sep="")
bins.size = 10000 
contacts <- cbind(contacts, bins.id)
colnames(contacts)[6] <- "frequency"
colnames(contacts)[7] <- "id"

################################################################################
####################### Parameters #############################################

num.time.steps = 600 # Length of simulation in time steps
graph.resolution = 1 #save occupancy data at every nth time step. Plots will have this resolution at the x-axis 

test.replicates = 1 # How many times to simulate, replicates
kon.group<-c(0.5) #binding probabilities for every binding try
koff1.group<-c(0.2) # dissociation probabilities for each bound particle
koff2.group<-c(0.05) #dissociation probabilities for each zipped fragments
ke1.group<-c(1e-2)
ke2.group<-c(2e-3)
m.group = c(2) #bindings allowed to occur per tethering
search.window.group = c(250) #the genomic distance of the tethering effect (per side)
rad54.group <- c(12) #proportional to the length of invading strand
rdh54.group <- c(2) #proportional to the number of rad54
misalignments.cutoff <- 5 #How many mismatches are allowed before break the zipping phase for the current donor
crosslink.density <- 500 #minimum density to get a probability of detection eguals to 1
additional.donors <- 2 # Additional donors ( without 'real' donor(s))

# Since the data needs to be outputted to files with human-readable names,we have to label the parameters with strings.
# For example 0005 is really 0.005
kon.group.names<- gsub("\\.", "", as.character(kon.group))
koff1.group.names<- gsub("\\.", "", as.character(koff1.group))
koff2.group.names<- gsub("\\.", "", as.character(koff2.group))
ke1.group.names<- gsub("\\.", "", as.character(ke1.group))
ke2.group.names<- gsub("\\.", "", as.character(ke2.group))
rad54.group.names<-gsub("\\.", "", as.character(rad54.group))
rdh54.group.names<-gsub("\\.", "", as.character(rdh54.group))

# print(kon.group.names)
# print(koff1.group.names)
# print(koff2.group.names)

################################################################################

################################################################################
################### Creation of the chromosome contact frequency matrix ########
# We have to check that the bins are the same between the 2 tables (sequences.bins and contacts) ;
# For example, for LY sequences.bins we have 2 more bins than in the contacts dataframe, so we remove them ;
# The comparison is made with the chromosome id and the start position for each bin 
chr_pos_occurences = c()
for (i in 2:ncol(sequences.bins)){
  chr_pos_occurences= c(chr_pos_occurences, 
                        paste(str_split(colnames(sequences.bins[i]), "_")[[1]][1], str_split(colnames(sequences.bins[i]), "_")[[1]][2], sep="_"))
}
chr_pos_contacts = c()
for (i in 1:length(bins.id)){
  chr_pos_contacts= c(chr_pos_contacts, 
                      paste(str_split(bins.id[i], "_")[[1]][1], str_split(bins.id[i], "_")[[1]][2], sep="_"))
}

remove = which(chr_pos_occurences %!in% chr_pos_contacts)+1
sequences.bins <- subset(sequences.bins, select=-remove)
contacts <- subset(contacts, select=-remove)
sequences.contacts.bins = mapply("*", sequences.bins, contacts$frequency)
rm(sequences.bins, contacts, chr_pos_occurences, chr_pos_contacts, remove)

################################################################################
############################ Single run simulation #############################
#profvis({
kon = 1; koff = 1; m = 1; sw = 1; koff2 = 1; rad54 = 1; rdh54 = 1; ke1 = 1; ke2 = 1; #for single Job run

kon.prob=kon.group[kon]
koff1.prob=koff1.group[koff]
koff2.prob=koff2.group[koff2]
ke1.prob=ke1.group[ke1]
ke2.prob=ke2.group[ke2]
nb.rad54=rad54.group[rad54]
nb.rdh54=rdh54.group[rdh54]

bindings.per.tethering = m.group[m]
search.window = search.window.group[sw]

kon.name=kon.group.names[kon]
koff1.name=koff1.group.names[koff]
koff2.name=koff2.group.names[koff2]
ke1.name=ke1.group.names[ke1]
ke2.name=ke2.group.names[ke2]
rad54.name=rad54.group.names[rad54]
rdh54.name=rdh54.group.names[rdh54]

# Initialize the occupied.rad51 list that contains :
#   the state of the invading fragment : bound or unbound ,
#   the nature of the invading strand : negative ,
#   the respective genome bins were heterologies and homologies are found,
#   the nature of the donor :
#     "H" if heterology ;
#     "donor#" if it's a donor (in the donors.list$id)
#   the position of rad51 in invading fragment for each bound microhomologies ;
occupied.rad51 = list(bound = "unbound",strand = "negative", genome.bins = c(), donor.invasions = "", pos.microhomology = 1)

# Generate N additional donors (other than the expected / experimental donor),
#   with a random number of mutations (between 10% and 40% snp ),
#   mutation occurs in random position of the sequence ,
#   the template sequence (sequence of the current inviding sequence).
# Change N to change the number of potential donors ;
donors.list = donors.generator(template = LY, realdonor.id = real.id, realdonor.location = real.bin, 
                               bins = bins.id, N=additional.donors)

saver = 0  #keeps track of how many individual simulations you want to save ;
# Bigtracker : tracker to know the in which replicate we are for which fragment, 
#   ex : for L1000 (2nd fragment of the list), for the 3th replicate, bigtracker = 8
bigtracker = 0 

####################" Initialize the population time series ####################

# Population time series are dataframes that allows us to see the provability to detect an interaction between 
#   the invading strand and the donor(s). This probability is proportional to the number of test replicates (inner runs) we make.

# In the case we introduce additional donors (not only expected or theoretical donor) in the simulation,
#   we need to create 2 more population time series in order to distinguish 
#   the general behavior of the simulation (with inclusion of all homologies for all donors) 
#   from the behavior of theoretical  donor associations ;

pop.time.series = as.data.frame(matrix(0.0,num.time.steps*3,4))
names(pop.time.series) = c("time.step","length", "homologies", "zip")
pop.time.series$time.step = rep(seq(1,num.time.steps,1),3)
pop.time.series$length = rep(invading.fragments$names, each = num.time.steps)

for (i in 1:(additional.donors+1)){
  col = paste("prob.detect", as.character(donors.list$id[i]), sep=".")
  pop.time.series[col] = rep(0, 3*num.time.steps)
}

################################################################################
############# Some other statistics dataframes #################################

# Saves the time step of each fragment sizes, at each replicates for 
#   the first homology with the real (expected) donor,
#   the two hundredth (also called twoh) homology with the real (expected) donor, 
#   the number of time steps between the first and the two hundredth ;
#   the first zipped fragment between real donor and invading strand ;
#   the probability of detect (based on the experimental crosslink density) half-detect (~ 250/500 nts) ;

occupancy.firsts = as.data.frame(matrix(-1, 3*test.replicates, 6))
names(occupancy.firsts) = c("length", "first.bound", "twoh.bound", "first.twoh.time.diff", "first.zip", "half.detect")
occupancy.firsts$length = rep(invading.fragments$names, times = test.replicates)

extensions.stats =as.data.frame(matrix(-1, 3*test.replicates, 4))
names(extensions.stats) = c("length", "time.step", "ke", "clipping.pos")
extensions.stats$length = rep(invading.fragments$names, times = test.replicates)

# Dataframe with the number of time each bins for each chromosome is contacted during the searching phase 
chromosome.contacts <- as.data.frame(matrix(0,num.time.steps*3, length(bins.id)+2))
colnames(chromosome.contacts) = c("time.step", "length", bins.id)
chromosome.contacts$time.step = rep(seq(1,num.time.steps,1),3)
chromosome.contacts$length = rep(invading.fragments$names, each = num.time.steps)

################################################################################
########################### Directory settings #################################

dirname=paste(num.time.steps, kon.name, koff1.name, koff2.name, ke1.name, ke2.name,
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
cat("Ke1 : ", ke1.prob, "\n")
cat("Ke2 : ", ke2.prob, "\n")
cat("Tethering per windows : ", bindings.per.tethering, "\n")
cat("Search windows : ", search.window, "\n")
cat("Number of rad54 : ", nb.rad54, "\n")
cat("Number of rdh54 : ", nb.rdh54, "\n\n")
for(i in 1:length(donors.list$id)){
  cat("Donor #",i-1,", on bin",donors.list$bins[i]," with", round(donors.list$mutations[i], 1)," % of SNPs", "\n")
}
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
    names(binding.ts) = c('time.step', "length", "bound", "heterologies", "homologies")
    binding.ts$time.step = rep(seq(1,num.time.steps, graph.resolution),3)
    binding.ts$length = rep(invading.fragments$names, each = (num.time.steps / graph.resolution))
  }
  
  # exonucleases are involved in the resection process of broken strands before starting the homologies search via rad51
  # The time for this operation is quite random, therefore we simulate it with a normal law :
  exonuclease.job <- as.integer(rnorm(n=1, mean = 25, sd = 10))
  if(exonuclease.job <=1){exonuclease.job = 1}
  
  for (fragment in 1:length(invading.fragments$names)){
    # initialize the ID and state of the current invading strand 
    bigtracker = bigtracker +1
    fragment.type = invading.fragments$names[fragment]
    invading.sequence = invading.fragments$sequences[fragment]
    
    current.donor = ""
    donors.blacklist = c()
    donors.list$invasion = rep("no", additional.donors+1)
    
    SEI.binding.tries = floor((nchar(invading.sequence)-7)/8)
    
    occupied.rad51$bound = "unbound"
    occupied.rad51$genome.bins = c()
    occupied.rad51$donor.invasions = c()
    occupied.rad51$pos.microhomology = c()
    
    # State of the invading fragment (occupancy) with donor(s)
    donors.occupancy = as.data.frame(matrix(0, nchar(invading.sequence),6))
    names(donors.occupancy) = c('bp', 'bound', "bound.id", "donor.id", "zipped", "bins")
    donors.occupancy$bp = 1:nchar(invading.sequence)
    donors.occupancy$bound = "no"
    donors.occupancy$bound.id = "unbound"
    donors.occupancy$donor.id = "unknown"
    donors.occupancy$zipped = "no"
    donors.occupancy$bins = "unknown"
    
    first = 0 #the first homology bound to real donor
    twoh = 0 # the first two hundred homologies bound to the real donor
    
    first.zip <- 0 #the first zipped fragment to the real donor
    half.detect <- 0 #when the probability detection is equal to 0.5 for zipped fragment to the real donor
    
    # We have to place randomly some rad54 and rdh54 in the invading fragment to induce the zipping ;
    # The number of rad54 depends of the length of the fragment,
    # and the number of rdh54 depends of the number of rad54 ;
    prop.rad54 <- floor((as.integer(fragment.type)/max(as.integer(invading.fragments$names)))*nb.rad54) #number of rad54 to be placed into the invading strand ;
    if(prop.rad54<=1){prop.rad54=1}
    prop.rdh54 <- floor((as.integer(fragment.type)/max(as.integer(invading.fragments$names)))*nb.rdh54) # number of rdh54 to be placed into the invading strand;
    if(prop.rdh54<=1){prop.rdh54=1}
    
    rad54.rdh54.locations <- rad54.rdh54.placement(number.rad54 = prop.rad54, number.rdh54 = prop.rdh54, invading.sequence = invading.sequence) 
    pos.rad54 <- rad54.rdh54.locations[[1]] #positions of rad54 in the invading strand;
    pos.rdh54 <- rad54.rdh54.locations[[2]] #positions of rdh54 in the invading strand;
    zipped.fragments.list <- as.data.frame(matrix(0,0,3)) #all the macrohomologies after zipping with start/end positions
    names(zipped.fragments.list ) = c("start", "end", "sequences")
    unzipped.rad54 <- pos.rad54 #positions of non-overlapped rad54
    
    #probability of detection proportional to the length of invading strand :
    crosslink.density <- 500
    
    # Loop through the time-steps
    for (time.step in 1:num.time.steps){
      if(kon.prob == 0){
        next
      }
      # Seach homologies in the binding tethering window : new.microhomologizer 
      if (occupied.rad51$bound != "unbound" & time.step > exonuclease.job){
        if (length(occupied.rad51$donor.invasions) != sum(occupied.rad51$donor.invasions == "H")){
          new.bindings = new.microhomologizer(occupied.rad51, search.window, bindings.per.tethering, kon.prob=kon.prob)
          occupied.rad51$genome.bins = c(occupied.rad51$genome.bins, new.bindings$genome.bins)
          occupied.rad51$donor.invasions = c(occupied.rad51$donor.invasions,new.bindings$donor.invasions)
          occupied.rad51$pos.microhomology = c(occupied.rad51$pos.microhomology, new.bindings$pos.microhomology)
        }
      }
      
      # Check if exonuclease opening DSB is done ;
      if(time.step > exonuclease.job){
        # Search new homologies for the free sites on the invading fragment
        new.bindings = genome.wide.sei(SEI.binding.tries)
        
        if (occupied.rad51$bound == "unbound"){
          occupied.rad51 = new.bindings
          if(length(new.bindings$pos.microhomology) > 0){
            occupied.rad51$bound = "bound"
          }
          
        }else{
          # print("bound and adding")
          occupied.rad51$genome.bins = c(occupied.rad51$genome.bins, new.bindings$genome.bins)
          occupied.rad51$donor.invasions = c(occupied.rad51$donor.invasions, new.bindings$donor.invasions)
          occupied.rad51$pos.microhomology = c(occupied.rad51$pos.microhomology, new.bindings$pos.microhomology)
        }
      }
      
      ############# Dissociate large number of heterologies #############################
      
      #Kind of zipping but for heterologies :
      # In the case where 200 or more hétérologies from the same bin are bound to the invading strand,
      # as there are heterologous, we simulate an SEI that will failed because of mismatches during zipping
      # and lead to the dissociation of these heterologies :
      
      heterologies.stats = as.data.frame(table(donors.occupancy$bins[which(donors.occupancy$bound.id == "heterology")]))
      if(nrow(heterologies.stats) > 0){
        names(heterologies.stats) = c("bins", "freq")
        heterologies.stats = heterologies.stats %>% filter(freq >= 200)
        remove.rad51 = which(occupied.rad51$genome.bins %in% as.character(heterologies.stats$bins))
        if(length(remove.rad51)>0){
          occupied.rad51$genome.bins = occupied.rad51$genome.bins[-remove.rad51]
          occupied.rad51$donor.invasions = occupied.rad51$donor.invasions[-remove.rad51]
          occupied.rad51$pos.microhomology = occupied.rad51$pos.microhomology[-remove.rad51]
        }
      }
      
      ############################################################################
      ###################### KOFF1 ###############################################
      #simulate random dissociation(s)
      num.bound = length(occupied.rad51$donor.invasions)
      # print(num.bound)
      preserved = sample(c(FALSE,TRUE), num.bound, replace = TRUE, prob = c(koff1.prob,1-koff1.prob)) #dissociate if FALSE
      preserved[which(occupied.rad51$pos.microhomology[preserved] %in% donors.occupancy$zipped=="yes")] = TRUE #koff1 can't dissociate zipped homologies
      
      occupied.rad51$genome.bins = occupied.rad51$genome.bins[preserved]
      occupied.rad51$donor.invasions = occupied.rad51$donor.invasions[preserved]
      occupied.rad51$pos.microhomology = occupied.rad51$pos.microhomology[preserved]
      
      if (sum(!preserved)==num.bound){
        occupied.rad51$bound = "unbound"
      }
      
      ############################################################################
      ############################## Occupancy ###################################
      donors.occupancy$bound = "no"
      donors.occupancy$bound.id = "unbound"
      donors.occupancy$donor.id = "unknown"
      donors.occupancy$bins = "unknown"
      
      if(occupied.rad51$bound != "unbound"){
        
        rad51.cover.index = occupied.rad51$pos.microhomology
        for(i in 1:7){rad51.cover.index = c(rad51.cover.index, occupied.rad51$pos.microhomology+i)}
        
        donors.occupancy$bound[rad51.cover.index]  = "yes"
        donors.occupancy$bins[rad51.cover.index] = c(sapply(occupied.rad51$genome.bins, function(x){rep(x, each = 8)}))
        donors.occupancy$donor.id[rad51.cover.index] = c(sapply(occupied.rad51$donor.invasions, function(x){rep(x, each = 8)}))
        donors.occupancy$bound.id[which(donors.occupancy$donor.id == "H")] = "heterology"
        donors.occupancy$bound.id[which(donors.occupancy$donor.id != "H" & donors.occupancy$donor.id != "unknown")] = "homology"
        
        donors.occupancy$bound[which(donors.occupancy$zipped == "yes")] = "yes"
        donors.occupancy$bound.id[which(donors.occupancy$zipped == "yes") ] = "homology"
        donors.occupancy$donor.id[which(donors.occupancy$zipped == "yes") ] = current.donor
        donors.occupancy$bins[which(donors.occupancy$zipped == "yes")] = donors.list$bins[which(donors.list$id == current.donor)] 
      }
      ############################################################################
      ################################# Zipping ##################################
      # When the twoh microhomology state is enable, the zipping occurs until all rad54 are zipped;
      if(length(unzipped.rad54 > 0) & current.donor != "" && donors.list$invasion[which(donors.list$id == current.donor)] != "failed"){
        
        if (donors.list$invasion[which(donors.list$id == current.donor)] =="no"){donors.list$invasion[which(donors.list$id == current.donor)] ="yes"}
        
        for (pos in unzipped.rad54){
          new.zip = zipping2.0(pos, zipped.fragments.list, donor= current.donor, limit = misalignments.cutoff)
          
          if(length(new.zip) > 1){
            #i.e new.zip is a vector,
            #i.e the zipping successes
            
            unzipped.rad54 = unzipped.rad54[which(unzipped.rad54 != pos)] #remove the current overlapped rad54 from the list
            zipped.fragments.list = rbind(zipped.fragments.list, new.zip) # add the zipped fragment to list of all the zip
            names(zipped.fragments.list) = c("start", "end", "sequences")
            current.zip.start <- as.integer(new.zip[1])
            current.zip.end <- as.integer(new.zip[2])
            
            donors.occupancy$zipped[current.zip.start : current.zip.end] = "yes" #set the state of zipped nucleotides as "yes
            donors.occupancy$bound[which(donors.occupancy$zipped == "yes")] = "yes"
            donors.occupancy$bound.id[which(donors.occupancy$zipped == "yes") ] = "homology"
            donors.occupancy$donor.id[which(donors.occupancy$zipped == "yes") ] = current.donor
            donors.occupancy$bins[which(donors.occupancy$zipped== "yes")] = donors.list$bins[which(donors.list$id == current.donor)]
            
          }else if (new.zip == -1){
            #i.e zipping failed because the current donor as too much differences with the invading strand
            # Therefore we know that the current donor is not good enough to lead to homologous recombination,
            # We have to search for another potential donor, and remove the current donor from the list;
            
            #Dissociate all rad51 bound to this wrong donor
            remove.rad51 = which(occupied.rad51$donor.invasions == current.donor)
            occupied.rad51$genome.bins = occupied.rad51$genome.bins[-remove.rad51]
            occupied.rad51$donor.invasions = occupied.rad51$donor.invasions[-remove.rad51]
            occupied.rad51$pos.microhomology = occupied.rad51$pos.microhomology[-remove.rad51]
            
            donors.occupancy$bound[which(donors.occupancy$donor.id==current.donor)] = "no"
            donors.occupancy$bound.id[which(donors.occupancy$donor.id==current.donor)] = "unbound"
            donors.occupancy$zipped = "no"
            donors.occupancy$bins[which(donors.occupancy$donor.id==current.donor)] = "unknown"
            donors.occupancy$donor.id[which(donors.occupancy$donor.id==current.donor)] = "unknown"
            
            zipped.fragments.list <- as.data.frame(matrix(0,0,3))
            names(zipped.fragments.list ) = c("start", "end", "sequences")
            unzipped.rad54 = pos.rad54
            donors.list$invasion[which(donors.list$id == current.donor)] = "failed"
            donors.blacklist = c(donors.blacklist, current.donor)
            current.donor = ""
            
            if(length(occupied.rad51$donor.invasions) == 0){
              occupied.rad51$bound = "unbound"
              break
            }
          }
        }
      }
      
      ############################################################################
      ######################## KOFF2 #############################################
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
            donors.occupancy$bins[current.zip.start : current.zip.end] = "unknown"
            
            unzipped.rad54 = c(unzipped.rad54, current.zip.start) #the rad54 into the sequence are no more overlapped by any microhomology
            remove.rad51 <- which(occupied.rad51$pos.microhomology %in% (current.zip.start : current.zip.end))
            
            #remove binding sites from the donor
            occupied.rad51$genome.bins = occupied.rad51$genome.bins[-remove.rad51]
            occupied.rad51$pos.microhomology = occupied.rad51$pos.microhomology[-remove.rad51]
            occupied.rad51$donor.invasions = occupied.rad51$donor.invasions[-remove.rad51]
            
            if(length(occupied.rad51$donor.invasions) == 0 | length(occupied.rad51$pos.microhomology) == 0){
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
      
      if(length(which(donors.occupancy$zipped == "yes")) == 0){
        current.donor = ""
      }
      
      ############################################################################
      
      #first homology to the real donor
      if(length(which(donors.occupancy$bound.id == "homology" & donors.occupancy$donor.id == real.id)) > 0 && first == 0){
        first = 1 
        occupancy.firsts$first.bound[bigtracker] = time.step
      }
      #two hundredth homology to the real donor
      if(length(which(donors.occupancy$bound.id == "homology" & donors.occupancy$donor.id == real.id)) >= 200){
        if(twoh == 0){
          twoh = 1
          occupancy.firsts$twoh.bound[bigtracker] = time.step
          occupancy.firsts$first.twoh.time.diff[bigtracker] = time.step - occupancy.firsts$first.bound[bigtracker]
        }
      }
      
      if(length(which(donors.occupancy$zipped == "yes" & donors.occupancy$donor.id == real.id)) >= 200 && first.zip == 0){
        first.zip <- 1
        occupancy.firsts$first.zip[bigtracker] = time.step
      }
      
      if(current.donor == ""){
        for (candidate.donor in unique(donors.occupancy$donor.id)){
          if(candidate.donor != "unknown" && candidate.donor != "H" && length(which(donors.occupancy$donor.id == candidate.donor)) >= 200 && candidate.donor %!in% donors.blacklist){
            
            current.donor = candidate.donor
          }
        }
      }
      
      ############################################################################
      ####################### Prob detection #####################################
      
      prob.detection.donors = rep(0, additional.donors+1)
      for (i in 1:length(donors.list$id)){
        prob.detection.donors[i] = length(which(donors.occupancy$donor.id == donors.list$id[i]))/ crosslink.density
        if (prob.detection.donors[i] >= 1){
          prob.detection.donors[i] = 1
        }
      }
      
      prob.detection.homo = length(which(donors.occupancy$bound.id=="homology"))/ crosslink.density
      if (prob.detection.homo >= 1){prob.detection.homo = 1}
      prob.detection.zip = length(which(donors.occupancy$zipped=="yes"))/ crosslink.density
      if (prob.detection.zip >= 1){prob.detection.zip = 1}
      if(current.donor == real.id & prob.detection.zip >= 0.5){
        half.detect <- 1
        occupancy.firsts$half.detect[bigtracker] = time.step
      }
      
      ############################################################################
      ######################## KE1 & KE2 #########################################
      if(length(unzipped.rad54)<prop.rad54){
        yy = runif(1)
        ## KE1 :
        #KE1 is effective only if the last rad54 is overlapped and zipped,
        # and if the total number of zipped nts is larger than 20% of the sequence length ;
        if (tail(pos.rad54,1) %!in% unzipped.rad54 & length(which(donors.occupancy$zipped=="yes"))>0.2*nchar(invading.sequence)){
          if(yy < ke1.prob){
            extensions.stats$time.step[bigtracker] = time.step
            extensions.stats$ke[bigtracker] = 1
            break
          }
          
        }else{
          ## KE2 :
          #zipping.window : distance between the first and last zipped nucleotids ; 
          #KE2 is effective only if the portion of zipped nts into the zipping window is larger than 20% of the sequence length;
          zipping.window <- max(as.integer(zipped.fragments.list$end)) - min(as.integer(zipped.fragments.list$start))+1
          if (zipping.window - sum(nchar(zipped.fragments.list$sequences)) > nchar(invading.sequence)*0.2 ){
            if(yy < ke2.prob){
              extensions.stats$time.step[bigtracker] = time.step
              extensions.stats$ke[bigtracker] = 2
              extensions.stats$clipping.pos[bigtracker] = max(as.integer(zipped.fragments.list$end))
              break
            }
          }
        }
      }
      
      ############################################################################
      ########################## Single runs #####################################
      
      if(saver < 3 ){
        # tabulate occupancies vs. time step and length
        binding.ts$bound[binding.ts$time.step == time.step & binding.ts$length == fragment.type] = length(which(donors.occupancy$bound == "yes"))
        binding.ts$heterologies[binding.ts$time.step == time.step & binding.ts$length == fragment.type] = length(which(donors.occupancy$bound.id == "heterology"))
        binding.ts$homologies[binding.ts$time.step == time.step & binding.ts$length == fragment.type] = length(which(donors.occupancy$bound.id == "homology"))
      }
      
      ############################################################################
      ################ Population time series and chromosomes bins ###############
      
      pop.time.series$homologies[pop.time.series$time.step == time.step & pop.time.series$length == fragment.type] = 
        pop.time.series$homologies[pop.time.series$time.step == time.step & pop.time.series$length == fragment.type] + prob.detection.homo
      
      pop.time.series$zip[pop.time.series$time.step == time.step &  pop.time.series$length == fragment.type] = 
        pop.time.series$zip[pop.time.series$time.step == time.step & pop.time.series$length == fragment.type] + prob.detection.zip
      
      ############################################################################
      
      for (i in 1:length(donors.list$id)){
        pop.time.series[pop.time.series$time.step == time.step & pop.time.series$length == fragment.type, i+4] = 
          pop.time.series[pop.time.series$time.step == time.step & pop.time.series$length == fragment.type, i+4] + prob.detection.donors[i]
      }
      
      if(occupied.rad51$bound != "unbound"){
        occupied.bins = as.data.frame(table(occupied.rad51$genome.bins))
        names(occupied.bins) = c("bins", "freq")

        chromosome.contacts[chromosome.contacts$time.step == time.step & chromosome.contacts$length == fragment.type, as.character(occupied.bins$bins)] =
          chromosome.contacts[chromosome.contacts$time.step == time.step & chromosome.contacts$length == fragment.type, as.character(occupied.bins$bins)] + occupied.bins$freq
      }
      
      #print(c(fragment.type, time.step, trial))
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
#})

write.csv(chromosome.contacts, file=paste(dirnew_data,"/chromosomes_contacts.csv",sep=""))
population.time.series(dirnew_data = dirnew_data, dirnew_plots = dirnew_pop, donors.list = donors.list, pop.time.series = pop.time.series)
stats.plots(dirnew_plots = dirnew_contacts, occupancy.firsts = occupancy.firsts)
