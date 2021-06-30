#! /usr/bin/env Rscript

options(bitmapType = "cairo") #fix some graphical display issues with X11 (PSMN)
rm(list=ls()) #clean global environment

###Set working directory
#setwd("/home/nicolas/Documents/INSA/Stage4BiM/DSB_homologous_recombination_simulation/")
setwd("/mnt/5EA60736A6070E69/Documents/INSA/Stage4BiM/DSB_homologous_recombination_simulation/")

# Directory where you want to save timeseries and plots. Need the slash at the end if you want sub-directories underneath. 
rootdir = paste(getwd(), "/datas/", sep="")

###################### Import librairies #######################################

library(ggplot2)
library(stringr)
library(dplyr)
library(Rcpp)
library(text.alignment)

################################################################################
############################## Import the datas ################################

#Information concerning the 'real' or expected/experiental donor :
real.id = "LYS"
real.bin = "chr2_470001_480001"
LY = (tolower("ATGACTAACGAAAAGGTCTGGATAGAGAAGTTGGATAATCCAACTCTTTCAGTGTTACCACATGACTTTTTACGCCCACAACAAGAACCTTATACGAAACAAGCTACATATTCGTTACAGCTACCTCAGCTCGATGTGCCTCATGATAGTTTTTCTAACAAATACGCTGTCGCTTTGAGTGTATGGGCTGCATTGATATATAGAGTAACCGGTGACGATGATATTGTTCTTTATATTGCGAATAACAAAATCTTAAGATTCAATATTCAACCAACGTGGTCATTTAATGAGCTGTATTCTACAATTAACAATGAGTTGAACAAGCTCAATTCTATTGAGGCCAATTTTTCCTTTGACGAGCTAGCTGAAAAAATTCAAAGTTGCCAAGATCTGGAAAGGACCCCTCAGTTGTTCCGTTTGGCCTTTTTGGAAAACCAAGATTTCAAATTAGACGAGTTCAAGCATCATTTAGTGGACTTTGCTTTGAATTTGGATACCAGTAATAATGCGCATGTTTTGAACTTAATTTATAACAGCTTACTGTATTCGAATGAAAGAGTAACCATTGTTGCGGACCAATTTACTCAATATTTGACTGCTGCGCTAAGCGATCCATCCAATTGCATAACTAAAATCTCTCTGATCACCGCATCATCCAAGGATAGTTTACCTGATCCAACTAAGAACTTGGGCTGGTGCGATTTCGTGGGGTGTATTCACGACATTTTCCAGGACAATGCTGAAGCCTTCCCAGAGAGAACCTGTGTTGTGGAGACTCCAACACTAAATTCCGACAAGTCCCGTTCTTTCACTTATCGCGACATCAACCGCACTTCTAACATAGTTGCCCATTATTTGATTAAAACAGGTATCAAAAGAGGTGATGTAGTGATGATCTATTCTTCTAGGGGTGTGGATTTGATGGTATGTGTGATGGGTGTCTTGAAAGCCGGCGCAACCTTTTCAGTTATCGACCCTGCATATCCCCCAGCCAGACAAACCATTTACTTAGGTGTTGCTAAACCACGTGGGTTGATTGTTATTAGAGCTGCTGGACAATTGGATCAACTAGTAGAAGATTACATCAATGATGAATTGGAGATTGTTTCAAGAATCAATTCCATCGCTATTCAAGAAAATGGTACCATTGAAGGTGGCAAATTGGACAATGGCGAGGATGTTTTGGCTCCATATGATCACTACAAAGACACCAGAACAGGTGTTGTAGTTGGACCAGATTCCAACCCAACCCTATCTTTCACATCTGGTTCCGAAGGTATTCCTAAGGGTGTTCTTGGTAGACATTTTTCCTTGGCTTATTATTTCAATTGGATGTCCAAAAGGTTCAACTTAACAGAAAATGATAAATTCACAATGCTGAGCGGTATTGCACATGATCCAATTCAAAGAGATATGTTTACACCATTATTTTTAGGTGCCCAATTGTATGTCCCTACTCAAGATGATATTGGTACACCGGGCCGTTTAGCGGAATGGATGAGTAAGTATGGTTGCACAGTTACCCATTTAACACCTGCCATGGGTCAATTACTTACTGCCCAAGCTACTACACCATTCCCTAAGTTACATCATGCGTTCTTTGTGGGTGACATTTTAACAAAACGTGATTGTCTGAGGTTACAAACCTTGGCAGAAAATTGCCGTATTGTTAATATGTACGGTACCACTGAAACACAGCGTGCAGTTTCTTATTTCGAAGTTAAATCAAAAAATGACGATCCAAACTTTTTGAAAAAATTGAAAGATGTCATGCCTGCTGGTAAAGGTATGTTGAACGTTCAGCTACTAGTTGTTAACAGGAACGATCGTACTCAAATATGTGGTATTGGCGAAATAGGTGAGATTTATGTTCGTGCAGGTGGTTTGGCCGAAGGTTATAGAGGATTACCAGAATTGAATAAAGAAAAATTTGTGAACAACTGGTTTGTTGAAAAAGATCACTGGAATTATTTGGATAAGGATAATGGTGAACCTTGGAGACAATTCTGGTTAGGTCCAAGAGATAGATTGTACAGAACGGGTGATTTAGGTCGTTATCTACCAAACGG"))

invading.fragments <- list(names = c(), sequences = c(), type= c(), probs = c(), bins.probs = c())

contacts <- read.csv("./LYS2/leftDSB_contacts_100000_110000_10kb.csv")
bins.id <- paste(as.character(contacts$chrom), "_", as.character(contacts$start_pos), "_", as.character(contacts$end_pos), sep="")
contacts <- cbind(contacts, bins.id)
colnames(contacts)[6] <- "frequency"
colnames(contacts)[7] <- "id"

split.dirs<-list.dirs(path = "./LY_split_homologies", full.names = FALSE, recursive=FALSE)

for (pp in split.dirs){
  split.files <- list.files(paste("./LY_split_homologies/", pp, sep=""), full.names = TRUE, recursive = FALSE)
  sequence <- readLines(paste("./LY_split_homologies/", pp, "/", pp, sep=""))
  len <- nchar(sequence)
  forward.sequences <- read.table(paste(split.files[grep(pattern=".txt", split.files)], sep = ""), sep="", header = TRUE)
  probs = forward.sequences$total / sum(forward.sequences$total)
  
  sequences.bins <- read.csv(paste(split.files[grep(pattern=".csv", split.files)], sep = ""))[-1]
  
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
  
  remove = which(!(chr_pos_occurences %in% chr_pos_contacts))+1
  sequences.bins <- subset(sequences.bins, select=-remove)
  contacts <- subset(contacts, select=-remove)
  sequences.contacts.bins = mapply("*", sequences.bins, contacts$frequency)
  
  invading.fragments$names = c(invading.fragments$name, as.character(len))
  invading.fragments$type = c(invading.fragments$type, pp)
  invading.fragments$sequences = c(invading.fragments$sequences, sequence)
  invading.fragments$probs = c(invading.fragments$probs, list(probs))
  invading.fragments$bins.probs = c(invading.fragments$bins.probs, list(sequences.contacts.bins))
  
  rm(sequences.bins, chr_pos_occurences, chr_pos_contacts, remove, sequences.contacts.bins)
}

################################################################################
####################### Parameters #############################################

num.time.steps = 600 # Length of simulation in time steps
graph.resolution = 1 #save occupancy data at every nth time step. Plots will have this resolution at the x-axis 

test.replicates = 20 # How many times to simulate, replicates
kon.group<-c(0.6) #binding probabilities for every binding try
koff1.group<-c(0.2) # dissociation probabilities for each bound particle
koff2.group<-c(0.02) #dissociation probabilities for each zipped fragments
ke1.group<-c(1e-2)
ke2.group<-c(1e-3)
m.group = c(5) #bindings allowed to occur per tethering
search.window.group = c(500) #the genomic distance of the tethering effect (per side)
rad54.group <- c(15) #proportional to the length of invading strand
rdh54.group <- c(4) #proportional to the number of rad54
misalignments.cutoff <- 5 #How many mismatches are allowed before break the zipping phase for the current donor
crosslink.density <- 500 #minimum density to get a probability of detection equals to 1
additional.donors <- 0 # Additional donors ( without 'real' donor(s))


# Since the data needs to be outputted to files with human-readable names,we have to label the parameters with strings.
# For example 0005 is really 0.005
kon.group.names<- gsub("\\.", "", as.character(kon.group))
koff1.group.names<- gsub("\\.", "", as.character(koff1.group))
koff2.group.names<- gsub("\\.", "", as.character(koff2.group))
ke1.group.names<- gsub("\\.", "", as.character(ke1.group))
ke2.group.names<- gsub("\\.", "", as.character(ke2.group))
rad54.group.names<-gsub("\\.", "", as.character(rad54.group))
rdh54.group.names<-gsub("\\.", "", as.character(rdh54.group))

#########################################################################################################
#################################### FUNCTIONS ##########################################################

#negate operator for %in% :
'%!in%' <- function(x,y)!('%in%'(x,y)) 

#########################################################################################################
#########################################################################################################

cppFunction(
  "IntegerVector find_occupancies_remastered(const IntegerVector occupiedRAD51, const IntegerVector additional_removals, int lower_window = 1,  int upper_window = 0, int n = 2069){
  IntegerVector indices ;
  IntegerVector remove (n);
  int x1;
  int x2;
  int i;
  int j;
  
  for (i = 0; i < occupiedRAD51.size(); i++){
    x1 = occupiedRAD51[i]-7;
    x2 = occupiedRAD51[i]+7;
    if(x1 < 1){x1 = 1;}
    for(j = x1; j<= x2; j++){
      remove[j] = j;
    }
  }
  
  if(additional_removals[0]!=0){
    for (i = 0; i < additional_removals.size(); i++){
      x1 = additional_removals[i]-7;
      x2 = additional_removals[i]+7;
      if(x1 < 1){x1 = 1;}
      for(j = x1; j<= x2; j++){
        remove[j] = j;
      }
    }
  }
  
  if(upper_window == 0 || upper_window > n-7){upper_window = n-7;}
  for(int w = lower_window; w <= upper_window; w++){
    if(remove[w] == 0){
      indices.push_back(w);
    }
  }
  return(indices);
}")

#########################################################################################################
#########################################################################################################

genome.wide.sei = function(initial.binding.tries){
  
  # Choose region of LY weighted by available microhomologies (MHs);
  # Dont let those already occupied be chosen (by the use of the find.occupancies() function) ;
  
  if(occupied.rad51$bound=="unbound"){
    open.sites = 1:(nchar(invading.sequence) -7)
  }else{
    open.sites = find_occupancies_remastered(occupiedRAD51 = occupied.rad51$pos.microhomology, additional_removals = c(0), n = nchar(invading.sequence))
  }
  
  # open.sites = find.occupancies() #indexes of unoccupied sites
  
  if (length(open.sites)== 0){ #if all the sites are occupied
    return(list(bound = occupied.rad51$bound, strand = "negative", genome.bins = c(), donor.invasions = c(), pos.microhomology = c()))
  }
  
  #matches : vector of possible bounding sites for MHs
  matches = c()
  bins = c()
  
  for (i in 1:initial.binding.tries){
    if (length(open.sites) == 1){ 
      # If there is only one available binding site :
      matches[i] = open.sites
      contact.freq = sequences.contacts.bins[matches[i], ] #check for all bins the frequence of contact for the current matche ;
      possible.bins = bins.id[which(contact.freq > 0)] #Keep only the bins where the matche has a non-null contact frequence ;
      bins = c(bins, sample(x=possible.bins, size=1, prob = contact.freq[contact.freq > 0])) #select a bin among the possible.bins by using frequence of contact as probability ;
      
    }else{ 
      # Draw a site among those available which will be matched with a MH according to the respective weighted probabilities :
      matches[i] = sample(x=open.sites, size=1, prob = microhomology.probs[open.sites])
      contact.freq = sequences.contacts.bins[matches[i], ] #check for all bins the frequence of contact for the current matche ;
      possible.bins = bins.id[which(contact.freq > 0)] #Keep only the bins where the matche has a non-null contact frequence ;
      bins = c(bins, sample(x=possible.bins, size=1, prob = contact.freq[contact.freq > 0])) #select a bin among the possible.bins by using frequence of contact as probability ;
    }
    
    # Where there is a match with a MH, we consider that the site concerned is no longer available ;
    # We remove it from the index open.sites;
    # To be sure that the next matches will not overlap with the previous one, 
    #  we have to remove from the index the 7 positions upstream and downstream of the match :
    
    #open.sites = open.sites[sequences.contacts.bins[open.sites,current.bin] > 0]
    open.sites = open.sites[-which(open.sites %in% (matches[i]-7):(matches[i] + 7))]
    
    if (length(open.sites) < 1){ 
      break
    }
  }
  
  # Choose if/which binds;
  # successes : list of sample positions among the matches where the probability of association is good enough to have a bond;
  # The probability of association is called 'Kon'.
  # We keep only the matches where a bond will occur;
  successes=sample( c(TRUE, FALSE), length(matches), replace = TRUE, prob = c(kon.prob, (1-kon.prob)) )
  matches = matches[successes]
  bins = bins[successes]
  
  if (length(matches)<1){
    return(list(bound = occupied.rad51$bound,strand = "negative", genome.bins = c(), donor.invasions = c(), pos.microhomology = c()))
  }
  
  
  # Set IDs of each bound (Heterology (H) vs LYS); 
  # If LYS, set as genomic position in + strand notation
  identities = c()
  for (b in 1:length(matches)){
    if(bins[b] %in% donors.list$bins){ #if the bin where the current microhomology is bound countains a donor (in the donor.list)
      donor = donors.list$id[which(donors.list$bins == bins[b])] #list of the donor(s) contained into the current.bin
      # If we found a microhomology in a bin that contains a potential donor ,
      #   we consider a probability of 1/2 for this microhomology to homologous, and thus 1/2 to be heterologous in the other case.
      yy = runif(1)
      if(yy <= 1/2){ #probability to be a donor 
        if(length(donor)>1){
          donor = sample(donor, size = 1) #rare case where we have more than one donor into a bin
          identities = c(identities, donor) #homology, bound to a potential donor
        }else{
          identities = c(identities, donor) #homology, bound to a potential donor
        }
      }else{
        identities = c(identities, "H") #heterology
      }
    }else{
      identities = c(identities, "H") #heterology
    }
  }
  
  if (occupied.rad51$bound != "unbound"){
    remove = which(matches %in% occupied.rad51$pos.microhomology)
    if (length(remove)>0){
      identities = identities[-remove]
      matches = matches[-remove]
      bins = bins[-remove]
    }
  } 
  
  
  # LYS alignment or misalignments :
  return(list(bound=occupied.rad51$bound, strand = "negative", genome.bins = bins, donor.invasions = identities, pos.microhomology = matches))
}  
#########################################################################################################
#########################################################################################################

new.microhomologizer = function(occupied.rad51, window, bindings.per.tethering, kon.prob){
  
  genome.bindings = which(occupied.rad51$donor.invasions %!in% donors.blacklist) #list of all bound MHs with the whole genome (genome.wide.sei step)
  new.bindings = list(bound=occupied.rad51$bound, strand = "negative", genome.bins = c(), donor.invasions = c(), pos.microhomology = c()) #initialize the list that will be return (copy of an empty occupied.rad51 list)
  
  # Check for unbound sites.
  #   If not, don't go futher and just return an empty list.
  #   All the rad51 are matched somewhere in the genome. 
  
  if (length(find_occupancies_remastered(occupiedRAD51 = occupied.rad51$pos.microhomology, additional_removals = c(0), n = nchar(invading.sequence))) == 0){
    return(new.bindings)
  }
  
  # bindings : list of sites occupied by another MHs into the search window around the current micros locus ;
  bindings = c()
  bins = c()
  identities= c()
  
  for (binding.index in genome.bindings){
    if (length(bindings) > 0){
      if (length(find_occupancies_remastered(occupiedRAD51 = occupied.rad51$pos.microhomology, additional_removals = bindings, n = nchar(invading.sequence))) == 0){
        break
      }
    }else{ 
      if (length(find_occupancies_remastered(occupiedRAD51 = occupied.rad51$pos.microhomology, additional_removals = c(0), n = nchar(invading.sequence))) == 0){
        break
      } 
    }
    
    #current.selocus : index of the MH we are currently looking around it (search window) to place another MHs;
    current.selocus = occupied.rad51$pos.microhomology[binding.index]
    current.bin = occupied.rad51$genome.bins[binding.index]
    current.id = occupied.rad51$donor.invasions[binding.index]
    
    # "additionals" will be past as argument during the next calls of the find.occupancies() function.
    if (length(bindings) <=0){
      additionals = "none"
    }else{
      additionals = bindings
    }
    
    # We look if we have available open sites around the current MH locus, i.e into the search windows around it;
    # open.sites = find.occupancies(lower.window = current.selocus - window,
    #                               upper.window = current.selocus + window, 
    #                               additional.removals = additionals)
    
    
    lower.window = ifelse(current.selocus - window <1, 1, current.selocus - window)
    upper.window = ifelse(current.selocus + window >nchar(invading.sequence)-7, nchar(invading.sequence)-7, current.selocus + window)
    if (additionals[1]== "none"){additionals=c(0)}
    
    open.sites = find_occupancies_remastered(occupiedRAD51 = occupied.rad51$pos.microhomology, additional_removals = additionals, 
                                             lower_window = lower.window, upper_window = upper.window, n = nchar(invading.sequence))
    
    if (length(open.sites) <= 0){next}
    current.bindings = c()
    
    # Here, for each binding Mh, we want to bind j other MHs in a window, this is the tethering phase. 
    # Of course, we can only bind new microhomologies onto free sites of the invading fragment (i.e only the open.sites can be matched).
    
    # In the case where the current MH is heterologous :
    #   We look at other MHs in the window distance that are in the same bin that the current heterology.
    
    if(current.id == "H"){
      for (j in 1:bindings.per.tethering){
        if (length(open.sites)==1){
          # if only one open.site remains free, i.e. the invading fragment is almost fully bound
          contact.freq = sequences.contacts.bins[open.sites, ] #check for all bins the frequence of contact for remaining open.site ;
          possible.bins = bins.id[which(contact.freq > 0)] #Keep only the bins where the matche has a non-null contact frequence ;
          if(current.bin %in% possible.bins){ #check if we can bind an MH that is in the same bin as the current hetetology
            yy = runif(1)
            if(yy <= kon.prob){ #kon : probability of association for each binding
              candidate = open.sites
              current.bindings = c(current.bindings,candidate)
              bins = c(bins, current.bin)
              identities = c(identities, current.id)
              open.sites = open.sites[-which(open.sites %in% (candidate-7):(candidate+7))]
              if (length(open.sites) <1){break}
            }
          }
        }else{
          #Most often case where we have more than one open.sites.
          # That the same process as above (with one remaining open.site), but we just sample the possible open.sites
          candidate = sample(open.sites, size = 1) #select one sites among the open.sites
          contact.freq = sequences.contacts.bins[candidate, ]
          possible.bins = bins.id[which(contact.freq > 0)]
          if(current.bin %in% possible.bins){
            yy = runif(1)
            if(yy <= kon.prob){ #kon
              current.bindings = c(current.bindings,candidate)
              bins = c(bins, current.bin)
              identities = c(identities, current.id)
              open.sites = open.sites[-which(open.sites %in% (candidate-7):(candidate+7))]
              if (length(open.sites) <1){break}
            }
          }
        }
      }
    }else{
      # In the case where the current MH is homologous, i.e. bound to a donor (in the donors.list) :
      #   We make also a tethering but for homologous MHs
      for (j in 1:bindings.per.tethering){
        if (length(open.sites)==1){
          yy = runif(1)
          if(yy <= kon.prob){
            candidate = open.sites
            current.bindings = c(current.bindings,candidate)
            bins = c(bins, current.bin)
            identities = c(identities, current.id)
            open.sites = open.sites[-which(open.sites %in% (candidate-7):(candidate+7))]
            if (length(open.sites) <1){break}
          }
        }else{
          candidate = sample(open.sites, size = 1)
          yy = runif(1)
          if(yy <= kon.prob){
            current.bindings = c(current.bindings,candidate)
            bins = c(bins, current.bin)
            identities = c(identities, current.id)
            open.sites = open.sites[-which(open.sites %in% (candidate-7):(candidate+7))]
            if (length(open.sites) <1){break}
          }
        }
      }
    }
    bindings = c(bindings, current.bindings)
  }
  
  new.bindings$genome.bins = c(new.bindings$genome.bins, bins)
  new.bindings$pos.microhomology = c(new.bindings$pos.microhomology, bindings)
  new.bindings$donor.invasions  = c(new.bindings$donor.invasions, identities)
  
  if (occupied.rad51$bound != "unbound"){
    #remove MHs ids in new.bindings that have already been counted as donor in occupied.rad51 :
    remove = which(bindings %in% occupied.rad51$pos.microhomology)
    if (length(remove) > 0){
      new.bindings$genome.bins = new.bindings$genome.bins[-remove]
      new.bindings$donor.invasions = new.bindings$donor.invasions[-remove]
      new.bindings$pos.microhomology = new.bindings$pos.microhomology[-remove]}
  }
  return(new.bindings)
}

#########################################################################################################
#########################################################################################################

donors.generator <- function(template,realdonor.id, realdonor.location, bins, N = 0){
  new.donors.list<-list(sequence = c(template), bins = c(realdonor.location), id = c(realdonor.id), 
                        invasion = c("no", rep("no", times = N)), mutations = c(0))
  
  bases <- c("a", "t", "g", "c")
  if(N >= 1){
    for (n in 1:N){
      new.donor <- template
      new.donor.id = paste("donor", as.character(n), sep="")
      lower.limit <- floor(0.05*nchar(template))
      upper.limit <- floor(0.4*nchar(template))
      nb.snp <-sample(lower.limit:upper.limit, size = 1)
      snp.location <- sample(1:nchar(template), size = nb.snp, replace = FALSE)
      
      for (i in snp.location){
        snp <- sample(bases[-which(bases == substr(template, i, i))], size = 1)
        substr(new.donor, i, i) <- snp
      }
      new.donors.list$sequence = c(new.donors.list$sequence, new.donor)
      new.donors.list$bins = c(new.donors.list$bins, sample(bins, size = 1))
      new.donors.list$id = c(new.donors.list$id, new.donor.id)
      new.donors.list$mutations = c(new.donors.list$mutations, (nb.snp/nchar(template))*100 )
    }
  }
  return(new.donors.list)
}

#########################################################################################################
#########################################################################################################

rev.comp<-function(x,rev=TRUE){ 
  #Compute the reverse complement of a sequence ;
  
  x<-toupper(x)
  y<-rep("N",nchar(x))
  xx<-unlist(strsplit(x,NULL))
  for (bbb in 1:nchar(x)){
    if(xx[bbb]=="A") y[bbb]<-"T"        
    if(xx[bbb]=="C") y[bbb]<-"G"        
    if(xx[bbb]=="G") y[bbb]<-"C"        
    if(xx[bbb]=="T") y[bbb]<-"A"
  }
  if(rev==FALSE) {
    for(ccc in (1:nchar(x))){
      if(ccc==1) yy<-y[ccc] else yy<-paste(yy,y[ccc],sep="")
    }
  }
  if(rev==T){
    zz<-rep(NA,nchar(x))
    for(ccc in (1:nchar(x))){
      zz[ccc]<-y[nchar(x)+1-ccc]
      if(ccc==1) yy<-zz[ccc] else yy<-paste(yy,zz[ccc],sep="")
    }
  }
  return(tolower(yy))
}

#########################################################################################################
#########################################################################################################

rad54.rdh54.placement <- function(number.rad54, number.rdh54, invading.sequence){
  
  #the last rad54 is the most important, it will start the extension (recombination) step once it will be zipped ;
  pos.last.rad54 <- nchar(invading.sequence) - as.integer(runif(1, min = 8, max=24)) #min/max chosen arbitrary
  location.rad54 <- c(pos.last.rad54)
  number.rad54 = number.rad54 -1
  location.rdh54 <- c()
  
  while (number.rad54 > 0){ #place the required rad54 randomly (according uniform distro) over the invading fragment ;
    new.location <- 0
    while (new.location == 0 | new.location %in% location.rad54){
      new.location <- floor(runif(1, min = 1, max=pos.last.rad54))
    }
    location.rad54 = c(location.rad54, new.location)
    number.rad54 = number.rad54 - 1
  }
  
  while (number.rdh54 > 0){
    new.location <-  0
    while (new.location == 0 | new.location %in% location.rad54 | new.location %in% location.rdh54){
      new.location <- floor(runif(1, min = 1, max=pos.last.rad54))
    }
    location.rdh54 = c(location.rdh54, new.location) #location.rad54 becomes location.rdh54
    number.rdh54 = number.rdh54 - 1
  }
  
  return(list(sort(location.rad54), sort(location.rdh54)))
}

#########################################################################################################
#########################################################################################################
zipping2.0 <- function(rad54, zipping.list, donor, limit){
  
  #return code :
  # new.zip : vector containing position stat and stop for the zipped fragment and its nucleotids sequence ;
  # 0 : if the zipping can't occur ;
  # -1 : the alignment is too bad, means that the current donor is not good enough and we have to change it (return to homology search step) ;
  
  #Check if the current rad54 is overlapped by an homologous microhomology ;
  if(donors.occupancy$bound.id[rad54] != "homology"){
    return(0)
  }
  
  #initialize the portion of nucleotids to zip ;
  # start : current rad54 ;
  # stop : nearest rad54 or rdh54 from the start position ;
  # fragment.to.zip : sequence of nts between start and stop ;
  
  start <- rad54
  if(start == max(pos.rad54)){
    stop <- nchar(invading.sequence)
  }else{
    stop = min(pos.rad54[pos.rad54>start], pos.rdh54[pos.rdh54>start])-1
  }
  
  fragment.to.zip <- substr(invading.sequence, start = start, stop = stop)
  #Minimum requiered length to be zipped : 16 nts ;
  if(nchar(fragment.to.zip) < 16){
    return(0)
  }
  
  donor.seq = donors.list$sequence[which(donors.list$id == donor)] #sequence of the current  donor
  
  #Here we use the algorithm of Smith Waterman (sw) to make a local alignment of 2 strings (a & b) with different length;
  # This algorithm give us a score which takes account the matches, the mismatches or the gaps between the 2 strings ;
  # This score is called "similary" and is normalized by the length of the string we want to align (b) ;
  sw <- as.data.frame(smith_waterman(a=donor.seq, b=fragment.to.zip, edit_mark = "*"))
  
  if(sw$similarity >= 5/8){
    #If the similarity score is good enough, we check the number of consecutive misalignments ;
    # We decide arbitrary that if there are more than 5 CONSECUTIVE misalignments, the fragment can't be zipped because of its instability ;
    miss <- strsplit(sw$b_aligned, split = "")[[1]]
    consecutive.miss <- ifelse(length(which(miss == "*"))>0, max(rle(miss)$length[which(rle(miss)$value=="*")]), 0)
    if(consecutive.miss <= limit){
      return(c(start, stop, fragment.to.zip))
      
    }else{
      return(-1)
    }
  }else{
    return(-1)
  }
}

#########################################################################################################
#########################################################################################################


#########################################################################################################
################################Outputs functions #######################################################

single.runs <-function(dirnew_singles, binding.ts, saver, w=14, h=8){
  
  outname=paste(dirnew_singles,"/Total_Occupancy_",saver,".png",sep="")
  occ_plot<-
    ggplot(data = binding.ts) + geom_step(aes(x = time.step, y = bound, color = length)) +
    labs(x = "time step", y = "Total Occupancy (bp)") + theme(text = element_text(size = 16))+
    scale_y_continuous(limits = c(0, max(binding.ts$bound)+1))
  ggsave(outname,plot=occ_plot, width = w, height = h)
  
  outname=paste(dirnew_singles,"/Occupancy_Heterologies_",saver,".png",sep="")
  het_plot<-
    ggplot(data = binding.ts) + geom_step(aes(x = time.step, y = heterologies, color = length)) +
    labs(x = "time step", y = "Occupancy at Heterologies (bp)") + theme(text = element_text(size = 16))+
    scale_y_continuous(limits = c(0, max(binding.ts$heterologies)+1))
  ggsave(outname,plot=het_plot, width = w, height = h)
  
  outname=paste(dirnew_singles,"/Occupancy_Homologies_",saver,".png",sep="")
  het_plot<-
    ggplot(data = binding.ts) + geom_step(aes(x = time.step, y = homologies, color = length)) +
    labs(x = "time step", y = "Occupancy at Homologies (bp)") + theme(text = element_text(size = 16))+
    scale_y_continuous(limits = c(0, max(binding.ts$homologies)+1))
  ggsave(outname,plot=het_plot, width = w, height = h)
  
}

#########################################################################################################
#########################################################################################################


population.time.series <- function(dirnew_data, dirnew_plots, donors.list, pop.time.series, w=14, h=8){
  
  # pop time series for all homologies
  df <- subset(x=pop.time.series, select=c(1, 2, 3))
  colnames(df) = c("time.step", "length", "homologies")
  df.name <- paste("pop_timeseries", "homologies", sep = "_")
  write.table(df, file=paste(dirnew_data, "/", df.name ,".txt", sep=""))
  outname=paste(dirnew_plots, "/", df.name, ".png",sep="")
  
  pop.plot<-
    ggplot(data = df) + geom_step(aes(x = time.step, y = homologies, color = length)) +
    labs(x = "time step", y = "Probability of detection for homologies") +  theme(text = element_text(size = 14))+
    scale_y_continuous(limits = c(0, max(df$homologies)+1))
  ggsave(outname, plot=pop.plot, width = w, height = h)
  
  #pop time series for all zipped homologies 
  df <- subset(x=pop.time.series, select=c(1, 2, 4))
  colnames(df) = c("time.step", "length", "zip")
  df.name <- paste("pop_timeseries", "zipped_homologies", sep = "_")
  write.table(df, file=paste(dirnew_data, "/", df.name ,".txt", sep=""))
  outname=paste(dirnew_plots, "/", df.name, ".png",sep="")
  
  pop.plot<-
    ggplot(data = df) + geom_step(aes(x = time.step, y = zip, color = length)) +
    labs(x = "time step", y = "Probability of detection for zips") + theme(text = element_text(size = 14))+
    scale_y_continuous(limits = c(0, max(df$zip)+1))
  ggsave(outname, plot=pop.plot, width = w, height = h)
  
  
  for (i in 1:length(donors.list$id)){
    df <- subset(x=pop.time.series, select=c(1, 2, i+4))
    colnames(df) = c("time.step", "length", "prob.detect")
    df.name <- paste("pop_timeseries", as.character(donors.list$id[i]),as.character(donors.list$bins[i]), sep = "_")
    write.table(df, file=paste(dirnew_data, "/", df.name ,".txt", sep=""))
    outname=paste(dirnew_plots, "/", df.name, ".png",sep="")
    
    pop.plot<-
      ggplot(data = df) + geom_step(aes(x = time.step, y = prob.detect, color = length)) +
      labs(x = "time step", y = paste("Probability of detection", as.character(donors.list$id[i]), sep = " ")) + 
      theme(text = element_text(size = 14))+
      scale_y_continuous(limits = c(0, max(df$prob.detect)+1))
    ggsave(outname, plot=pop.plot, width = w, height = h)
    
  }
}

#########################################################################################################
#########################################################################################################

stats.plots <- function(dirnew_plots, occupancy.firsts, w=10, h=8){
  
  write.table(occupancy.firsts,file=paste(dirnew_data,"/", "occupancy_first.txt", sep = ""))
  write.table(dloop.stats,file=paste(dirnew_data,"/", "dloop_invasions.txt", sep = ""))
  write.table(extensions.stats,file=paste(dirnew_data,"/", "extensions_stats.txt", sep = ""))
  
  file = paste(dirnew_plots,"/dloop_invasion_count.png",sep="")
  dloop.hist1 <- 
    ggplot(dloop.stats, aes(x = as.character(time.step), y = count))+
    geom_bar(
      aes(fill = length), stat = "identity", color = "white",
      position = position_dodge(0.6))
  ggsave(file,plot=dloop.hist1, width = w, height = h)
  
  file = paste(dirnew_plots,"/dloop_invasion_average_size.png",sep="")
  dloop.hist2 <- 
    ggplot(dloop.stats, aes(x = as.character(time.step), y = average.size))+
    geom_bar(
      aes(fill = length), stat = "identity", color = "white",
      position = position_dodge(0.6))
  ggsave(file,plot=dloop.hist2, width = w, height = h)
  
  file = paste(dirnew_plots,"/first_contact_time_hist.png",sep="")
  first.hist<-
    ggplot(occupancy.firsts[c(which(occupancy.firsts$first.bound != -1)),], 
           aes(x=first.bound, fill=length)) +
    geom_histogram(binwidth = 0.5, alpha = 0.5, position="identity")
  ggsave(file,plot=first.hist, width = w, height = h)
  
  file = paste(dirnew_plots,"/first_contact_time_boxplot.png",sep="")
  first.boxplot<-
    ggplot(occupancy.firsts[c(which(occupancy.firsts$first.bound!= -1)),], aes(x=length, y=first.bound, fill=length)) +
    geom_boxplot(outlier.colour ="red", position = position_dodge(1)) +
    stat_summary(fun = mean, geom = "point", shape = 8, size = 4)
  ggsave(file,plot=first.boxplot, width = w, height = h)
  
  file = paste(dirnew_plots,"/200_contact_time_hist.png",sep="")
  first.hist<-
    ggplot(occupancy.firsts[c(which(occupancy.firsts$twoh.bound != -1)),], 
           aes(x=twoh.bound, fill=length)) +
    geom_histogram(binwidth = 0.5, alpha = 0.5, position="identity")
  ggsave(file,plot=first.hist, width = w, height = h)
  
  file = paste(dirnew_plots,"/200_contact_time_boxplot.png",sep="")
  first.boxplot<-
    ggplot(occupancy.firsts[c(which(occupancy.firsts$first.twoh.time.diff != -1)),], 
           aes(x=length, y=twoh.bound, fill=length)) +
    geom_boxplot(outlier.colour ="red", position = position_dodge(1)) +
    stat_summary(fun = mean, geom = "point", shape = 8, size = 4)
  ggsave(file,plot=first.boxplot, width = w, height = h)
  
  file = paste(dirnew_plots,"/1st_to_200_contact_timediff_hist.png",sep="")
  first.hist<-
    ggplot(occupancy.firsts[c(which(occupancy.firsts$first.twoh.time.diff!= -1)),], 
           aes(x=first.twoh.time.diff, fill=length)) +
    geom_histogram(binwidth = 0.5, alpha = 0.5, position="identity")
  ggsave(file,plot=first.hist, width = w, height = h)
  
  file = paste(dirnew_plots,"/1st_to_200_contact_timediff_boxplot.png",sep="")
  first.boxplot<-
    ggplot(occupancy.firsts[c(which(occupancy.firsts$first.twoh.time.diff != -1)),],
           aes(x=length, y=first.twoh.time.diff, fill=length)) +
    geom_boxplot(outlier.colour ="red", position = position_dodge(1)) +
    stat_summary(fun = mean, geom = "point", shape = 8, size = 4)
  ggsave(file,plot=first.boxplot, width = w, height = h)
  
  file = paste(dirnew_plots,"/first_zip_contact_time_boxplot.png",sep="")
  first.boxplot<-
    ggplot(occupancy.firsts[c(which(occupancy.firsts$first.zip != -1)),],
           aes(x=length, y=first.zip, fill=length)) +
    geom_boxplot(outlier.colour ="red", position = position_dodge(1)) +
    stat_summary(fun = mean, geom = "point", shape = 8, size = 4)
  ggsave(file,plot=first.boxplot, width = w, height = h)
  
  file = paste(dirnew_plots,"/half_detect_boxplot.png",sep="")
  first.boxplot<-
    ggplot(occupancy.firsts[c(which(occupancy.firsts$half.detect != -1)),],
           aes(x=length, y=half.detect, fill=length)) +
    geom_boxplot(outlier.colour ="red", position = position_dodge(1)) +
    stat_summary(fun = mean, geom = "point", shape = 8, size = 4)
  ggsave(file,plot=first.boxplot, width = w, height = h)
  
  file = paste(dirnew_plots,"/start_extensions.png",sep="")
  extensions.boxplot<-
    ggplot(extensions.stats[c(which(extensions.stats$time.step!= -1)),],
           aes(x=length, y=time.step, fill=length)) +
    geom_boxplot(outlier.colour ="red", position = position_dodge(1)) +
    stat_summary(fun = mean, geom = "point", shape = 8, size = 4)
  ggsave(file,plot=extensions.boxplot, width = w, height = h)
  
  file = paste(dirnew_plots,"/ke_occurences.png",sep="")
  ke.hist<-
    ggplot(extensions.stats[c(which(extensions.stats$ke!= -1)),],
           aes(x=ke)) + geom_histogram(binwidth = 0.5, alpha = 0.5, position="identity", color="black", fill="lightblue")
  ggsave(file,plot=ke.hist, width = w, height = h)
  
}

#########################################################################################################
#########################################################################################################

######################################### Simulation Start ##############################################
#########################################################################################################


# DIRECTORY NAMING CONVENTIONS
# directory name convention : num.time.steps _ kon.prob _ koff.prob _ bindings.per.tethering _ search.window
# no spaces, just did spaces to make it clearer

# FILE NAMING CONVENTIONS
# timeseries_time.steps_kon_koff1_tethering.window_trial#.txt

# Simulation for each combination kon/koff/m/sw
# i.e len(kon) * len(koff) * len(m) * len(sw) simulations 
for(kon in 1:length(kon.group)){
  for(koff in 1:length(koff1.group)){
    for(koff2 in 1:length(koff2.group)){
      for(ke1 in 1:length(ke1.group)){
        for(ke2 in 1:length(ke2.group)){
          for(m in 1:length(m.group)){
            for(sw in 1:length(search.window.group)){
              for(rad54 in 1:length(rad54.group)){
                for (rdh54 in 1:length(rdh54.group)){
    
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
                  
                  pop.time.series = as.data.frame(matrix(0.0,num.time.steps*length(invading.fragments$names),4))
                  names(pop.time.series) = c("time.step","length", "homologies", "zip")
                  pop.time.series$time.step = rep(seq(1,num.time.steps,1),length(invading.fragments$names))
                  pop.time.series$length = rep(invading.fragments$names, each = num.time.steps)
                  
                  for (i in 1:(additional.donors+1)){
                    col = paste("prob.detect", as.character(donors.list$id[i]), sep=".")
                    pop.time.series[col] = rep(0, length(invading.fragments$names)*num.time.steps)
                  }
                  
                  ################################################################################
                  ############# Some other statistics dataframes #################################
                  
                  # Saves the time step of each fragment sizes, at each replicates for 
                  #   the first homology with the real (expected) donor,
                  #   the two hundredth (also called twoh) homology with the real (expected) donor, 
                  #   the number of time steps between the first and the two hundredth ;
                  #   the first zipped fragment between real donor and invading strand ;
                  #   the probability of detect (based on the experimental crosslink density) half-detect (~ 250/500 nts) ;
                  
                  occupancy.firsts = as.data.frame(matrix(-1, length(invading.fragments$names)*test.replicates, 6))
                  names(occupancy.firsts) = c("length", "first.bound", "twoh.bound", "first.twoh.time.diff", "first.zip", "half.detect")
                  occupancy.firsts$length = rep(invading.fragments$names, times = test.replicates)
                  
                  # Saves the time step where the extension is started, store if it is ke1 or ke1 :
                  extensions.stats =as.data.frame(matrix(-1, length(invading.fragments$names)*test.replicates, 4))
                  names(extensions.stats) = c("length", "time.step", "ke", "clipping.pos")
                  extensions.stats$length = rep(invading.fragments$names, times = test.replicates)
                  
                  # Captures the number of dloops per single end fragment and their average size at time step 200, 400, 600:
                  dloop.stats = as.data.frame(matrix(0, length(invading.fragments$names)*3, 4))
                  names(dloop.stats) = c("length", "time.step", "count", "average.size")
                  dloop.stats$length = rep(invading.fragments$names, times = 3)
                  dloop.stats$time.step = rep(c(200, 400, 600), each= length(invading.fragments$names))
                  
                  # Dataframe with the number of time each bins for each chromosome is contacted during the searching phase 
                  chromosome.contacts <- as.data.frame(matrix(0,num.time.steps*length(invading.fragments$names), length(bins.id)+2))
                  colnames(chromosome.contacts) = c("time.step", "length", bins.id)
                  chromosome.contacts$time.step = rep(seq(1,num.time.steps,1),length(invading.fragments$names))
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
                      binding.ts = as.data.frame(matrix(0, (num.time.steps/graph.resolution)*length(invading.fragments$names),5))
                      names(binding.ts) = c('time.step', "length", "bound", "heterologies", "homologies")
                      binding.ts$time.step = rep(seq(1,num.time.steps, graph.resolution), length(invading.fragments$names))
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
                      microhomology.probs = invading.fragments$probs[[fragment]]
                      sequences.contacts.bins = invading.fragments$bins.probs[[fragment]]
                      
                      
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
                      
                      rad54.rdh54.locations <- rad54.rdh54.placement(number.rad54 = nb.rad54, number.rdh54 = nb.rdh54, invading.sequence = invading.sequence) 
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
                        if(length(unzipped.rad54) > 0 & current.donor != "" && current.donor %!in% donors.blacklist){
                          
                          if (donors.list$invasion[which(donors.list$id == current.donor)] =="no"){
                            donors.list$invasion[which(donors.list$id == current.donor)] ="yes"
                          }
                          
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
                              
                              break
                            }
                          }
                        }
                        
                        if(length(occupied.rad51$donor.invasions) == 0){
                          occupied.rad51$bound = "unbound"
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
                        if(length(unzipped.rad54)<nb.rad54){
                          ## KE1 :
                          #KE1 is effective only if the last rad54 is overlapped and zipped,
                          # and if the total number of zipped nts is larger than 20% of the sequence length ;
                          if ( max(pos.rad54) %!in% unzipped.rad54 & length(which(donors.occupancy$zipped=="yes"))>0.2*nchar(invading.sequence)){
                            yy = runif(1)
                            if(yy < ke1.prob){
                              extensions.stats$time.step[bigtracker] = time.step
                              extensions.stats$ke[bigtracker] = 1
                              
                              dloop.stats$count[which(dloop.stats$time.step>time.step & dloop.stats$length == fragment.type)[1]] = 
                                dloop.stats$count[which(dloop.stats$time.step>time.step & dloop.stats$length == fragment.type)[1]] + nrow(zipped.fragments.list)
                              
                              dloop.stats$average.size[which(dloop.stats$time.step>time.step & dloop.stats$length == fragment.type)[1]] = 
                                dloop.stats$average.size[which(dloop.stats$time.step>time.step & dloop.stats$length == fragment.type)[1]] + sum(nchar(zipped.fragments.list$sequences))
                              
                              break
                            }
                            
                          }else{
                            ## KE2 :
                            #For each zipped fragment larger than 32 nts, try KE2 probability, if it pass, start the extension at the end of the i-th zipped fragment ,
                            # and break the fragment loop (go to the next fragment)
                            stop <- FALSE
                            for(i in 1:nrow(zipped.fragments.list)){
                              if (nchar(zipped.fragments.list[i, ]$sequences) >= 32){
                                yy = runif(1)
                                if(yy < ke2.prob){
                                  extensions.stats$time.step[bigtracker] = time.step
                                  extensions.stats$ke[bigtracker] = 2
                                  extensions.stats$clipping.pos[bigtracker] = as.integer(zipped.fragments.list[i, ]$end)
                                  stop <- TRUE
                                  
                                  dloop.stats$count[which(dloop.stats$time.step>time.step & dloop.stats$length == fragment.type)[1]] = 
                                    dloop.stats$count[which(dloop.stats$time.step>time.step & dloop.stats$length == fragment.type)[1]] + nrow(zipped.fragments.list)
                                  
                                  dloop.stats$average.size[which(dloop.stats$time.step>time.step & dloop.stats$length == fragment.type)[1]] = 
                                    dloop.stats$average.size[which(dloop.stats$time.step>time.step & dloop.stats$length == fragment.type)[1]] + sum(nchar(zipped.fragments.list$sequences))
                                  
                                  break
                                }
                              }
                            }
                            if (stop){break}
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
                        
                        ############################################################################
                        
                        if(occupied.rad51$bound != "unbound"){
                          occupied.bins = as.data.frame(table(occupied.rad51$genome.bins))
                          names(occupied.bins) = c("bins", "freq")
                          
                          chromosome.contacts[chromosome.contacts$time.step == time.step & chromosome.contacts$length == fragment.type, as.character(occupied.bins$bins)] =
                            chromosome.contacts[chromosome.contacts$time.step == time.step & chromosome.contacts$length == fragment.type, as.character(occupied.bins$bins)] + occupied.bins$freq
                        }
                       
                        ############################################################################
                        if (time.step == 200 || time.step == 400 || time.step == 600){
                          if (nrow(zipped.fragments.list) >0){
                            dloop.stats$count[dloop.stats$time.step == time.step & dloop.stats$length == fragment.type] = 
                              dloop.stats$count[dloop.stats$time.step == time.step & dloop.stats$length == fragment.type] + nrow(zipped.fragments.list)
                            
                            dloop.stats$average.size[dloop.stats$time.step == time.step & dloop.stats$length == fragment.type] = 
                              dloop.stats$average.size[dloop.stats$time.step == time.step & dloop.stats$length == fragment.type] + sum(nchar(zipped.fragments.list$sequences))
                          }
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
                  dloop.stats$average.size[dloop.stats$count > 0] = dloop.stats$average.size[dloop.stats$count > 0] / dloop.stats$count[dloop.stats$count > 0]
                  
                  write.csv(chromosome.contacts, file=paste(dirnew_data,"/chromosomes_contacts.csv",sep=""))
                  population.time.series(dirnew_data = dirnew_data, dirnew_plots = dirnew_pop, donors.list = donors.list, pop.time.series = pop.time.series)
                  stats.plots(dirnew_plots = dirnew_contacts, occupancy.firsts = occupancy.firsts)
                  
                }
              }
            }
          }
        }
      }
    }  
  }
}

