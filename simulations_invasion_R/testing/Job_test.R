rm(list=ls()) #clean global environment

###Set working directory
setwd("/home/nicolas/Documents/INSA/Stage4BiM/DSB_homologous_recombination_simulation/")

# Directory where you want to save timeseries and plots. Need the slash at the end if you want sub-directories underneath. 
rootdir = "/home/nicolas/Documents/INSA/Stage4BiM/DSB_homologous_recombination_simulation/datas/"


###################### Import librairies #######################################

library(ggplot2)
library(stringr)
library(dplyr)
library(Rcpp)

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
sequences.bins <- read.csv("./LYS2/LY_occurences_per_8bp_(for_rev_donor)_with_bins.csv")[-1] #doesn't take the first column containing sequences

# Import the experimental contacts of the left DSB 10kb with the genome wide :
contacts <- read.csv("./LYS2/leftDSB_contacts_100000_110000_10kb.csv")
bins.id <- paste(as.character(contacts$chrom), "_", as.character(contacts$start_pos), "_", as.character(contacts$end_pos), sep="")
contacts <- cbind(contacts, bins.id)
colnames(contacts)[6] <- "frequency"
colnames(contacts)[7] <- "id"

#########################################################################################################
#################################### FUNCTIONS ##########################################################

#negate operator for %in% :
'%!in%' <- function(x,y)!('%in%'(x,y)) 

#########################################################################################################
#########################################################################################################

find.occupancies = function(lower.window ="none", upper.window = "none", additional.removals = "none"){
  # Find the positions where there is no MH bounded ;
  # Takes into account distance search window, if they are specified ;
  # Additional positions not to be retained can be specified as parameters ;
  # Return a vector of indices of free binding sites ;
  
  indices = 1:(nchar(lys2.fragment) -7)
  if (occupied.rad51$bound=="unbound"){
    return(indices)
  }
  
  remove = c()
  for (i in 0:7){
    remove = c(remove, (occupied.rad51$lys2.microhomology - i), (occupied.rad51$lys2.microhomology + i))
    if(length(additional.removals) > 1){
      remove = c(remove, (additional.removals - i), (additional.removals + i))
    }
  }
  
  if (lower.window != "none"){remove=c(remove, 0:lower.window)} #remove downstream sites of the search window
  if (upper.window != "none"){remove=c(remove, upper.window:nchar(lys2.fragment))} # same, but for the upstream sites
  
  remove = remove[which(remove > 0)] #remove the negative and null indices
  return(indices[-remove])
}

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
    open.sites = 1:(nchar(lys2.fragment) -7)
  }else{
    open.sites = find_occupancies_remastered(occupiedRAD51 = occupied.rad51$lys2.microhomology, additional_removals = c(0), n = nchar(lys2.fragment))
  }
  
  # open.sites = find.occupancies() #indexes of unoccupied sites
  
  if (length(open.sites)== 0){ #if all the sites are occupied
    return(list(bound = occupied.rad51$bound, strand = "negative", genome.bins = c(), donor.invasions = c(), lys2.microhomology = c()))
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
    return(list(bound = occupied.rad51$bound,strand = "negative", genome.bins = c(), donor.invasions = c(), lys2.microhomology = c()))
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
      if(yy <= 0.5){ #probability to be a donor 
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
    remove = which(matches %in% occupied.rad51$lys2.microhomology)
    if (length(remove)>0){
      identities = identities[-remove]
      matches = matches[-remove]
      bins = bins[-remove]
    }
  } 
  
  
  # LYS alignment or misalignments :
  return(list(bound=occupied.rad51$bound, strand = "negative", genome.bins = bins, donor.invasions = identities, lys2.microhomology = matches))
}  
#########################################################################################################
#########################################################################################################

new.microhomologizer = function(occupied.rad51, window, bindings.per.tethering, kon.prob){
  
  genome.bindings = which(occupied.rad51$donor.invasions %!in% donors.blacklist) #list of all bound MHs with the whole genome (genome.wide.sei step)
  new.bindings = list(bound=occupied.rad51$bound, strand = "negative", genome.bins = c(), donor.invasions = c(), lys2.microhomology = c()) #initialize the list that will be return (copy of an empty occupied.rad51 list)
  
  # Check for unbound sites.
  #   If not, don't go futher and just return an empty list.
  #   All the rad51 are matched somewhere in the genome. 
  
  if (length(find_occupancies_remastered(occupiedRAD51 = occupied.rad51$lys2.microhomology, additional_removals = c(0), n = nchar(lys2.fragment))) == 0){
    return(new.bindings)
  }
  
  # bindings : list of sites occupied by another MHs into the search window around the current micros locus ;
  bindings = c()
  bins = c()
  identities= c()
  
  for (binding.index in genome.bindings){
    if (length(bindings) > 0){
      if (length(find_occupancies_remastered(occupiedRAD51 = occupied.rad51$lys2.microhomology, additional_removals = bindings, n = nchar(lys2.fragment))) == 0){
        break
      }
    }else{ 
      if (length(find_occupancies_remastered(occupiedRAD51 = occupied.rad51$lys2.microhomology, additional_removals = c(0), n = nchar(lys2.fragment))) == 0){
        break
      } 
    }
    
    #current.selocus : index of the MH we are currently looking around it (search window) to place another MHs;
    current.selocus = occupied.rad51$lys2.microhomology[binding.index]
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
    upper.window = ifelse(current.selocus + window >nchar(lys2.fragment)-7, nchar(lys2.fragment)-7, current.selocus + window)
    if (additionals[1]== "none"){additionals=c(0)}
    
    open.sites = find_occupancies_remastered(occupiedRAD51 = occupied.rad51$lys2.microhomology, additional_removals = additionals, 
                                             lower_window = lower.window, upper_window = upper.window, n = nchar(lys2.fragment))
    
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
  new.bindings$lys2.microhomology = c(new.bindings$lys2.microhomology, bindings)
  new.bindings$donor.invasions  = c(new.bindings$donor.invasions, identities)
  
  if (occupied.rad51$bound != "unbound"){
    #remove MHs ids in new.bindings that have already been counted as donor in occupied.rad51 :
    remove = which(bindings %in% occupied.rad51$lys2.microhomology)
    if (length(remove) > 0){
      new.bindings$genome.bins = new.bindings$genome.bins[-remove]
      new.bindings$donor.invasions = new.bindings$donor.invasions[-remove]
      new.bindings$lys2.microhomology = new.bindings$lys2.microhomology[-remove]}
  }
  return(new.bindings)
}

#########################################################################################################
#########################################################################################################

donors.generator <- function(template, bins, N = 0){
  new.donors.list<-list(sequence = c(template), bins = c("chr2_470001_480001"), 
                        id = c("LYS"), invasion = c("no", rep("no", times = N)))
  bases <- c("a", "t", "g", "c")
  
  if(N >= 1){
    for (n in 1:N){
      new.donor <- template
      new.donor.id = paste("donor", as.character(n), sep="")
      lower.limit <- floor(0.1*nchar(template))
      upper.limit <- floor(0.4*nchar(template))
      snp.location <- sample(1:nchar(template), size = (sample(lower.limit:upper.limit, size = 1)), replace = FALSE)
      
      for (i in snp.location){
        snp <- sample(bases[-which(bases == substr(template, i, i))], size = 1)
        substr(new.donor, i, i) <- snp
      }
      new.donors.list$sequence = c(new.donors.list$sequence, new.donor)
      new.donors.list$bins = c(new.donors.list$bins, sample(bins, size = 1))
      new.donors.list$id = c(new.donors.list$id, new.donor.id)
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

rad54.rdh54.placement <- function(nb.rad54, nb.rdh54, lys2.fragment){
  
  location.rad54 <- c()
  location.rdh54 <- c()
  while (nb.rad54 > 0){ #place the required rad54 randomly (according uniform distro) over the invading fragment ;
    new.location <- 0
    while (new.location == 0 | new.location %in% location.rad54){
      new.location <- floor(runif(1, min = 1, max=str_length(lys2.fragment)))
    }
    location.rad54 = c(location.rad54, new.location)
    nb.rad54 = nb.rad54 - 1
  }
  
  while (nb.rdh54 > 0){
    new.location <-  0
    while (new.location == 0 | new.location %in% location.rad54 | new.location %in% location.rdh54){
      new.location <- floor(runif(1, min = 1, max=str_length(lys2.fragment)))
    }
    location.rdh54 = c(location.rdh54, new.location) #location.rad54 becomes location.rdh54
    nb.rdh54 = nb.rdh54 - 1
  }
  
  return(list(sort(location.rad54), sort(location.rdh54)))
}

#########################################################################################################
#########################################################################################################
check.before.zipping <- function(current.rad54, donor){
  
  microhomologies.left <- 0
  microhomologies.right <- 0
  left <- 1
  right <- 1
  
  if (donors.occupancy$bound[current.rad54] == "yes" && 
      donors.occupancy$bound.id[current.rad54] == "homology" &&
      donors.occupancy$donor.id[current.rad54] == donor){
    
    while(left !=0 && right != 0){
      while(left != 0){
        if (donors.occupancy$bound[current.rad54 - left] == "yes" && 
            donors.occupancy$bound.id[current.rad54 - left] == "homology" && 
            donors.occupancy$donor.id[current.rad54 - left] == donor &&
            left < current.rad54){
          
          microhomologies.left = microhomologies.left +1
          left = left+1
        }else{
          left = 0
        }
      }
      
      while(right != 0){
        if(donors.occupancy$bound[current.rad54 + right] == "yes" && 
           donors.occupancy$bound.id[current.rad54 + right] == "homology" && 
           donors.occupancy$donor.id[current.rad54 + right] == donor &&
           (current.rad54 + right) < str_length(lys2.fragment)){
          
          microhomologies.right = microhomologies.right +1
          right = right+1
        }else{
          right = 0
        }
      }
    }
  } 
  return(sum(microhomologies.left + 1 + microhomologies.right))
} 

#########################################################################################################
#########################################################################################################
zipping <- function(rad54, zipping.list, donor, limit){
  
  pos <- rad54
  zip.indexe <- c()
  zip.fragment <-"" 
  donor.seq = donors.list$sequence[which(donors.list$id == donor)]
  consecutive.mismatches <- 0
  
  while(pos %!in% pos.rdh54 && 
        pos %!in% pos.rad54[which(pos.rad54 != rad54)] && 
        pos < nchar(lys2.fragment)){
    
    if(donors.occupancy$bound.id[pos] == "homology" && donors.occupancy$donor.id[pos] == donor){
      new.nt <- substr(lys2.fragment, pos, pos)
      zip.indexe = c(zip.indexe, pos)
      zip.fragment = paste(zip.fragment, new.nt, sep="")
      pos = pos + 1
      consecutive.mismatches = 0
      
    }else{
      if (str_sub(string = lys2.fragment, start = pos, end = pos) ==  str_sub(string = donor.seq, start = pos, end = pos)){
        new.nt <- substr(lys2.fragment, pos, pos)
        zip.indexe = c(zip.indexe, pos)
        zip.fragment = paste(zip.fragment, new.nt, sep="")
        pos = pos + 1
        consecutive.mismatches = 0
        
      }else{
        consecutive.mismatches = consecutive.mismatches + 1
        if (consecutive.mismatches >= limit){
          return(-1) # too much misalignment, this donor isn't good enough to lead an HR
          
        }else if (consecutive.mismatches < limit){
          new.nt <- substr(lys2.fragment, pos, pos)
          zip.indexe = c(zip.indexe, pos)
          zip.fragment = paste(zip.fragment, new.nt, sep="")
          pos = pos + 1
        }
      }
    }
  }
  
  if(nchar(zip.fragment) > 16){
    return(c(as.integer(zip.indexe[1]),  as.integer(tail(zip.indexe,1)), zip.fragment) )
    
  }else{
    return(0) #zipped sequence is too short 
  }
}

#########################################################################################################
#########################################################################################################

#########################################################################################################
################################Outputs functions #######################################################

single.runs <-function(dirnew_singles, binding.ts, saver){
  
  outname=paste(dirnew_singles,"/Total_Occupancy_",saver,".png",sep="")
  occ_plot<-
    ggplot(data = binding.ts) + geom_step(aes(x = time.step, y = bound, color = length)) +
    labs(x = "time step", y = "Total Occupancy (bp)") + theme_minimal() + theme(text = element_text(size = 16))+
    scale_y_continuous(limits = c(0, max(binding.ts$bound)+1))
  ggsave(outname,plot=occ_plot)
  
  outname=paste(dirnew_singles,"/Occupancy_Heterologies_",saver,".png",sep="")
  het_plot<-
    ggplot(data = binding.ts) + geom_step(aes(x = time.step, y = heterologies, color = length)) +
    labs(x = "time step", y = "Occupancy at Heterologies (bp)") + theme_minimal()+ theme(text = element_text(size = 16))+
    scale_y_continuous(limits = c(0, max(binding.ts$heterologies)+1))
  ggsave(outname,plot=het_plot)
  
  outname=paste(dirnew_singles,"/Occupancy_Homologies_",saver,".png",sep="")
  het_plot<-
    ggplot(data = binding.ts) + geom_step(aes(x = time.step, y = homologies, color = length)) +
    labs(x = "time step", y = "Occupancy at Homologies (bp)") + theme_minimal()+ theme(text = element_text(size = 16))+
    scale_y_continuous(limits = c(0, max(binding.ts$homologies)+1))
  ggsave(outname,plot=het_plot)
  
}

#########################################################################################################
#########################################################################################################


population.time.series <- function(dirnew_data, dirnew_plots, donors.list, pop.time.series){
  
  # pop time series for all homologies
  df <- subset(x=pop.time.series, select=c(1, 2, 3))
  colnames(df) = c("time.step", "length", "homologies")
  df.name <- paste("pop_timeseries", "homologies", sep = "_")
  write.table(df, file=paste(dirnew_data, "/", df.name ,".txt", sep=""))
  outname=paste(dirnew_plots, "/", df.name, ".png",sep="")
  
  pop.plot<-
    ggplot(data = df) + geom_step(aes(x = time.step, y = homologies, color = length)) +
    labs(x = "time step", y = "Probability of detection for homologies") + theme_minimal()+ theme(text = element_text(size = 14))+
    scale_y_continuous(limits = c(0, max(df$homologies)+1))
  ggsave(outname, plot=pop.plot)
  
  #pop time series for all zipped homologies 
  df <- subset(x=pop.time.series, select=c(1, 2, 4))
  colnames(df) = c("time.step", "length", "zip")
  df.name <- paste("pop_timeseries", "zipped_homologies", sep = "_")
  write.table(df, file=paste(dirnew_data, "/", df.name ,".txt", sep=""))
  outname=paste(dirnew_plots, "/", df.name, ".png",sep="")
  
  pop.plot<-
    ggplot(data = df) + geom_step(aes(x = time.step, y = zip, color = length)) +
    labs(x = "time step", y = "Probability of detection for zips") + theme_minimal()+ theme(text = element_text(size = 14))+
    scale_y_continuous(limits = c(0, max(df$zip)+1))
  ggsave(outname, plot=pop.plot)
  
  
  for (i in 1:length(donors.list$id)){
    df <- subset(x=pop.time.series, select=c(1, 2, i+4))
    colnames(df) = c("time.step", "length", "prob.detect")
    df.name <- paste("pop_timeseries", as.character(donors.list$id[i]),as.character(donors.list$bins[i]), sep = "_")
    write.table(df, file=paste(dirnew_data, "/", df.name ,".txt", sep=""))
    outname=paste(dirnew_plots, "/", df.name, ".png",sep="")
    
    pop.plot<-
      ggplot(data = df) + geom_step(aes(x = time.step, y = prob.detect, color = length)) +
      labs(x = "time step", y = paste("Probability of detection", as.character(donors.list$id[i]), sep = " ")) + 
      theme_minimal()+ theme(text = element_text(size = 14))+
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
####################### Parameters #############################################

num.time.steps = 800 # Length of simulation in time steps
graph.resolution = 1 #save occupancy data at every nth time step. Plots will have this resolution at the x-axis 

test.replicates = 25 # How many times to simulate, replicates
kon.group<-c(0.5) #binding probabilities for every binding try
koff1.group<-c(0.2) # dissociation probabilities for each bound particle
koff2.group<-c(0.02) #dissociation probabilities for each zipped fragments
m.group = c(2) #bindings allowed to occur per tethering
search.window.group = c(250) #the genomic distance of the tethering effect (per side)
rad54.group <- c(1/200) #proportional to the length of invading strand
rdh54.group <- c(1/10) #proportional to the number of rad54
misalignments.cutoff <- 6 #How many mismatches are allowed before break the zipping phase for the current donor 
additional.donors <- 2 # Additional donors ( without LYS2)

# Since the data needs to be outputted to files with human-readable names,we have to label the parameters with strings.
# For example 0005 is really 0.005
kon.group.names<- gsub("\\.", "", as.character(kon.group))
koff1.group.names<- gsub("\\.", "", as.character(koff1.group))
koff2.group.names<- gsub("\\.", "", as.character(koff2.group))
rad54.group.names<-gsub("\\.", "", as.character(round(rad54.group, 4)))

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
rdh54.name=gsub("\\.", "", as.character(round(rad54.prop*rdh54.prop, 4)))

# Initialize the occupied.rad51 list that contains :
#   the state of the invading fragment : bound or unbound ,
#   the nature of the invading strand : negative ,
#   the respective genome bins were heterologies and homologies are found,
#   the nature of the donor :
#     "H" if heterology ;
#     "donor#" if it's a donor (in the donors.list$id)
#     "LYS" if the donor is LYS2
#   the position of rad51 in invading fragment for each bound microhomologies ;

occupied.rad51 = list(bound = "unbound",strand = "negative", genome.bins = c("chr2_470001_480001"), donor.invasions = "LYS", lys2.microhomology = 368)

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
    names(binding.ts) = c('time.step', "length", "bound", "heterologies", "homologies")
    binding.ts$time.step = rep(seq(1,num.time.steps, graph.resolution),3)
    binding.ts$length = rep(ly.names, each = (num.time.steps / graph.resolution))
  }
  
  # exonucleases are involved in the resection process of broken strands before starting the homologies search via rad51
  # The time for this operation is quite random, therefore we simulate it with a normal law :
  exonuclease.job <- as.integer(rnorm(n=1, mean = 25, sd = 10))
  if(exonuclease.job <=1){exonuclease.job = 1}
  
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
    donors.occupancy = as.data.frame(matrix(0, nchar(lys2.fragment),6))
    names(donors.occupancy) = c('bp', 'bound', "bound.id", "donor.id", "zipped", "bins")
    donors.occupancy$bp = 1:nchar(lys2.fragment)
    donors.occupancy$bound = "no"
    donors.occupancy$bound.id = "unbound"
    donors.occupancy$donor.id = "unknown"
    donors.occupancy$zipped = "no"
    donors.occupancy$bins = "unknown"
    
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
    
    #probability of detection proportional to the length of invading strand :
    crosslink.density <- 500 * 1.5 * (nchar(lys2.fragment) / as.integer(max(ly.type))) 
    
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
          occupied.rad51$lys2.microhomology = c(occupied.rad51$lys2.microhomology, new.bindings$lys2.microhomology)
        }
      }
      
      # Check if exonuclease opening DSB is done ;
      if(time.step > exonuclease.job){
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
          occupied.rad51$lys2.microhomology = occupied.rad51$lys2.microhomology[-remove.rad51]
        }
      }
      
      ##########################################################################
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
      ############################## Occupancy #################################
      donors.occupancy$bound = "no"
      donors.occupancy$bound.id = "unbound"
      donors.occupancy$donor.id = "unknown"
      donors.occupancy$bins = "unknown"
      
      if(occupied.rad51$bound != "unbound"){
        
        rad51.cover.index = occupied.rad51$lys2.microhomology
        for(i in 1:7){rad51.cover.index = c(rad51.cover.index, occupied.rad51$lys2.microhomology+i)}
        
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
      ##########################################################################
      ################################# Zipping ################################
      # When the twoh microhomology state is enable, the zipping occurs until all rad54 are zipped;
      if(length(unzipped.rad54 > 0) && current.donor != "" &&
         donors.list$invasion[which(donors.list$id == current.donor)] != "failed"){
        
        if (donors.list$invasion[which(donors.list$id == current.donor)] =="no"){
          donors.list$invasion[which(donors.list$id == current.donor)] ="yes"
        }
        
        for (pos in unzipped.rad54){
          
          # Check if the sequence to zip is big enough ;
          #   We decided >= 16 (2*8 nts) arbitrary (could be more or less)
          if(donors.occupancy$zipped[pos] != "yes" & check.before.zipping(pos, donor = current.donor) >= 16){
            new.zip = zipping(pos, zipped.fragments.list, donor= current.donor, limit = misalignments.cutoff)
            
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
              occupied.rad51$lys2.microhomology = occupied.rad51$lys2.microhomology[-remove.rad51]
              
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
            donors.occupancy$bins[current.zip.start : current.zip.end] = "unknown"
            
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
          if(candidate.donor != "unknown" && candidate.donor != "H" && length(which(donors.occupancy$donor.id == candidate.donor)) >= 200 && 
             donors.list$invasion[which(donors.list$id == candidate.donor)] == "no"){
            
            current.donor = candidate.donor
          }
        }
      }
      
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
      
      
      
      if(saver < 3 ){
        # tabulate occupancies vs. time step and length
        binding.ts$bound[binding.ts$time.step == time.step & 
                           binding.ts$length == ly.type] = length(which(donors.occupancy$bound == "yes"))
        
        binding.ts$heterologies[binding.ts$time.step == time.step & 
                                  binding.ts$length == ly.type] = length(which(donors.occupancy$bound.id == "heterology"))
        
        binding.ts$homologies[binding.ts$time.step == time.step & 
                                binding.ts$length == ly.type] = length(which(donors.occupancy$bound.id == "homology"))
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
        occupied.bins = as.data.frame(table(occupied.rad51$genome.bins))
        names(occupied.bins) = c("bins", "freq")
        chromosome.contacts[chromosome.contacts$time.step == time.step &
                              chromosome.contacts$length == ly.type, as.character(occupied.bins$bins)] = 
          chromosome.contacts[chromosome.contacts$time.step == time.step &
                                chromosome.contacts$length == ly.type, as.character(occupied.bins$bins)] + occupied.bins$freq
      }
      
      #print(c(ly.type, time.step, trial))
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