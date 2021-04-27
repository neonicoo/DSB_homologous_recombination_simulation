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

'%!in%' <- function(x,y)!('%in%'(x,y)) #negate operator for %in%

# Directory where you want to save timeseries and plots. Need the slash at the end if you want sub-directories underneath. 
rootdir = "/home/nicolas/Documents/INSA/Stage4BiM/DSB_homologous_recombination_simulation/datas/";

# genome-wide microhomology counts
forward.sequences <- read.table("./Occurences_per_8bp_motif(for+rev_donor).txt", sep="\t", header = TRUE)
forward.sequences = forward.sequences[,c("start", "sequence", "total")]
row.names(forward.sequences) = 1:nrow(forward.sequences)
microhomology.probs = forward.sequences$total / sum(forward.sequences$total)

# within-lys microhomologies (misalignments)
L500.self.micros = as.data.frame(matrix(c("aacaagct","aacaagct",98,319,319,98),2,3),stringsAsFactors = F); names(L500.self.micros) = c("L500", "position1", "position2"); L500.self.micros$position3 = NA
L1000.selfmicros <- read.delim("L1000_self-microhomologies.txt", stringsAsFactors=FALSE); L1000.selfmicros$position3 = NA
LY2000.selfmicros <- read.csv("LY2000_self-microhomologies.txt", sep="", stringsAsFactors=FALSE)

# Name the DNA sequences of the invading strands

LY = (tolower("ATGACTAACGAAAAGGTCTGGATAGAGAAGTTGGATAATCCAACTCTTTCAGTGTTACCACATGACTTTTTACGCCCACAACAAGAACCTTATACGAAACAAGCTACATATTCGTTACAGCTACCTCAGCTCGATGTGCCTCATGATAGTTTTTCTAACAAATACGCTGTCGCTTTGAGTGTATGGGCTGCATTGATATATAGAGTAACCGGTGACGATGATATTGTTCTTTATATTGCGAATAACAAAATCTTAAGATTCAATATTCAACCAACGTGGTCATTTAATGAGCTGTATTCTACAATTAACAATGAGTTGAACAAGCTCAATTCTATTGAGGCCAATTTTTCCTTTGACGAGCTAGCTGAAAAAATTCAAAGTTGCCAAGATCTGGAAAGGACCCCTCAGTTGTTCCGTTTGGCCTTTTTGGAAAACCAAGATTTCAAATTAGACGAGTTCAAGCATCATTTAGTGGACTTTGCTTTGAATTTGGATACCAGTAATAATGCGCATGTTTTGAACTTAATTTATAACAGCTTACTGTATTCGAATGAAAGAGTAACCATTGTTGCGGACCAATTTACTCAATATTTGACTGCTGCGCTAAGCGATCCATCCAATTGCATAACTAAAATCTCTCTGATCACCGCATCATCCAAGGATAGTTTACCTGATCCAACTAAGAACTTGGGCTGGTGCGATTTCGTGGGGTGTATTCACGACATTTTCCAGGACAATGCTGAAGCCTTCCCAGAGAGAACCTGTGTTGTGGAGACTCCAACACTAAATTCCGACAAGTCCCGTTCTTTCACTTATCGCGACATCAACCGCACTTCTAACATAGTTGCCCATTATTTGATTAAAACAGGTATCAAAAGAGGTGATGTAGTGATGATCTATTCTTCTAGGGGTGTGGATTTGATGGTATGTGTGATGGGTGTCTTGAAAGCCGGCGCAACCTTTTCAGTTATCGACCCTGCATATCCCCCAGCCAGACAAACCATTTACTTAGGTGTTGCTAAACCACGTGGGTTGATTGTTATTAGAGCTGCTGGACAATTGGATCAACTAGTAGAAGATTACATCAATGATGAATTGGAGATTGTTTCAAGAATCAATTCCATCGCTATTCAAGAAAATGGTACCATTGAAGGTGGCAAATTGGACAATGGCGAGGATGTTTTGGCTCCATATGATCACTACAAAGACACCAGAACAGGTGTTGTAGTTGGACCAGATTCCAACCCAACCCTATCTTTCACATCTGGTTCCGAAGGTATTCCTAAGGGTGTTCTTGGTAGACATTTTTCCTTGGCTTATTATTTCAATTGGATGTCCAAAAGGTTCAACTTAACAGAAAATGATAAATTCACAATGCTGAGCGGTATTGCACATGATCCAATTCAAAGAGATATGTTTACACCATTATTTTTAGGTGCCCAATTGTATGTCCCTACTCAAGATGATATTGGTACACCGGGCCGTTTAGCGGAATGGATGAGTAAGTATGGTTGCACAGTTACCCATTTAACACCTGCCATGGGTCAATTACTTACTGCCCAAGCTACTACACCATTCCCTAAGTTACATCATGCGTTCTTTGTGGGTGACATTTTAACAAAACGTGATTGTCTGAGGTTACAAACCTTGGCAGAAAATTGCCGTATTGTTAATATGTACGGTACCACTGAAACACAGCGTGCAGTTTCTTATTTCGAAGTTAAATCAAAAAATGACGATCCAAACTTTTTGAAAAAATTGAAAGATGTCATGCCTGCTGGTAAAGGTATGTTGAACGTTCAGCTACTAGTTGTTAACAGGAACGATCGTACTCAAATATGTGGTATTGGCGAAATAGGTGAGATTTATGTTCGTGCAGGTGGTTTGGCCGAAGGTTATAGAGGATTACCAGAATTGAATAAAGAAAAATTTGTGAACAACTGGTTTGTTGAAAAAGATCACTGGAATTATTTGGATAAGGATAATGGTGAACCTTGGAGACAATTCTGGTTAGGTCCAAGAGATAGATTGTACAGAACGGGTGATTTAGGTCGTTATCTACCAAACGG"))
L = (tolower("ATGACTAACGAAAAGGTCTGGATAGAGAAGTTGGATAATCCAACTCTTTCAGTGTTACCACATGACTTTTTACGCCCACAACAAGAACCTTATACGAAACAAGCTACATATTCGTTACAGCTACCTCAGCTCGATGTGCCTCATGATAGTTTTTCTAACAAATACGCTGTCGCTTTGAGTGTATGGGCTGCATTGATATATAGAGTAACCGGTGACGATGATATTGTTCTTTATATTGCGAATAACAAAATCTTAAGATTCAATATTCAACCAACGTGGTCATTTAATGAGCTGTATTCTACAATTAACAATGAGTTGAACAAGCTCAATTCTATTGAGGCCAATTTTTCCTTTGACGAGCTAGCTGAAAAAATTCAAAGTTGCCAAGATCTGGAAAGGACCCCTCAGTTGTTCCGTTTGGCCTTTTTGGAAAACCAAGATTTCAAATTAGACGAGTTCAAGCATCATTTAGTGGACTTTGCTTTGAATTTGGATACCAGTAATAATGCGCATGTTTTGAACTTAATTTATAACAGCTTACTGTATTCGAATGAAAGAGTAACCATTGTTGCGGACCAATTTACTCAATATTTGACTGCTGCGCTAAGCGATCCATCCAATTGCATAACTAAAATCTCTCTGATCACCGCATCATCCAAGGATAGTTTACCTGATCCAACTAAGAACTTGGGCTGGTGCGATTTCGTGGGGTGTATTCACGACATTTTCCAGGACAATGCTGAAGCCTTCCCAGAGAGAACCTGTGTTGTGGAGACTCCAACACTAAATTCCGACAAGTCCCGTTCTTTCACTTATCGCGACATCAACCGCACTTCTAACATAGTTGCCCATTATTTGATTAAAACAGGTATCAAAAGAGGTGATGTAGTGATGATCTATTCTTCTAGGGGTGTGGATTTGATGGTATGTGTGATGGGTGTCTTGAAAGCCGGCGCAACCTTTTCAGTTATCGACCCTGCATATCCCCCAGCCAGACAAACCATTTACTTAGGTGTTGCTAAACCACGTGGGTTGATTGTTATTA"))
L500 = (tolower("ATGACTAACGAAAAGGTCTGGATAGAGAAGTTGGATAATCCAACTCTTTCAGTGTTACCACATGACTTTTTACGCCCACAACAAGAACCTTATACGAAACAAGCTACATATTCGTTACAGCTACCTCAGCTCGATGTGCCTCATGATAGTTTTTCTAACAAATACGCTGTCGCTTTGAGTGTATGGGCTGCATTGATATATAGAGTAACCGGTGACGATGATATTGTTCTTTATATTGCGAATAACAAAATCTTAAGATTCAATATTCAACCAACGTGGTCATTTAATGAGCTGTATTCTACAATTAACAATGAGTTGAACAAGCTCAATTCTATTGAGGCCAATTTTTCCTTTGACGAGCTAGCTGAAAAAATTCAAAGTTGCCAAGATCTGGAAAGGACCCCTCAGTTGTTCCGTTTGGCCTTTTTGGAAAACCAAGATTTCAAATTAGACGAGTTCAAGCATCATTTAGTGGACTTTGCTTTGAATTTGGATACCAG"))

ly.names = c("500", "1000", "2000")
ly.sequences = c(L500, L, LY)


#Import of the chr2.fa sequence file from the yeast genome (S288) :
yeast.genome.chr2 <- read.fasta("./yeast-genome/S288c-R64-2-1-v2014/chr2.fa" ,
                                seqtype = 'DNA', as.string = TRUE, 
                                forceDNAtolower  = TRUE, set.attributes = FALSE)

yeast.genome.chr2 <- yeast.genome.chr2[[1]] #select just the nucleotides sequence

num.time.steps = 600 # Length of simulation in time steps
test.replicates = 50 # How many times to simulate, replicates
graph.resolution = 1 #save occupancy data at every nth time step. Plots will have this resolution at the x-axis 
kon.group<-c(0.005,0.05,0.1,0.4,0.7,0.9) #binding probabilities for every binding try
koff1.group<-c(0,0.0001,0.05,0.6) # dissociation probabilities for each bound particle
koff2.group<-c(0, 0.0005, 0.05, 0.5) #dissociation probabilities for each zipped fragments
m.group = c(2,5) #bindings allowed to occur per tethering
search.window.group = c(250,500) #the genomic distance of the tethering effect (per side)

# Since the data needs to be outputted to files with human-readable names,we have to label the parameters with strings.
# For example 0005 is really 0.005
kon.group.names<- gsub("\\.", "", as.character(kon.group))
koff1.group.names<- gsub("\\.", "", as.character(koff1.group))
koff2.group.names<- gsub("\\.", "", as.character(koff2.group))

print(kon.group.names)
print(koff1.group.names)
print(koff2.group.names)

#########################################################################################################
#################################### FUNCTIONS ##########################################################

find.occupancies = function(lower.window ="none", upper.window = "none", additional.removals = "none"){
  # Find the positions where there is no MH bounded ;
  # Takes into account distance search window, if they are specified ;
  # Additional positions not to be retained can be specified as parameters ;
  # Return a vector of indices of free binding sites ;
  
  indices = 1:(nchar(lys2.fragment) -7)
  if (occupied.rad51$bound=="unbound"){
    return(indices)
  }
  
  remove = c(occupied.rad51$lys2.microhomology)
  for (i in 1:7){
    remove = c(remove, (occupied.rad51$lys2.microhomology - i), 
               (occupied.rad51$lys2.microhomology + i))
  }
  
  if (lower.window != "none"){remove=c(remove, 0:lower.window)} #remove downstream sites of the search window
  if (upper.window != "none"){remove=c(remove, upper.window:nchar(lys2.fragment))} # same, but for the upstream sites
  
  if (additional.removals[1] != "none"){
    for (i in 0:7){
      remove = c(remove, (additional.removals - i), (additional.removals + i))
    }
  }
  remove = remove[which(remove > 0)] #remove the negative and null indices 
  return(indices[-remove])
}

#########################################################################################################
#########################################################################################################


genome.wide.sei = function(initial.binding.tries){
  
  # Choose region of LY weighted by available microhomologies (MHs);
  # Dont let those already occupied be chosen (by the use of the find.occupancies() function) ;
  
  
  open.sites = find.occupancies() #indexes of unoccupied sites
  
  if (length(open.sites)== 0){ #if all the sites are occuped
    return(list(bound = occupied.rad51$bound,strand = "negative", donor.invasions = c(), lys2.microhomology = c()))
  }
  
  #matches : vector of possible bounding sites for MHs
  matches = c()
  for (i in 1:initial.binding.tries){
    if (length(open.sites) == 1){ 
      # If there is only one available binding site :
      matches[i] = open.sites
    }else{ 
      # Draw a site among those available which will be matched with a MH according to the respective weighted probabilities :
      matches[i] = sample(x=open.sites, size=1, prob = microhomology.probs[open.sites]) 
    }
    
    # Where there is a match with a MH, we consider that the site concerned is no longer available ;
    # We remove it from the index open.sites;
    # To be sure that the next matches will not overlap with the previous one, 
    #  we have to remove from the index the 7 positions upstream and downstream of the match :
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
  
  if (length(matches)<1){
    return(list(bound = occupied.rad51$bound,strand = "negative",donor.invasions = c(), lys2.microhomology = c()))
  }
  
  # Set IDs of each bound (Heterology (H) vs LYS); 
  # If LYS, set as genomic position in + strand notation
  
  # id.probs : probability for a MH to be homologous or heterologous, according it's occurrence in the genome wide (and self-occurrences) 
  # identities : list of id (H or LYS) according the above id.probs for each match ;
  id.probs = sapply(matches, function(x){(1 + ifelse(x %in% self.micros$position1, 1,0))/forward.sequences$total[x]})
  identities = sapply(1:length(matches), function(x) {sample(c("H", "LYS"), 1, prob = c(sapply(id.probs[x], function(x) max(0, 1-x)), id.probs[x]))})
  
  # donor.ids : vector of homologous MHs  ;
  donor.ids = matches[which(identities == "LYS")]
  
  # If the MH is also a self-micro,
  #  Choose randomly a binding site  between all the self micros occurrences (position1, position2, position3 etc...:
  for (index in which(donor.ids %in% self.micros$position1)){
    y = which(self.micros$position1 == donor.ids[index])[1]
    sampling.micros = c(self.micros$position1[y], self.micros$position2[y])
    if(!is.na(self.micros$position3[y])){
      sampling.micros = c(sampling.micros,self.micros$position3[y])
    }
    donor.ids[index] = sample(as.numeric(sampling.micros), size = 1)
  }
  
  # 473927 - donor.ids : position of the last nt for the ly.sequence in the chr II - postion of the donor where the MH is "LYS"
  identities[which(identities == "LYS")] = 473927 - donor.ids
  
  if (occupied.rad51$bound != "unbound"){
    remove = which((identities !="H") & (as.character(identities) %in% occupied.rad51$donor.invasions) )
    if (length(remove)>0){
      identities = identities[-remove]; matches = matches[-remove]
    }
  } 
  
  # LYS alignment or misalignments :
  return(list(bound=occupied.rad51$bound, strand = "negative", donor.invasions = identities, lys2.microhomology = matches))
}  
#########################################################################################################
#########################################################################################################

new.microhomologizer = function(occupied.rad51, window, bindings.per.tethering){
  # When a binding is homologous, we have to search others MHs in a distance search windows ;
  # The number of bindings per search window is the bindings.per.tethering variable ;
  
  
  # correct.binding : vector of indexes for the micro-homologies (LYS) donors ;
  # new.bindings : deep copy of an empty occupied.rad51 ;
  correct.bindings = as.numeric(which(occupied.rad51$donor.invasions != "H"))
  new.bindings = list(bound=occupied.rad51$bound, strand = "negative", donor.invasions = c(), lys2.microhomology = c())
  
  # Check for unbound sites :
  if (length(find.occupancies()) == 0){
    return(new.bindings)
  }
  
  # bindings : list of sites occupied by another MHs into the search window around the current micros locus ;
  bindings = c()
  
  for (binding.index in correct.bindings){
    if (length(bindings) > 0){
      if (length(find.occupancies(additional.removals = bindings)) == 0){break}
    }else{ 
      if (length(find.occupancies()) == 0){break} 
    }
    #current.selocus : index of the MH we are currently looking around it (search window) to place another MHs;
    current.selocus = occupied.rad51$lys2.microhomology[binding.index]
    
    if (length(bindings) <=0){
      additionals = "none"
    }else{
      additionals = bindings
    }
    
    
    # We look if we have available open sites around the current MH locus, i.e into the search windows around it;
    open.sites = find.occupancies(lower.window = current.selocus - window,
                                  upper.window = current.selocus + window, 
                                  additional.removals = additionals)
    
    if (length(open.sites) <= 0){next}
    current.bindings = c()
    
    
    # Bind an unoccupied site located in the search window around our current locus ;
    for (j in 1:bindings.per.tethering){
      if (length(open.sites)==1){
        current.bindings[j] = open.sites
      }else{
        candidate = sample(open.sites, size = 1)
        yy = runif(1)
        if(yy <= kon){
          current.bindings = c(current.bindings,candidate)
        }
      }
      # Remove the candidates and the 7 positions upstream and downstream it from the free sites index ;
      open.sites = open.sites[-which(open.sites %in% (current.bindings[j]-7):(current.bindings[j] + 7))]
      if (length(open.sites) <1){break}
    }
    bindings = c(bindings, current.bindings)
  }
  
  # donor.ids : same as in the genome.wide.sei function  ;
  donor.ids = bindings
  for (index in which(donor.ids %in% self.micros$position1)){
    y = which(self.micros$position1 == donor.ids[index])[1]
    sampling.micros = c(self.micros$position1[y], self.micros$position2[y])
    if(!is.na(self.micros$position3[y])){
      sampling.micros =   c(sampling.micros,self.micros$position3[y])
    }
    donor.ids[index] = sample(as.numeric(sampling.micros), size = 1)
  }
  identities = 473927 - donor.ids
  new.bindings$lys2.microhomology = c(new.bindings$lys2.microhomology, bindings)
  new.bindings$donor.invasions    = c(new.bindings$donor.invasions, identities)
  
  
  if (occupied.rad51$bound != "unbound"){
    #remove MHs ids in new.bindings that have already been counted as donor in occupied.rad51 :
    remove = which((new.bindings$donor.invasions !="H") & (as.character(new.bindings$donor.invasions) %in% occupied.rad51$donor.invasions) )
    if (length(remove) > 0){
      new.bindings$donor.invasions = new.bindings$donor.invasions[-remove]
      new.bindings$lys2.microhomology = new.bindings$lys2.microhomology[-remove]}
  }
  return(new.bindings)
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

rad54.rdh54.placement <- function(nb.rad54, nb.rdh54){
  
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
check.before.zipping <- function(current.rad54){
  
  microhomologies.left <- 0
  microhomologies.right <- 0
  left <- 1
  right <- 1
  
  if (lys2.occupancy$bound[current.rad54] == "yes" && lys2.occupancy$id[current.rad54] == "homology"){
    while(left !=0 && right != 0){
      while(left != 0){
        if (lys2.occupancy$bound[current.rad54 - left] == "yes" && lys2.occupancy$id[current.rad54 - left] == "homology" && left < current.rad54){
          microhomologies.left = microhomologies.left +1
          left = left+1
        }else{
          left = 0
        }
      }
      
      while(right != 0){
        if(lys2.occupancy$bound[current.rad54 + right] == "yes" && lys2.occupancy$id[current.rad54 + right] == "homology" && (current.rad54 + right) < str_length(lys2.fragment)){
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
zipping <- function(rad54, zipping.list){
  
  pos <- rad54
  zip.indexe <- c()
  zip.fragment <-"" 
  zip.junction <- 0
  new.zipping.list <- zipping.list
  
  while(pos %!in% pos.rdh54 && lys2.occupancy$id[pos] == "homology" && pos < nchar(lys2.fragment)+1){
    if (pos %in% zipping.list$start){
      zip.junction = 1
      zip.fragment = paste(zip.fragment, zipping.list$sequences[which(zipping.list$start == pos)], sep = "")
      new.zipping.list = new.zipping.list[-c(which(new.zipping.list$start == pos)),]
      new.zipping.list = rbind(new.zipping.list, c(rad54, rad54+nchar(zip.fragment)-1, zip.fragment))
      break
      
    }else{
      new.nt <- substr(lys2.fragment, pos, pos)
      zip.indexe = c(zip.indexe, pos)
      zip.fragment = paste(zip.fragment, new.nt, sep="")
      pos = pos+1
    }
  }
  
  if(zip.junction == 0){
    new.zipping.list  = rbind(zipping.list , c(as.integer(zip.indexe[1]), as.integer(tail(zip.indexe,1)), zip.fragment))
  }
  
  names(new.zipping.list ) = c("start", "end", "sequences")
  return(new.zipping.list)
}


#########################################################################################################
######################################### Temporary simulation ##########################################


kon = 2; koff = 3; m = 2; sw = 2; koff2 = 3
kon.prob=kon.group[kon]
koff1.prob=koff1.group[koff]
koff2.prob=koff2.group[koff2]

bindings.per.tethering = m.group[m]
search.window = search.window.group[sw]

kon.name=kon.group.names[kon]
koff1.name=koff1.group.names[koff]
koff2.name=koff2.group.names[koff2]

# Initialize the occupied.rad51 vector, genomic (start) position of RAD51 particles (bp / 8) of invaded strand ;
occupied.rad51 = list(bound = "unbound",strand = "negative", donor.invasions = 473927 - 368, lys2.microhomology = 368)

pop.time.series = as.data.frame(matrix(0,num.time.steps*3,3))
names(pop.time.series) = c("time.step","prob.detect","length")
pop.time.series$length = rep(ly.names, each = num.time.steps)
pop.time.series$time.step = rep(seq(1,num.time.steps,1),3)

# Create a variable, saver, that keeps track of how many individual simulations you want to save. 
saver = 0
bigtracker = 0

occupancy.firsts = as.data.frame(matrix(-1, 3*test.replicates, 4))
names(occupancy.firsts) = c("length", "first.bound", "twoh.bound", "first.twoh.time.diff")
occupancy.firsts$length = rep(ly.names, times = test.replicates)

stats.zipping = as.data.frame(matrix(-1, 3*test.replicates, 3))
names(stats.zipping) = c("length", "first.zip", "half.detect")
stats.zipping$length = rep(ly.names, times = test.replicates)

dirname=paste(num.time.steps, kon.name, koff1.name, koff2.name, bindings.per.tethering, search.window, sep="_")

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

# Now make all the replicates
for (trial in 1:test.replicates){ 
  # initialize tabulation of bound microhomologies/heterologies
  # print(trial)
  
  if(saver < 3){
    ly.binding.ts = as.data.frame(matrix(0, (num.time.steps/graph.resolution)*3,4))
    names(ly.binding.ts) = c('time.step', 'bound', "length", "heterologies")
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
    
    SEI.binding.tries = floor((nchar(lys2.fragment)-7)/8)
    
    occupied.rad51$bound = "unbound"
    occupied.rad51$donor.invasions = c()
    occupied.rad51$lys2.microhomology = c()
    
    first = 0 
    twoh = 0
    
    first.zip <- 0
    half.detect <- 0
    
    lys2.occupancy = as.data.frame(matrix(0, nchar(lys2.fragment),4))
    names(lys2.occupancy) = c('bp', 'bound', "id", "zipped")
    lys2.occupancy$bp = 1:nchar(lys2.fragment)
    lys2.occupancy$bound = "no"
    lys2.occupancy$id = "heterology"
    lys2.occupancy$zipped = "no"
    
    rad54 <- floor(0.005*str_length(lys2.fragment)) #number of rad54 to be placed into the invading strand ;
    rdh54 <- floor(0.1*rad54)+1 # number of rdh54 to be placed into the invading strand;
    rad54.rdh54.locations <- rad54.rdh54.placement(rad54, rdh54) 
    pos.rad54 <- rad54.rdh54.locations[[1]] #positions of rad54 in the invading strand;
    pos.rdh54 <- rad54.rdh54.locations[[2]] #positions of rdh54 in the invading strand;
    
    zipped.fragments.list <- as.data.frame(matrix(0,0,3))
    names(zipped.fragments.list ) = c("start", "end", "sequences")
    
    unzipped.rad54 <- pos.rad54[which(lys2.occupancy$zipped[pos.rad54]=="no")]
    
    # Loop through the time-steps
    for (time.step in 1:num.time.steps){
      if(kon.prob == 0){
        next
      }
      
      if (occupied.rad51$bound != "unbound"){
        if (length(occupied.rad51$donor.invasions) != sum(occupied.rad51$donor.invasions == "H")){
          new.bindings = new.microhomologizer(occupied.rad51, search.window, bindings.per.tethering)
          occupied.rad51$donor.invasions = c(occupied.rad51$donor.invasions,new.bindings$donor.invasions)
          occupied.rad51$lys2.microhomology = c(occupied.rad51$lys2.microhomology, new.bindings$lys2.microhomology)
        }
      }
      

      if(twoh == 1 && length(unzipped.rad54 > 0)){
        for (pos in unzipped.rad54){
          if(lys2.occupancy$zipped[pos] != "yes" && check.before.zipping(pos) >= 16){
            zipped.fragments.list = zipping(pos, zipped.fragments.list)

            if(nrow(zipped.fragments.list) > 0){
              if(nchar(tail(zipped.fragments.list$sequences,1)) < 16){
                zipped.fragments.list = zipped.fragments.list[-c(dim(zipped.fragments.list)[1]), ]
                break

              }else{
                current.zip.start <- as.integer(tail(zipped.fragments.list,1)$start)
                current.zip.end <- as.integer(tail(zipped.fragments.list,1)$end)
                lys2.occupancy$zipped[current.zip.start : current.zip.end] = "yes"
                unzipped.rad54 = pos.rad54[which(lys2.occupancy$zipped[pos.rad54]=="no")]
                break
              }
            }
          }
        }
      }

      if(nrow(zipped.fragments.list) > 0){
        for(i in 1:nrow(zipped.fragments.list)){
          preserved.zip <- sample(c(FALSE, TRUE), size =1, replace = TRUE, prob = c(koff2.prob,1-koff2.prob))
          if(!preserved.zip){
            current.zip.start <- as.integer(zipped.fragments.list[i, ]$start)
            current.zip.end <- as.integer(zipped.fragments.list[i, ]$end)
            lys2.occupancy$zipped[current.zip.start : current.zip.end] = "no"
            unzipped.rad54 = pos.rad54[which(lys2.occupancy$zipped[pos.rad54]=="no")]
          }
        }
      }
      
      new.bindings = genome.wide.sei(SEI.binding.tries)
      
      if (occupied.rad51$bound == "unbound"){
        occupied.rad51 = new.bindings;
        if(length(new.bindings$lys2.microhomology) > 0){
          occupied.rad51$bound = "bound"
        }
      }
      else{
        # print("bound and adding")
        occupied.rad51$donor.invasions = c(occupied.rad51$donor.invasions, new.bindings$donor.invasions)
        occupied.rad51$lys2.microhomology = c(occupied.rad51$lys2.microhomology, new.bindings$lys2.microhomology)
      }
      
      if(occupied.rad51$bound != "unbound"){
        for (i in 1:length(occupied.rad51$lys2.microhomology)){
          lys2.occupancy$bound[occupied.rad51$lys2.microhomology[i]:(occupied.rad51$lys2.microhomology[i] + 7)] = "yes"
          if(occupied.rad51$donor.invasions[i] != "H"){
            lys2.occupancy$id[occupied.rad51$lys2.microhomology[i]:(occupied.rad51$lys2.microhomology[i] + 7)]  = "homology"
          }
        }
        prob.detection = length(which(lys2.occupancy$zipped == "yes"))
        prob.detection = prob.detection/500
        
      }else{
        prob.detection = 0
      }
      
      if(prob.detection >= 1){
        prob.detection = 1
      }

      
      lys2.occupancy$id[lys2.occupancy$bound == "no"] = "unbound"
      
      if(length(which(lys2.occupancy$id == "homology")) > 0 && first == 0){
        first = 1;
        occupancy.firsts$first.bound[bigtracker] = time.step
      }
      
      if(length(which(lys2.occupancy$id == "homology")) >= 200 && twoh == 0){
        twoh = 1;
        occupancy.firsts$twoh.bound[bigtracker] = time.step
        occupancy.firsts$first.twoh.time.diff[bigtracker] = time.step - occupancy.firsts$first.bound[bigtracker]
      }
      
      if(length(which(lys2.occupancy$zipped == "yes")) > 0 && first.zip == 0){
        first.zip = 1
        stats.zipping$first.zip[bigtracker] = time.step
      }
      
      if(prob.detection > 0.5 & half.detect == 0){
        half.detect = 1
        stats.zipping$half.detect[bigtracker] = time.step
        }
      
      if(saver <3 ){
        # tabulate occupancies vs. time step and length
        ly.binding.ts$bound[ly.binding.ts$time.step == time.step & ly.binding.ts$length == ly.type] = length(which(lys2.occupancy$bound == "yes"))
        ly.binding.ts$heterologies[ly.binding.ts$time.step == time.step & ly.binding.ts$length == ly.type] = length(which(lys2.occupancy$id == "heterology"))
      }
      
      pop.time.series$prob.detect[pop.time.series$time.step == time.step & pop.time.series$length == ly.type] = pop.time.series$prob.detect[pop.time.series$time.step == time.step & pop.time.series$length == ly.type] + prob.detection
      
      #simulate random dissociation
      num.bound = length(occupied.rad51$donor.invasions)
      #  print(num.bound)
      preserved = sample(c(FALSE,TRUE), num.bound, replace = TRUE, prob = c(koff1.prob,1-koff1.prob)) #dissociate if FALSE
      
      occupied.rad51$donor.invasions = occupied.rad51$donor.invasions[preserved]
      occupied.rad51$lys2.microhomology = occupied.rad51$lys2.microhomology[preserved]
      
      if (sum(!preserved)==num.bound){
        occupied.rad51$bound = "unbound"
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
    
    occ_plot<-ggplot(data = ly.binding.ts) + geom_step(aes(x = time.step, y = bound, color = length)) +
      labs(x = "time step", y = "Total Occupancy (bp)") + theme_minimal() + theme(text = element_text(size = 16))+
      scale_y_continuous(limits = c(0, 2070))
    ggsave(outname,plot=occ_plot)
    
    outname=paste(dirnew_singles,"/Occupancy_Heterologies_",saver,".png",sep="")
    het_plot<-ggplot(data = ly.binding.ts) + geom_step(aes(x = time.step, y = heterologies, color = length)) +
      labs(x = "time step", y = "Occupancy at Heterologies (bp)") + theme_minimal()+ theme(text = element_text(size = 16))+           scale_y_continuous(limits = c(0, 2070))
    ggsave(outname,plot=het_plot)
    
    outname=paste(dirnew_singles,"/Occupancy_Lys2_",saver,".png",sep="")
    ly.binding.ts$homologies = ly.binding.ts$bound - ly.binding.ts$heterologies
    lys2_plot<-ggplot(data = ly.binding.ts) + geom_step(aes(x = time.step, y = homologies, color = length)) +
      labs(x = "time step", y = "Occupancy at Lys2 (bp)") + theme_minimal()+ theme(text = element_text(size = 16))+
      scale_y_continuous(limits = c(0, 2070))
    ggsave(outname,plot=lys2_plot)
  }
  
  saver=saver+1
}#end process

# population timeseries
write.table(pop.time.series, file=paste(dirnew_data,"/population_timeseries.txt",sep=""))

outname=paste(dirnew_plots,"/population_time_series.png",sep="")
pop.plot<-ggplot(data = pop.time.series) + geom_step(aes(x = time.step, y = prob.detect, color = length)) +
  labs(x = "time step", y = "Probability of Detection") + theme_minimal()+ theme(text = element_text(size = 16))+
  scale_y_continuous(limits = c(0, test.replicates))
ggsave(outname,plot=pop.plot)


#### Histograms+ boxplot for the first and twoh MH
occupancy.firsts2 <- occupancy.firsts[-c(which(occupancy.firsts$first.bound == -1)),]

final.firsts = as.data.frame(matrix(-1,test.replicates,3))
names(final.firsts) = c("500","1000","2000")
final.firsts$`500` = occupancy.firsts$first.bound[which(occupancy.firsts$length == 500)]
final.firsts$`1000` = occupancy.firsts$first.bound[which(occupancy.firsts$length == 1000)]
final.firsts$`2000` = occupancy.firsts$first.bound[which(occupancy.firsts$length == 2000)]

fname = "first_contact_time.txt";
write.table(final.firsts,file=paste(dirnew_data,"/", fname, sep = ""))

file = paste(dirnew_plots,"/first_contact_time_hist.png",sep="")
first.hist<-ggplot(occupancy.firsts2, aes(x=first.bound, fill=length)) + 
  geom_histogram(binwidth = 0.5, alpha = 0.5, position="identity") 
ggsave(file,plot=first.hist)

file = paste(dirnew_plots,"/first_contact_time_boxplot.png",sep="")
first.boxplot<- ggplot(occupancy.firsts2, aes(x=length, y=first.bound, fill=length)) + 
  geom_boxplot(outlier.colour ="red", position = position_dodge(1)) +
  stat_summary(fun = mean, geom = "point", shape = 8, size = 4)
ggsave(file,plot=first.boxplot)

final.firsts$`500` = occupancy.firsts$twoh.bound[which(occupancy.firsts$length == 500)]
final.firsts$`1000` = occupancy.firsts$twoh.bound[which(occupancy.firsts$length == 1000)]
final.firsts$`2000` = occupancy.firsts$twoh.bound[which(occupancy.firsts$length == 2000)]

fname = "200_contact_time.txt";
write.table(final.firsts,file=paste(dirnew_data,"/", fname, sep = ""))

file = paste(dirnew_plots,"/200_contact_time_hist.png",sep="")
first.hist<-ggplot(occupancy.firsts2, aes(x=twoh.bound, fill=length)) + 
  geom_histogram(binwidth = 0.5, alpha = 0.5, position="identity")
ggsave(file,plot=first.hist)

file = paste(dirnew_plots,"/200_contact_time_boxplot.png",sep="")
first.boxplot<- ggplot(occupancy.firsts2, aes(x=length, y=twoh.bound, fill=length)) + 
  geom_boxplot(outlier.colour ="red", position = position_dodge(1)) +
  stat_summary(fun = mean, geom = "point", shape = 8, size = 4)
ggsave(file,plot=first.boxplot)

final.firsts$`500` = occupancy.firsts$first.twoh.time.diff[which(occupancy.firsts$length == 500)]
final.firsts$`1000` = occupancy.firsts$first.twoh.time.diff[which(occupancy.firsts$length == 1000)]
final.firsts$`2000` = occupancy.firsts$first.twoh.time.diff[which(occupancy.firsts$length == 2000)]

fname = "first_200_contact_time_diff.txt";
write.table(final.firsts,file=paste(dirnew_data,"/", fname, sep = ""))

file = paste(dirnew_plots,"/1st_to_200_contact_timediff_hist.png",sep="")
first.hist<-ggplot(occupancy.firsts2, aes(x=first.twoh.time.diff, fill=length)) + 
  geom_histogram(binwidth = 0.5, alpha = 0.5, position="identity")
ggsave(file,plot=first.hist)

file = paste(dirnew_plots,"/1st_to_200_contact_timediff_boxplot.png",sep="")
first.boxplot<- ggplot(occupancy.firsts2, aes(x=length, y=first.twoh.time.diff, fill=length)) + 
  geom_boxplot(outlier.colour ="red", position = position_dodge(1)) +
  stat_summary(fun = mean, geom = "point", shape = 8, size = 4)
ggsave(file,plot=first.boxplot)


#### Zipping detection 
stats.zipping2 <- stats.zipping[-c(which(stats.zipping$first.zip== -1)),]

file = paste(dirnew_plots,"/first_zip_boxplot.png",sep="")
first.zip.boxplot <- ggplot(stats.zipping, aes(x=length, y=first.zip, color=length)) + 
  geom_boxplot(fill = "white", position = position_dodge(1), size = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 8, size = 3)
ggsave(file,plot=first.zip.boxplot)

file = paste(dirnew_plots,"/half_detection_boxplot.png",sep="")
first.zip.boxplot <- ggplot(stats.zipping, aes(x=length, y=half.detect, color=length)) + 
  geom_boxplot(fill = "white", position = position_dodge(1), size = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 8, size = 3)
ggsave(file,plot=first.zip.boxplot)



