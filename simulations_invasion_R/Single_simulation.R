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

#Import of the chr2.fa sequence file from the yeast genome (S288) :
yeast.genome.chr2 <- read.fasta("./yeast-genome/S288c-R64-2-1-v2014/chr2.fa" ,
                                seqtype = 'DNA', as.string = TRUE, 
                                forceDNAtolower  = TRUE, set.attributes = FALSE)

yeast.genome.chr2 <- yeast.genome.chr2[[1]] #select just the nucleotides sequence

# genome-wide microhomology counts
forward.sequences <- read.table("./Occurences_per_8bp_motif(for+rev_donor).txt", sep="\t", header = TRUE)
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


num.time.steps = 600 # Length of simulation in time steps
test.replicates = 5 # How many times to simulate, replicates
graph.resolution = 1 #save occupancy data at every nth time step. Plots will have this resolution at the x-axis 
kon.group<-c(0.005,0.05,0.1,0.4,0.7,0.9) #binding probabilities for every binding try
koff1.group<-c(0,0.0001,0.05,0.6) # dissociation probabilities for each bound particle
m.group = c(2,5) #bindings allowed to occur per tethering
search.window.group = c(250,500) #the genomic distance of the tethering effect (per side)


# Since the data needs to be outputted to files with human-readable names,we have to label the parameters with strings.
# For example 0005 is really 0.005
kon.group.names<- gsub("\\.", "", as.character(kon.group))
koff1.group.names<- gsub("\\.", "", as.character(koff1.group))

print(kon.group.names)
print(koff1.group.names)

#########################################################################################################
#################################### FUNCTIONS ##########################################################

find.occupancies = function(lower.window ="none", upper.window = "none", additional.removals = "none"){
  #' Find the positions where there is no MH bounded ;
  #' 
  # the indices of upper.window and lower.window are included in the removal
  # all function parameters are indices in 1:nchar(lys2.fragment)
  indices = 1:(nchar(lys2.fragment) -7)
  if (occupied.rad51$bound=="unbound"){
    return(indices)
  }
  #remove  occupied 8-bp sites along single end
  remove = c(occupied.rad51$lys2.microhomology)
  for (i in 1:7){
    remove = c(remove, (occupied.rad51$lys2.microhomology - i), 
               (occupied.rad51$lys2.microhomology + i))}
  if (lower.window != "none"){remove=c(remove, 0:lower.window)}
  if (upper.window != "none"){remove=c(remove, upper.window:nchar(lys2.fragment))}
  if (additional.removals[1] != "none"){  for (i in 0:7){
    remove = c(remove, (additional.removals - i), (additional.removals + i))}}
  remove = remove[which(remove > 0)]
  return(indices[-remove])
}
#########################################################################################################
#########################################################################################################


genome.wide.sei = function(initial.binding.tries){
  #choose region of LY weighted by available microhomologies
  #dont let those already occupied be chosen
  open.sites = find.occupancies() #indexes of unoccupied sites
  
  if (length(open.sites)== 0){ #if all the sites are occuped
    return(list(bound = occupied.rad51$bound,strand = "negative", donor.invasions = c(), lys2.microhomology = c()))
  }
  
  #matches : vector of possible bounding sites for MHs
  matches = c()
  for (i in 1:initial.binding.tries){
    if (length(open.sites) == 1){ 
      #if there is only one available site,the last ramaining site appends the vector matches
      matches[i] = open.sites
    }else{ 
      #We draw a site among those available which will be matched with a MH according to the respective probabilities of the MHs :
      matches[i] = sample(x=open.sites, size=1, prob = microhomology.probs[open.sites]) 
    }
    #Where there is a match with a MH, we consider that the site concerned is no longer available, we remove it from the index open.sites;
    #To be sure that the next matches will not overlap with the previous one, we have to remove from the index the 7 positions upstream and downstream of the match.
    open.sites = open.sites[-which(open.sites %in% (matches[i]-7):(matches[i] + 7))]
    if (length(open.sites) < 1){ 
      #if there is no more available sites where we could bound another MH
      break
    }
  }
  
  #Choose if/which binds;
  #successes : list of sample position among the matches where the probability of association is good enough to have a bond;
  successes=sample( c(TRUE, FALSE), length(matches), replace = TRUE, prob = c(kon.prob, (1-kon.prob)) )
  matches = matches[successes] #We keep only the matches where a bond will occur;
  
  if (length(matches)<1){
    #if no one match has a good enough probability Kon to be au bounding site ;
    return(list(bound = occupied.rad51$bound,strand = "negative",donor.invasions = c(), lys2.microhomology = c()))
  }
  
  # Set IDs of each bound (Heterology (H) vs LYS); 
  # If LYS, set as genomic position in + strand notation
  
  # id.probs : vector the matches which are among the self-MHs positions, divided by the total number of apparitions in the genome for the MH itself;
  # identities : list of id (H or LYS) according the above id.probs for each match ;
  id.probs = sapply(matches, function(x){(1 + ifelse(x %in% self.micros$position1, 1,0))/forward.sequences$total[x]})
  identities = sapply(1:length(matches), function(x) {sample(c("H", "LYS"), 1, prob = c(sapply(id.probs[x], function(x) max(0, 1-x)), id.probs[x]))})
  
  # donor.ids : vector of extrated SE that invaded lys2 (and not heterology (H) ) ;
  donor.ids = matches[which(identities == "LYS")]
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
    #remove MH ids that have already been counted as donor
    remove = which((identities !="H") & (as.character(identities) %in% occupied.rad51$donor.invasions) )
    if (length(remove)>0){
      identities = identities[-remove]; matches = matches[-remove]
    }
  } 
  #LYS alignment or misalignments
  #for LYS, match to correct by heterology within lys already allowed
  return(list(bound=occupied.rad51$bound, strand = "negative", donor.invasions = identities, lys2.microhomology = matches))
}  
#########################################################################################################
#########################################################################################################

new.microhomologizer = function(occupied.rad51, window, bindings.per.tethering){
  # correct.binding : vector of indexes for the micro-homologies (LYS) donors :
  correct.bindings = as.numeric(which(occupied.rad51$donor.invasions != "H"))
  # new.bindings : deep copy of an empty occupied.rad51 ;
  new.bindings = list(bound=occupied.rad51$bound, strand = "negative", donor.invasions = c(), lys2.microhomology = c())
  
  if (length(find.occupancies()) == 0){
    #if all the sites are occuped ;
    return(new.bindings)
  }
  
  # bindings : list of sites occupied by another MHs into a search window around some MHs locus ;
  bindings = c()
  
  for (binding.index in correct.bindings){
    # if all the sites are occuped once we add the binding has additional removals :
    if (length(bindings) > 0){
      if (length(find.occupancies(additional.removals = bindings)) == 0){break}
    }else{ 
      if (length(find.occupancies()) == 0){break} 
    }
    #current.selocus : index of a lys2.microhomology;
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
    
    
    # Here we want for each tethering, bound an open site loated in the search window of our current locus ;
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
      # We have to remove the candidate from the open sites and the 7 positions upstream and downstream it ;
      open.sites = open.sites[-which(open.sites %in% (current.bindings[j]-7):(current.bindings[j] + 7))]
      if (length(open.sites) <1){break}
    }
    bindings = c(bindings, current.bindings)
  } #end binding.index loop
  
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
    #remove MH ids in new.dindings that have already been counted as donor in occupied.rad51 ;
    remove = which((new.bindings$donor.invasions !="H") & (as.character(new.bindings$donor.invasions) %in% occupied.rad51$donor.invasions) )
    if (length(remove) > 0){
      new.bindings$donor.invasions = new.bindings$donor.invasions[-remove]
      new.bindings$lys2.microhomology = new.bindings$lys2.microhomology[-remove]}
  }
  return(new.bindings)
}

#########################################################################################################
#########################################################################################################

rev.comp<-function(x,rev=TRUE){ #Compute the reverse complement of a seq
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
count.consecutive.micros <- function(detect.rad54){
  # We have to know before zipping if the 200 MHs are consecutive or not ,
  # If not, it means that we have a "gap" between 2 homologies, and it's not a good thing if we want to zipp using rad54...
  
  microhomologies.left <- 0
  microhomologies.right <- 0
  
  if(detect.rad54 != 0){
    
    left <- 1
    while(lys2.occupancy$bound[detect.rad54 - left] == "yes" && lys2.occupancy$id[detect.rad54 - left] == "homology" && left < detect.rad54){
      microhomologies.left = microhomologies.left +1
      left = left+1
    }
    
    right <- 1
    while(lys2.occupancy$bound[detect.rad54 + right] == "yes" && lys2.occupancy$id[detect.rad54 + right] == "homology" && (detect.rad54 + right) < str_length(lys2.fragment)){
      microhomologies.right = microhomologies.right +1
      right = right+1
    }
  }
  return(c(microhomologies.left, microhomologies.right))
}

#########################################################################################################
#########################################################################################################

zipping <- function(detect.rad54, l, r){
  #overlapped.rad54 : all the others rad54 overlapped by the macrohomology
  overlapped.rad54 <- c(detect.rad54)
  for (pos in pos.rad54){
    if(pos %in% (detect.rad54 - l) : r){
      overlapped.rad54 = c(overlapped.rad54, pos)
    }
  }
  
  last.rad54 <- tail(sort(overlapped.rad54),1) #last rad54 position in the macrohomology
  if (length(which(pos.rdh54 < last.rad54) != 0)){
    nearby.rdh54 <- tail(sort(pos.rdh54[which(pos.rdh54 < last.rad54)]),1) #the nearest rdh54 from the last rad54 of the macrohomology
  }else{
    nearby.rdh54 <- 1
  }
  
  #When the macrohomology is ready to be zipped,
  #it is rad54 that will 'pump' each nucleotide by removing its associated rad51
  #in order to send it to invade the donor strand, until rad54 meets an rdh54
  
  zipped.indexes <- nearby.rdh54 : last.rad54 #start and end indexes of the zipping
  zipped.fragment <- str_sub(lys2.fragment, zipped.indexes[1], tail(zipped.indexes,1)) #the zipped fragment
  
  #We have to remove all the associated RAD51 from the zipped nucléotides :
  
  if(length(which(occupied.rad51$lys2.microhomology %in% zipped.indexes)) != 0){
    rad51.indexes2remove <- which(occupied.rad51$lys2.microhomology %in% zipped.indexes)
    occupied.rad51$donor.invasions <- occupied.rad51$donor.invasions[-rad51.indexes2remove]
    occupied.rad51$lys2.microhomology <- occupied.rad51$lys2.microhomology[-rad51.indexes2remove]
  }
  
  return(list(zipped.indexes, zipped.fragment))
}

#########################################################################################################
#########################################################################################################

template.copying <- function(zipped.indexes, zipped.fragment, start, end){
  ###Recombination with template start :
  
  revcomp.invading.fragment <- rev.comp(zipped.fragment)
  
  # If we recombine the left side of the rev.comp zipped fragment (and thus the right side of the initial fragment):
  if(str_locate_all(string = yeast.genome.chr2, pattern = rev.comp(lys2.fragment))[[1]][2] >= 473926){
    
    while(start > 469748){
      preserved2 = sample(c(TRUE, FALSE), size = nchar(revcomp.invading.fragment), replace = TRUE, prob = c(koff2,1-koff2))
      
      if (length(which(preserved2 == TRUE)) > floor(0.05*nchar(revcomp.invading.fragment))){ #koff2 risk of dissociation
        break
        
      }else{
        start = start-1
        new.nt <- str_sub(yeast.genome.chr2, start, start)
        revcomp.invading.fragment  = paste(new.nt, revcomp.invading.fragment, sep="")
      }
    }
    
    
    invading.fragment <- rev.comp(revcomp.invading.fragment)
    new.lys2.fragment <- paste(str_sub(lys2.fragment, start = 1, end = zipped.indexes[1]-1), invading.fragment, sep="") #the repaired DSB
    
    # If we recombine the right side of the rev.comp zipped fragment (and thus the left side of the initial fragment) :
  }else if (str_locate_all(string = yeast.genome.chr2, pattern = rev.comp(lys2.fragment))[[1]][1] <= 469748){
    while(end < 473926){
      preserved2 = sample(c(TRUE, FALSE), size = nchar(revcomp.invading.fragment), replace = TRUE, prob = c(koff2,1-koff2))
      
      if (length(which(preserved2 == TRUE)) > floor(0.05*nchar(revcomp.invading.fragment))){
        break
        
      }else{
        end = end +1
        new.nt <- str_sub(yeast.genome.chr2, end, end)
        revcomp.invading.fragment  = paste(revcomp.invading.fragment, new.nt, sep="")
      }
    }
    
    invading.fragment <- rev.comp(revcomp.invading.fragment)
    new.lys2.fragment <- paste(invading.fragment, str_sub(lys2.fragment, start = 1, end = zipped.indexes[1]-1), sep="")
  }
  
  return(new.lys2.fragment)
}

#########################################################################################################
#########################################################################################################


########################################################################################################
######################################### Single simulation ############################################


kon = 3; koff = 3; m = 2; sw = 2; 
kon.prob=kon.group[kon]
koff1.prob=koff1.group[koff]
bindings.per.tethering = m.group[m]
search.window = search.window.group[sw]

kon.name=kon.group.names[kon]
koff1.name=koff1.group.names[koff]

# Initialize the occupied.rad51 vector, genomic (start) position of RAD51 particles (bp / 8) of invaded strand ;
occupied.rad51 = list(bound = "unbound",strand = "negative", donor.invasions = 473927 - 368, lys2.microhomology = 368)

pop.time.series = as.data.frame(matrix(0,num.time.steps*3,3))
names(pop.time.series) = c("time.step","prob.detect","length")
pop.time.series$length = rep(ly.names, each = num.time.steps)
pop.time.series$time.step = rep(seq(1,num.time.steps,1),3)


invasion.stats <- as.data.frame(matrix(0,3*test.replicates,5))
names(invasion.stats) <- c('replicate.trial', 'length', 'time.step', 'invasion.trials', 'recombine.success')

saver = 0
bigtracker = 0
summary.stats = as.data.frame(matrix(0,3*test.replicates,3))
names(summary.stats) = c('length','first.time','twoh.time')

for (trial in 1:test.replicates){ 
  # initialize tabulation of bound microhomologies/heterologies
  
  if(saver < 3){
    # ly.binding.ts : over the course of an invasion round, 
    # it tabulates the number of bound microhomologies and heterologies for each of the L500, L and LY
    ly.binding.ts = as.data.frame(matrix(0, (num.time.steps/graph.resolution)*3,4))
    names(ly.binding.ts) = c('time.step', 'bound', "length", "heterologies")
    ly.binding.ts$time.step = rep(seq(1,num.time.steps, graph.resolution),3)
    ly.binding.ts$length = rep(ly.names, each = (num.time.steps / graph.resolution))
  }
  for (fragment in 1:3){
    # initialize the ID and state of the current invading strand 
    bigtracker=bigtracker+1
    ly.type = ly.names[fragment]
    lys2.fragment = ly.sequences[fragment]
    recombined.lys2.fragment <- ""
    if (fragment == 1){
      self.micros =L500.self.micros
    }else{
      if(fragment == 2){
        self.micros = L1000.selfmicros
      }else{
        self.micros = LY2000.selfmicros
      }
    }
    
    nb.rad54 <- floor(0.025*str_length(lys2.fragment)) #number of rad54 to be placed into the invading strand ;
    nb.rdh54 <- floor(0.25 * nb.rad54) # number of rdh54 to be placed into the invading strand;
    pos.rad54 <- c() #positions of rad54 in the invading strand;
    pos.rdh54 <- c() #positions of rdh54 in the invading strand;
    detect.rad54 <- 0 #the pos where a rad54 is overlapped by a rad51-MH complex
    
    start.dloop <- 0 #statement variable to engage a dloop invasion
    invasion.trials <- 1
    
    koff2 <- 0.00075 #probability for a SEI to be dissociated during the D-LOOP
    
    while (nb.rad54 > 0){ #place the requiered rad54 randomly (according uniform distro) over the invading fragment ;
      new.pos <- 0
      while (new.pos == 0 || new.pos %in% pos.rad54){
        new.pos <- floor(runif(1, min = 0, max=str_length(lys2.fragment)))
      }
      pos.rad54 = c(pos.rad54, new.pos)
      nb.rad54 = nb.rad54 - 1
    }
    
    while (nb.rdh54 > 0){ #Consider that a proportion of rad54 becomes rdh54
      new.pos <- sample(pos.rad54, size = 1)
      pos.rdh54 = c(pos.rdh54, new.pos) #pos.rad54 becomes pos.rdh54
      pos.rad54 = pos.rad54[-which(pos.rad54 == new.pos)] #remove the new pos.rdh54 from the pos.rad54 vector
      nb.rdh54 = nb.rdh54 - 1
    }
    
    # -7 because we made n-7 microhomolies of 8nt for a fragment of n nt ;
    # /8 because we just want 1 nt per MH of 8 nt ;
    SEI.binding.tries = floor((nchar(lys2.fragment)-7)/8)
    
    occupied.rad51$bound = "unbound"
    occupied.rad51$donor.invasions = c()
    occupied.rad51$lys2.microhomology = c()
    
    firster = 0 
    twoh = 0
    
    lys2.occupancy = as.data.frame(matrix(0, nchar(lys2.fragment),3));names(lys2.occupancy) = c('bp', 'bound', "id")
    lys2.occupancy$bp = 1:nchar(lys2.fragment)
    lys2.occupancy$bound = "no"
    lys2.occupancy$id = "heterology"
    # Loop through the time-steps
    
    for (time.step in 1:num.time.steps){
      
      if (occupied.rad51$bound != "unbound"){
        if (length(occupied.rad51$donor.invasions) != sum(occupied.rad51$donor.invasions == "H")){
          new.bindings = new.microhomologizer(occupied.rad51, search.window, bindings.per.tethering)
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
        occupied.rad51$donor.invasions = c(occupied.rad51$donor.invasions, new.bindings$donor.invasions)
        occupied.rad51$lys2.microhomology = c(occupied.rad51$lys2.microhomology, new.bindings$lys2.microhomology)
      }
      
      if(length(which(occupied.rad51$donor.invasions != "H")) != 0 && firster == 0){
        summary.stats$length[bigtracker] = ly.type
        summary.stats$first.time[bigtracker] = time.step
        firster = 1
      }
      
      if(occupied.rad51$bound != "unbound"){
        for (i in 1:length(occupied.rad51$lys2.microhomology)){
          lys2.occupancy$bound[occupied.rad51$lys2.microhomology[i]:(occupied.rad51$lys2.microhomology[i] + 7)]="yes"
          if(occupied.rad51$donor.invasions[i] != "H"){
            lys2.occupancy$id[occupied.rad51$lys2.microhomology[i]:(occupied.rad51$lys2.microhomology[i] + 7)]  = "homology"
          }
          
        }
        prob.detection = length(which(lys2.occupancy$id == "homology"))
        prob.detection = prob.detection/500 #take into a count crosslink density 1/500 ;
      }
      else{
        prob.detection = 0;
      }
      
      if(prob.detection >= 1){
        prob.detection = 1
      }
      if(length(which(lys2.occupancy$id == "homology")) >= 200 & twoh == 0){
        summary.stats$twoh.time[bigtracker] = time.step
        twoh = 1
      }
      if(saver < 3){
        # tabulate occupancies vs. time step and length
        ly.binding.ts$bound[ly.binding.ts$time.step == time.step & ly.binding.ts$length == ly.type] = length(which(lys2.occupancy$bound == "yes"))
        ly.binding.ts$heterologies[ly.binding.ts$time.step == time.step & ly.binding.ts$length == ly.type] = length(which(lys2.occupancy$bound == "yes" & lys2.occupancy$id == "heterology"))
      }
      pop.time.series$prob.detect[pop.time.series$time.step == time.step & pop.time.series$length == ly.type] = pop.time.series$prob.detect[pop.time.series$time.step == time.step & pop.time.series$length == ly.type] + prob.detection
      
      ### simulate random dissociation
      num.bound = length(occupied.rad51$donor.invasions)
      # print(num.bound)
      preserved = sample(c(FALSE,TRUE), num.bound, replace = TRUE, prob = c(koff1.prob,1-koff1.prob)) #dissociate if FALSE
      
      occupied.rad51$donor.invasions = occupied.rad51$donor.invasions[preserved]
      occupied.rad51$lys2.microhomology = occupied.rad51$lys2.microhomology[preserved]
      
      
      if (sum(!preserved)==num.bound){
        occupied.rad51$bound = "unbound"
      }
      
      ########  DLOOP invasion #########################
      
      # If a protein rad54 is overlaped by a micro-homology's donor AND we somewhere in the invading strand more than 200bp homologies :
      
      
      if (start.dloop == 0 && twoh ==1){
        consecutive.micros <-c()  #list of consecutive MHs around an overlapped rad54, when it occurs
        for (pos in pos.rad54){
          if (lys2.occupancy$bound[pos] == "yes" && lys2.occupancy$id[pos] == "homology"){
            consecutive.micros <- count.consecutive.micros(pos)
          }
        }
        
        #limit: consecutive bp necessary to enable zipping
        #As we know this cut-off is between 200 - 250 bp, we define it randomly using a draw from normal distribution that we will add to 225 bp
        #N(mean = 0, std = 20)
        
        
        limit <- 225 + rnorm(1, 0, 20)
        if (sum(consecutive.micros) + 1 > limit){
          start.dloop <- 1
          zip <- zipping(pos,l=consecutive.micros[1], r=consecutive.micros[2])
          # LY/L/L500 are in fact the reverse complements of the corresponding fragment in lys2 gene from the chr2 :
          
          if (length(zip[[1]]) > 16){
            #Get the first and the last position (on the genome) of the alignement between the rev-comp-zipped fragment
            #and the genome/chr the during the D-LOOP invasion step ;
            #rev.comp(zip[[2]]) : reverse complement of the zipped fragment
            donors.locations <- as.data.frame(str_locate_all(pattern = rev.comp(zip[[2]]), str = yeast.genome.chr2))
            first.match.invasion <- donors.locations$start[which(donors.locations$start %in% 469748:473927  & donors.locations$end %in% 469748:473927)]
            last.match.invasion <- donors.locations$end[which(donors.locations$start %in% 469748:473927  & donors.locations$end %in% 469748:473927)]
            #print(c("pos :", pos, "start :", first.match.invasion, "end : ", last.match.invasion))
            
            recombined.lys2.fragment <- template.copying(zipped.indexes = zip[[1]], zipped.fragment = zip[[2]], 
                                                         start = first.match.invasion, end = last.match.invasion)
            if(recombined.lys2.fragment != 0){
              #print(nchar(recombined.lys2.fragment))
              break #if the recombination successes, get out the time-set search homologies loop 
              
              
            }else{
              start.dloop <- 0
              invasion.trials = invasion.trials+1
            }
            
          }else{
            start.dloop <- 0
            invasion.trials = invasion.trials+1
          }

        }else{
          start.dloop <- 0
          invasion.trials = invasion.trials+1
        }
      }
      ###################################################
      
      print(c(ly.type, trial, time.step))
      
    }#next time step
    
    invasion.stats$length[bigtracker] = ly.type
    invasion.stats$replicate.trial[bigtracker] = trial
    invasion.stats$time.step[bigtracker] = time.step
    invasion.stats$invasion.trials[bigtracker] = invasion.trials
    invasion.stats$recombine.success[bigtracker] = ifelse(invasion.trials < num.time.steps, "yes", "no")
    
  }#next fragment
  
  if(saver < 3){
    saver = saver + 1;
  }
  
} #end process
