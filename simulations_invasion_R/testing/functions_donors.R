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
      
      contact.freq = sequences.contacts.bins[matches[i], ]
      possible.bins = bins.id[which(contact.freq > 0)]
      bins[i] = sample(x=possible.bins, size=1, prob = contact.freq[contact.freq > 0])
      
    }else{ 
      # Draw a site among those available which will be matched with a MH according to the respective weighted probabilities :
      matches[i] = sample(x=open.sites, size=1, prob = microhomology.probs[open.sites])
      
      contact.freq = sequences.contacts.bins[matches[i], ]
      possible.bins = bins.id[which(contact.freq > 0)]
      bins[i] = sample(x=possible.bins, size=1, prob = contact.freq[contact.freq > 0])
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
    if(bins[b] %in% donors.list$bins){
      this.donor = donors.list$id[which(donors.list$bins == bins[b])]
      if(length(this.donor)>1){
        this.donor = sample(this.donor, size = 1)
      }
      identities = c(identities, this.donor)
    }else{
      identities = c(identities, "H")
    }
  }

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
    #remove = which((identities !="H") & (as.character(identities) %in% occupied.rad51$donor.invasions) )
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

new.microhomologizer = function(occupied.rad51, window, bindings.per.tethering){
  # When a binding is homologous, we have to search others MHs in a distance search windows ;
  # The number of bindings per search window is the bindings.per.tethering variable ;
  
  # correct.binding : vector of indexes for the micro-homologies (LYS) donors ;
  # new.bindings : deep copy of an empty occupied.rad51 ;
  correct.bindings = which(occupied.rad51$donor.invasions != "H" & occupied.rad51$donor.invasions %!in% donors.blacklist)
  new.bindings = list(bound=occupied.rad51$bound, strand = "negative", genome.bins = c(), donor.invasions = c(), lys2.microhomology = c())
  
  # Check for unbound sites :
  if (length(find.occupancies()) == 0){
    return(new.bindings)
  }
  
  # bindings : list of sites occupied by another MHs into the search window around the current micros locus ;
  bindings = c()
  bins = c()
  identities= c()
  
  for (binding.index in correct.bindings){
    if (length(bindings) > 0){
      if (length(find.occupancies(additional.removals = bindings)) == 0){break}
    }else{ 
      if (length(find.occupancies()) == 0){break} 
    }
    #current.selocus : index of the MH we are currently looking around it (search window) to place another MHs;
    current.selocus = occupied.rad51$lys2.microhomology[binding.index]
    current.bin = occupied.rad51$genome.bins[binding.index]
    this.donor = occupied.rad51$donor.invasions[binding.index]
    
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
        bins = c(bins, current.bin)
        identities = c(identities, this.donor)
        
      }else{
        candidate = sample(open.sites, size = 1)
        yy = runif(1)
        if(yy <= kon){
          current.bindings = c(current.bindings,candidate)
          bins = c(bins, current.bin)
          identities = c(identities, this.donor)
        }
      }
      # Remove the candidates and the 7 positions upstream and downstream it from the free sites index ;
      open.sites = open.sites[-which(open.sites %in% (current.bindings[j]-7):(current.bindings[j] + 7))]
      if (length(open.sites) <1){break}
    }
    bindings = c(bindings, current.bindings)
  }
  
  # donor.ids : same as in the genome.wide.sei function  ;
  donor.ids = bindings[which(identities == "LYS")]
  
  for (index in which(donor.ids %in% self.micros$position1)){
    y = which(self.micros$position1 == donor.ids[index])[1]
    sampling.micros = c(self.micros$position1[y], self.micros$position2[y])
    if(!is.na(self.micros$position3[y])){
      sampling.micros =   c(sampling.micros,self.micros$position3[y])
    }
    donor.ids[index] = sample(as.numeric(sampling.micros), size = 1)
  }
  
  identities[which(identities == "LYS")] = 473927 - donor.ids
  new.bindings$genome.bins = c(new.bindings$genome.bins, bins)
  new.bindings$lys2.microhomology = c(new.bindings$lys2.microhomology, bindings)
  new.bindings$donor.invasions    = c(new.bindings$donor.invasions, identities)
  
  if (occupied.rad51$bound != "unbound"){
    #remove MHs ids in new.bindings that have already been counted as donor in occupied.rad51 :
    #remove = which((new.bindings$donor.invasions !="H") & (as.character(new.bindings$donor.invasions) %in% occupied.rad51$donor.invasions) )
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

donors.generator <- function(template, bins, N = 10){
  new.donors.list<-list(sequence = c(LY), bins = c("chr2_470001_480001"), id = c("LYS"), invasion = c("no", rep("no", times = N)))
  bases <- c("a", "t", "g", "c")
  
  if(N >= 1){
    for (n in 1:N){
      new.donor <- template
      new.donor.id = paste("!LYS", as.character(n), sep="")
      lower.limit <- floor(0.05*nchar(template))
      upper.limit <- floor(0.3*nchar(template))
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
zipping <- function(rad54, zipping.list, donor, limit = 10){
  
  pos <- rad54
  zip.indexe <- c()
  zip.fragment <-"" 
  counter <- 0
  donor.seq = donors.list$sequence[which(donors.list$id == donor)]
  
  while(pos %!in% pos.rdh54 && 
        pos %!in% pos.rad54[which(pos.rad54 != rad54)] && 
        pos < nchar(lys2.fragment)){
    
    if(donors.occupancy$bound.id[pos] == "homology" && donors.occupancy$donor.id[pos] == donor){
      new.nt <- substr(lys2.fragment, pos, pos)
      zip.indexe = c(zip.indexe, pos)
      zip.fragment = paste(zip.fragment, new.nt, sep="")
      pos = pos + 1
      
    }else{
      if (str_sub(string = lys2.fragment, start = pos, end = pos) ==  str_sub(string = donor.seq, start = pos, end = pos)){
        new.nt <- substr(lys2.fragment, pos, pos)
        zip.indexe = c(zip.indexe, pos)
        zip.fragment = paste(zip.fragment, new.nt, sep="")
        pos = pos + 1
        
      }else{
        counter = counter + 1
        if (counter > limit){
          return(-1) #wrong donor
          
        }else if (counter <= limit){
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
