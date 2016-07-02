MDA <- function (path, title, tauGuess, YFlag){
    
  requireNamespace("graphics")
  
  tau5 <- 0
  
  if (YFlag == 1){ # to Choose Yasara output format
    
    print("YASARA MD Analysis")
    
    filePath = paste (path, title, sep="")        
    tab <- loadStats(filePath) # load ASCII file in table format
    
    # converts columns to vectors
    time <- as.numeric(paste(tab[2:nrow(tab)-1,1])) # time COL 1 - on standard md_analyze output
    rmsd <- as.numeric(paste(tab[2:nrow(tab)-1,11])) # backbone - rmsd COL 10 
    energy <- as.numeric(paste(tab[2:nrow(tab)-1,9])) # Total - energy COL 2 
    
  }
  else if (YFlag == 0) { # to Choose NAMD output format the flag have to be 0
    
    print("NAMD MD Analysis")
    
        
    filePath = paste (path, title, sep="")        
    tab <- loadStats(filePath) # carica il file ASCII in una tabella    
    
    # converts columns to vectors
    time <- as.numeric(paste(tab[,1])) # time
    rmsd <- as.numeric(paste(tab[,2])) # CA - rmsd
    energy <- as.numeric(paste(tab[,3])) # Total - energy  
      
  }
  
  D<-data.frame("time"=time,"rmsd"=rmsd,"energy"=energy) # 
    
  if (tauGuess == 0)
  {
    costantidbExp <- doubExponentialFitting(D, title) # starts fitting with a double exponential function
    
    if (is.na(costantidbExp)){ # if double exponential fitting is not doable then applies an exponential fitting
      costantisingExp <- exponentialFitting(D, title)
      tau5<- (1/costantisingExp[2]) * 5 # choose time constant
    }    
    else {
      tau5<- ceiling((1/costantidbExp[3]) * 5) # choose the slowest of the time constants and then multiplies to 5
    }  
    
    if (is.na(tau5) || tau5 > 10000) { # if guessed tau is null or too much high
      print ("############## Error - Tau Guess will be used!!! ")
      tau5<-tauGuess 
    } 
    
  } 
  else{
    print ("############## Tau Guess will be used!!! ")
    tau5<-tauGuess 
  }
    
  stats <- statistics(D, tau5)
    
  print (paste ("Analysis of Energy is done from (ps): ", tau5))   
  
  tmp <- paste(filePath,"_rmsd.svg", sep="")
  dev.copy(svg,tmp)
  dev.off()
  
  tmp <- paste(filePath,"_energy.svg", sep="")
  plot(time[time>tau5], energy[time>tau5], type="l", xlab="time (ps)", ylab="Energy (Kcal/mol)", main="Energy vs Time", sub=title)
  dev.copy(svg,tmp)
  dev.off()
  
  tmp <- paste(filePath,"_hist.svg", sep="")
  histogram <- hist(energy[time>tau5],freq=FALSE, xlab="Energy (Kcal/mol)", main="Energies Histogram",sub=title, density=5, nclass=15)
  dev.copy(svg,tmp)
  dev.off()  
    
  tmp <- paste(filePath,"_stats.csv", sep="")
  write.csv(cbind(stats, tau5), file = tmp)
  
  tmp <- paste(filePath,"_histo.csv", sep="")
  write.csv(rbind(histogram$breaks, histogram$counts, histogram$density, deparse.level=1), file = tmp)
  
  #RMSD histogram
  tmp <- paste(filePath,"_rmsd_hist.svg", sep="")
  rmsd_histogram <- hist(rmsd[time>tau5],freq=FALSE, xlab="RMSD (A)", main="RMSD Histogram",sub=title, density=5, nclass=15)
  dev.copy(svg,tmp)
  dev.off()  
  
  tmp <- paste(filePath,"_rmsd_hist.csv", sep="")
  write.csv(rbind(rmsd_histogram$breaks, rmsd_histogram$counts, rmsd_histogram$density, deparse.level=1), file = tmp)
  
    
  histogram
    
  return(stats)
  
}