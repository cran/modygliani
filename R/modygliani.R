#' modygliani: MOlecular DYnamics GLobal ANalysis
#' 
#' @description Modygliani performs a fitting of the RMSD trajectories. Internal energies are analyzed on the basis of the RMSD fitting. This process is iterated as many times as the number of files provided. Modygliani allows for the comparison of different trajectories. See the references for further information. Example input files are in 'inst/extdata/'.Example output files are stored in 'inst/extdata/output'.
#' @param path A string containing the path to the properly formatted ASCII files. NAMD file style have to be composed of 3 columns including time, RMSD and Energy. Yasara styles are as the standard output.
#' @param names A vector of strings containing the names of the peptides/proteins to compare. The number of files to analyze have to correspond to the length of this vector and vice versa. File names cannot contain an extension.
#' @param tauGuess A numerical value containing the time constant that is guessed by the user. tayGuess is a useful parameter when fitting fails. If tauGuess is to 0, Modygliani guesses tau in order to attempt fittings.
#' @param miniRMSD A numerical value corresponding to the minimum value of RMSD
#' @param maxiRMSD A numerical value corresponding to the maximum value of RMSD
#' @param miniEne A numerical value corresponding to the minimum value of Energy
#' @param maxiEne A numerical value corresponding to the maximum value of Energy
#' @param YFlag A 1|0 numerical value. The values are 1 for Yasara dynamics and 0 for NAMD
#' @param title A string containing the title of the charts
#' @param colors A vector of strings containing color shades, see the example for further details
#' @author Luca Belmonte, Sheref S. Mansy
#' @references Belmonte, L. Rossetto, D. Forlin, M. Scintilla, S. Bonfio, C. Mansy, S. S. "Cysteine containing dipeptides show a metal specificity that matches the composition of seawater" Phys. Chem. Chem. Phys., 2016, DOI: 10.1039/C6CP00608F
#' @examples
#' ## fitting and comparison of four different MD trajectories
#' names <- c("fe_CG","co_CG","ni_CG","zn_CG") # file names
#' colors<- c("gray0","gray25","gray50","gray75")
#' path <-"inst/extdata/"
#' modygliani(path, names, 0, 0, 8, -200, 250, 0, "Modygliani example on NAMD trajectories", colors)
#' @export
#' @importFrom grDevices dev.copy dev.off svg
#' @importFrom graphics arrows barplot hist lines par plot
#' @importFrom stats nls predict sd
#' @importFrom utils read.csv write.csv
#' @importFrom graphics par plot rect text
modygliani <- function(path, names, tauGuess, miniRMSD, maxiRMSD, miniEne, maxiEne, YFlag, title, colors){
  
  requireNamespace("graphics")
      
  print("#######################################################")
  print("")
  print("MoDyGliAni: Molecular Dynamics Global Analysis")
  print("")
  print("If you are using Modygliani for your own research please Cite: ")
  print("Belmonte L, Rossetto D, Forlin M, Scintilla S, Bonfio C and Mansy SS ")
  print("'Cysteine containing dipeptides show metal specificity that matches the composition of seawater'")
  print("Physical Chemistry Chemical Physics, 2016")
  print("#######################################################")
  print("")  
  
  writeOnFS <- 0
  
  writeOnFS <- as.numeric(readline(paste(paste("modygliani needs to write results on", path,sep=" "), "Do you agree? [1=Yes, 0=No] ", sep=" ")))
  
  if((!is.na(writeOnFS==1))&&(writeOnFS==1)){
    i <- 1
    
    lunghezzaCiclo <- length(names)
    
    Energy_AVG<-NULL
    Energy_SD<-NULL
    
    RMSD_AVG<-NULL
    RMSD_SD<-NULL
    
    for (i in 1:lunghezzaCiclo) {
      
      stats <- MDA(path, names[i], tauGuess, YFlag)
      
      RMSD_AVG<-cbind(RMSD_AVG,stats[,1])
      RMSD_SD<-cbind(RMSD_SD,stats[,2])
      
      Energy_AVG<-cbind(Energy_AVG,stats[,5])
      Energy_SD<-cbind(Energy_SD,stats[,6])    
      
    }
    colnames(Energy_AVG)<-names
    
    old.par <- par(mfrow=c(2, 1))
    
    barRMSD<- barplot(RMSD_AVG, ylim=c(miniRMSD, maxiRMSD), col=colors, beside=TRUE, ylab="Average RMSD (Angstrom)", main=title, las=1)
    #legend("top",legend=c("Fe","Co","Ni","Cu","Zn"), fill=colors, bty = "n", ncol=5)
    
    arrows(barRMSD,RMSD_AVG+RMSD_SD,barRMSD,RMSD_AVG, angle=90, code=1)
    arrows(barRMSD,RMSD_AVG-RMSD_SD,barRMSD,RMSD_AVG, angle=90, code=1)
    
    barEne<- barplot(Energy_AVG, ylim=c(miniEne, maxiEne), col=colors, beside=TRUE, xlab="Chains", ylab="Average Energy (kcal/mol)", las=3)
    
    arrows(barEne,Energy_AVG+Energy_SD,barEne,Energy_AVG, angle=90, code=1)
    arrows(barEne,Energy_AVG-Energy_SD,barEne,Energy_AVG, angle=90, code=1)
    
    par(old.par)
    
    dev.copy(svg,paste(path,"aggregates.svg"))
    dev.off()
    
    
    print("#######################################################")
    print("Molecular Dynamics Analysis - FINISHED")
    print( paste ("Check your output in ", path, "*.dat & *.png"))
    print("#######################################################")
    
    summary<-rbind(Energy_AVG, Energy_SD, RMSD_AVG, RMSD_SD, deparse.level = 1)
    orderedSummary<- summary[, order(summary[1,], decreasing=FALSE)]
    rownames(orderedSummary) <- c("Energy_avg", "Energy_rmsd", "RMSD_avg", "RMSD_rmsd")
    
    tmp <- paste(path,"summary.csv", sep="")
    write.csv(orderedSummary, file = tmp)  
  }
  else {
    print("#######################################################")
    print("Molecular Dynamics Analysis - ABORTED!!!")
    print("#######################################################")
  }  
    
}
