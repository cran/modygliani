statistics <- function(D, limit){
  
  time = D[,1]
  RMSD = D[,2]
  Energy = D[,3]
  
  RMSD_avg = mean(RMSD, na.rm = TRUE)
  RMSD_sdv = sd(RMSD, na.rm = TRUE)
  RMSD_min = min(RMSD, na.rm = TRUE)
  RMSD_max = max(RMSD, na.rm = T)
  
  Ene_avg = mean(Energy[time > limit], na.rm = TRUE)
  Ene_sdv = sd(Energy[time > limit], na.rm = T)
  Ene_min = min(Energy[time > limit], na.rm = T)
  Ene_max = max(Energy[time > limit], na.rm = T)
  
  aggregates = data.frame("RMSD_avg"=RMSD_avg, "RMSD_sdv"=RMSD_sdv, "RMSD_min"=RMSD_min, "RMSD_max"=RMSD_max, "Ene_avg"=Ene_avg, "Ene_sdv"=Ene_sdv, "Ene_min"=Ene_min, "Ene_max"=Ene_max)
  
  return(aggregates)
  
}