exponentialFitting <- function (D, title){
  
  time=D[,1]
  rmsd=D[,2]
  
  plot(time, rmsd,type="l", main="RMSD and Exp Fitting", sub=title,xlab="Time (ps)", ylab="RMSD (A)")
  
    out <- tryCatch(
    {
      print("Try to apply an Exponential Fitting...")
      nlPrediction <- nls(rmsd ~  Const - exp(-A * time),data=D, start=c("Const"=2,"A"=0.01))
      print (summary(nlPrediction))
      lines(time, predict(nlPrediction), col = 2)
      return(summary(nlPrediction)$coefficients[,1])
    },
    error=function(cond){
      print("ERROR, Exponential Fitting cannot be applied!!!")
      return(NA)
    },
    warning=function(cond){
      print("WARNING, probably it is better you use a NEW Fitting Function")
      print (summary(nlPrediction))
      lines(time, predict(nlPrediction), col = 2)
      return(summary(nlPrediction)$coefficients[,1])
    }
  ) #end tryCatch
  return(out)
    
}