doubExponentialFitting <- function (D, title){
  
  time=D[,1]
  rmsd=D[,2]
  
  plot(time, rmsd,type="l",main="RMSD and Doub Exp Fitting",sub=title,xlab="Time (ps)", ylab="RMSD (A)")
  
  out <- tryCatch(
    {
      print("Try to apply a Double Exponential Fitting...")
      nlPrediction <- nls(rmsd ~  Const - exp(-A1 * time) - exp(-A2 * time),data=D, start=c("Const"=2,"A1"=0.12,"A2"=0.012),lower=c(0.1,0.004,0.003), trace=TRUE,algorithm="port")
      print (summary(nlPrediction))
      lines(time, predict(nlPrediction), col = 2)
      return(summary(nlPrediction)$coefficients[,1])
    },
    error=function(cond){
      print("ERROR, Double Exponential Fitting cannot be applied!!!")
      return(NA)
    },
    warning=function(cond){
      print("WARNING, probably it is better you use Exponentinal Fitting")
      print (summary(nlPrediction))
      lines(time, predict(nlPrediction), col = 2)
      return(summary(nlPrediction)$coefficients[,1])
    }
  ) #end tryCatch
  return(out)
  
}