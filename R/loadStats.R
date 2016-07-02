loadStats <- function (filePath){
  
  print (paste ("Reading analysis output file", filePath))
    
  table <- read.csv(filePath, header = TRUE, sep = "\t")
  
  table
  
}