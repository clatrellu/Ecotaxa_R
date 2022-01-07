setwd("/Users/claratrellu/Documents/année sab/plankton planet/Rscripts")
require(data.table)
require(ggplot2)



# object.info should have the column biovol
# NBSS is calculated for one sample only, calculate several NBSS spectra to compare samples

NBSS <- function(object.info,sample.info,sample){ 
  
  data <- sample.info[sample_id==sample,]
  objects <- object.info[sample_id==sample,]
  print(data)
  
  Bvmin <- 4/3*pi*((data$min_mesh/2)*10**-3)**3
  intervals <- c(Bvmin)
  top <- (data$max_mesh*10**-3/2)**3*4/3*pi
  while (Bvmin < top){
    add <- 2**0.25*Bvmin # Bvmax = Bvmin*2^0.25
    intervals <- append(intervals,add)
    Bvmin <- add
  }
  print(intervals)
  y <- list()
  x <- list()
  for (i in 1:(length(intervals)-1)){
    a <- intervals[i]
    b <- intervals[i+1]
    Bvtot <- b-a
    add_ <- objects[(biovol>a)&(biovol<b),sum(biovol)]
    add_ <- add_/Bvtot
    y <- append(y,add_)
    x <- append(x,b)
    
  }
  nbss.plot <- data.table(Spectra=x,NBSS=y)
  
  return(nbss.plot)

}
