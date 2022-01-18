setwd("/Users/claratrellu/Documents/année sab/plankton planet/Rscripts/Ecotaxa_R")
require(data.table)
require(ggplot2)


#' NBSS : Normalized Biovolume Size Spectra calculation for a given sample, 
#' calculate several NBSS spectra to compare samples
#'
#' @param object.info A data table containing objects in rows (not necessarily 
#' from the same sample) and details about the objects in the columns. 
#' @param sample.info A data table containing samples in rows and details about 
#' the latter in columns.
#' @param sample a string containing the name of the sample for which to calculate the NBSS
#' @return a data table containing two columns, the first with the size spectra, 
#' given in mm^3, the second is the NBSS for each interval of size.
#' @note The variables used are called biovolumes but to simplify we actually use 
#' volumes. Since biovolumes are the volumes normalized by the same total volume 
#' for all objects in a sample, the result is the same, except that size spectra 
#' are given in mm^3 instead of mm^3/m^3

NBSS <- function(object.info,sample.info,sample){ 
  
  data <- sample.info[sample_id==sample,] # extract the row of interest
  objects <- object.info[sample_id==sample,] # take only the objects of this sample
  
  Bvmin <- 4/3*pi*((data$min_mesh/2)*10**-3)**3 # lower bound
  intervals <- c(Bvmin)
  top <- (data$max_mesh*10**-3/2)**3*4/3*pi # upper bound
  while (Bvmin < top){
    add <- 2**0.25*Bvmin # Bvmax = Bvmin*2^0.25
    intervals <- append(intervals,add) # interval serving for x axis of nbss plot
    Bvmin <- add
  }
  y <- list()
  x <- list()
  for (i in 1:(length(intervals)-1)){
    a <- intervals[i]
    b <- intervals[i+1]
    Bvtot <- b-a
    add_ <- objects[(volume>a)&(volume<b),sum(biovol)] # sum biovolumes of a certain range of size
    add_ <- add_/Bvtot
    y <- append(y,add_)
    x <- append(x,b)
    
  }
  nbss.plot <- data.table(Spectra=as.numeric(x),NBSS=as.numeric(y))
  return(nbss.plot)
}

#' Plotting of the NBSS - Normalized Biovolume Size Spectra - of a given sample
#' 
#' @param objects A data table containing objects as rows and details on these objects in the 
#' columns. This argument will be passed in the NBSS function
#' @param samples A data table containing samples in rows and details about 
#' the latter in columns.
#' @param sample_name A string containing the name of the sample from which to extract objects 
#' and do an NBSS analysis.
#' @param ESD If ESD is true, by default, the x axis is expressed as ESD (Eqivalent spherical 
#' diameter for a given volume) in [µm]. Otherwise, volumes will be on the x-axis in [mm^3]
#' @return an object of the type ggplot which can then be plotted.
#' 
#' @examples 
#' p <- NBSS.plot(objects=object.info,samples=sample.info,sample_name=sample_id,ESD=TRUE)
#' p 
#' 
NBSS.plot <- function(objects,samples,sample_name,ESD=TRUE){
  
  data <- NBSS(objects,samples,sample_name)
  #convert back to ESD :
  data[,Spectra:=2*10**3*(Spectra*3/(4*pi))**(1/3)] 
  p <- ggplot(data,aes(x=Spectra,y=NBSS)) + 
    geom_point() +
    scale_x_log10() +
    labs(x="Equivalent Spherical Diameter",y="NBSS [mm^3/mm^3/m^3]",title = paste("Normalized Biovolume Size Spectra for the sample",sample_name))
  
  return(p)
}
