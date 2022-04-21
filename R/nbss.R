require(data.table)
require(ggplot2)


#' NBSS : Normalized Biovolume Size Spectra calculation for a given sample, 
#' calculate several NBSS spectra to compare samples
#'
#' @param objects A data table containing objects in rows (not necessarily 
#' from the same sample) and details about the objects in the columns. 
#' @param samples A data table containing samples in rows and details about 
#' the latter in columns.
#' @param sample a string containing the name of the sample for which to calculate the NBSS
#' @return a data table containing two columns, the first with the size spectra, 
#' given in mm^3, the second is the NBSS for each interval of size.
#' @note The variables used are called biovolumes but to simplify we actually use 
#' volumes. Since biovolumes are the volumes normalized by the same total volume 
#' for all objects in a sample, the result is the same, except that size spectra 
#' are given in mm^3 instead of mm^3/m^3

NBSS <- function(planktotable,N=TRUE){ 
  
  if (nrow(planktotable)<1000){
    print("Warning : this sample contains less than 1000 object, therefore the NBSS plot might not be representative enough for such small datasets")
  }
  
  Bvmin <- 4/3*pi*((unique(planktotable$acq_minimum_mesh)/2)*10**-3)**3 # lower bound
  intervals <- c(Bvmin)
  top <- (unique(planktotable$acq_maximum_mesh)*10**-3/2)**3*4/3*pi # upper bound
  while (Bvmin < top){
    add <- 2**0.25*Bvmin # Bvmax = Bvmin*2^0.25
    intervals <- c(intervals,add) # interval serving for x axis of nbss plot
    Bvmin <- add
  }
  y <- list()
  x <- list()
  for (i in 1:(length(intervals)-1)){
    a <- intervals[i]
    b <- intervals[i+1]
    Bvtot <- b-a
    add_ <- planktotable[(volume>a)&(volume<b),sum(biovolume)] # sum biovolumes of a certain range of size
    if (N){
      add_ <- add_/Bvtot
    }
    y <- append(y,add_)
    x <- append(x,b)
    
  }
  nbss.plot <- data.table(Spectra=as.numeric(x),NBSS=as.numeric(y))
  nbss.plot <- nbss.plot[which(NBSS!=0),] # remove size intervals where there are no objects + log transformation
  nbss.plot[,NBSS:=log10(NBSS)]
  nbss.plot[]
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
#' p <- NBSS.plot(objects=objects,samples=samples,sample_name=sample_id,ESD=TRUE)
#' p 
#' 
NBSS.plot <- function(planktotable,ESD=TRUE){
  
  sample_name <- planktotable[,unique(sample_id)]
  data <- NBSS(planktotable,N=TRUE)
  
  #convert back to ESD :
  
  data[,Spectra:=2*10**3*(Spectra*3/(4*pi))**(1/3)]
  
  min_mesh <- unique(planktotable$acq_minimum_mesh)
  max_mesh <- unique(planktotable$acq_maximum_mesh)
  
  p <- ggplot(data,aes(x=Spectra,y=NBSS)) + 
    geom_point()  + scale_x_log10(limits=c(min_mesh,max_mesh))+
    labs(x="Equivalent Spherical Diameter [µm]",y="NBSS [mm^3/mm^3/m^3]",title = paste("NBSS:",sample_name))
  
  return(p)
}

BSS.plot <- function(planktotable,ESD=TRUE){
  
  sample_name <- planktotable[,unique(sample_id)]
  data <- NBSS(planktotable,N=FALSE)
  
  #convert back to ESD :
  
  data[,Spectra:=2*10**3*(Spectra*3/(4*pi))**(1/3)]
  
  min_mesh <- unique(planktotable$acq_minimum_mesh)
  max_mesh <- unique(planktotable$acq_maximum_mesh)
  
  p <- ggplot(data,aes(x=Spectra,y=NBSS)) + 
    geom_point()  + scale_x_log10(limits=c(min_mesh,max_mesh))+
    labs(x="Equivalent Spherical Diameter [µm]",y="BSS [mm^3/m^3]",title = paste("BSS:",sample_name))
  
  return(p)
}
  