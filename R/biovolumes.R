require(data.table)
require(ggplot2)

#' This function calculates the volume of a particle based on an ellipsoid 
#' approximation. It uses the major and minor axis of the image as well as the 
#' pixel size. Units are in [mm^3] 
#'
#'@param sample.info A data table containing samples in rows and details about 
#' the latter in columns.
#'@param object.info A data table containing objects in rows (not necessarily 
#' from the same sample) and details about the objects in the columns. 
#'@note This function does not return anything, it adds directly a column 
#'containing the volume of each object in the data table object.info


ellipsoid.vol <- function(sample.info,object.info){
  for (id in sample.info$sample_id){
    pixel <- sample.info[sample_id==id,pixel_size]*10**-3 # transform to mm vs in µm in ecotaxa
    object.info[sample_id==id,':='(major=major*pixel,minor=minor*pixel)] #cf piQv tutorial
    object.info[sample_id==id,volume:=4/3*pi*major/2*(minor/2)**2] # mm^-3
  }
}

#' This function calculates the biovolume of a particle. This is the percentage 
#' volume a particle [mm^3] takes up in a unit of filtered water[m^3].
#'@param sample.info A data table containing samples in rows and details about 
#' the latter in columns.
#'@param object.info A data table containing objects in rows (not necessarily 
#' from the same sample) and details about the objects in the columns. 
#'@note This function does not return anything, it adds directly a column 
#'containing the biovolume of each object in the data table object.info
#'@note It is relevant to talk about biovolume when looking about at a certain 
#'catgeory of plankton, when we sum individual biovolumes to show what part
#'certain species take up in a unit of water.
#'
biovolume <- function(sample.info,object.info) {
  sample.info[is.na(dilution_factor),dilution_factor:=1] # default dilution factor is 1
  sample.info[,norm_vol:=concentrated_sample_volume/
                (dilution_factor*imaged_volume*filtered_volume)]
  for (id in sample.info$sample_id){
    norm_vol_ <- sample.info[sample_id==id,norm_vol]
    object.info[sample_id==id,biovol:=volume*norm_vol_]
  }
}

#' This function computes the sum of biovolumes per category within each sample
#' @param object.info A data table containing objects in rows (not necessarily 
#' from the same sample) and details about the objects in the columns. 
#' @return the result is a data table with : in the rows the different combinations 
#' of samples and categories, and containing an additional column: the summed 
#' biovolumes of each category within their sample.

summed.biovol <- function(object.info){
  if ("biovol" %in% colnames(object.info)){
    sample.names <- unique(object.info[,sample_id])
    if (length(sample.names)==1){ # if all objects are from the same sample 
      summed <- object.info[,sum(biovol),by=category]
      setnames(summed,c("sample_id","summed_biovol"))
    }
    else {
      summed <- object.info[,sum(biovol),by=list(sample_id,category)]
      setnames(summed,c("sample_id","category","summed_biovol"))
    }
    replace_na(summed,list(rep(0,ncol(summed))))
    return(summed)
  }
  else{
    print("You need to calculate the biovolumes first")
  } 
}
