require(data.table)
require(ggplot2)

#' This function calculates the biovolume based on the ellipsoid approximation 
#' of the volume of a particle. It uses the major and minor axis of the image as 
#' well as the pixel size. units are in [mm^3] 
#'
#'@param sample.info A data table containing samples in rows and details about 
#' the latter in columns.
#'@param object.info A data table containing objects in rows (not necessarily 
#' from the same sample) and details about the objects in the columns. 
#'@note This function does not return anything, it adds directly a column 
#'containing the volume of each object in the data table object.info
#'
#roxygen2::roxygenise()
ellipsoid.vol <- function(sample.info,object.info){
  for (id in sample.info$sample_id){ # for each sample
    pixel <- sample.info[sample_id==id,pixel_size]*10**-3 # we want it in mm vs in µm on ecotaxa
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
  sample.info[is.na(dilution_factor),dilution_factor:=1] # when no dilution factor is entered, it means that the sample has not been concentrated -> dilution of '1'
  sample.info[,norm_vol:=concentrated_sample_volume/(dilution_factor*imaged_volume*filtered_volume)]
  for (id in sample.info$sample_id){ # for each sample
    norm_vol_ <- sample.info[sample_id==id,norm_vol]
    object.info[sample_id==id,biovol:=volume*norm_vol_]
  }
}

#' This function computes the sum of biovolumes per category within each sample
#' @param object.info A data table containing objects in rows (not necessarily 
#' from the same sample) and details about the objects in the columns. 
#' @return the result is a data table with in the rows the different combinations 
#' of samples and categories, and containing an additional column: the summed 
#' biovolumes of each category within their sample.

summed.biovol <- function(object.info){ # object.info should have the column biovol
  if ("biovol" %in% colnames(object.info)){
    sample.names <- unique(object.info[,sample_id])
    if (length(sample.names)==1){ # if all objects are from the same sample 
      summed <- as.numeric(object.info[,sum(biovol),by=category])
    }
    else {
    summed <- object.info[,sum(biovol),by=list(sample_id,category)]
    }
    setnames(summed,c("sample_id","category","summed_biovol"))
    return(summed)
  }
  else{
    print("You must compute the biovolume first, for example using ellipsoid.biovol")
  } 
}

