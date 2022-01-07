setwd("/Users/claratrellu/Documents/année sab/plankton planet/Rscripts")
require(data.table)
require(ggplot2)


ellipsoid.biovol <- function(sample.info,object.info){
  
  for (id in sample.info$sample_id){ # for each sample
    pixel <- sample.info[sample_id==id,pixel_size]*10**-3 # we want it in mm vs in µm on ecotaxa
    object.info[sample_id==id,':='(major=major*pixel,minor=minor*pixel)] #cf piQv tutorial
    object.info[sample_id==id,biovol:=4/3*pi*major/2*(minor/2)**2] # mm^-3
  }
#return(object.info) # par référence ou pas? 
}

norm_biovol <- function(sample.info,object.info) {
  sample.info[is.na(dilution_factor),dilution_factor:=1] # when no dilution factor is entered, it means that the sample has not been concentrated -> dilution of '1'
  sample.info[,norm_vol:=concentrated_sample_volume/(dilution_factor*imaged_volume*filtered_volume)]
  for (id in sample.info$sample_id){ # for each sample
    norm_vol_ <- sample.info[sample_id==id,norm_vol]
    object.info[sample_id==id,n_biovol:=biovol*norm_vol_]
  }
}

# returns the sum of the biovolumes per sample

summed.biovol <- function(object.info_){ # object.info should have the column biovol
  if ("biovol" %in% colnames(object.info_)){
    sample.names <- unique(object.info_[,sample_id])
    if (length(sample.names)==1){ # if all objects are from the same sample 
      summed <- as.numeric(object.info_[,sum(biovol)])
    }
    else {
    summed <- object.info_[,sum(biovol),by=sample_id]
    setnames(summed,c("sample_id","summed_biovol"))
    }
  }
  else{
    print("You must compute the biovolume first, for example using ellipsoid.biovol")
  } 
  return(summed)
}

