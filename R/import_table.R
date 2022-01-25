require(data.table)

#' This function is mandatory to use the functions provided in this package. 
#' It will allow you to process the original tsv file and copy the data in 
#' different 'data.table' to be then analysed by specific functions. Within this 
#' function, the functions from biovolume and relative abundance are used to add 
#' these informations to the tables
#' 
#' @param file a .tsv file from a Planktoscope aquisition
#' @param volumes If true, the functions ellipsoid.vol and biovolume will be 
#' called to add two additional columns to the returned data table object.info
#' containing the volume and biovolume of each object
#' @return the result is a list of four data tables. $for.veg to be used by the 
#' vegan package, sample in rows and categories in columns, the entries being 
#' the count of objects being in the corresponding sample and category. $object.info 
#' contains all objects and their specifications (eg size, color...). $sample.info
#' contains all samples and their specifications (eg speed, volume filtered, 
#' concentration factor). $counts is the summary of the absolute and relative 
#' abundancies of categories in each sample.

import.table <- function(file,volumes=TRUE,unwanted=c("not-living")){

  data <- fread(file) # file is in tsv format
  
  # get rid of the rows containing unwanted objects
  data<-data[!grepl(paste(unwanted,collapse="|"),object_annotation_hierarchy),]
  
  # summary of absolute abundances
  counts <- data[,.N, by=list(sample_id,object_annotation_category)]
  setnames(counts,c("sample_id","category","count"))
  
  # calculation of relative abundance and composition
  relative.abundance(counts)
  # abundancies grouped to be used by the vegan package, contingency table
  for.veg <- counts
  for.veg<-for.veg %>% dcast(sample_id~category,fill=0,value.var = "count")

  # information unique to every object
  object.info <- data[, .(object_id,sample_id,object_area,object_width,object_height,
                        object_annotation_category,object_annotation_hierarchy, 
                        object_equivalent_diameter,
                        object_major,object_minor, object_stdsaturation)]
  setnames(object.info,c("object_id","sample_id","area","width","height","category","hierarchy","ESD",
                       "major","minor","std_saturation"))
  
  # information in common to every object within a sample
  sample.info <- data[, .(sample_id,sample_project,process_id,process_pixel, 
                        sample_concentrated_sample_volume,sample_dilution_factor,
                        sample_speed_through_water,sample_total_volume,
                        acq_id,acq_imaged_volume,acq_minimum_mesh,acq_maximum_mesh)]
  sample.info <- unique(sample.info)
  setnames(sample.info,c("sample_id","project","process_id","pixel_size",
                         "concentrated_sample_volume","dilution_factor","speed","filtered_volume",
                         "acq_id","imaged_volume","min_mesh","max_mesh"))
  
  if (volumes) {
    ellipsoid.vol(sample.info,object.info)
    biovolume(sample.info,object.info)
  }
  
  result <- list("for.veg"=for.veg,"object.info"=object.info,"sample.info"=sample.info,"counts"=counts)
  return (result)
}
