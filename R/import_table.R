setwd("/Users/claratrellu/Documents/année sab/plankton planet/Rscripts/Ecotaxa_R")
require(data.table)

import.table <- function(file){

  data <- fread(file) # file is in tsv format
  # kind of contingency table, useful for vegan
  plot.data <- data[!grepl("^not-living|duplicate|multiple$|t001",object_annotation_hierarchy),
                    .N,
                    by=list(sample_id,object_annotation_category)]
  setnames(plot.data,c("sample","taxo","count"))
  
  # abundancies grouped to be used the vegan package
  for.veg <- plot.data
  for.veg<-for.veg %>% dcast(sample~taxo,fill=0,value.var = "count")

  # information unique to every object
  object.info <- data[!grepl("^not-living|duplicate|multiple$|t001",object_annotation_hierarchy),
                      .(object_id,sample_id,object_area,object_width,object_height,
                        object_annotation_category,object_annotation_hierarchy, 
                        object_equivalent_diameter,
                        object_major,object_minor, object_stdsaturation)]
  setnames(object.info,c("object_id","sample_id","area","width","height","category","hierarchy","ESD",
                       "major","minor","std_saturation"))
  
  # information in common to every object within a sample
  sample.info <- data[!grepl("^not-living|duplicate|multiple$|t001",object_annotation_hierarchy),
                      .(sample_id,sample_project,process_id,process_pixel, 
                        sample_concentrated_sample_volume,sample_dilution_factor,
                        sample_speed_through_water,sample_total_volume,
                        acq_id,acq_imaged_volume,acq_minimum_mesh,acq_maximum_mesh)]
  sample.info <- unique(sample.info)
  setnames(sample.info,c("sample_id","project","process_id","pixel_size",
                         "concentrated_sample_volume","dilution_factor","speed","filtered_volume",
                         "acq_id","imaged_volume","min_mesh","max_mesh"))
  
  result <- list("for.veg"=for.veg,"object.info"=object.info,"sample.info"=sample.info,"plot.data"=plot.data)
  return (result)
}
