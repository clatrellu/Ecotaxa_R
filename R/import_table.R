require(data.table)
require(tidyr)

#' Import data from tsv file
#' @description This function is mandatory to use the functions provided in this package. 
#' It will allow you to process the original tsv file and copy the data in 
#' different 'data.table' to be then analysed by specific functions. Within this 
#' function, the functions from biovolume and relative abundance are used to add 
#' these parameters in new columns to the tables
#' @param file one or several .tsv file(s) from a Planktoscope aquisition classified 
#' with Ecotaxa. If multiple files are inserted, make sure they are grouped in a list.
#' @param volumes If true, the functions ellipsoid.vol and biovolume will be 
#' called to add two additional columns to the returned data table object.info
#' containing the volume, calculated based on the ellipsoid approximation and 
#' biovolume of each object
#' @param unwanted a vector or list of strings containing the names of the categories 
#' you wish to discard. The function will search for these keywords in the field
#' object_annotation_hierarchy, therefore you can target large categories that include
#' several species. For example if category A contains categories a,b and c, and 
#' you want to omit a,b and c you can simply put A as argument. By default, the 
#' function will keep all categories.
#' @return the result is a list of four data tables. $for.veg to be used by the 
#' vegan package, sample in rows and categories in columns, the entries being 
#' the count of objects being in the corresponding sample and category. $object.info 
#' contains all objects (rows) and their specifications (eg size, color...). $sample.info
#' contains all samples (rows) and their specifications (eg speed, volume filtered, 
#' concentration factor). $counts is the summary of the absolute and relative 
#' abundances of categories in each sample.

import.table <- function(files,volumes=TRUE,unwanted=" "){
  
  # load files
  if (length(files)>1){ # multiple files as input
    data <- data.table()
    for (file in files){
      data <- rbind(data,fread(file))
      print("Note: You have given more than one file in input therefore the two tsv files will be combined")
    }
  } else { # only one file
    data <- fread(files)
  }

  # get rid of the rows containing unwanted objects
  data<-data[!grepl(paste(unwanted,collapse="|"),object_annotation_hierarchy),]
  
  # vectors with names of the columns of interest, je trouvais plus propre de déclarer les vecteurs avant
  for.object <- c("object_id","sample_id","object_area","object_width","object_height",
              "object_annotation_category","object_annotation_hierarchy", 
              "object_equivalent_diameter",
              "object_major","object_minor", "object_stdsaturation")
  for.sample <- c("sample_id","sample_project","process_id","process_pixel", 
              "sample_concentrated_sample_volume","sample_dilution_factor",
              "sample_speed_through_water","sample_total_volume",
              "acq_id","acq_imaged_volume","acq_minimum_mesh","acq_maximum_mesh")
  
  
  # rows with empty entries are omitted  ? e.g. in the Tara data, several entries were empty which caused errors downstream
  #data <- na.omit(data,unique(c(for.object,for.sample)))
  
  # summary of absolute abundances
  counts <- data[,.N, by=list(sample_id,object_annotation_category)]
  setnames(counts,c("sample_id","category","count"))

  # calculation of relative abundance and composition + addition of columns named 
  # "rel_abundance" and "composition"
  relative.abundance(counts)
  
  # abundances grouped to be used by the vegan package, contingency table
  for.veg <- counts
  for.veg <- for.veg %>% dcast(sample_id~category,fill=0,value.var = "rel_abundance")

  # information unique to every object
  object.info <- data[, ..for.object]
  setnames(object.info,c("object_id","sample_id","area","width","height","category",
                         "hierarchy","ESD", "major","minor","std_saturation"))
  
  # information in common to every object within a sample
  sample.info <- data[, ..for.sample]
  sample.info <- unique(sample.info)
  setnames(sample.info,c("sample_id","project","process_id","pixel_size",
                         "concentrated_sample_volume","dilution_factor","speed",
                         "filtered_volume", "acq_id","imaged_volume","min_mesh","max_mesh"))
  
  if (volumes) {
    ellipsoid.vol(sample.info,object.info) # adds additional columns
    biovolume(sample.info,object.info)
  }
  
  result <- list("for.veg"=for.veg,"object.info"=object.info,
                 "sample.info"=sample.info,"counts"=counts)
  return (result)
}
