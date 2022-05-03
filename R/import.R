#' Import data from tsv file downloaded from Ecotaxa
#' @import data.table
#' @description This function imports one of several Ecotaxa tsv files into  
#' one data.table object with computed volumes.
#' @param files one or several tsv file names from a Planktoscope aquisition classified 
#' with Ecotaxa. If multiple files are inserted, make sure they are grouped in a list.
#' @param unwanted a vector or list of strings containing the names of the categories 
#' you wish to discard. 
#' @return A data table with all the fields from the original data.
#' @details 
#' Biovolume and volume information is added.
#' 
#' Date and time are removed from sample ids.
#' @export import.data

import.data <- function(files=NULL, unwanted = NA, taxo_column = "object_annotation_hierarchy"){
  if(length(files)>1){
    lapply(X = files, FUN = fread) |> 
      rbindlist() |>
      rm_unwanted(taxo_column = taxo_column,unwanted=unwanted) |>
      add.volumes() |>
      simplify_ids()
  }else if(length(files)==1){
    fread(files) |> 
      rm_unwanted(taxo_column = taxo_column,unwanted=unwanted) |>
      add.volumes() |>
      simplify_ids()
  } else{
    print("Please provide at least one file name")
  }
}


#' Remove unwanted objects from your data.table
#'@description When vignettes are classified on Ecotaxa, some of them are put 
#'in categories containing objects that have been segmented but that are not
#'interesting to study a sample. This can be for example artefacts of the 
#'acquisition, like bubbles, or detritus.
#'
#'@param input data.table from which we want to remove the unwanted objects
#'@param taxo_column to precise in which column to search for the unwanted
#'keywords
#'@unwanted a list containing the classification of the organisms to ignore
#'@return a data.table with specified lines removed
#'@details The function will search for these keywords in the field specified by 
#' the argument taxo_column, therefore you can target large categories that include
#' several species. For example if category A contains categories a,b and c, and 
#' you want to omit a,b and c you can simply put A as argument. By default, the 
#' function will keep all categories.
rm_unwanted <-  function(input,taxo_column ="object_annotation_hierarchy",unwanted=NA){
  if(sum(is.na(unwanted))==0){
    pattern <- paste(unwanted,collapse = "|")
    input[!grepl(pattern = pattern,get(taxo_column))]
  }else{
    input
  }
}


#' Add columns with the computed volume and biovolume of each particle
#' @description To compare samples, it is sometimes useful to look at biovolumes. 
#' Here we propose a way of computing the volume, and further the biovolume, of 
#' the objects based on the elipsoid approximation
#' @param planktotable data.table to which the columns will be added
#' @return A data.table with columns containing individual volumes and biovolumes
#' @details to compute biovolumes, we need to normalize the individual volumes
#' of the particle with the total filtered volume of the sample. 
add.volumes <- function(planktotable){
  planktotable[,volume:={
    pixel <- unique(process_pixel)*10**-3
    major_pix <- object_major * pixel
    minor_pix <- object_minor * pixel
    4/3*pi*major_pix/2*(minor_pix/2)**2
  },by=sample_id]
  
  planktotable[is.na(sample_dilution_factor),sample_dilution_factor:=1]
  
  planktotable[,biovolume:={
    norm_vol <- sample_concentrated_sample_volume/(sample_dilution_factor*acq_imaged_volume*sample_total_volume)
    volume*norm_vol
  },by=sample_id]
  planktotable[]
}

#' Makes the sample IDs neater
#' 
#' 
simplify_ids <- function(planktotable){
  planktotable[,sample_id:=sub("^\\d+/\\d+/\\d+_\\d+:\\d+_","",sample_id)]
  planktotable[]
}