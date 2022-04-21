#' Import data from tsv file downloaded from Ecotaxa
#' @import data.table
#' @description This function imports one of several Ecotaxa tsv files into  
#' one data.table object with computed volumes.
#' @param files one or several tsv file names from a Planktoscope aquisition classified 
#' with Ecotaxa. If multiple files are inserted, make sure they are grouped in a list.
#' @param unwanted a vector or list of strings containing the names of the categories 
#' you wish to discard. The function will search for these keywords in the field
#' specified by the argument taxo_column, therefore you can target large categories that include
#' several species. For example if category A contains categories a,b and c, and 
#' you want to omit a,b and c you can simply put A as argument. By default, the 
#' function will keep all categories.
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

rm_unwanted <-  function(input,taxo_column ="object_annotation_hierarchy",unwanted=NA){
  if(sum(is.na(unwanted))==0){
    pattern <- paste(unwanted,collapse = "|")
    input[!grepl(pattern = pattern,get(taxo_column))]
  }else{
    input
  }
}

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

simplify_ids <- function(planktotable){
  planktotable[,sample_id:=sub("^\\d+/\\d+/\\d+_\\d+:\\d+_","",sample_id)]
  planktotable[]
}