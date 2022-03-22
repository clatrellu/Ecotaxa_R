#' @import data.table
rm_unwanted <-  function(input,taxo_column ="object_annotation_hierarchy",unwanted=NA){
  if(sum(is.na(unwanted))==0){
    pattern <- paste(unwanted,collapse = "|")
    input[!grepl(pattern = pattern,get(taxo_column))]
  }else{
    input
  }
}

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

import.data <- function(files=NULL, unwanted = NA, taxo_column = "object_annotation_hierarchy"){
  if(length(files)>1){
    lapply(X = files, FUN = fread) |> 
      rbindlist() |>
      rm_unwanted(taxo_column = taxo_column,unwanted=unwanted)
  }else if(length(files)==1){
    fread(files) |> rm_unwanted(taxo_column = taxo_column,unwanted=unwanted)
  } else{
    print("Please provide at least one file name")
  }
}


get.counts <- function(planktotable,transfo="none"){
  counts <- planktotable[,.N, by=list(sample_id,object_annotation_category)]
  setnames(counts,c("sample_id","category","count"))
  if(transfo=="none"){
    list(wide=counts,
         long=dcast(sample_id~category,data=counts,fill=0,value.var = "count"))
  }else if(transfo=="total"){
    counts[,prop:=count/sum(count),by=sample_id]
    list(wide=counts,
         long=dcast(sample_id~category,data=counts,fill=0,value.var = "prop"))
  }

}

get.info <- function(planktotable,type="sample"){
  tmp_regex <- paste0("^",type,"_")
  cols <- grep(tmp_regex,names(planktotable),value=T)
  if(type!="sample"){
    cols <- c("sample_id",cols)
  }
  output <- planktotable[,.SD,.SDcols=cols]
  setnames(output,cols,sub(tmp_regex,"",cols))
  unique(output)
}

