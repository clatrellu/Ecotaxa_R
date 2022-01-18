# calculation of relative abundance according to the Hellinger method
require(data.table)

#' The relative abundance is calculated simply by dividing the absolute count
#' of organisms from a certain category and sample from the total number of organisms
#' in this sample
#' A column containing the composition is also added, using the formula 
#' composition = sqrt(relative_abundance)
#' 
#' @param counts a data table with three columns, the sample, the category and 
#' the number of occurences of this category within the sample.
#' @note this function does not return anything as it modifies directly the data
#' table passed as argument.
#' 

relative.abundance <- function(counts){
  totals <- counts[,sum(count),by=sample_id]
  for (samp in totals[,sample_id]){
    counts[sample_id==samp,rel_abundance:=100*count/totals[sample_id==samp,V1]]
    counts[sample_id==samp,composition:=sqrt(rel_abundance)]
  }
}
