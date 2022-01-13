# calculation of relative abundance according to the Hellinger method

setwd("/Users/claratrellu/Documents/année sab/plankton planet/Rscripts/Ecotaxa_R")
require(data.table)


relative.abundance <- function(for_veg){
  for.veg<-for_veg
  totals <- for.veg[,sum(count),by=sample]
    
  for (samp in totals[,sample]){
    for.veg[sample==samp,rel_abundance:=100*count/totals[sample==samp,V1]]
    for.veg[sample==samp,composition:=sqrt(rel_abundance)]
  }
  return(for.veg)
}
