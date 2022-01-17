# calculation of relative abundance according to the Hellinger method
require(data.table)

relative.abundance <- function(counts){
  totals <- counts[,sum(count),by=sample_id]
  for (samp in totals[,sample_id]){
    counts[sample_id==samp,rel_abundance:=100*count/totals[sample_id==samp,V1]]
    counts[sample_id==samp,composition:=sqrt(rel_abundance)]
  }
}
