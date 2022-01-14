setwd("/Users/claratrellu/Documents/année sab/plankton planet/Rscripts")
require(data.table)
require(ggplot2)

basic.div <- function(plot.data){
  
  #rename sample names
  plot.data[,sample:=sub("^[^A-Za-z]+_","",sample)]
  
  taxo_order <- plot.data[,sum(count),by = taxo][order(V1,decreasing = TRUE),taxo]
  
  plot.data[,taxo:=factor(taxo,levels = taxo_order)]
  
  # group all small abundances
  plot.data[as.numeric(factor(taxo))%in%c(12:26),taxo:= "other"]
  
  setnames(plot.data,"taxo","categories")
  
  ggplot(plot.data, aes(y = sample, x = count, fill = categories)) + 
    geom_col(position="fill") +
    labs(title = "Différence de compositions entre les filets à plancton", y = "Net/speed")+
    theme(plot.title = element_text(hjust=0.5,face="bold"))
  
}