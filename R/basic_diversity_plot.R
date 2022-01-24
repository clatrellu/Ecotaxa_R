require(data.table)
require(ggplot2)

# doesn't work for Tara samples

basic.div <- function(plot.data){
  
  #rename sample names
  plot.data[,sample:=sub("^[^A-Za-z]+_","",sample)]
  
  category_order <- plot.data[,sum(count),by = category][order(V1,decreasing = TRUE),category]
  
  plot.data[,category:=factor(category,levels = category_order)]
  
  # group all small abundances
  plot.data[as.numeric(factor(category))%in%c(12:26),category:= "other"]
  
  setnames(plot.data,"category","categories")
  
  ggplot(plot.data, aes(y = sample, x = count, fill = categories)) + 
    geom_col(position="fill") +
    labs(title = "Différence de compositions entre les filets à plancton", y = "Net/speed")+
    theme(plot.title = element_text(hjust=0.5,face="bold"))
  
}
