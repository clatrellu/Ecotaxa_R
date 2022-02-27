require(data.table)
require(ggplot2)

#' This function returns a plot to show the composition and diversity of samples
#' @param data a data table with each row corresponding to a pair of category and
#' sample. The table must also have columns for relative and/or absolute abundances
#' @param threshold_rare some categories might contain very few organisms which 
#' can be grouped together in a new category : 'other'. By default this threshold 
#' is set to 0.03 (3%), which corresponds to the proportion of organisms represented
#' by this category. If you do not want to group the least abundant categories, set
#' this threshold to 0.
#' @return an object of type ggplot
#' @example  
#' plot<-basic.div(data.to.plot,0.05)
#' plot



basic.div <- function(data, threshold_rare=0.03){
  
  category_order <- data[,sum(count),by = category][order(V1,decreasing = TRUE),category]
  
  data[,category:=factor(category,levels = category_order)]
  
  #group all small abundances
  data<-data[rel_abundance < threshold_rare*100,category := "other"]

  ggplot(data, aes(y = sample_id, x = count, fill = category)) + 
    geom_col(position="fill") +
    labs(title = "Basic Diversity Plot", y = "Sample") +
    theme(plot.title = element_text(hjust=0.5,face="bold"))
  
}
