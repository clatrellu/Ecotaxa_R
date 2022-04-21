#' Extract total counts and volumes
#' @description This function counts how many objects has been observed 
#' per category (object_annotation_category) for each sample. It also 
#' sums the volumes and biovolumes for each category.
#' @param planktotable a data table obtained using the function import.data()
#' @return A data table with 5 columns. The sample id, object category, counts,
#' volume and biovolumes.
#' @export get.counts.and.vol

get.counts.and.vol <- function(planktotable){
  counts <- planktotable[,.N, by=list(sample_id,object_annotation_category)]
  setnames(counts,c("sample_id","category","count"))
  counts[,count_proportion:=count/sum(count),by=sample_id]
  
  volumes <- planktotable[,.(sum(volume),sum(biovolume)),by=list(sample_id,object_annotation_category)]
  setnames(volumes,c("sample_id","category","volume","biovolume"))
  
  output <- merge(counts,volumes,by=c("sample_id","category"))
  output[]
}

#' Get description of samples or objects
#' @description This function extract the fields related to the objects or the samples
#' @param planktotable a data table obtained using the function import.data()
#' @param type type of fields to extract. Possible options are "sample" and "object"
#' @return A data tale with information related to samples or objects
#' @export get.info

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

#' Get a community data matrix
#' @description This function convert a long format table to a community
#' data matrix that can be directly used by functions from the library vegan 
#' @param x a long format table obtained using the function get.counts.and.vol()
#' @param value.var value column to used. Possible options are "count", "volume"
#' and "biovolume"
#' @return A community data matrix (data frame) with categories as columns and samples as rows.
#' @export get.df.vegan

get.df.vegan <- function(x,value.var="count"){
  dcast(sample_id~category,data=x,fill=0,value.var = value.var) |>
    data.frame(row.names="sample_id")
}


#' Taxa selection for plots
#' @description This function add an extra column containing the most abundant categories
#' @param x a long format table obtained using the function get.counts.and.vol()
#' @param categ column in the input table containing the category information
#' @param quantity column in the input table to quantifuy each category
#' @param ngroups maximum number of categories to keep
#' @return Add an extra column to the input table. This column contains the most abundant
#' categories, the number of categories being defined by ngroups. The other categories are
#' names "others"
#' @export plot_categ

plot_categ <- function(x,categ,quantity,ngroups=9){
  y <- x[,.(total=sum(get(quantity))),by=categ]
  if(nrow(y)>ngroups){
    y <- y[order(total,decreasing=TRUE),get(categ)][1:(ngroups-1)]
    x[get(categ)%in%y,category_plot:=get(categ)]
    x[!get(categ)%in%y,category_plot:="others"]
    x[,category_plot:=factor(category_plot,levels=c(y,"others"))]
  }else{
    x[,category_plot:=get(categ)]
    print("Less groups in the table than the number requested")
  }
}
