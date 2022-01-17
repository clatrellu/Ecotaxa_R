source(file.path(getwd(), "R", "import_table.R"))
source(file.path(getwd(),"R","nbss.R"))
#source("cleanplot.pca.R")
source("PCA.R")
require(vegan)

# import
tables <- import.table("../global.tsv")
veg <- tables$for.veg
info.object <- tables$object.info
info.sample <- tables$sample.info
data.plot <- tables$counts


# PCA from the relative abundancies

forpca<-data.plot %>% dcast(sample_id~category,fill=0,value.var = "rel_abundance") # sample in rows, species in columns
forpca[,sample_id:=sub("^[^A-Za-z]+_","",sample_id)] # have only the right net names, not ordered according to date and time
forpca<-forpca[order(sample_id)]
env <- rda(X=forpca[,-1],scale=TRUE)
biplot(env,display = c("species","sites"),scaling=1,main="PCA-scaling1-from biovolumes") # rename sites=forpca$sample_id


# PCA from biovolumes

sumspca<-summed.biovol(info.object)
forpca2<-sumspca %>% dcast(sample_id~category,fill=0,value.var = "summed_biovol")
forpca2[,sample_id:=sub("^[^A-Za-z]+_","",sample_id)] # have only the right net names, not ordered according to date and time
forpca2<-forpca2[order(sample_id)]
env2 <- rda(X=forpca2[,-1],scale=TRUE)
biplot(env2,display = c("sites"),scaling=1,main="PCA-scaling1-from biovolumes")
legend(1,y=1,legend=forpca2$sample_id)


#nbss plot test
g<-NBSS.plot(info.object,info.sample,info.sample[1,1])
g
