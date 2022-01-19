source(file.path(getwd(), "R", "import_table.R"))
source(file.path(getwd(),"R","nbss.R"))
require(vegan)
library(magrittr)

vec=c("not-living","duplicate","multiple")

# import
tables <- import.table("../global.tsv",unwanted=vec)
veg <- tables$for.veg
info.object <- tables$object.info
info.sample <- tables$sample.info
data.plot <- tables$counts

data.plot<-data.plot[,sample_id:=sub("^[^A-Za-z]+_","",sample_id)] # to have only the right net names, not ordered according to date and time
# PCA from the relative abundancies

forpca<-data.plot %>% dcast(sample_id~category,fill=0,value.var = "rel_abundance") # sample in rows, species in columns
forpca<-forpca[order(sample_id)]
env <- rda(X=forpca[,-1],scale=TRUE)
biplot(env,display = c("species","sites"),scaling=1,main="PCA-scaling1-from relative abundances") # rename sites=forpca$sample_id


# PCA from biovolumes

sumspca<-summed.biovol(info.object)
forpca2<-sumspca %>% dcast(sample_id~category,fill=0,value.var = "summed_biovol")
forpca2<-forpca2[order(sample_id)]
env2 <- rda(X=forpca2[,-1],scale=TRUE)
biplot(env2,display = c("sites"),scaling=1,main="PCA-scaling1-from biovolumes")
legend(1,y=1,legend=forpca2$sample_id)


#nbss plot test
g<-NBSS.plot(info.object,info.sample,info.sample[7,1])
g


# compare total biovolumes in between samples

total_biovolumes <- sumspca[,sum(summed_biovol),by=sample_id]
total_biovolumes[,sample_id:=sub("^[^A-Za-z]+_","",sample_id)]
p<-ggplot(data=total_biovolumes,mapping=aes(y=sample_id,x=V1)) + geom_col()
p

# PCOA and Bray-Curtis

veg[,sample_id:=sub("^[^A-Za-z]+_","",sample_id)]
veg<-veg[order(sample_id)]
d.bray<-vegdist(veg[,-1],method="bray") # Bray-Curtis distance actually by default in vegan
dataa<-cmdscale(d.bray,k=2,eig=TRUE)
daata<-data.table(veg[,1],PCo1=dataa$points[,1],PCo2=dataa$points[,2])

my.colors <- c(rep("Coryphaena",8),rep("Decknet",3),rep("Fanon",3))
d<-ggplot(daata,aes(x=PCo1,y=PCo2,color=my.colors)) + geom_point() +
    labs(title="PCoA to differentiate the different nets") 
d



# alternatives
# d.bray<-vegdist(veg[,-1],method="bray") # Bray-Curtis distance actually by default in vegan
# forpcoa<-cmdscale(d.bray,k=2,eig = TRUE) # cmdscale function to do principal coordinate
# biplot(forpcoa$points,display="sites")

# veg[,sample_id:=sub("^[^A-Za-z]+_","",sample_id)]
# veg<-veg[order(sample_id)]
# d.bray<-vegdist(veg[,-1],method="bray") # Bray-Curtis distance actually by default in vegan
# forpcoa<-pcoa(d.bray)
# biplot(forpcoa)



