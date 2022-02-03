# to use on tara data

source(file.path(getwd(), "R", "import_table.R"))
source(file.path(getwd(),"R","nbss.R"))
source(file.path(getwd(),"R","basic_diversity_plot.R"))
require(vegan)
library(magrittr)
install.packages("gridExtra")
library(gridExtra)

vec=c("not-living","multiple","duplicate")

# import
tables <- import.table("../tara1.tsv",unwanted=vec)
veg <- tables$for.veg
info.object <- tables$object.info
info.sample <- tables$sample.info
data.plot <- tables$counts

# PCA from the relative abundancies
forpca<-data.plot %>% dcast(sample_id~category,fill=0,value.var = "rel_abundance") # sample in rows, species in columns
forpca<-forpca[order(sample_id)]
env <- rda(X=forpca[,-1],scale=TRUE)
biplot(env,display = c("sites"),scaling=1,main="PCA-scaling1-from relative abundances") # rename sites=forpca$sample_id



# using biovolumes : 
sums<-summed.biovol(info.object)
sumspca<-sums %>% dcast(sample_id~category,fill=0,value.var = "summed_biovol")
#sumspca<-sumspca[,sample_id:=sub("^[^A-Za-z]+_","",sample_id)]
sumspca<-sumspca[order(sample_id)]
f<-setnafill(sumspca[,-1],fill=0)

# PCA from biovolumes
env2 <- rda(X=f,scale=TRUE)
env2 <- summary(env2)
forgg<-data.table(sumspca[,1],PC1=env2$sites[,1],PC2=env2$sites[,2])
s <- ggplot(data=forgg,mapping=aes(x=PC1,y=PC2,color=sample_id)) + geom_point()
s

# PCoA on biovolumes
d.bray.biov<-vegdist(sumspca[,-1],na.rm = TRUE)
pcoa.biovol<-cmdscale(d.bray.biov,k=2,eig=TRUE)
pPCo1 <- round(100*pcoa.biovol$eig[1]/sum(pcoa.biovol$eig))
pPCo2 <- round(100*pcoa.biovol$eig[2]/sum(pcoa.biovol$eig))
pcoa.b.plot<-data.table(sumspca[,1],PCo1=pcoa.biovol$points[,1],PCo2=pcoa.biovol$points[,2])
t<-ggplot(pcoa.b.plot,aes(x=PCo1,y=PCo2,color=sample_id)) + geom_point() +
  labs(title="PCoA to differentiate the different stations using biovolumes",
       x = paste("PCo1 ",pPCo1," %"),y=paste("PCo2 ",pPCo2," %")) 
t

# PCoA on standardized biovolumes
std.sums<-decostand(sumspca[,-1],"total",na.rm=T)
d.bray.std<-vegdist(std.sums,na.rm=T)
pcoa.std<-cmdscale(d.bray.std,k=2,eig=TRUE)
pcoa.std.plot<-data.table(sumspca[,1],PCo1=pcoa.std$points[,1],PCo2=pcoa.std$points[,2])
pPCo1 <- round(100*pcoa.std$eig[1]/sum(pcoa.std$eig))
pPCo2 <- round(100*pcoa.std$eig[2]/sum(pcoa.std$eig))
u<-ggplot(pcoa.std.plot,aes(x=PCo1,y=PCo2,color=sample_id)) + geom_point(show.legend=F) +   
  labs(title="PCoA to differentiate the different stations using standardized biovolumes", 
       x = paste("PCo1 ",pPCo1," %"), y=paste("PCo2 ",pPCo2," %")) 
u

#compare both
grid.arrange(t, u, ncol = 2, nrow = 1)
#nbss plot test
g<-BSS.plot(info.object,info.sample,info.sample[3,1])
g


# compare total biovolumes in between samples
total_biovolumes <- sumspca[,sum(summed_biovol),by=sample_id]
#total_biovolumes[,sample_id:=sub("^[^A-Za-z]+_","",sample_id)]
p<-ggplot(data=total_biovolumes,mapping=aes(y=sample_id,x=V1)) + geom_col()
p

# PCOA on absolute abundances
#veg[,sample_id:=sub("^[^A-Za-z]+_","",sample_id)]
veg<-veg[order(sample_id)]
d.bray<-vegdist(veg[,-1],method="bray") # Bray-Curtis distance actually by default in vegan
pcoa<-cmdscale(d.bray,k=2,eig=TRUE)
pcoa.plot<-data.table(veg[,1],PCo1=pcoa$points[,1],PCo2=pcoa$points[,2])
d<-ggplot(pcoa.plot,aes(x=PCo1,y=PCo2,color=my.colors)) + geom_point() +
  labs(title="PCoA on absolute abundances",
       x = paste("PCo1 ",round(100*pcoa$eig[1]/sum(pcoa$eig))," %"),
       y = paste("PCo2 ",round(100*pcoa$eig[2]/sum(pcoa$eig))," %"))  +
  theme(plot.title.position="panel")
d

# PCOA on absolute abundances with log transformation
veg.log <- log1p(veg[,-1])
d.bray.log<-vegdist(veg.log[,-1],method="bray") # Bray-Curtis distance actually by default in vegan
pcoa.l<-cmdscale(d.bray.log,k=2,eig=TRUE)
pcoa.log<-data.table(veg[,1],PCo1=pcoa.l$points[,1],PCo2=pcoa.l$points[,2])
l<-ggplot(pcoa.log,aes(x=PCo1,y=PCo2,color=my.colors)) + geom_point() +
  labs(title="PCoA on absolute abundances with log transformation",
       x = paste("PCo1 ",round(100*pcoa.l$eig[1]/sum(pcoa.l$eig))," %"),
       y = paste("PCo2 ",round(100*pcoa.l$eig[2]/sum(pcoa.l$eig))," %"))  +
  theme(plot.title.position="panel")
l

grid.arrange(d, l, ncol = 2, nrow = 1)

