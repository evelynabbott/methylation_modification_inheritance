########################################
##### make standard naming system ######
########################################
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(viridis)
library(vegan)
library(ggpmisc)
library(gridExtra)
library(limma)
library(adegenet)

rm(list=ls())

#for mdRAD -------------------
counts=read.table(file="feature_counts_out.tsv",sep = "\t",header = T,row.names='Geneid')
counts = counts[,6:ncol(counts)]

#remove genes with few counts
cut=3
cc=counts
means=apply(cc,1,mean)
table(means>cut) #6055 genes below cutoff
counts=cc[means>cut,]

#remove tagseq samples
# counts = counts %>%
#   select(-contains('tag'))

#fix names
n<-gsub("Lamarck.mdRAD.","",colnames(counts))
n<-gsub("Lamararck.mdRAD.","",n)
colnames(counts)<-n
n<-gsub(".bam","",colnames(counts))
colnames(counts)<-n


coldata=read.table(file="coldata.csv",sep=",",header = T)

#are samids unique?
test=unique(coldata$samid)
length(test) #yes

#match up and name adults
col_adult = coldata %>%
  filter(str_detect(samtype, "adult"))
counts_ad=counts[, grep("adult", names(counts))]
names(counts_ad)=gsub("_.*","",names(counts_ad))
names(counts_ad)=gsub("*adult.","",names(counts_ad))
#names(counts_ad)=as.factor(counts_ad)
ord=names(counts_ad)
col_adult=col_adult[match(ord,col_adult$samid),]
ord1=col_adult$samid
identical(ord,ord1)
#test=counts_ad
names(counts_ad)=c(1:24)
#coltest=col_adult
col_adult$samid=as.factor(col_adult$samid)
col_adult$samplenumber = c(1:24)

#match up and name sperm
col_sperm = coldata %>%
  filter(str_detect(samtype, "sperm"))
counts_sperm=counts[, grep("sperm", names(counts))]
names(counts_sperm)=gsub("_.*","",names(counts_sperm))
names(counts_sperm)=gsub("*sperm.","",names(counts_sperm))
ord=names(counts_sperm)
col_sperm=col_sperm[match(ord,col_sperm$samid),]
ord1=col_sperm$samid
identical(ord,ord1)
names(counts_sperm)=c(25:36)
col_sperm$samplenumber = c(25:36)

#match up and name single larvae
col_singlelarvae = coldata %>%
  filter(str_detect(samtype, "singlelarvae"))
counts_sl=counts[, grep("singlelarvae", names(counts))]
names(counts_sl)=gsub("_.*","",names(counts_sl))
names(counts_sl)=gsub("*singlelarvae.","",names(counts_sl))
col_singlelarvae$samid=gsub("-", ".", col_singlelarvae$samid)
ord=names(counts_sl)
col_singlelarvae=col_singlelarvae[match(ord,col_singlelarvae$samid),]
ord1=col_singlelarvae$samid
identical(ord,ord1)
names(counts_sl)=c(37:68)
col_singlelarvae$samplenumber = c(37:68)

#match up and name pooled larvae
col_pooledlarvae = coldata %>%
  filter(str_detect(samtype, "pooledlarvae"))
counts_pooledlarvae=counts[, grep("pooledlarvae", names(counts))]
names(counts_pooledlarvae)=gsub("_.*","",names(counts_pooledlarvae))
names(counts_pooledlarvae)=gsub("*pooledlarvae.","",names(counts_pooledlarvae))
ord=names(counts_pooledlarvae)
col_pooledlarvae=col_pooledlarvae[match(ord,col_pooledlarvae$samid),]
ord1=col_pooledlarvae$samid
identical(ord,ord1)
names(counts_pooledlarvae)=c(69:125)
col_pooledlarvae$samplenumber = c(69:125)


#recombine metadata with sample numbers
coldata=rbind(col_adult,col_sperm,col_singlelarvae,col_pooledlarvae)
coldata <- coldata[order(coldata$samplenumber, decreasing = FALSE),]
rownames(coldata)=coldata$samplenumber

counts=cbind(counts_ad,counts_sperm,counts_sl,counts_pooledlarvae)

save(counts,coldata,file = "all_samps_counts_coldata.Rdata")
# 

########################################
#### make columns work for deseq #######
########################################

#as is, model matrix is not full rank
#make a single column for genotype

rm(list=ls())
load("all_samps_counts_coldata.Rdata")

#make separate columns for parent geno and parent temp ------
#how many of each tissue type?
s=coldata[(coldata$samtype=="sperm"),]
a=coldata[(coldata$samtype=="adult"),]
sl=coldata[(coldata$samtype=="singlelarvae"),]
pl=coldata[(coldata$samtype=="pooledlarvae"),]

#larvae
coldata1=coldata
coldata1=coldata[!(coldata$samtype=="pooledlarvae"|coldata$samtype=="singlelarvae"),]
coldata=coldata[(coldata$samtype=="pooledlarvae"|coldata$samtype=="singlelarvae"),]
coldata$sire_geno <- as.factor(str_extract(coldata$sire, "[0-9]+"))
coldata$sire_temp <- as.factor(str_extract(coldata$sire, "[aA-zZ]+"))
coldata$dam_geno <- as.factor(str_extract(coldata$dam, "[0-9]+"))
coldata$dam_temp <- as.factor(str_extract(coldata$dam, "[aA-zZ]+"))
coldata=dplyr::bind_rows(coldata1,coldata)

#adults and sperm
coldata1=coldata
coldata1=coldata[!(coldata$samtype=="adult"|coldata$samtype=="sperm"),]
coldata=coldata[(coldata$samtype=="adult"|coldata$samtype=="sperm"),]
coldata$sire_geno <- as.factor(str_extract(coldata$samid, "[0-9]+"))
coldata$sire_temp <- as.factor(str_extract(coldata$samid, "[aA-zZ]+"))
coldata$dam_geno <- as.factor(str_extract(coldata$samid, "[0-9]+"))
coldata$dam_temp <- as.factor(str_extract(coldata$samid, "[aA-zZ]+"))
coldata=dplyr::bind_rows(coldata1,coldata)

rownames(coldata)=as.numeric(coldata$samplenumber)
coldata <- coldata[order(coldata$samplenumber, decreasing = FALSE),]

save(counts,coldata,file="all_samps_counts_coldata.Rdata")

#keep "moderate" samples ---------------------------------------
rm(list=ls())
load("all_samps_counts_coldata.Rdata")
#rownames(coldata) = coldata$samplenumber
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata,
                             design=~1)

dds=DESeq(dds)

resultsNames(dds)
res = results(dds)
rld = vst(dds)
rld.df=assay(rld)
vsd=rld.df

#sample clustering to detect outliers
gen=t(vsd)
par(cex=0.5)
sampleTree = hclust(dist(gen), method = "average");
plot(sampleTree, main = "", sub="", xlab="", cex.lab = 0.1, cex.axis = 0.7, cex.main = 1)
abline(h=50,col="red")

#color by log counts
coldata$countsum=log(colSums(counts))
coldata$countnum=colSums(counts)

#color according to life stage
# coldata$color = ifelse(grepl("adult",coldata$samtype,ignore.case = T),"red",
#                        ifelse(grepl("pooledlarvae",coldata$samtype,ignore.case = T),"orange",
#                               ifelse(grepl("singlelarvae",coldata$samtype,ignore.case = T),"green","blue")))



coldata$color = ifelse(coldata$samplenumber == "27" |
                         coldata$samplenumber == "58" |
                         coldata$samplenumber == "13" |
                         coldata$samplenumber == "50" |
                         coldata$samplenumber == "1" |
                         coldata$samplenumber == "15" |
                         coldata$samplenumber == "9", "red",
                       ifelse(coldata$samplenumber == "18" |
                              coldata$samplenumber == "7" |
                              coldata$samplenumber == "2" | 
                              coldata$samplenumber == "26" | 
                              coldata$samplenumber == "78","blue","black"))

#PCA
dds.pcoa=capscale(dist(t(vsd),method="manhattan")~1)
plot(dds.pcoa$CA$eig/sum(dds.pcoa$CA$eig))
scores=dds.pcoa$CA$u

#
lab=coldata$samplenumber
#lab=coldata$countsum
col=coldata$color
#
axes2plot=c(1,2)
plot(scores[,axes2plot],col=col,pch=19,main="",sub="",cex=1.5)
ordispider(scores[,axes2plot],lab,label=T,col="grey70",cex=1.5)


#plot with ggplot
sc1=data.frame(scores(dds.pcoa,scaling=1,choices=c(1:2))$sites)

sc1=cbind(sc1,coldata)

a=ggplot(sc1,aes(MDS1,MDS2,color=countsum))+
  geom_point(aes(color=countsum),size=3)+
  theme_classic()+
  scale_color_viridis(direction = -1)

#plot counts as histogram
b=ggplot(coldata,aes(countsum))+
  geom_histogram(bins=200,color="black",fill="white")+
  theme_classic()+
  ylab("number of samples")+
  xlab("log(counts)")+
  geom_vline(aes(xintercept=mean(countsum)),
             color="red", linetype="dashed", size=0.5)


grid.arrange(a,b,nrow=1)


#remove outliers (samples with low counts)
coldata=coldata[coldata$countsum > 12,]
coldata = coldata[!(coldata$color == "red" | coldata$color == "blue"),]
mean(coldata$countsum) #14.3
keeps = coldata$samplenumber
counts=counts[,keeps]

save(counts,coldata,file="testcounts_coldata.Rdata")

#test with outliers removed ----------------------------------
rm(list=ls())
load("testcounts_coldata.Rdata")
#did this fix the outliers?
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata,
                             design=~1)

dds=DESeq(dds)

resultsNames(dds)
res = results(dds)
rld = vst(dds)
rld.df=assay(rld)
vsd=rld.df

#sample clustering to detect outliers
gen=t(vsd)
par(cex=0.5)
sampleTree = hclust(dist(gen), method = "average");
plot(sampleTree, main = "", sub="", xlab="", cex.lab = 0.1, cex.axis = 0.7, cex.main = 1)
abline(h=50,col="red")

#PCA
dds.pcoa=capscale(dist(t(vsd),method="manhattan")~1)
plot(dds.pcoa$CA$eig/sum(dds.pcoa$CA$eig))
scores=dds.pcoa$CA$u


lab=coldata$samplenumber
#lab=coldata$countsum
#col=coldata$color

axes2plot=c(1,2)
plot(scores[,axes2plot],pch=19,main="",sub="",cex=1.5)
ordispider(scores[,axes2plot],lab,label=T,col="grey70",cex=1.5)


#plot with ggplot
sc1=data.frame(scores(dds.pcoa,scaling=1,choices=c(1:2))$sites)

sc1=cbind(sc1,coldata)

a=ggplot(sc1,aes(MDS1,MDS2,color=countsum))+
  geom_point(aes(color=countsum),size=3)+
  theme_classic()+
  scale_color_viridis(direction = -1)
  
#plot counts as histogram
b=ggplot(coldata,aes(countsum))+
  geom_histogram(bins=200,color="black",fill="white")+
  theme_classic()+
  ylab("number of samples")+
  xlab("log(counts)")+
  geom_vline(aes(xintercept=mean(countsum)),
             color="red", linetype="dashed", size=0.5)


grid.arrange(a,b,nrow=1)

#remove color column
coldata = coldata[, -which(names(coldata) %in% c("color"))]

save(counts,coldata,file="all_samps_counts_coldata.Rdata")

#-------------------------
# MAKE DATAFRAME FOR JUST ADULTS AND LARVAE #
rm(list=ls())
load("all_samps_counts_coldata.Rdata")

coldata=coldata[(coldata$samtype=="adult" | coldata$samtype=="pooledlarvae"),]

coldatap=coldata[(coldata$samtype=="pooledlarvae"),]
coldataa=coldata[(coldata$samtype=="adult"),]

coldatap$geno=paste0(coldatap$sire_geno,coldatap$dam_geno)
coldataa$geno=coldataa$sire_geno

coldata=rbind(coldataa,coldatap)

keeps = coldata$samplenumber
counts=counts[,colnames(counts) %in% keeps]

#add "all treat" to represent parent and larvae treatments
coldata=coldata %>%
  mutate(
    all_treat = case_when(
      sire_temp == "H" & dam_temp == "H" ~ "1",
      sire_temp == "H" & dam_temp == "M" ~ "0.75",
      sire_temp == "H" & dam_temp == "A" ~ "0.5",
      sire_temp == "M" & dam_temp == "M" ~ "0.5",
      sire_temp == "M" & dam_temp == "H" ~ "0.75",
      sire_temp == "M" & dam_temp == "A" ~ "0.25",
      sire_temp == "A" & dam_temp == "A" ~ "0",
      sire_temp == "A" & dam_temp == "M" ~ "0.25",
      sire_temp == "A" & dam_temp == "H" ~ "0.5",
      TRUE ~ ""
    ))


save(counts,coldata,file="counts_coldata.Rdata")




