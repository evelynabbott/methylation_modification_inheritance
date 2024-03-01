#PCAs

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(viridis)
library(vegan)
library(ggpmisc)
library(gridExtra)
library(limma)
library(adegenet)
library(pheatmap)

#####################
# ADULT x PL M INCL #
#####################
rm(list=ls())
load("M_counts_coldata.Rdata")

coldata=coldata[(coldata$samtype=="adult" | coldata$samtype=="pooledlarvae"),]

coldatap=coldata[(coldata$samtype=="pooledlarvae"),]
coldataa=coldata[(coldata$samtype=="adult"),]

coldatap$geno=paste0(coldatap$sire_geno,coldatap$dam_geno)
coldataa$geno=coldataa$sire_geno

coldata=rbind(coldataa,coldatap)

rownames(coldata)=coldata$samplenumber
keeps=rownames(coldata)
counts=counts[,keeps]
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, 
                             design=~1)
dds=DESeq(dds)

resultsNames(dds)
res = results(dds)
rld = vst(dds)
rld.df=assay(rld)
vsd=rld.df

#vegan
dd=dist(t(vsd),method="manhattan")/nrow(vsd)
dds.pcoa=capscale(dd~1)
plot(dds.pcoa$CA$eig/sum(dds.pcoa$CA$eig))
sc1=data.frame(scores(dds.pcoa,scaling=1,choices=c(1:2))$sites)

coldata$geno2=ifelse((coldata$geno=="2926"|coldata$geno=="2629"),"26x29",
                     ifelse((coldata$geno=="2630"|coldata$geno=="3026"),"30x26",
                            ifelse((coldata$geno=="3126"|coldata$geno=="2631"),"26x31",
                                   ifelse((coldata$geno=="3029"|coldata$geno=="2930"),"29x30",
                                          ifelse((coldata$geno=="2931"|coldata$geno=="3129"),"31x29",
                                                 ifelse((coldata$geno=="3031"|coldata$geno=="3130"),"3130",
                                                        ifelse((coldata$geno=="26"),"26",
                                                               ifelse((coldata$geno=="29"),"29",
                                                                      ifelse((coldata$geno=="30"),"30","31")))))))))

coldata$geno2=as.factor(coldata$geno2)
coldata$geno2=factor(coldata$geno2,levels=c("26","29","30","31","26x29","30x26","26x31","29x30","31x29","3130"))
coldata1=coldata
coldata1$samtype = factor(coldata$samtype, levels = c("adult","pooledlarvae"))
coldata1$shapes=ifelse(coldata$geno2 == "26","a",
                      ifelse(coldata$geno2 == "29", "b",
                             ifelse(coldata$geno2 == "30", "c",
                                    ifelse(coldata$geno2 == "31", "d", "e"))))

coldata1$shapes = factor(coldata1$shapes, levels = c("e","a","b","c","d"))

shapes1 = c(19,15,17,25,23)

all=ggplot(sc1,aes(MDS1,MDS2,color=coldata1$geno2))+
  geom_point(aes(color=coldata1$geno2,shape=coldata1$shapes),alpha=0.5,size=3)+
  scale_shape_manual(values=shapes1)+
  scale_color_manual(values=c("26"="black",
                              "29"="darkgray",
                              "30"="azure4",
                              "31"="azure3",
                              "26x29"="#00798c",
                              "30x26"="#d1495b",
                              "26x31"="#edae49",
                              "29x30"="#66a182",
                              "31x29"="royalblue3",
                              "3130"="purple"))+
  theme_classic()+
  guides(color=guide_legend("Genotype"))+
  stat_ellipse(aes(color=coldata1$samtype))+
  theme(legend.position="none")

#Heatmap
pheatmap(cor(vsd),fontsize = 7, 
         labels_row= coldata$geno2,
         labels_col = coldata$samtype,
         cluster_rows = T,
         cluster_cols=T,
         legend = F,
         treeheight_row = 0,
         main="")

#adonis
adonis2(t(vsd)~samtype,data=coldata,method="manhattan")

###################
###### PL ######## 
###################
#rm(list=ls())
load("M_counts_coldata.Rdata")

coldata=coldata[(coldata$samtype=="pooledlarvae"),]

coldata$geno=paste0(coldata$sire_geno,coldata$dam_geno)

rownames(coldata)=coldata$samplenumber
keeps=rownames(coldata)
counts=counts[,keeps]
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, 
                             design=~1)
dds=DESeq(dds)

resultsNames(dds)
res = results(dds)
rld = vst(dds)
rld.df=assay(rld)
vsd=rld.df

#vegan
dd=dist(t(vsd),method="manhattan")/nrow(vsd)
dds.pcoa=capscale(dd~1)
plot(dds.pcoa$CA$eig/sum(dds.pcoa$CA$eig))
sc=data.frame(scores(dds.pcoa,scaling=1,choices=c(1:2))$sites)

coldata$geno2=ifelse((coldata$geno=="2926"|coldata$geno=="2629"),"26x29",
                     ifelse((coldata$geno=="2630"|coldata$geno=="3026"),"30x26",
                            ifelse((coldata$geno=="3126"|coldata$geno=="2631"),"26x31",
                                   ifelse((coldata$geno=="3029"|coldata$geno=="2930"),"29x30",
                                          ifelse((coldata$geno=="2931"|coldata$geno=="3129"),"31x29","31x30")))))

coldata$geno2=as.factor(coldata$geno2)
coldata$geno2=factor(coldata$geno2,levels=c("26x29","30x26","26x31","29x30","31x29","31x30"))
coldata2=coldata
sc2=sc

larv=ggplot(sc2,aes(MDS1,MDS2,color=coldata2$geno2))+
  geom_point(aes(color=coldata2$geno2),alpha=0.5,size=3)+
  scale_color_manual(values=c("26x29"="#00798c",
                              "30x26"="#d1495b",
                              "26x31"="#edae49",
                              "29x30"="#66a182",
                              "31x29"="royalblue3",
                              "31x30"="purple"))+
  theme_classic()+
  stat_ellipse()+
  guides(color=guide_legend("Genotype"))+
  theme(legend.position="none")


pheatmap(cor(vsd),fontsize = 7, 
         labels_row= coldata$geno2,
         labels_col = coldata$geno2,
         cluster_rows = T,
         cluster_cols=T,
         legend = F,
         treeheight_row = 0,
         main="")


#Adults only --------------------------------
load("M_counts_coldata.Rdata")

coldata=coldata[(coldata$samtype=="adult"),]

coldata$geno=coldata$sire_geno

rownames(coldata)=coldata$samplenumber
keeps=rownames(coldata)
counts=counts[,keeps]
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, 
                             design=~1)
dds=DESeq(dds)

resultsNames(dds)
res = results(dds)
rld = vst(dds)
rld.df=assay(rld)
vsd=rld.df

#vegan
dd=dist(t(vsd),method="manhattan")/nrow(vsd)
dds.pcoa=capscale(dd~1)
plot(dds.pcoa$CA$eig/sum(dds.pcoa$CA$eig))
sc=data.frame(scores(dds.pcoa,scaling=1,choices=c(1:2))$sites)

#coldata$geno=factor(coldata$geno,levels=c("26","29","30","31"))
coldata3=coldata
sc3=sc

coldata3$shapes=ifelse(coldata3$geno == "26","a",
                       ifelse(coldata3$geno == "29", "b",
                              ifelse(coldata3$geno == "30", "c","d")))

coldata3$shapes = factor(coldata3$shapes, levels = c("a","b","c","d"))

shapes = c(15,17,25,23)

adult=ggplot(sc3,aes(MDS1,MDS2))+
  geom_point(aes(fill=coldata3$geno,shape=coldata3$geno),alpha=0.5,size=3)+
  scale_shape_manual(values=shapes)+
  scale_fill_manual(values=c("26"="black",
                             "29"="darkgray",
                             "30"="azure4",
                             "31"="azure3"))+
  theme_classic()+
  stat_ellipse(aes(group=coldata3$geno))+
  guides(color=guide_legend("Genotype"))
  theme(legend.position="none")


grid.arrange(all,larv,adult,nrow=1)

