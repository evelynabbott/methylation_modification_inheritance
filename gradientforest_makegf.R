#TEST gradient forest
#does sire treatment matter more than dam treatment?

#gradient forest

rm(list=ls())
library(RDAforest)
library(gradientForest)
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
library(RDAforest)

#source("/Users/evelynabbott/Dropbox/project/Projects/lamk/lamk_new/RDAforest/R/RDAforest_functions.R")
source("RDAforest_functions.R")
load("M_counts_coldata.Rdata")

coldata=coldata[(coldata$samtype=="adult" | coldata$samtype=="pooledlarvae"),]

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

#vegan ----------------
dd=dist(t(vsd),method="manhattan")/nrow(vsd)
dds.pcoa=capscale(dd~1)
plot(dds.pcoa$CA$eig/sum(dds.pcoa$CA$eig))
#sc1=scores(dds.pcoa,scaling=1,choices=c(1:2)$sites)
sc=data.frame(scores(dds.pcoa,scaling=1,choices=c(1:2))$sites)
#sc3=scores(dds.pcoa$CA$u)

colnames(coldata)
coldata$all_treat=as.numeric(coldata$all_treat)
ggplot(sc,aes(MDS1,MDS2,color=coldata$all_treat))+
  geom_point()+
  theme_void()+
  coord_equal()

#make X --------------------
#load("M_counts_coldata.Rdata")

rownames(coldata)=coldata$samplenumber

coldata$adult=ifelse(grepl("adult", coldata$samtype),"1","0")
coldata$pl=ifelse(grepl("pooledlarvae",coldata$samtype),"1","0")

#pooled larvae
pl=coldata[(coldata$samtype=="pooledlarvae"),]

pl=pl %>% 
  mutate(
    all_treat = case_when(
      sire_temp == "H" & dam_temp == "H" ~ "1",
      sire_temp == "H" & dam_temp == "A" ~ "0.5",
      sire_temp == "M" & dam_temp == "M" ~ "0.5",
      sire_temp == "A" & dam_temp == "A" ~ "0",
      sire_temp == "A" & dam_temp == "H" ~ "0.5"
    ))

pl$t6=ifelse(pl$sire_geno=="26"|pl$dam_geno=="26",0.5,0)
pl$t9=ifelse(pl$sire_geno=="29"|pl$dam_geno=="29",0.5,0)
pl$th0=ifelse(pl$sire_geno=="30"|pl$dam_geno=="30",0.5,0)
pl$th1=ifelse(pl$sire_geno=="31"|pl$dam_geno=="31",0.5,0)


pl$t6s=ifelse(pl$sire_geno=="26",1,0)
pl$t9s=ifelse(pl$sire_geno=="29",1,0)
pl$th0s=ifelse(pl$sire_geno=="30",1,0)
pl$th1s=ifelse(pl$sire_geno=="31",1,0)

pl$t6d=ifelse(pl$dam_geno=="26",1,0)
pl$t9d=ifelse(pl$dam_geno=="29",1,0)
pl$th0d=ifelse(pl$dam_geno=="30",1,0)
pl$th1d=ifelse(pl$dam_geno=="31",1,0)

pl$s_h=ifelse(pl$sire_temp=="H",1,0)
pl$s_m=ifelse(pl$sire_temp=="M",1,0)
pl$s_a=ifelse(pl$sire_temp=="A",1,0)

pl$d_h=ifelse(pl$dam_temp=="H",1,0)
pl$d_m=ifelse(pl$dam_temp=="M",1,0)
pl$d_a=ifelse(pl$dam_temp=="A",1,0)

pl$h = ifelse(pl$sire_temp == "H" & pl$dam_temp == "H",1,
              ifelse(pl$sire_temp == "H" & pl$dam_temp != "H",0.5,
                     ifelse(pl$sire_temp != "H" & pl$dam_temp == "H",0.5,0)))
pl$m = ifelse(pl$sire_temp == "M" & pl$dam_temp == "M",1,
              ifelse(pl$sire_temp == "H" & pl$dam_temp != "H",0.5,
                     ifelse(pl$sire_temp != "M" & pl$dam_temp == "M",0.5,0)))
pl$a = ifelse(pl$sire_temp == "A" & pl$dam_temp == "A",1,
              ifelse(pl$sire_temp == "A" & pl$dam_temp != "A",0.5,
                     ifelse(pl$sire_temp != "A" & pl$dam_temp == "A",0.5,0)))

#adult - column averages for irrelevant columns
#sire or dam genotype
a=coldata[(coldata$samtype=="adult"),]
plmeans=as.vector(colMeans(pl[,c(23:30)]))

#adults with sire/dam genotype as the average of those columns in the larvae
a$t6s=plmeans[c(1)]
a$t9s=plmeans[c(2)]
a$th0s=plmeans[c(3)]
a$th1s=plmeans[c(4)]

a$t6d=plmeans[c(5)]
a$t9d=plmeans[c(6)]
a$th0d=plmeans[c(7)]
a$th1d=plmeans[c(8)]

#sire or dam genotype
plmeans=as.vector(colMeans(pl[,c(31:36)]))
a$s_h=plmeans[c(1)]
a$s_m=plmeans[c(2)]
a$s_a=plmeans[c(3)]
a$d_h=plmeans[c(4)]
a$d_m=plmeans[c(5)]
a$d_a=plmeans[c(6)] 

a$h = ifelse(a$sire_temp == "H",1,0)
a$m = ifelse(a$sire_temp == "M",1,0)
a$a = ifelse(a$sire_temp == "A",1,0)

#adults with sire/dam genotype as the average of those columns in the larvae
# a$all_treat= ifelse(grepl("H", a$sire_temp),"1",
#                     ifelse(grepl("M", a$sire_temp),"0.5","0"))

a$t6 = ifelse(grepl("26", a$sire_geno),"1","0")
a$t9 = ifelse(grepl("29", a$sire_geno),"1","0")
a$th0 = ifelse(grepl("30", a$sire_geno),"1","0")
a$th1 = ifelse(grepl("31", a$sire_geno),"1","0")

#combine, subset, and make numeric
coldata=rbind(pl,a)
which(colnames(coldata) == "adult") #17
coldata=coldata[,-c(1:16)]
coldata <- coldata %>% mutate_if(is.character, as.numeric)

#makegf
plot(dds.pcoa$CA$eig)
sb0=makeGF(dds.pcoa,coldata,keep=c(1:10))
plot(sb0)
a=data.frame(importance(sb0))
df=data.frame("importance"=a$importance.sb0.,"var"=rownames(a))
df1=df[order(df$importance),]

varlabs=rev(c("31","29","30","larvae","26","adult","C","30 sire","30 dam","31 dam","H","29 dam","26 dam","31 sire","29 sire","26 sire","H sire","M sire","M dam","H dam","C sire","C dam","M"))
varcols=rev(c("#84C0AA","#84C0AA","#84C0AA","#D17836","#84C0AA","#D17836","#E1AD49","#84C0AA","#84C0AA","#84C0AA","#E1AD49","#84C0AA","#84C0AA","#84C0AA","#84C0AA","#84C0AA","#E1AD49","#E1AD49","#E1AD49","#E1AD49","#E1AD49","#E1AD49","#E1AD49"))

df1$varcols=varcols


#colors:
#geno: #84C0AA
#lifestage: #D17836
#treat:#E1AD49


#pretty plot
gg1=ggplot(data=df1,aes(x=reorder(var,importance),y=importance,fill=varcols))+
  geom_bar(stat="identity")+
  theme_classic()+
  xlab("")+
  coord_flip()+
  scale_fill_manual(values=c("#84C0AA","#D17836","#E1AD49"))+
  #scale_fill_manual(values=c("#84C0AA","#84C0AA","#84C0AA","#D17836","#84C0AA","#D17836","#E1AD49","#84C0AA","#84C0AA","#84C0AA","#E1AD49","#84C0AA","#84C0AA","#84C0AA","#84C0AA","#84C0AA","#E1AD49","#E1AD49","#E1AD49","#E1AD49","#E1AD49","#E1AD49","#E1AD49"))+
  scale_x_discrete(labels=varlabs)+
  theme(legend.position = "none")

#export df1 as table
df1 = df1[order(df1$importance),]
df1$variable_names = varlabs
dftable = data.frame("importance"=df1$importance,"variable_names"=df1$variable_names)
dftable=t(dftable)
write.csv(dftable,file="gradforest_table.csv",quote = F)

#sum like categories
sumgeno= df1[df1$varcols == "#84C0AA",]
sumls = df1[df1$varcols == "#D17836",]
sumtreat = df1[df1$varcols == "#E1AD49",]

geno = sum(sumgeno$importance)
ls = sum(sumls$importance)
treat = sum(sumtreat$importance)

dfsums = data.frame("importance"=c(geno,ls,treat),"variable"=c("genotype","lifestage","treatment"))

#plot pretty
gg2=ggplot(data=dfsums,aes(x=reorder(variable,importance),y=importance,fill=variable))+
  geom_bar(stat="identity")+
  theme_classic()+
  xlab("")+
  coord_flip()+
  scale_fill_manual(values = c("#84C0AA","#D17836","#E1AD49"))

grid.arrange(gg1,gg2,nrow=2)


