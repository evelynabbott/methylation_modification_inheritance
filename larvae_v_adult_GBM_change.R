library(DESeq2)
library(tidyverse)
library(ggplot2)
library(viridis)
library(vegan)
library(ggpmisc)
library(gridExtra)
library(limma)
library(adegenet)
library(ggpubr)
library(lsmeans)

rm(list=ls())

rm(list=ls())
load("M_counts_coldata.Rdata")

coldata=coldata[(coldata$samtype=="adult"),]

coldata$geno=paste0(coldata$sire_geno,coldata$dam_geno)
#coldata$all_treat=paste0(coldata$sire_temp,coldata$dam_temp)

rownames(coldata)=coldata$samplenumber
keeps=rownames(coldata)
counts=counts[,keeps]
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata,
                             design=~sire_geno+sire_temp)
dds=DESeq(dds)
resultsNames(dds)

ah=results(dds,name="sire_temp_H_vs_A")
ah=data.frame("sire_temp_H_vs_A"=ah$log2FoldChange,"gene"=rownames(ah))

am=results(dds,name="sire_temp_M_vs_A")
am=data.frame("sire_temp_M_vs_A"=am$log2FoldChange)

adults=cbind(ah,am)

save(adults,file="adult_HM_blobs.Rdata")

#larvae ------------------------------
rm(list=ls())
load("M_counts_coldata.Rdata")

coldata=coldata[(coldata$samtype=="pooledlarvae"),]

coldata$geno=paste0(coldata$sire_geno,coldata$dam_geno)
coldata$all_treat=paste0(coldata$sire_temp,coldata$dam_temp)

#any parent temp
coldata$numtreat = ifelse(coldata$all_treat == "HH","H_two",
                          ifelse(coldata$all_treat == "AA","zero",
                                 ifelse(coldata$all_treat == "AH" | coldata$all_treat == "HA","H_one",
                                        ifelse(coldata$all_treat == "MM","M_two","none"))))

coldata$numtreat = factor(coldata$numtreat,levels=c("zero","H_one","H_two","M_two"))

rownames(coldata)=coldata$samplenumber
keeps=rownames(coldata)
counts=counts[,keeps]
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata,
                             design=~geno+numtreat)
dds=DESeq(dds)
resultsNames(dds)

lh1=results(dds,name="numtreat_H_one_vs_zero")
lh1=data.frame("numtreat_H_one_vs_zero"=lh1$log2FoldChange,"gene"=rownames(lh1))

lh2=results(dds,name="numtreat_H_two_vs_zero")
lh2=data.frame("numtreat_H_two_vs_zero"=lh2$log2FoldChange)

lm2=results(dds,name="numtreat_M_two_vs_zero")
lm2=data.frame("numtreat_M_two_vs_zero"=lm2$log2FoldChange)

larvae=cbind(lh1,lh2,lm2)

save(larvae,file="larvae_HM_blobs.Rdata")

#time to plot --------------------------------
rm(list=ls())
load("adult_HM_blobs.Rdata")
load("larvae_HM_blobs.Rdata")

all=full_join(adults,larvae,by="gene")

# cor(y=all$anytreat_1_vs_0,x=all$adult_H_v_A,method = "pearson")
x=all$sire_temp_H_vs_A
y=all$numtreat_H_one_vs_zero
summary(lm(y ~ x))

summary(lm(formula= numtreat_H_one_vs_zero ~ sire_temp_H_vs_A,data=all))


a=ggplot(all, aes(sire_temp_H_vs_A,numtreat_H_one_vs_zero))+
  geom_hex()+
  theme_classic()+
  geom_smooth(method='lm', formula= y~x)+
  scale_fill_viridis(trans="log")+
  ylim(-4,4)+
  xlim(-4,4)+
  #coord_cartesian(ylim = c(-3, 2))+
  xlab("adult H v C")+
  ylab("larvae 1 H v 0")+
  theme(legend.position="none")
# stat_cor(
#   aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
#   label.x = 0.8,label.y = -2.8,digits=5)

summary(lm(formula= numtreat_H_one_vs_zero ~ sire_temp_H_vs_A,data=all))
#Adjusted R-squared:  0.05355 ,  p-value: < 2.2e-16
#b = 0.12261

b=ggplot(all, aes(sire_temp_H_vs_A,numtreat_H_two_vs_zero))+
  geom_hex()+
  theme_classic()+
  geom_smooth(method='lm', formula= y~x)+
  scale_fill_viridis(trans="log")+
  ylim(-4,4)+
  xlim(-4,4)+
  #coord_cartesian(ylim = c(-3, 2))+
  xlab("adult H v C")+
  ylab("larvae 2 H v 0")+
  theme(legend.position="none")
# stat_cor(
#   aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
#   label.x = 0.8,label.y = -2.8,digits=3)

summary(lm(formula= numtreat_H_two_vs_zero ~ sire_temp_H_vs_A,data=all))
#Adjusted R-squared:  0.05238 , p-value: < 2.2e-16
#b = 0.15450  

c=ggplot(all, aes(sire_temp_M_vs_A,numtreat_M_two_vs_zero))+
  geom_hex()+
  theme_classic()+
  geom_smooth(method='lm', formula= y~x)+
  scale_fill_viridis(trans="log")+
  #coord_cartesian(ylim = c(-3, 2))+
  ylim(-4,4)+
  xlim(-4,4)+
  xlab("adult M v C")+
  ylab("larvae 2 M v 0")+
  theme(legend.position="none")
# stat_cor(
#   aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), #get r2
#   label.x = 0.8,label.y = -2.8)

summary(lm(formula= numtreat_M_two_vs_zero ~ sire_temp_M_vs_A,data=all))
#Adjusted R-squared: 0.00251,  p-value: 9.919e-08
#b = 0.02095  

grid.arrange(a,b,c,nrow=1)

