#volcano plot all: 1 v 2 parents
#2 v 0; 1 v 0
#H & M

##################
### MAKE VSDS ####
##################

rm(list=ls())
#adults ---------------------------------------
load("M_counts_coldata.Rdata")

coldata=coldata[(coldata$samtype=="adult"),]

rownames(coldata)=coldata$samplenumber
keeps=rownames(coldata)
counts=counts[,keeps]
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, 
                             design=~sire_geno+sire_temp)
dds=DESeq(dds)
resultsNames(dds)

#volcano plot
adh = results(dds, name ="sire_temp_H_vs_A")
adm = results(dds, name ="sire_temp_M_vs_A")

save(adh,adm,file="adult_HM_volcano.Rdata")

#larvae --------------------------
rm(list=ls())
load("M_counts_coldata.Rdata")

coldata=coldata[(coldata$samtype=="pooledlarvae"),]

coldata$geno=paste0(coldata$sire_geno,coldata$dam_geno)
coldata$all_treat=paste0(coldata$sire_temp,coldata$dam_temp)

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

#"numtreat_H_one_vs_zero" "numtreat_H_two_vs_zero" "numtreat_M_two_vs_zero"
plh2 = results(dds, name ="numtreat_H_two_vs_zero")
plh1 = results(dds, name ="numtreat_H_one_vs_zero")
plm = results(dds, name ="numtreat_M_two_vs_zero")

save(plh2,plh1,plm,file="larvae_HM_volcano.Rdata")

################# -------------------------
## MAKE PLOTS ###
#################
rm(list=ls())
load("larvae_HM_volcano.Rdata")
load("adult_HM_volcano.Rdata")

#adh, adm, plh1, plh2, plm

#ADULT HEAT ---------------
adh = data.frame(adh) %>% 
  arrange(pvalue) %>% 
  mutate(sig = !is.na(padj) & padj < 0.1) %>% 
  as_tibble()
adh
nsig2 = sum(adh$sig)
tot2 = nrow(adh)
adh=na.omit(adh)
a=ggplot(data=adh, aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
  geom_point(alpha=0.1,size=2) +
  theme_classic()+
  ylim(0,5)+
  xlim(-2,2)+
  scale_color_manual(values=c('black', 'red')) +
  #ggtitle("adults H v C")+
  labs(x=bquote(log[2]~'fold difference'),
       y=bquote("-"*log[10]~'p-value')) +
  #subtitle = "adults H v C") +
  theme(legend.position = "none")

#ADULT MED ---------------
adm = data.frame(adm) %>% 
  arrange(pvalue) %>% 
  mutate(sig = !is.na(padj) & padj < 0.1) %>% 
  as_tibble()
adm
nsig2 = sum(adm$sig)
tot2 = nrow(adm)
adm=na.omit(adm)
b=ggplot(data=adm, aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
  geom_point(alpha=0.1,size=2) +
  theme_classic()+
  ylim(0,5)+
  xlim(-2,2)+
  #ggtitle("adults M v C")+
  scale_color_manual(values=c('black', 'red')) + 
  labs(x=bquote(log[2]~'fold difference'),
       y=bquote("-"*log[10]~'p-value'))+
  #subtitle="adults M v C") +
  theme(legend.position = "none")

#LARVAE 1 H -----------------------
plh1 = data.frame(plh1) %>% 
  arrange(pvalue) %>% 
  mutate(sig = !is.na(padj) & padj < 0.1) %>% 
  as_tibble()
plh1
nsig2 = sum(plh1$sig)
tot2 = nrow(plh1)
plh1=na.omit(plh1)
c=ggplot(data=plh1, aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
  geom_point(alpha=0.1,size=2) +
  theme_classic()+
  ylim(0,5)+
  xlim(-2,2)+
  #ggtitle("larvae 1 H v 0 H")+
  scale_color_manual(values=c('black', 'red')) + 
  labs(x=bquote(log[2]~'fold difference'),
       y=bquote("-"*log[10]~'p-value'))+
  #subtitle = "larvae 1 H v 0 H") +
  theme(legend.position = "none")


#LARVAE 2 H --------------
plh2 = data.frame(plh2) %>% 
  arrange(pvalue) %>% 
  mutate(sig = !is.na(padj) & padj < 0.1) %>% 
  as_tibble()
plh2
nsig2 = sum(plh2$sig)
tot2 = nrow(plh2)
plh2=na.omit(plh2)
d=ggplot(data=plh2, aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
  geom_point(alpha=0.1,size=2) +
  theme_classic()+
  ylim(0,5)+
  xlim(-2,2)+
  scale_color_manual(values=c('black', 'red')) + 
  #ggtitle("larvae 2 H v 0 H")+
  labs(x=bquote(log[2]~'fold difference'),
       y=bquote("-"*log[10]~'p-value'))+
  #subtitle = "larvae 2 H v 0 H") +
  theme(legend.position = "none")


#LARVAE 2 M --------------
plm = data.frame(plm) %>% 
  arrange(pvalue) %>% 
  mutate(sig = !is.na(padj) & padj < 0.1) %>% 
  as_tibble()
plm
nsig2 = sum(plm$sig)
tot2 = nrow(plm)
plm=na.omit(plm)
e=ggplot(data=plm, aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
  geom_point(alpha=0.1,size=2) +
  theme_classic()+
  ylim(0,5)+
  xlim(-2,2)+
  #ggtitle("larvae 2 M v 0 M")+
  scale_color_manual(values=c('black', 'red')) + 
  labs(x=bquote(log[2]~'fold difference'),
       y=bquote("-"*log[10]~'p-value'))+
  #subtitle = "larvae 2 M v 0 M") +
  theme(legend.position = "none")


grid.arrange(a,b,c,d,e,nrow=2)
grid.arrange(a,c,d,b,e, nrow=2)






