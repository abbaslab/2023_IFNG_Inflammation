library(survminer);library(survival)
library(survivalAnalysis) ##for multivariable cox model
library(cowplot)  ##for theme_cowplot
library(limma)  ## for lmfit in t statistic
library(tidyverse)
library(ggrepel)

############################################################################
#Figure 4A
############################################################################
nes.score=data.frame(cell=c("AML","CD4","CD8","NK"), nes=c(-1.533232,1.096463,1.63423,1.573604),
                     error=c(0.3807304,0.2042948,0.4550599,0.4070179))
ggplot(data=nes.score, aes(x=cell, y=nes)) +
  geom_bar(stat="identity",fill=as.character(ColAssign(letters[1:20])[c(1,3:5)])) +
  ylab("NES ") + theme_classic() +
  geom_hline(yintercept = 0,linetype="solid",size=0.75) + 
  geom_errorbar(aes(ymin=nes-error, ymax=nes+error), width=.1) +
  theme(legend.position="none",
        axis.title.y = element_text(size = 12),axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12,angle = 0,hjust = 0.5),axis.text.y = element_text(size = 12))
ggsave("Research/manuscript/cell_of_origin/paper/Figure4/fig4_AML_T_NK_IFNG_Production_NES.pdf",width=8,height =10,unit="cm")


############################################################################
#Figure 4B-E code adapted from tutorial of cellchat and multichnet
############################################################################

############################################################################
#Figure 4H
############################################################################
distance=read.csv("Research/manuscript/cell_of_origin/paper/Figure4/CD34toCD3dist.csv")
distance=reshape2::melt(distance[,1:3], id.vars = c("Name"), measure.vars = c("CD34HLAEtoCD3","CD34toCD3"),
                        variable.name = "group", value.name = "dist")

ggplot(distance,aes(x=group,y=dist))+#scale_y_continuous(trans=squish_trans(0,6),breaks = seq(6,15))+
  geom_boxplot(aes(fill=group),width=0.5,color='gray40',outlier.shape=NA,size=0.2)+
  ylab("Distance to T cells") + theme_classic() + geom_line(aes(group=Name)) + geom_point() +
  scale_fill_manual(values=c("bisque2","darkseagreen3","indianred2"))+
  theme(legend.position="none",
        axis.title.y = element_text(size = 12),axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12,angle = 30,hjust = 1),axis.text.y = element_text(size = 12))

distance %>% pairwise_wilcox_test(dist ~ group, p.adjust.method = "BH",paired=TRUE)
t.test(dist ~ group,data=distance, paired=TRUE)
wilcox.test(distance$dist[1:10],distance$dist[11:20],paired = TRUE)


