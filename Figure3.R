library(SCP)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(purrr) # for map function
library(ggpubr)  # for ggarrange 
library(patchwork)  #for plot_layout
library(readxl)
library(ggcorrplot)
library(ggsignif)  #for stat_signif

##############################
#Fig 3A
##############################
df=allcell@meta.data[,c("orig.ident","class2","cyto4","HALLMARK_IFNG")]
df=df[df$class2 %in% c("T","AML"),]
#temp is the T cells and HSC in class1 of normal456
temp=normal456@meta.data[,c("orig.ident", "class2","class1","HALLMARK_INTERFERON_GAMMA_RESPONSE")]
temp$class2=as.character(temp$class2);temp$class2[temp$class1=="T"]="T"
temp=temp[temp$class2 %in% c("T","HSC","CMP.LMPP", "GMP"),]
temp$cyto4="Healthy";temp[temp$class2 %in% c("HSC","CMP.LMPP", "GMP"),"class2"]='AML'
temp=temp[,c(1,2,4,5)];colnames(temp)=c("orig.ident", "class2", "HALLMARK_IFNG", "cyto4")
temp=temp[,c("orig.ident", "class2", "cyto4", "HALLMARK_IFNG")]
df=rbind(df,temp)


df$cyto4=factor(df$cyto4,levels = c("Healthy","Diploid-Nonmono","Del5.5q","Double del","Del7.7q", "Diploid-mono"))
gm=df %>% group_by(cyto4,class2) %>% summarise(HALLMARK=mean(HALLMARK_IFNG),.groups = 'drop')
gm=as.data.frame(reshape(as.data.frame(gm[,c(1,2,3)]), timevar = "class2", idvar = "cyto4", direction = "wide"))
gm$ratio=gm$HALLMARK.AML/gm$HALLMARK.T
rownames(gm)=gm[,1];gm=gm[,2:4];gm=t(gm)
gm1=rbind(rep(0.04,5),gm);gm1=rbind(rep(0.12,5),gm1)
gm1=as.data.frame(gm1);rownames(gm1)[1:2]=c("max","min")

radarchart(gm1[c(1:4),], axistype = 1,
           # Customize the polygon
           pcol = c("#00AFBB","#E7B800"),  plwd = 3.5, plty = 1, caxislabels = c(0.04, 0.06, 0.08, 0.10, 0.12),
           # Customize the grid
           cglcol = "grey", cglty = 1, cglwd = 1.5, calcex=1.8,
           # Customize the axis   Variable labels
           axislabcol = "grey", vlcex = 2)
legend(x = "bottom", legend = rownames(gm1[c(3,4),]), horiz = TRUE,
       bty = "n", pch = 20 , col = c("#00AFBB","#E7B800"),
       text.col = "black", cex = 1.5, pt.cex = 2)



##############################
#Fig 3B
##############################
df=amlnew@meta.data[,c("cyto4","HALLMARK_INTERFERON_GAMMA_RESPONSE")]  ##amlnew is the AML object
df$cyto4=factor(df$cyto4,levels = c("Diploid-Nonmonocytic","Diploid-monocytic","Del5.5q","Del7.7q", "Double deletion"))
p=ggplot(df,aes_string(x="cyto4", y="HALLMARK_INTERFERON_GAMMA_RESPONSE",fill="cyto4")) + ylim(0,0.22) +
  geom_boxplot(width=0.5, color="black", alpha=0.9,outlier.shape = NA) + theme_classic() + #geom_violin() +
  scale_fill_manual(values=c("darkseagreen3","indianred2", "darkslategray", "gold2","deepskyblue"),guide=guide_legend(c('Diploid-Nonmono', 'Diploid-mono',"Del5.5q","Del7.7q", "Double del"))) + ylab("HALLMARK IFNG RESPONSE") +
  scale_x_discrete(labels=c("Diploid-Nonmono","Diploid-mono","Del5.5q","Del7.7q","Double del")) +
  theme(legend.position = "none", axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 22,angle = 45,hjust = 1),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 26),axis.ticks = element_blank()) 
kruskal.test(df$HALLMARK_INTERFERON_GAMMA_RESPONSE~df$cyto4)



##############################
#Fig 3C
##############################
df=amlnew@meta.data[,c("class2","Andy.prediction","HALLMARK_INTERFERON_GAMMA_RESPONSE","hla1","hla2")]
df$Andy.prediction=factor(df$Andy.prediction,levels = c("LSPC-Quiescent","LSPC-Primed","LSPC-Cycle","GMP-like","ProMono-like", "Mono-like","cDC-like"))
gm=groupMeans(t(df[,3:5]),df$Andy.prediction)
gm=as.data.frame(t(gm))
gm=gm[c("LSPC-Quiescent","LSPC-Cycle","LSPC-Primed","GMP-like","ProMono-like", "Mono-like","cDC-like"),]
gm$x=as.character(seq(1,7))
pheatmap(gm[1], color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), border_color = "black", cluster_cols = FALSE, cluster_rows = FALSE,angle_col=90, scale = "column",fontsize_row=12,fontsize_col = 15,cellheight=23, cellwidth = 20,main = "Regulon activity")


##############################
#Fig 3E
##################################regulon activity of 19 TFs from Eagle 2022 BCD paper
cr.tf=c("IRF8","MEF2C","MYB","MEIS1","SPI1","RUNX2","SREBF1","GFI1","ARID2","ZEB2","STAT5B","CEBPA","FOSL2","GFI1B","MEF2D","RUNX1","FLI1","E2F3","LYL1")
cr.tf=cr.tf[cr.tf %in% rownames(regulonActivity_byCellType_Scaled)]
pheatmap(regulonActivity_byCellType_Scaled[cr.tf,], color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), border_color = "black", cluster_cols = FALSE, cluster_rows = TRUE,angle_col=90, scale = "none",fontsize_row=12,fontsize_col = 15,cellheight=20, cellwidth = 20,main = "Regulon activity")






