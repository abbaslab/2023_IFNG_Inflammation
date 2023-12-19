library(survminer);library(survival)
library(survivalAnalysis) ##for multivariable cox model
library(cowplot)  ##for theme_cowplot
library(limma)  ## for lmfit in t statistic
library(tidyverse)
library(ggrepel)
#############################################################################
#figure 5A: scatter plot of ranked correlation
#############################################################################
ifng.cor=read.table("Research/manuscript/cell_of_origin/paper/Figure4/IFNG_gene_cor_SCcells.txt",header = TRUE)
ifng.cor=ifng.cor[2:nrow(ifng.cor),];ifng.cor$rank=c(1:nrow(ifng.cor)) #Exclude IFNG itself, assign the rank
ifng.cor$group=ifelse(ifng.cor$correlation>0,"positive","negative")
ggplot(ifng.cor, aes(x = rank, y = correlation, colour = group)) + geom_point(size=0.5)  + 
  theme_classic() + labs(x="Rank of gene correlation with IFNG score",y="Correlation") +
  theme(panel.grid=element_blank(),legend.position="none",
        axis.title = element_text(size=26),axis.text=element_text(size=28)) +
  scale_color_manual(values = c("positive"="red","negative"="blue")) +
  geom_text(aes(label=ifelse(gene %in% c("CD74","HLA-E","IFITM3","CIITA","IFITM1","IFITM2"),as.character(gene),'')),hjust=-0.1,fontface=3,size=4) +
  geom_text(aes(label=ifelse(gene %in% c("CFD","AZU1","SOX4","SPINK2"),as.character(gene),'')),hjust=1,fontface=3,size=4) +
  #geom_vline(xintercept = 16971,linetype="dashed",size=1) + 
  geom_hline(yintercept = 0,linetype="dashed",size=1) + 
  scale_x_continuous(limits=c(0, 26000),breaks = c(0,5000,10000,15000,20000,25000))



### bottom heatmap
ifng.cor=ifng.cor[!is.na(ifng.cor$correlation),]
ifng <- as.matrix(t(ifng.cor$correlation))
row.names(ifng) <- "ifng.cor"
names(ifng) <- ifng.cor$gene
#colour_bk <- c(#"#F0F0F0",
#colorRampPalette(c("red","#D9EF8B"))(length(ifng.cor[ifng.cor$group=="positive",4])-1),
#"gray90",
#colorRampPalette(c("#D9EF8B","blue"))(dim(ifng.cor)[1]-length(ifng.cor[ifng.cor$group=="positive",4])))
show_gene <- c("IFITM3")
col_anno <-  columnAnnotation(
  show=anno_mark(at= which(names(ifng) %in% show_gene),
                 labels = names(ifng)[which(names(ifng) %in% show_gene)],
                 link_width = unit(1, "mm"),link_height=unit(5, "mm"),labels_gp = gpar(fontsize = 8)))
p=Heatmap(ifng,border=T,show_row_names=F,cluster_columns = F,show_heatmap_legend=F,
          top_annotation =col_anno, width = unit(10, "cm"), height = unit(1, "cm"))#col = colour_bk,



##########################################
#Figure 5d
##########################################
survival=bulk.meta
survival=merge(survival,data.frame(IFITM3=t(bulk.tpm)[,"IFITM3"],Sample=rownames(t(bulk.tpm))),by="Sample")
survival[(!is.na(survival$OS_days) & survival$OS_days> 5*365 ),"Dead_Alive"]="Alive"
survival[(!is.na(survival$OS_days) & survival$OS_days> 5*365 ),"OS_days"]= 5*365
fit=survfit(Surv(OS_days,Dead_Alive == "Dead") ~ IFITM3y>quantile(IFITM3,0.5), survival)
p=ggsurvplot(fit,data=survival,pval=T,palette = c("blue","red"),xlab="Survival Time (Years)",
             xscale=365,xlim=c(0,5*365),break.x.by = 1*365,surv.median.line = "hv",
             legend.title="IFITM3", risk.table = TRUE,#title="BEAT2, >0.5",
             legend.labs = c("Below Median", "Above Median"),pval.coord = c(1000,0.8),pval.size =10) 

##########################################
#Figure 5f
##########################################
############################################################################
#Crispr data show gene effect on specific genes
############################################################################
depmap=read.csv("Research/manuscript/cell_of_origin/paper/Figure4/depmap/CRISPR_gene_effect.csv")
colnames(depmap) <- gsub("\\..*","",colnames(depmap))
depmap1=depmap[,c("DepMap_ID","IFITM3")]
depmap.meta=read.csv("Research/manuscript/cell_of_origin/paper/Figure4/depmap/sample_info.csv")
depmap.meta=depmap.meta[depmap.meta$lineage_subtype=="AML",]
depmap.meta=depmap.meta[,-c(2,4:29)]
depmap1=depmap1[depmap1$DepMap_ID %in% depmap.meta$DepMap_ID,]
depmap1$y="AML"

######celline figure show each cell line as a barplot, rank from low(more death) to high
depmap1=merge(depmap1,depmap.meta[,c("DepMap_ID","stripped_cell_line_name")],by="DepMap_ID")
depmap1=depmap1[order(depmap1$IFITM3..10410.,decreasing = T),]
depmap1$stripped_cell_line_name=factor(depmap1$stripped_cell_line_name)
ggplot(depmap1, aes(y=IFITM3..10410., x=reorder(stripped_cell_line_name,+IFITM3..10410.),fill=y)) + 
  theme_classic() + xlab("") +
  geom_bar(position="stack",stat="identity",alpha=0.8) + #labs(fill='Cell type') +
  scale_fill_manual(values = "steelblue")+
  theme(axis.title.y = element_text(size = 12),axis.text.y =element_text(size=10),
        legend.position = 'none', axis.text.x =element_text(angle = 90,vjust = 0.5,hjust=1,size=10))



