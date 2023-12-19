library(Seurat);library(ComplexHeatmap);library(symphony);library(pheatmap)
library(dplyr);library(tidydr);library(readxl)
library(data.table);library(stringr)
library(ggplot2);library(circlize)
library(patchwork);library(corrplot);library(SCP)
library(ggcorrplot)
library(rstatix)  # add pvalue position
############################################
#fig1A
#################################
bulk.hallmark=read.csv("Research/manuscript/cell_of_origin/Figure1_new/scoring/TCGA_Abbas_BEAT_ssGSEA_tpm_Bofei.csv")
rownames(bulk.hallmark)=bulk.hallmark$X; bulk.hallmark=bulk.hallmark[,2:12071]
corr.ifng=cor(bulk.hallmark[ ,colnames(bulk.hallmark) != "HALLMARK_INTERFERON_GAMMA_RESPONSE"], bulk.hallmark$HALLMARK_INTERFERON_GAMMA_RESPONSE)
corr.ifng=data.frame(pathway=rownames(corr.ifng),corr=corr.ifng[,1])
path.name=corr.ifng[corr.ifng$corr > 0.7 | corr.ifng$corr < -0.7,"pathway"]
path.name=c(path.name,"HALLMARK_INTERFERON_GAMMA_RESPONSE")
corr1=cor(bulk.hallmark[ ,path.name])

path.name=c("REACTOME_INTERFERON_GAMMA_SIGNALING","REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM", "GOBP_REGULATION_OF_RESPONSE_TO_BIOTIC_STIMULUS","GOBP_RESPONSE_TO_INTERFERON_GAMMA","GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS","GOBP_IMMUNE_RESPONSE","REACTOME_ADAPTIVE_IMMUNE_SYSTEM","HALLMARK_IL6_JAK_STAT3_SIGNALING","GOBP_CELL_ACTIVATION","GOBP_ACTIVATION_OF_IMMUNE_RESPONSE","GOBP_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE", "GOBP_DENDRITIC_CELL_ANTIGEN_PROCESSING_AND_PRESENTATION","GOBP_INFLAMMATORY_RESPONSE", "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION", "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I","REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION","REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION", "GOBP_REGULATION_OF_ANTIGEN_PROCESSING_AND_PRESENTATION","GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY", "REACTOME_ORGANELLE_BIOGENESIS_AND_MAINTENANCE","REACTOME_CILIUM_ASSEMBLY","GOBP_MANNOSYLATION","GOBP_GLUTAMINE_METABOLIC_PROCESS","GOBP_METHYLATION","GOBP_INTRACILIARY_TRANSPORT") 


corr.ifng=corr.ifng[corr.ifng$pathway %in% path.name,]
corr.ifng=corr.ifng[order(corr.ifng$corr,decreasing = FALSE),]
corr.ifng$pathway=factor(corr.ifng$pathway,levels = corr.ifng$pathway)
ggplot(corr.ifng,aes(x=pathway, y=corr,fill=pathway)) + labs(x=NULL, y=NULL) +
  geom_bar(aes(y = corr,x = pathway),position="dodge",stat="identity") + #ylim(0,10.5) + 
  #scale_size(name="NES",breaks=c(1.5,1.75,2,2.25,2.5),limits=c(1,4))+
  #scale_fill_gradientn(colours =colorRampPalette(rev(brewer.pal(n = 7,name ="RdYlBu")))(100), breaks=c(-1,0,1), limits=c(-1,1)) + 
  scale_fill_manual(values = c(rep("steelblue1",6),rep("hotpink1",19))) + guides(fill=FALSE) + theme_bw() +
  theme(axis.text.x = element_text(size=8, angle=0),axis.text.y = element_text(size=10)) + coord_flip()



############################################
# Fig 1B
############################################
ggscatter(bulk.meta, x="HALLMARK.IFNG.RESPONSE", y="Percent_Blasts_in_BM", add="reg.line",conf.int = FALSE,shape=1,
          cor.coef = TRUE,cor.coef.size = 3,cor.method = "pearson", cor.coef.coord =c(0.22,1),
          ylab = "Blast frequency", xlab = "HALLMARK IFNG RESPONSE",size=1) + 
  theme(legend.position = "right", axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))


############################################
# Fig 1c
############################################
bulk.meta$CG_group1=as.character(bulk.meta$CG_group)
bulk.meta$CG_group1[bulk.meta$CG_group1=="Diploid-Nonmono" & bulk.meta$FAB_Simple=="M4"]="Diploid-mono"
bulk.meta$CG_group1[bulk.meta$CG_group1=="Diploid-Nonmono" & !(bulk.meta$FAB_Simple %in% c("M0","M1","M2"))]="Diploid-ONS"
bulk.meta$CG_group2=as.character(bulk.meta$CG_group1)
bulk.meta$CG_group2[bulk.meta$CG_group2=="Other" & bulk.meta$inv16==TRUE]="Inv(16)"
bulk.meta$CG_group2[bulk.meta$CG_group2=="Other" & bulk.meta$t8_21==TRUE]="t(8;21)"
trial=bulk.meta[,c("HALLMARK.IFNG.RESPONSE","CG_group","CG_group1","CG_group2")]
trial=trial[trial$CG_group %in% c("Diploid-mono","Diploid-Nonmono"),]
trial$CG_group1=factor(trial$CG_group1,levels = c("Diploid-ONS","Diploid-Nonmono","Diploid-mono"))
ggplot(trial,aes(x=CG_group1,y=HALLMARK.IFNG.RESPONSE)) + 
  geom_boxplot(aes(fill=CG_group1),width=0.5,color='gray40',outlier.shape=NA,size=0.2)+
  geom_jitter(fill='gray80',color='black',shape=21,width =0.2,size=1.25,stroke = 0.2)+
  ylab("HALLMARK IFNG RESPONSE") + theme_classic() + #ylim(0,0.2) +
  scale_fill_manual(values=c("bisque2","darkseagreen3","indianred2"))+
  theme(legend.position="none",
        axis.title.y = element_text(size = 12),axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10,angle = 30,hjust = 1),axis.text.y = element_text(size = 12))


############################################
# Fig 1F
############################################
#pie chart showing the percentage of patients with abnormal IFNG level at diagnosis
pie(c(29,14),labels = c("Abnormal","Normal"), col=c("#1f77b4","#ff7f0e"))


############################################
# Fig 1G
############################################
corrplot(cor(bulk.meta[,c("HALLMARK.IFNG.RESPONSE","hla1","HLA2","exhaustion_markers","cd8_dysfunction","senescence")]),method = 'circle',type="upper",col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),tl.col = "black") +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 10),legend.text = element_text(size = 10)) 
#colorRampPalette(c("blue","white","red"))(200)

