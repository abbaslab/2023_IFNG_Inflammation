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
ggsave("Research/manuscript/cell_of_origin/paper/Figure1/fig1_bulk_IFNG_pathway_corr_barplot.pdf",height=15,width=20, units="cm")




############################################
# Fig 1B
############################################
ggscatter(bulk.meta, x="HALLMARK.IFNG.RESPONSE", y="Percent_Blasts_in_BM", add="reg.line",conf.int = FALSE,shape=1,
          cor.coef = TRUE,cor.coef.size = 3,cor.method = "pearson", cor.coef.coord =c(0.22,1),
          ylab = "Blast frequency", xlab = "HALLMARK IFNG RESPONSE",size=1) + 
  theme(legend.position = "right", axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))
ggsave("Research/manuscript/cell_of_origin/paper/Figure1/fig1_bulk_ifng_percent_blast_cor.pdf",width=10,height =10,unit="cm")


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
ggsave("Research/manuscript/cell_of_origin/paper/Figure1/fig1_bulk_ifng_diploid_box.pdf",width=8,height =10,unit="cm")



############################################
# Fig 1F
############################################
#pie chart showing the percentage of patients with abnormal IFNG level at diagnosis
pdf("Research/manuscript/cell_of_origin/paper/Submission/NatComm/Revision/cytokine_IFNG_abnormal_pie1.pdf",width=8, height=8)
pie(c(29,14),labels = c("Abnormal","Normal"), col=c("#1f77b4","#ff7f0e"))
dev.off()

############################################
# Fig 1G
############################################
pdf("Research/manuscript/cell_of_origin/paper/Figure1/fig1_bulk_IFNG_corr_matrix.pdf", width = 10, height = 10)
corrplot(cor(bulk.meta[,c("HALLMARK.IFNG.RESPONSE","hla1","HLA2","exhaustion_markers","cd8_dysfunction","senescence")]),method = 'circle',type="upper",col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),tl.col = "black") +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 10),legend.text = element_text(size = 10)) 
dev.off()
#colorRampPalette(c("blue","white","red"))(200)










#correlation of HLA2 and IFNG score in bulk 
ggscatter(df, x="HALLMARK.IFNG.RESPONSE", y="hla1", color="CG_group",
          conf.int = TRUE,cor.coef = TRUE,cor.coef.size = 4,cor.method = "pearson",
          ylab = "HLA1 score", xlab = "HALLMARK IFNG RESPONSE",size=1.2) + 
  scale_color_manual(name="CG",values = c("darkseagreen3","indianred2", "darkslategray", "gold2","deepskyblue","mediumpurple1")) + guides(color = guide_legend(override.aes = list(size = 3))) + #size of legend items
  theme(legend.position = "right", axis.title = element_text(size = 9),
        legend.title = element_text(size = 9), legend.text = element_text(size = 8),
        axis.text = element_text(size = 9),axis.ticks = element_blank()) +
  geom_smooth(method = "lm", color = "black")
ggsave("Research/manuscript/cell_of_origin/paper/Figure1/fig1_bulk_IFNG_HLA1_cor.pdf",width=12,height =9,unit="cm")



############################################
#heatmaps of IFNG genes by CG_group
#################################
## ifng-related gene (common.gene) is defined in Figure3
#venn diagram for intersection
ggvenn(list(HALLMARK=hallmark, REACTOME=reactome), 
       fill_color = c("#0073C2FF", "#CD534CFF"),stroke_size = 1, set_name_size = 8,text_size = 8)
ggsave("Research/manuscript/cell_of_origin/paper/Figure1/fig1_IFNG_common_gene.pdf",width=18,height =18,unit="cm")


###heatmap input
mat=bulk.tpm[hla[c(1:6,8:21)],]
mat=log2(mat+1)
group.list=bulk.meta[,c("Sample","CG_group1","CEBPA_mut","DNMT3A_mut","IDH1_mut","IDH2_mut","NPM1_mut","RUNX1_mut","TP53_mut")]
#group.list=group.list[group.list$Sample %in% df$Sample,]
group.list$CG_group1=factor(group.list$CG_group1,levels = c("Diploid-ONS","Diploid-Nonmono","Diploid-mono","Del5.5q","Del7.7q","Double del","Other"))
group.list=group.list[order(match(group.list$CG_group1,c("Diploid-ONS","Diploid-Nonmono","Diploid-mono","Del5.5q","Del7.7q","Double del","Other"))),]
mat=as.matrix(mat[,group.list$Sample])
gm=groupMeans(mat,group.list$CG_group1)
mat=t(scale(t(mat)));gm=t(scale(t(gm)))


bk <- c(seq(-2,-0.1,by=0.02),seq(0,2,by=0.02))
colour_bk <- c(colorRampPalette(c("#2166ac","#d1e5f0"))(83),
               colorRampPalette(c("#d1e5f0","#f7f7f7"))(15),
               colorRampPalette(c("#f7f7f7","#fddbc7"))(15),
               colorRampPalette(c("#fddbc7","#b2182b"))(84))
pdf("Research/manuscript/cell_of_origin/paper/Figure1/fig1_bulk_HLAgene_CG_heat.pdf",width=6,height=7)
Heatmap(gm,name="Expression",col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
        cluster_columns=TRUE, cluster_rows = FALSE, row_names_side = "left",
        row_split = c(rep("",7),rep(" ",13)),row_gap = unit(1.5, "mm"),
        column_title=NULL,show_heatmap_legend = T, show_column_names = TRUE,column_names_rot = 45,
        row_names_gp = gpar(fontsize = 12.5),column_names_gp = gpar(fontsize = 15),
        width = ncol(gm)*unit(10, "mm"), height = nrow(gm)*unit(6, "mm"))
dev.off()


pdf("Research/manuscript/cell_of_origin/paper/Figure1/fig1_bulk_IFNGgene_CG_heat.pdf",width=6,height=10)
pheatmap(gm, cluster_rows=TRUE, cluster_cols=TRUE,cellwidth = 22, cellheight = 18,
         scale="none",show_colnames=T,show_rownames=T, fontsize_row = 12, fontsize_col = 15,
         breaks = bk,color = colour_bk,angle_col=45)
dev.off()
##############
#IFNG and HLA genes by patients gathered by cyto4 
annotation=data.frame(CG=group.list[,2],CEBPA=group.list[,3],DNMT3A=group.list[,4],IDH1=group.list[,5],IDH2=group.list[,6],NPM1=group.list[,7],RUNX1=group.list[,8],TP53=group.list[,9])
rownames(annotation)=group.list[,1]
annotation[annotation==""]="NP"
var=c("darkseagreen3","indianred2", "darkslategray", "gold2","deepskyblue","mediumpurple1")
names(var)=c("Diploid-Nonmono","Diploid-mono","Del5.5q","Del7.7q","Double del","Other")
var1=c("hotpink1","steelblue1","grey")
names(var1)=c("Positive","Negative","NP")
anno_colors=list(CG=var,CEBPA=var1,DNMT3A=var1,IDH1=var1,IDH2=var1,NPM1=var1,RUNX1=var1,TP53=var1)
pdf("Research/manuscript/cell_of_origin/paper/Figure1/fig1e_bulk_hla_CG_patient.pdf",width=12,height=8)
pheatmap(mat, cluster_rows=TRUE, cluster_cols=FALSE, annotation = annotation, annotation_colors=anno_colors,
         scale="none",show_colnames=F,show_rownames=T, fontsize_row = 12, fontsize_col = 15,
         breaks = bk,color = colour_bk,angle_col=45)
dev.off()
rm(mat,group.list,bk,colour_bk,var1,var,anno_colors,annotation)




####pca on IFNG genes
mat=bulk.tpm[,common.gene]
mat=mat[rownames(mat) %in% bulk.meta$Sample,]
pca=prcomp(mat,center = TRUE,scale. = TRUE)
fviz_nbclust(pca$x, kmeans, method = 'wss')
km <- kmeans(as.data.frame(pca$x), 4, nstart = 25)
result=as.data.frame(pca$x);result$Sample=rownames(result)
result$kmeansCluster=km$cluster
result=merge(result,bulk.meta[,c("Sample","CG_group","Cohort")])

ggplot(result, aes_string(x = "PC1", y = "PC2", colour = "CG_group")) + geom_point(size=2)  + 
  xlab("PC1") + ylab("PC 2") +theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=6))) + 
  theme(legend.text=element_text(size=28), legend.title=element_text(size=26),
        axis.title=element_text(size=27),axis.text=element_text(size=28)) +
  scale_colour_manual(values = c("darkseagreen3","indianred2", "darkslategray", "gold2","deepskyblue","turquoise4"),name="Cluster")


table(result$kmeans,result$CG_group)



############################################
# ifng score by CG group, add healthy, test run
#################################
bulk.meta$CG_group[bulk.meta$Healthy_Disease !="D"]="Healthy"
temp=as.data.frame(t(bulk.tpm[c("CD74","IFITM3"),]))
temp$Sample=rownames(temp)
df=merge(bulk.meta[,c("Sample","CG_group","HALLMARK.IFNG.RESPONSE","HLA2")],temp,by="Sample")
df$CD74=log2(df$CD74+1);df$IFITM3=log2(df$IFITM3+1)



df$CG_group=factor(df$CG_group,levels = c("Healthy",'Diploid-Nonmono', 'Diploid-mono',"Del5.5q","Del7.7q","Double del","Other"))
kruskal.test(df$HALLMARK.IFNG.RESPONSE~df$CG_group)
ggplot(df,aes(x=CG_group,y=HALLMARK.IFNG.RESPONSE))+#scale_y_continuous(trans=squish_trans(0,6),breaks = seq(6,15))+
  geom_boxplot(aes(fill=CG_group),width=0.5,color='gray40',outlier.shape=NA,size=0.2)+
  geom_jitter(fill='gray80',color='black',shape=21,width =0.2,size=1.25,stroke = 0.2)+
  ylab("HALLMARK IFNG RESPONSE") + theme_classic() + #ylim(0,0.2) +
  geom_vline(xintercept = 2.5,linetype="dashed",size=1) + 
  #ggtitle("BULK exclude Other in CG_group, and NA/CMML-T/RAEB in FAB of diploid pts") + #facet_grid(.~cg2) +
  #scale_fill_manual(values=c("darkseagreen3","indianred2", "darkslategray", "gold2","deepskyblue","mediumpurple1"))+
  theme(legend.position="none",
        axis.title.y = element_text(size = 12),axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12,angle = 30,hjust = 1),axis.text.y = element_text(size = 12))
ggsave("Research/manuscript/cell_of_origin/paper/Figure1/fig1_bulk_cg_IFNG_box_healthy.pdf",height=15,width=12, units="cm")






############################################
#A heatmap of bulk patients ordered by increasing IFNG score and show clinical characteristics below
#################################
bulk.cg=bulk.meta[order(bulk.meta$HALLMARK.IFNG.RESPONSE,decreasing = FALSE),]
#bulk.cg=bulk.cg[bulk.cg$CG_group!="Other",]
bulk.cg$del17.17p="Negative"
bulk.cg$del17.17p[bulk.cg$del17p==TRUE | bulk.cg$loss17==TRUE]="Positive"
bulk.cg$median_IFNG=ifelse(bulk.cg$HALLMARK.IFNG.RESPONSE>median(bulk.cg$HALLMARK.IFNG.RESPONSE),"Above","Below")

plot.data=bulk.cg[,c("CG_group","FAB_Simple","ELN2017_incomplete","FLT3_mut","NPM1_mut","RUNX1_mut","TP53_mut","del17.17p","HALLMARK.IFNG.RESPONSE")] #"Is_DeNovo","Age","MLL_F","Percent_Blasts_in_BM"
#plot.data[plot.data$Age > 60 & !is.na(plot.data$Age),'Age']="Elder"
#plot.data[plot.data$Age <= 60 & !is.na(plot.data$Age),'Age']="Young"
plot.data[plot.data==""]="NP";plot.data[is.na(plot.data)]="NP"
plot.data$FAB_Simple[plot.data$FAB_Simple %in% c("RAEB","CMML-T","NP")]="Other"
plot.data$ELN2017_incomplete[plot.data$ELN2017_incomplete %in% c("IntermediateOrAdverse","NonAML","FavorableOrIntermediate","NotInitial","NP")]="Other"
plot.data$FAB_Simple=factor(plot.data$FAB_Simple,levels = c("M0","M1","M2","M4","M5","Other"))
plot.data$CG_group=factor(plot.data$CG_group,levels = c('Diploid-Nonmono', 'Diploid-mono','Del5.5q','Del7.7q',"Double del","Other"))
plot.data$ELN2017_incomplete=factor(plot.data$ELN2017_incomplete,levels = c('Favorable', 'Intermediate','Adverse',"Other"))

colnames(plot.data)[1:7]=unlist(lapply(colnames(plot.data[1:7]), function(x){str_split(x,"_")[[1]][1]})) 
colnames(plot.data)[c(8,9)]=c("Del17.17p","HALLMARK IFNG RESPONSE")

lgd1 = Legend(labels = c("Mutation", "No Mutation", "Missing"), legend_gp = gpar(fill = c("deepskyblue", "white", "lightgray")), title = "Mutation")
lgd2 = Legend(labels=c("M0","M1","M2","M4","M5","Other"), legend_gp = gpar(fill=c(brewer.pal(n = 10, name = 'Set3')[c(1:3,5,4)],"lightgray")), title="FAB")
lgd3=Legend(labels=c('Diploid-Nonmono', 'Diploid-mono','Del5.5q','Del7.7q',"Double del","Other"), legend_gp = gpar(fill=c("deeppink","gray28","chartreuse1", "gold2","deepskyblue","lightgray")), title="CG")
#lgd4=Legend(labels=c('DeNovo', 'Secondary'), legend_gp = gpar(fill=c("chartreuse","brown2")), title="Is DeNovo")
#lgd5=Legend(labels=c('Young', 'Elder'), legend_gp = gpar(fill=c("red","blue")), title="Age")
lgd6=Legend(labels=c("Adverse","Favorable","Intermediate","Other"),legend_gp = gpar(fill=c("#8DD3C7","#BEBADA","#80B1D3","lightgray")), title="ELN2017")
pd = packLegend(lgd3,lgd2,lgd6, lgd1, direction = "vertical", gap = unit("1", "cm"))
#column_ha <- HeatmapAnnotation(IFNG_score = anno_barplot(as.numeric(bulk.cg$HALLMARK_INTERFERON_GAMMA_RESPONSE),
#                                bar_width = 0.7, border = FALSE, axis = TRUE,
#                                axis_param = list(gp = gpar(fontsize = 20, lwd = 0.5)),
#                                gp = gpar(col = NA, fill = "grey"), show_annotation_name = F, height = unit(4, "cm")))
column_ha=HeatmapAnnotation(df=data.frame(
  IFNG_level=factor(c(rep("Low IFNG",336),rep("High IFNG",336)),levels = c("Low IFNG","High IFNG"))), annotation_name_side = "left",annotation_name_gp= gpar(fontsize = 10),show_legend = c(FALSE),
  col = list(IFNG_level=c("Low IFNG"="#8DD3C7","High IFNG"="#FB8072")),height=unit(600, "mm"))


pdf("Research/manuscript/cell_of_origin/paper/Figure1/fig1c_bulk_meta2.pdf", width = 12, height = 8)
Heatmap(t(plot.data[,1:8]), col = c("Positive" = "deepskyblue", "Negative"="white", "NP" = "lightgray",
                                    "M0"="#8DD3C7","M1"="#FFFFB3","M2"="#BEBADA","M4"="#80B1D3","M5"="#FB8072","Other"="lightgray",
                                    'Diploid-Nonmono'="deeppink", 'Diploid-mono'="gray28",'Del5.5q'="chartreuse1",'Del7.7q'="gold2","Double del"="deepskyblue",
                                    "Adverse"="#8DD3C7","Favorable"="#BEBADA","Intermediate"="#80B1D3"),
        row_title=" ",row_names_side = "left", border=FALSE, width=unit(7, "in"),column_title="bulk 672",
        column_split = c(rep("",336),rep(" ",336)),column_gap = unit(1, "mm"),show_column_names = FALSE,
        height=unit(4,"in"),  show_heatmap_legend = FALSE,top_annotation = column_ha)
draw(pd, just=c("right"),x=unit(30, "cm"))
dev.off()
rm(bulk.cg,plot.data,lgd1,column_ha,lgd2,pd,lgd3,lgd6)







############################################
# validate IFNG in nanostring (STM 2020)
#################################
nano.expri=read_excel("Research/AML/info_files/Nanostring_Sergio/Nanastring_Vadak2020_STM/nanostring_expression.xlsx")
nano.meta=read_excel("Research/AML/info_files/Nanostring_Sergio/Nanastring_Vadak2020_STM/nanostring_meta.xlsx")
nano.cyto=read.csv("Research/AML/info_files/Nanostring_Sergio/Nanastring_Vadak2020_STM/nanostring_cytogenetics.csv")
nano.expri=as.data.frame(nano.expri);rownames(nano.expri)=nano.expri$Patient_Id
nano.expri=nano.expri[,2:731]
nano.expri=as.data.frame(t(nano.expri))
##score HLA2 and hallmark ifng


##code to generate CG_group is the same as bulk.meta
nano.cyto$CG_group=nano.cyto$diploid
nano.cyto$CG_group[nano.cyto$CG_group==TRUE]="Diploid"
nano.cyto$CG_group[nano.cyto$del7p==TRUE | nano.cyto$del7q==TRUE | nano.cyto$loss7==TRUE]="Del7.7q"
nano.cyto$CG_group[nano.cyto$del5p==TRUE | nano.cyto$del5q==TRUE | nano.cyto$loss5==TRUE]="Del5.5q"
nano.cyto$CG_group[(nano.cyto$del5p==TRUE | nano.cyto$del5q==TRUE | nano.cyto$loss5==TRUE) & (nano.cyto$del7p==TRUE | nano.cyto$del7q==TRUE | nano.cyto$loss7==TRUE)]='Double deletion'
nano.cyto[nano.cyto$CG_group==FALSE,"CG_group"]="Other"
nano.cyto=nano.cyto[,-c(1,3:5)]
###merge meta data with CG group
nano.meta=merge(nano.meta,nano.cyto[,c("Patient_Id","CG_group")],by="Patient_Id")
nano.meta$CG_group[nano.meta$CG_group=="Diploid" & nano.meta$FAB %in% c("M4/M5","M5","M5a","M5b")]="Diploid-mono"
nano.meta$CG_group[nano.meta$CG_group=="Diploid"]="Diploid-Nonmono"
nano.meta=nano.meta[nano.meta$Age>=18,]
nano.meta=nano.meta[nano.meta$Relapse != "Yes",]

kruskal.test(nano.meta$HALLMARK_IFNG~nano.meta$CG_group)
ggplot(nano.meta,aes(x=CG_group,y=HALLMARK_IFNG)) + 
  geom_boxplot(aes(fill=CG_group),width=0.5,color='gray40',outlier.shape=NA,size=0.2)+
  geom_jitter(fill='gray80',color='black',shape=21,width =0.2,size=1.25,stroke = 0.2)+
  ylab("HALLMARK IFNG RESPONSE") + theme_classic() + #ylim(0,0.2) +
  geom_vline(xintercept = 3.5,linetype="dashed",size=1) + 
  scale_fill_manual(values=c("turquoise4","bisque2","indianred2","chartreuse1","gold2","deepskyblue"))+
  theme(legend.position="none",
        axis.title.y = element_text(size = 24),axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20,angle = 30,hjust = 1),axis.text.y = element_text(size = 25))





############################################
# Fig1D percent of T, B amd AML cells, correlation with IFNG, HLA-E, Age, Phase score at patient level
#################################
allcell.new=readRDS("/Users/bwang8/Research/AML/AML_subsets/allcell_new_noharmony.RDS")
tcell=readRDS("/Users/bwang8/Research/manuscript/poonam_CD8Tcell_paper/Tcells.rds")
allcell.new$orig.ident=factor(allcell.new$orig.ident,levels = levels(amlnew$orig.ident))
tcell=subset(tcell,subset=orig.ident %in% unique(amlnew$orig.ident))
tcell$orig.ident=factor(tcell$orig.ident,levels = levels(amlnew$orig.ident))
aml.frac=table(amlnew$orig.ident)/table(allcell.new$orig.ident)
composition=data.frame(PTs=names(aml.frac),aml_frac=as.numeric(aml.frac))
composition$T_frac=as.numeric(table(tcell$orig.ident)/table(allcell.new$orig.ident))
#composition$B_frac=as.numeric(table(allcell.new$orig.ident,allcell.new$class1)[,"B"]/table(allcell.new$orig.ident))
composition$B_frac=as.numeric(c(339,539,295,73,61,86,0,86,174,231,139,229,98,389,168,454,126,213,132,33)/table(allcell.new$orig.ident))
composition$Age=round(metadata$Age_at_Colln)
score=amlnew@meta.data[,c("hla2","HALLMARK_INTERFERON_GAMMA_RESPONSE","S.Score","G2M.Score")]
score.mean=groupMeans(t(score),amlnew$orig.ident)
composition=cbind(composition,as.data.frame(t(score.mean)))
colnames(composition)=c("PTs","AML_frac","Tcell_frac","Bcell_frac","Age","HLA2 score","HALLMARK_IFNG","S.Score","G2M.Score")
p=ggcorrplot(cor(composition[,2:ncol(composition)]),method = 'square',type="full",tl.cex = 15) +
  theme(axis.text.x = element_text(size = 22,angle = 45,hjust = 0.9),
        legend.title = element_text(size = 23),legend.text = element_text(size = 22),
        axis.text.y = element_text(size = 22),axis.ticks = element_blank()) 
pdf("Research/manuscript/cell_of_origin/paper/Figure1/fig1d_cor_martrix.pdf", width = 12, height = 12)
print(p)
dev.off()
rm(metadata,aml.cell,score,score.mean,aml.frac,allcell.new,tcell,p)































