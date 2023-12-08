############################################
#process bulk metadata
bulk.meta=read.csv("/Users/bwang8/Downloads/MetaData_All_v4.csv")
bulk.meta=bulk.meta[bulk.meta$Cohort %in% c("TCGA","Abbas","BEAT1","BEAT2"),] ##choose dataset
bulk.meta=bulk.meta[bulk.meta$Adult %in% c(TRUE,NA),] ##choose adult
bulk.meta=bulk.meta[bulk.meta$Healthy_Disease=="D",] ##choose only AML patients
bulk.meta=bulk.meta[bulk.meta$At_Diagnosis %in% c(TRUE,NA),] ##choose only newly diagnosed patients
bulk.meta=bulk.meta[!bulk.meta$FAB_Simple %in% c("M3","M6","M7"),] ##choose only newly diagnosed patients

bulk.meta$CG_group=bulk.meta$diploid
bulk.meta$CG_group[bulk.meta$CG_group==TRUE]="Diploid"
bulk.meta$CG_group[bulk.meta$del7p==TRUE | bulk.meta$del7q==TRUE | bulk.meta$loss7==TRUE]="Del7.7q"
bulk.meta$CG_group[bulk.meta$del5p==TRUE | bulk.meta$del5q==TRUE | bulk.meta$loss5==TRUE]="Del5.5q"
bulk.meta$CG_group[(bulk.meta$del5p==TRUE | bulk.meta$del5q==TRUE | bulk.meta$loss5==TRUE) & (bulk.meta$del7p==TRUE | bulk.meta$del7q==TRUE | bulk.meta$loss7==TRUE)]='Double del'
bulk.meta[bulk.meta$CG_group==FALSE,"CG_group"]="Other"

#bulk.meta$Sample[448,600]=gsub("\\.","\\-",bulk.meta$Sample[448,600])
bulk.meta$Sample[448:600]=substr(bulk.meta$Sample[448:600],1,12)


############################################
# histogram showing the distribution of IFNG scores
############################################
ggplot(bulk.meta, aes(x = HALLMARK.IFNG.RESPONSE)) + theme_classic() + geom_density(linewidth=0.8) +
  geom_histogram(fill="darkgrey",aes(y = ..density..),binwidth = 0.003297811) + 
    geom_vline(xintercept = min(bulk.meta$HALLMARK.IFNG.RESPONSE)+IQR(bulk.meta$HALLMARK.IFNG.RESPONSE),color="red",linetype="dashed",size=0.8) +
  geom_vline(xintercept = max(bulk.meta$HALLMARK.IFNG.RESPONSE)-IQR(bulk.meta$HALLMARK.IFNG.RESPONSE),color="red",linetype="dashed",size=0.8) +
  theme(axis.title = element_text(size = 12),axis.text =element_text(size=10))
ggsave("Research/manuscript/cell_of_origin/paper/FigureS/figsupp_bulk_IFNG_histogram.pdf",height=10,width=12, units="cm")

ggplot(bulk.meta,aes(x=Counts_TPM, y=HALLMARK.IFNG.RESPONSE,fill=Counts_TPM))+
  geom_boxplot(width=0.5, color="black", alpha=0.9,outlier.stroke=0.2,lwd=0.2) +
  scale_fill_manual(values=c("indianred2")) + xlab("AML patients") +
  ylab("HALLMARK IFNG RESPONSE") + theme_classic() + 
  theme(legend.position = "none", axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),axis.title.x = element_blank(),axis.ticks = element_blank())
ggsave("Research/manuscript/cell_of_origin/paper/FigureS/figsupp_bulk_IFNG_box.pdf",height=10,width=8, units="cm")



############################################
# patient contribution in each CG group in AML cells
############################################
meta.by.cyto=amlnew@meta.data;meta.by.cyto$cyto4=as.character(meta.by.cyto$cyto4)
meta.by.cyto$cyto4=factor(meta.by.cyto$cyto4,levels=c("Diploid-Nonmonocytic",'Diploid-monocytic', 'Del5.5q', 'Del7.7q', "Double deletion"))
phenoPropotion_flip(meta.by.cyto,x_name="cyto4",y_name="orig.ident",out_path="Research/manuscript/cell_of_origin/paper/FigureS/",plot_width=9,plot_height=5.9,legend_name="Patient")

phenoPropotion_flip <- function(data=patient_obs,x_name,y_name,legend_name,out_path,plot_height=6.2,plot_width){
  df_stats <- data %>% 
    dplyr::group_by(get(x_name),get(y_name),.drop = FALSE) %>% 
    dplyr::summarise(n=n())
  names(df_stats)[1:2] <- c(x_name,y_name)
  df_stats <- na.omit(df_stats)
  df_stats <- as.data.frame(df_stats)
  df_stats[,y_name] <-  factor(as.character(df_stats[,y_name]),levels=levels(data[,y_name]))#c('Diploid', 'Del5.5q', 'Del7.7q', "Double del")
  
  ggplot(df_stats, aes(x = get(x_name), y = n, fill = get(y_name))) + theme_classic()+ 
    geom_bar(position = "fill",stat = "identity",na.rm=T,width=0.8) + #position="stack" gives numbers
    scale_fill_manual(values= as.character(ColAssign(letters[1:20]))) +
    theme(legend.position="top",
          legend.title=element_text(size = 8),legend.text=element_text(size = 6),
          legend.key.height = unit(0.25,"cm"),legend.key.width = unit(0.25,"cm"),
          axis.title.x = element_text(size = 8),axis.text.x = element_text(size = 8),
          axis.title.y = element_text(size = 10),axis.text.y = element_text(size = 8))+
    xlab("")+ylab("Patient cell Fraction (%)")+labs(fill=legend_name)+
    scale_y_continuous(labels = scales::percent)+
    coord_flip()
  ggsave(file.path(out_path,paste0(x_name,"_",y_name,"_percent.pdf")),height=plot_height,width=plot_width,unit="cm") 
  ggsave(file.path(out_path,paste0(x_name,"_",y_name,"_percent.png")),height=plot_height,width=plot_width,unit="cm") 
}


############################################
# cells in each patient, count not percentage
#################################
meta=allcell@meta.data[,c("orig.ident","cyto4","class2")]
fraction=meta %>% group_by(orig.ident,class2)  %>% summarise(count=n()) 
fraction=as.data.frame(fraction)
fraction$orig.ident=factor(fraction$orig.ident,levels = c('PT20A',"PT23A","PT25A","PT12A","PT19A","PT30A","PT32A","PT9A", "PT13A","PT17A","PT22A","PT26A","PT14A","PT15A","PT21A","PT27A","PT28A","PT10A","PT16A","PT29A"))
fraction$class2=factor(fraction$class2,levels = c("AML","Progenitor","T","B","NK","Monocytes","Erythroids","DC","Other"))
ggplot(fraction, aes(fill=class2, y=count, x=orig.ident)) + theme_classic() + xlab("") +
  geom_bar(position="stack",stat="identity") + labs(fill='Cell type') +
  scale_fill_manual(values = as.character(ColAssign(letters[1:10])))+
  theme(plot.title = element_text(hjust = 0.5),axis.title.y = element_text(size = 12),
        axis.text.y =element_text(size=10), axis.text.x =element_text(angle = 90,vjust = 0.5,size=10))
ggsave("Research/manuscript/cell_of_origin/paper/FigureS/figsupp_allcell_pt_count.pdf",height=10,width=12, units="cm")



############################################
# pie chart of cell types in single cell
############################################
meta=allcell@meta.data
fraction=meta %>% group_by(class2)  %>% summarise(count=n()) %>% mutate(Percent = round(count/sum(count),4))
fraction=as.data.frame(fraction)
fraction$class2=factor(fraction$class2,levels = c("AML","Progenitor","T","B","NK","Monocytes","Erythroids","DC","Other"))
df <- fraction %>% mutate(csum = rev(cumsum(rev(count))), 
                          pos = count/2 + lead(csum, 1),pos = if_else(is.na(pos), count/2, pos))

ggplot(fraction, aes(x="", y=count, fill=class2))  + theme_void() + 
  coord_polar("y", start=0,direction = -1) + scale_fill_manual(values = as.character(ColAssign(letters[1:10]))) +      
  geom_bar(stat="identity", width=1, color="white") +
  theme(legend.position="right",legend.text=element_text(size=15), legend.title=element_blank()) +   
  geom_label_repel(data = df,aes(y = pos, label = paste0(Percent *100, "%")),size = 3, nudge_x = c(1,0.7,0,1,0,0,0,0.8,0), show.legend=FALSE) 
ggsave("Research/manuscript/cell_of_origin/paper/FigureS/figsupp_allcell_percentage.pdf",height=10,width=12, units="cm")

rm(meta,fraction,df)


############################################
# ncount violin plot by cell types 
############################################
df=allcell@meta.data[,c('nCount_RNA',"nFeature_RNA","mitoRatio","riboRatio","S.Score","G2M.Score","class2")]
df$class2=factor(df$class2,levels = c("AML","Progenitor","T","B","NK","Monocytes","Erythroids","DC","Other"))
ggplot(df,aes(x=class2,y=G2M.Score,fill=class2)) + theme_classic()+
  geom_violin(size=1) + scale_fill_manual(values=as.character(ColAssign(letters[1:10]))) + 
  theme(legend.position="none",title=element_blank(),axis.title = element_blank(),
        axis.text.x = element_text(size = 10,angle = 45,hjust = 1),axis.text.y = element_text(size = 12))+ 
  stat_summary(fun.y=mean, geom="point", shape=18,size=2, color="red")
ggsave("Research/manuscript/cell_of_origin/paper/FigureS/figsupp_allcell_g2mscore.pdf",height=10,width=10, units="cm")

ggplot(df,aes(x=class2,y=mitoRatio,fill=class2)) + theme_classic()+
  ggtitle("Mitochondrial Ratio") + ylim(0,1) +
  geom_violin(size=1) + scale_fill_manual(values=as.character(ColAssign(letters[1:10]))) + 
  theme(legend.position="none",axis.title = element_blank(),
        axis.text.x = element_text(size = 10,angle = 45,hjust = 1),axis.text.y = element_text(size = 12))+ 
  stat_summary(fun.y=mean, geom="point", shape=18,size=2, color="red")
ggsave("Research/manuscript/cell_of_origin/paper/FigureS/figsupp_allcell_mitoRatio.pdf",height=10,width=10, units="cm")





############################################
# AML cell count in each cytogenetic group
############################################
meta=amlnew@meta.data[,c("orig.ident","cyto4","class2")]
fraction=meta %>% group_by(cyto4)  %>% summarise(count=n()) 
fraction=as.data.frame(fraction); fraction=fraction[c(2,1,3:5),]
ggplot(fraction, aes(fill=cyto4, y=count, x=cyto4)) + theme_classic() + 
  geom_bar(position="stack",stat="identity") + xlab("") +
  scale_fill_manual(values = c("darkseagreen3","indianred2", "darkslategray", "gold2","deepskyblue"))+
  theme(axis.title.y = element_text(size = 12),legend.position = "none",
        axis.text.y =element_text(size=10), axis.text.x =element_text(angle = 45,vjust = 0.5,size=10))
ggsave("Research/manuscript/cell_of_origin/paper/FigureS/figsupp_aml_count_cyto4.pdf",height=10,width=8, units="cm")




############################################
# Lasry nature cancer AML single cell data, validate interferon gamma
#Larsy has 20 diagnistic adult AML patients, 2 diploid-mono, 6 diploid-nonmono 
############################################
larsy=readRDS("Research/AML/info_files/public_scAML/Lasry2022_adult.rds")
#larsy=subset(larsy,subset=Cell_type_identity=="malignant")
df=larsy@meta.data[,c("donor_id","malignant","HALLMARK_INTERFERON_GAMMA_RESPONSE","hla1","hla2","Broad_cell_identity","Cell_type_identity")]
df=df[df$donor_id %in% c("AML0160","AML2910","AML2123","AML3133","AML0361","AML4340","AML4897","AML2451"),]
df$cyto4=ifelse(df$donor_id %in% c("AML0160","AML2910"), 'Diploid-mono',"Diploid-Nonmono")
df$cyto4=factor(df$cyto4,levels = c("Diploid-Nonmono",'Diploid-mono'))
ggplot(df,aes_string(x="cyto4", y="HALLMARK_INTERFERON_GAMMA_RESPONSE",fill="cyto4"))+#geom_violin() +
  geom_boxplot(width=0.6, color="black", alpha=0.9,outlier.shape = NA,outlier.stroke=1,lwd=0.2) +
  scale_fill_manual(values=c("darkseagreen3","indianred2")) +
  ylab("HALLMARK IFNG RESPONSE") + theme_classic() + ylim(0,0.3) +
  #scale_x_discrete(labels=c('Diploid', 'Del5.5q', 'Del7.7q', "Double del"))+
  theme(legend.position = "none", axis.title = element_text(size = 8),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 10,angle = 30,hjust = 1),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),axis.ticks = element_blank()) + facet_grid(.~malignant)
ggsave("Research/manuscript/cell_of_origin/paper/FigureS/figsupp_larsy__IFNG_box.pdf",height=10,width=8, units="cm")

gm=df %>% group_by(malignant,cyto4) %>% summarise(HALLMARK=mean(HALLMARK_INTERFERON_GAMMA_RESPONSE),.groups = 'drop')
gm=as.data.frame(gm)
ggplot(data=gm, aes(x=cyto4, y=HALLMARK, group=malignant)) + ylab("Average IFNG score") +
  geom_line(aes(color=malignant), size=0.5) + theme_classic() +
  geom_point(aes(shape=malignant), size=0.8)+scale_linetype_manual(values=c("red2","cornflowerblue")) +
  theme(legend.position = "top", axis.title = element_text(size = 8),
        legend.text = element_text(size=7), legend.title=element_blank(),
        axis.text.x = element_text(size = 7,angle = 0,hjust = 0.5),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8),axis.ticks = element_blank()) 
ggsave("Research/manuscript/cell_of_origin/paper/FigureS/figsupp_larsy_average_line.pdf",height=6,width=8, units="cm")





###############################petti diploid
petti.score=petti@meta.data[,c("orig.ident","hallmarkIFN")]
for (i in 1:nrow(petti.score)){
  if (startsWith(petti.score[i,"orig.ident"],"452198_P")){petti.score[i,"cyto2"]="Diploid-mono"}
  else if (petti.score[i,"orig.ident"] %in% c("220882_P_1","220882_P_2","823477_P_1","823477_P_2","869586_P_1","869586_P_2", "508084_P_1", "508084_P_2" )){petti.score[i,"cyto2"]="Diploid-Nonmono"}
  else {petti.score[i,"cyto2"]="Other"}
} 
petti.score=petti.score[petti.score$cyto2!="Other",]
petti.score$cyto2=factor(petti.score$cyto2,levels = c("Diploid-Nonmono","Diploid-mono"))
p=ggplot(petti.score, aes_string(x="hallmarkIFN",fill="cyto2")) + xlim(0,0.28)+
  geom_density(alpha=0.4)+ theme_classic() + ylab("Density")+ labs(fill="Group")+
  scale_fill_manual(values=c("darkseagreen3","indianred2"))+ xlab("HALLMARK IFNG RESPONSE")+
  theme(axis.text.x = element_text(size = 30),legend.position = 'top',axis.title.x=element_text(size=24),
        axis.text.y = element_text(size = 30, angle = 0),axis.title.y=element_text(size=28),legend.text = element_text(size = 24),legend.title = element_text(size = 28))




#############################################################
#Malani Cancer Discovery paper, score RNA-seq, correlate with venetoclax resistance
# and across FAB classification
############################################################
drug.response=read_excel("Research/AML/info_files/public_AMLdata/Malini2022/Functional_Precision_Medicine_Tumor_Board_AML/File_3_Drug_response_sDSS_164S_17Healthy.xlsx")
meta=read_excel("/Users/bwang8/Research/AML/info_files/public_AMLdata/Malini2022/Functional_Precision_Medicine_Tumor_Board_AML/File_0_Common_sample_annotation_252S.xlsx")
rna.seq=read.csv("/Users/bwang8/Research/AML/info_files/public_AMLdata/Malini2022/malani_tpm_Bofei.csv")
parsi.gene=read.table("Research/manuscript/cell_of_origin/paper/Figure6/HALLMARK_IFNG_gene_lasso_coef.txt",header=TRUE)
parsi.gene=parsi.gene[parsi.gene$Coefficient !=0,]
clinical=read_excel("/Users/bwang8/Research/AML/info_files/public_AMLdata/Malini2022/clinical_data.xlsx")

####clean the data
meta=meta[meta$`Disease status`=="Diagnosis",]
meta$FAB_subtype[meta$FAB_subtype %in% c("FAB M4","FAB M4 eos","FAB M5","FAB M4; FAB M5")]= "Monocytic"
meta$FAB_subtype[meta$FAB_subtype %in% c("FAB M0","FAB M1","FAB M2","FAB M1; FAB M0")]= "Primitive"
meta$FAB_subtype[is.na(meta$FAB_subtype)]= "NOS"


rna.seq=distinct(rna.seq,Genes,.keep_all = TRUE)
rownames(rna.seq)=rna.seq$Genes
rna.seq=rna.seq[,-1]
rna.seq=rna.seq[complete.cases(rna.seq), ]

IFITM3=as.data.frame(t(rna.seq[c("IFITM3"),]));IFITM3$Sample_ID=rownames(IFITM3)
IFITM3$IFITM3=log2(IFITM3$IFITM3 + 1)
# 
# gsva_score=gsva(as.matrix(rna.seq), list(HALLMARK_IFNG=hallmark$HALLMARK_INTERFERON_GAMMA_RESPONSE),method="ssgsea",abs.ranking=F)
# gsva_score=as.data.frame(t(gsva_score))
# gsva_score$Sample_ID=rownames(gsva_score)
gsva_score=read.csv("/Users/bwang8/Research/manuscript/cell_of_origin/paper/Figure6/Malani_IFNG_score_BofeiTPM.csv")
gsva_score=merge(gsva_score,IFITM3,by="Sample_ID")
###Drug resistance correlation
drug.response=as.data.frame(drug.response)
venetoclax=drug.response[drug.response$Drug_name=="Venetoclax",]
venetoclax=as.data.frame(t(venetoclax))
colnames(venetoclax)=venetoclax[2,]
venetoclax$Sample_ID=rownames(venetoclax)
venetoclax=venetoclax[3:183,]
temp=merge(gsva_score,venetoclax,by="Sample_ID")
temp$Venetoclax=as.numeric(temp$Venetoclax)
temp=temp[!is.na(temp$Venetoclax),]
temp=temp[temp$Sample_ID %in% meta$Sample_ID,]
ggscatter(temp, x = "HALLMARK_IFNG", y = "Venetoclax", add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "IFNG signaling score", ylab = "Venetoclax sensitivity",size=2,
          cor.coef.size=7) + #,cor.coef.coord=c(4000,300)
  theme(legend.position = "none", axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11),
        axis.ticks = element_blank()) 
ggsave("Research/manuscript/cell_of_origin/paper/FigureS/figsupp_IFNG_Venetoclax_corr_Malani.pdf",height=8,width=8)






ifng.fab=merge(gsva_score, meta,by="Sample_ID",all.x=TRUE)
ifng.fab$FAB_subtype[ifng.fab$Sample_ID %in% c("Healthy1_CD34P","Healthy2_CD34P","Healthy3_CD34P","Healthy4_CD34P")]="CD34+ Healthy"
ifng.fab=ifng.fab[!is.na(ifng.fab$FAB_subtype),]
ifng.fab$FAB_subtype=factor(ifng.fab$FAB_subtype,levels = c("CD34+ Healthy","Primitive","Monocytic","NOS"))
ggplot(ifng.fab,aes(x=FAB_subtype ,y=HALLMARK_IFNG)) + 
  geom_boxplot(aes(fill=FAB_subtype ),width=0.5,color='gray40',outlier.shape=NA,size=0.2)+
  geom_jitter(fill='gray80',color='black',shape=21,width =0.2,size=1.25,stroke = 0.2)+
  ylab("HALLMARK IFNG RESPONSE") + theme_classic() + #ylim(0,0.2) +
  scale_fill_manual(values=c("lightcyan","darkseagreen3","indianred2","mediumpurple1"))+
  theme(legend.position="none",
        axis.title.y = element_text(size = 12),axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10,angle = 30,hjust = 1),axis.text.y = element_text(size = 12))
ggsave("Research/manuscript/cell_of_origin/paper/FigureS/figsupp_IFNG_fab_box_Malani.pdf",width=8,height =10,unit="cm")


####evaluate the survival by ifng or parsimonious scores
gsva_score$Patient_ID=substr(gsva_score$Sample_ID,1,7)
gsva_score$timepoint=substr(gsva_score$Sample_ID,9,10)
survival=merge(gsva_score,clinical[,c("Patient_ID","os_event","os_time")],by="Patient_ID")
survival=survival[survival$Sample_ID %in% meta$Sample_ID,]
survival$os_event[survival$os_time/365 > 5]=0
survival$os_time=ifelse(survival$os_time/365 > 5,5*365,survival$os_time)
survival=survival[survival$Sample_ID %in% meta$Sample_ID,]

fit=survfit(Surv(os_time,os_event == 1) ~ HALLMARK_IFNG_parsi>quantile(HALLMARK_IFNG_parsi,0.5), survival)
p=ggsurvplot(fit,data=survival,pval=T,palette = c("blue","red"),xlab="Survival Time (Days)",
             xscale=365,xlim=c(0,5*365),break.x.by = 1*365,surv.median.line = "hv",
             legend.title="Parsimonious IFNG", risk.table = FALSE,#title="BEAT2, >0.5",
             legend.labs = c("Below Median", "Above Median"),pval.size =6) 

p$plot+theme(axis.text.x=element_text(size=10),axis.title.y = element_text(size = 12),
             axis.text.y=element_text(size = 10),axis.title.x = element_text(size = 12),
             legend.text=element_text(size = 11),legend.title=element_text(size = 12)) 
ggsave("Research/manuscript/cell_of_origin/paper/FigureS/figsupp_IFNG_parsi_survival_Malani_diagnosis.pdf",width=12, height=10,unit="cm")








############################################
# FigS1f  ifng pathway score in TCGA bulk data, select diploid compare M5 VS M0-M2,M4
#################################
mono_compare=bulk.meta[bulk.meta$FAB_Simple %in% c("M0","M1","M2","M4","M5"),c("Sample","FAB_Simple","CG_group")]
mono_compare=mono_compare[mono_compare$CG_group == "Diploid",]

mono_compare=merge(mono_compare,bulk.hallmark[,c("Sample","HALLMARK_INTERFERON_GAMMA_RESPONSE")],by="Sample")
mono_compare$FAB=ifelse(mono_compare$FAB_Simple=="M5","M5","M0/M1/M2/M4")
p=ggplot(mono_compare, aes_string(x="HALLMARK_INTERFERON_GAMMA_RESPONSE",fill="FAB")) + 
  geom_density(alpha=0.4)+ theme_classic() + ylab("Density")+ 
  scale_fill_manual(values=c("#2ca02c","steelblue2"))+ xlab("HALLMARK IFNG RESPONSE")+
  theme(axis.text.x = element_text(size = 30),legend.position = 'top',axis.title.x=element_text(size=24),
        axis.text.y = element_text(size = 30, angle = 0),axis.title.y=element_text(size=28),legend.text = element_text(size = 28),legend.title = element_text(size = 32))
png("Research/manuscript/cell_of_origin/paper/FigureS1/figs1f_TCGA_diploid_FAB_box.png",res=200,width=8, height=10,units="in")
print(p)
dev.off()


temp=merge(bulk.meta[,c("Sample","CG_group","FAB_Simple")],bulk.hallmark[,c("Sample","HALLMARK_INTERFERON_GAMMA_RESPONSE")],by="Sample")
temp[temp$CG_group=='Diploid' & temp$FAB_Simple=="M5","CG_group"]= "Diploid-mono"
temp[temp$CG_group=='Diploid' & temp$FAB_Simple!="M5","CG_group"]= "Diploid-Non_M5"
temp[temp$CG_group=='Other' & temp$FAB_Simple=="M5","CG_group"]="Other-mono"
temp[temp$CG_group=='Other' & temp$FAB_Simple!="M5","CG_group"]="Other-Non_M5"

ggplot(temp,aes(x=CG_group,y=HALLMARK_INTERFERON_GAMMA_RESPONSE))+
  geom_boxplot(aes(fill=CG_group),width=0.5,color='gray40',outlier.shape=NA,size=0.3)+
  geom_jitter(fill='gray80',color='black',shape=21,width =0.2,size=2,stroke = 0.3)+
  ylab("HALLMARK IFNG RESPONSE") + #scale_fill_manual(values=c("indianred2", "gold2"))+
  theme_classic() + #scale_y_continuous(labels = scales::percent) + 
  #stat_signif(comparisons = list(c("Diploid","Del7.7q")),textsize=2.5,y_position=c(0.08,0.28))+
  theme(legend.position="none",
        axis.title.y = element_text(size = 24),axis.title.x = element_blank(),
        axis.text.x = element_text(size = 22,angle = 45,hjust = 1),axis.text.y = element_text(size = 22))







############################################
#supplemental files
############################################
####Correlation of GO.KEGG, HALLMARK pathways with IFNG Signaling
#####bulk cohort information
temp=bulk.meta[,c("Sample","Cohort","Age","OS_days","CG_group1","CG_group2","FAB_Simple","Percent_Blasts_in_BM","HALLMARK.IFNG.RESPONSE","HALLMARK.IFNG.Parsimony","hla1","HLA2")]
write.xlsx(temp,file="Research/manuscript/cell_of_origin/paper/supplementary_file/bulk_info_scores.xlsx",row.names = FALSE)


####table S3, DEG between 9 cell types
deg.allcell=FindAllMarkers(allcell)
write.xlsx(deg.allcell,file="Research/manuscript/cell_of_origin/paper/supplementary_file/deg_allcell_class2.xlsx")
#deg.allcell[!(grepl(pattern="^(MT|AJ|AL|AC|LIN|AP)", deg.allcell$gene)),]
deg.allcell.top100 = deg.allcell %>% filter(p_val_adj<0.05) %>% arrange(desc(avg_log2FC)) %>% group_by(cluster) %>% slice(1:100)
write.csv(deg.allcell.top100,file="Research/manuscript/cell_of_origin/paper/supplementary_file/deg_allcell_class2_top100.csv",row.names = TRUE)






