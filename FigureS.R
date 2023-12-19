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
# FigureS1A:histogram showing the distribution of IFNG scores
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
# Figure S2A: cell counts in each patient
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
# Figure S2H: pie chart of cell types in single cell
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


############################################
# Figure S2 B-G:ncount violin plot by cell types 
############################################
df=allcell@meta.data[,c('nCount_RNA',"nFeature_RNA","mitoRatio","riboRatio","S.Score","G2M.Score","class2")]
df$class2=factor(df$class2,levels = c("AML","Progenitor","T","B","NK","Monocytes","Erythroids","DC","Other"))
ggplot(df,aes(x=class2,y=G2M.Score,fill=class2)) + theme_classic()+
  geom_violin(size=1) + scale_fill_manual(values=as.character(ColAssign(letters[1:10]))) + 
  theme(legend.position="none",title=element_blank(),axis.title = element_blank(),
        axis.text.x = element_text(size = 10,angle = 45,hjust = 1),axis.text.y = element_text(size = 12))+ 
  stat_summary(fun.y=mean, geom="point", shape=18,size=2, color="red")

ggplot(df,aes(x=class2,y=mitoRatio,fill=class2)) + theme_classic()+
  ggtitle("Mitochondrial Ratio") + ylim(0,1) +
  geom_violin(size=1) + scale_fill_manual(values=as.character(ColAssign(letters[1:10]))) + 
  theme(legend.position="none",axis.title = element_blank(),
        axis.text.x = element_text(size = 10,angle = 45,hjust = 1),axis.text.y = element_text(size = 12))+ 
  stat_summary(fun.y=mean, geom="point", shape=18,size=2, color="red")




############################################
# Figure S3A:AML cell count in each cytogenetic group
############################################
meta=amlnew@meta.data[,c("orig.ident","cyto4","class2")]
fraction=meta %>% group_by(cyto4)  %>% summarise(count=n()) 
fraction=as.data.frame(fraction); fraction=fraction[c(2,1,3:5),]
ggplot(fraction, aes(fill=cyto4, y=count, x=cyto4)) + theme_classic() + 
  geom_bar(position="stack",stat="identity") + xlab("") +
  scale_fill_manual(values = c("darkseagreen3","indianred2", "darkslategray", "gold2","deepskyblue"))+
  theme(axis.title.y = element_text(size = 12),legend.position = "none",
        axis.text.y =element_text(size=10), axis.text.x =element_text(angle = 45,vjust = 0.5,size=10))

############################################
# Figure S3H: patient contribution in each CG group in AML cells
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
# Figure S4A: Density plot reactome IFNG
############################################
df=amlnew@meta.data[,c("cyto4","REACTOME_INTERFERON_GAMMA_SIGNALING")]
p=ggplot(df, aes_string(x="REACTOME_INTERFERON_GAMMA_SIGNALING",fill="cyto4")) + 
  geom_density(alpha=0.9)+ theme_classic() + ylab("Density")+ labs(fill="Group")+
  scale_fill_manual(values=c("darkseagreen3","indianred2", "darkslategray", "gold2","deepskyblue"))+ 
  xlab("Reactome IFNG RESPONSE")+
  scale_fill_discrete(labels=c('Diploid-mono', 'Diploid-Nonmono')) +
  theme(axis.text.x = element_text(size = 30),legend.position = 'top',axis.title.x=element_text(size=24),
        axis.text.y = element_text(size = 30, angle = 0),axis.title.y=element_text(size=28),legend.text = element_text(size = 24),legend.title = element_text(size = 28))



############################################
# Figure S4B: petti diploid
############################################
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


############################################
# Figure S4C: Lasry nature cancer AML single cell data, validate interferon gamma
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


############################################
# Figure S5A: IFNG expression
############################################
expr=as.data.frame(t(as.matrix(combined@assays$RNA@data[c("IFNG","CD3E","CD34"),])))
expr$Cell=colnames(combined);expr$class=combined$class2
expr$class[expr$Cell %in% cells.to.keep[cells.to.keep$class=="CD4","cellID"]]='CD4'
expr$class[expr$Cell %in% cells.to.keep[cells.to.keep$class=="CD8","cellID"]]='CD8'
expr <- reshape2::melt(expr, id.vars = c("Cell","class"), measure.vars = c("IFNG","CD3E","CD34"),
                       variable.name = "Feat", value.name = "Expr")

ggplot(expr, aes(factor(class), Expr, fill = class)) + 
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0), position="left", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(Feat), scales = "free") +  #, switch = "y"
  theme_cowplot(font_size = 10) +
  theme(legend.position = "none", panel.spacing = unit(0.4, "lines"),axis.text.x = element_text(size = 14),
        panel.background = element_rect(fill = NA, color = "black"),axis.ticks = element_blank(),
        axis.text.y = element_text(size = 14),#strip.background = element_blank(),
        strip.text = element_text(face = "bold",size=14),strip.text.y.left = element_text(angle = 0)) +
  xlab("") + ylab("Expression") +
  scale_fill_manual(values = as.character(ColAssign(letters[1:20])[c(1,3:5)]))



############################################
# Figure S5 B-D: Cell proportion by groups
############################################
nkcell=readRDS("/Users/bwang8/Research/manuscript/poonam_CD8Tcell_paper/nk_cells.rds")
meta=nkcell@meta.data[,c("orig.ident","cyto","group")]

fraction=meta %>% group_by(cyto,group)  %>% summarise(count=n()) 
fraction=fraction %>% group_by(cyto) %>% mutate(Percent = round(count/sum(count),4))
fraction=as.data.frame(fraction)
fraction$cyto=factor(fraction$cyto,levels = c("Diploid-Nonmonocytic","Diploid-monocytic","Del5.5q","Del7.7q","Double deletion"))
fraction$group=factor(fraction$group,levels = c("CD16 high","CD56 high","NKG2A high"))
ggplot(fraction, aes(fill=group, y=Percent, x=cyto)) + theme_classic() + xlab("") +
  geom_bar(position="stack",stat="identity") + labs(fill='Cell type') +
  scale_fill_manual(values = as.character(ColAssign(letters[1:20])[13:18]))+
  theme(plot.title = element_text(hjust = 0.5),axis.title.y = element_text(size = 18),
        axis.text.y =element_text(size=16), axis.text.x =element_text(angle = 45,hjust = 1,size=15))
temp=as.data.frame(reshape(as.data.frame(fraction[,c(1,2,3)]), timevar = "group", idvar = "cyto", direction = "wide"))
chisq.test(temp[,2:4])


############################################
# Figure S6 A-D; 7F: HR plot
############################################
survival=bulk.meta[,c("Sample","OS_days","Dead_Alive","HALLMARK.IFNG.Parsimony","hallmark.parsi.beat","Cohort","hla1")]
survival[(!is.na(survival$OS_days) & survival$OS_days> 5*365 ),"Dead_Alive"]="Alive"
survival[(!is.na(survival$OS_days) & survival$OS_days> 5*365 ),"OS_days"]= 5*365

survival[survival$Dead_Alive=="Alive",'Dead_Alive']=0
survival[survival$Dead_Alive=="Dead",'Dead_Alive']=1
survival$Dead_Alive=as.numeric(survival$Dead_Alive)
smoothCoxph(survival$OS_days,survival$Dead_Alive,survival$HALLMARK.IFNG.Parsimony,xlab="IFNG Parsimony",main="IFNG Parsimony") 
abline(v = median(survival$HALLMARK.IFNG.Parsimony), col="grey", lwd=3, lty=2)
gene.cox <- coxph(Surv(OS_days, Dead_Alive==1) ~ HALLMARK.IFNG.Parsimony, survival)
summary(gene.cox)


############################################
# Figure S6E: HR plot
############################################
############multivariate forest plot
cancer.survival=bulk.meta
cancer.survival$Dead_Alive=ifelse(cancer.survival$Dead_Alive =="Dead", 1, ifelse(cancer.survival$Dead_Alive == "Alive",0,NA))
cancer.survival$Age <- ifelse(cancer.survival$Age > 60, "Greater", "Less")
cancer.survival$Age <- factor(cancer.survival$Age, levels=c("Less", "Greater"))
cancer.survival=cancer.survival[,c("Sample","Cohort","OS_days","Dead_Alive",'Age',"Percent_Blasts_in_BM", 
                                   "Cytogenetic_risk_incomplete","ELN2017_incomplete","HALLMARK.IFNG.Parsimony")]
cancer.survival=merge(cancer.survival,data.frame(IFITM3=t(bulk.tpm)[,"IFITM3"],Sample=rownames(t(bulk.tpm))),by="Sample")
#survdiff(Surv(OS_days,Dead_Alive == 1) ~ IFITM3>quantile(IFITM3,0.5), cancer.survival)

cancer.survival=cancer.survival[!(cancer.survival$Cytogenetic_risk_incomplete %in% c("","Unknown")),]
cancer.survival=cancer.survival[!is.na(cancer.survival$OS_days),]
#depends on the shape of bulk.tpm
cancer.survival$HALLMARK.IFNG.Parsimony=as.numeric(cancer.survival$HALLMARK.IFNG.Parsimony)
colnames(cancer.survival)[9]="IFNG_Parsimony"
covariate_names <- c(`Age:Greater` ="Age > 60", Percent_Blasts_in_BM = "Blast %", IFNG_Parsimony="IFNG_Parsimony",
                     `Cytogenetic_risk_incomplete:Intermediate` = "Cytogenetic risk Intermediate",
                     `Cytogenetic_risk_incomplete:Adverse` = "Cytogenetic risk Adverse")
cancer.survival %>%analyse_multivariate(vars(OS_days, Dead_Alive),
                                        covariates = vars(Age, Cytogenetic_risk_incomplete, Percent_Blasts_in_BM, IFNG_Parsimony),
                                        covariate_name_dict = covariate_names) -> result

forest_plot(result,factor_labeller = covariate_names, labels_displayed = c("factor"),
            endpoint_labeller = c(time="OS_days"),orderer = ~order(HR),
            ggtheme = ggplot2::theme_bw(base_size = 15),relative_widths = c(1, 1.5, 1),
            HR_x_breaks = c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4)) 







#############################################################
#Fig S7G:survival in malini data
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






