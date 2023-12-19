library(ggpubr)
library(ggplot2)
library(ggpmisc)
library(Seurat)


############################################
# Fig2A Oncoprint
#################################
#Data preparation
metadata <- read_excel("Research/AML/info_files/patients_info/Copy of Metadata_Hemepath 4.14.2022.xlsx",sheet = "Clean_metadata")
clinical=as.data.frame(read_excel("Research/AML/info_files/patients_info/Copy of Metadata_Hemepath 4.14.2022.xlsx",sheet = "Clinical"))
metadata=as.data.frame(metadata[metadata$Bofei_AML_cellAnalysis=="Include in Analysis",])
metadata=metadata[metadata$Disease=="AML",]
plot.data <- metadata[,c(1,25,29:32,36:41,43,44,46:48,50)]#!(metadata$SampleID %like% "^N")
plot.data$AML_status=clinical$`AML Status at BL`
plot.data$AML_status=ifelse(plot.data$AML_status=="Denovo",'Denovo','Secondary')
rownames(plot.data)=plot.data$SampleID
plot.data=plot.data[,2:ncol(plot.data)];plot.data=as.data.frame(t(plot.data))
rownames(plot.data)=c("Karyotype","Del5/5q","Del7/7q","Trisomy8","Del17/17p","EZH2","TP53","FLT3","NPM1","RAS","SF3B1","IDH1","IDH2","DNMT3A","TET2","ASXL1","ELN2017","Status")


mut=as.data.frame(str_split_fixed(metadata$MDL_Pos_Mutations, ',',Inf))
mut=as.data.frame(add_column(mut, Patient = metadata$SampleID, .before = 1))
mut[18,3]="TP53"
mut_new=mut[,1:2]
for (i in colnames(mut)[3:ncol(mut)]) {
  temp=mut[,c("Patient",i)]
  names(temp)=names(mut_new)
  mut_new=rbind(mut_new,temp)
}
mut_new=reshape2::dcast(mut_new,Patient~V1)
rownames(mut_new)=mut_new$Patient
mut_new=mut_new[c(20,1:19),-c(1,2,23)]##subset(mut_new, select=-c(None,Patient))
plot.data=rbind(plot.data[1:5,],as.data.frame(t(mut_new)),plot.data[c("ELN2017","Status"),])


lgd1 = Legend(labels = c("Mutation", "No Mutation"), legend_gp = gpar(fill = c("blue", "lightgray")), title = "Mutation")
lgd2 = Legend(labels=c("Diploid", "Complex", "Not Complex"), legend_gp = gpar(fill=c("rosybrown","palevioletred", "lightpink")), title="Karyotype")
lgd3 = Legend(labels=c("Alternation", "No alternation"), legend_gp=gpar(fill=c("steelblue3","lightgray")), title= "Chromosome Alternations")
lgd4= Legend(labels=c("Favorable", "Intermediate", "Adverse"), legend_gp = gpar(fill=c("darkseagreen3","khaki", "paleturquoise")), title="ELN 2017")
lgd5= Legend(labels=c("Secondary", "Denovo"), legend_gp = gpar(fill=c("wheat3", "seashell")), title="AML Status")
pd = packLegend(lgd2, lgd3, lgd1,lgd4,lgd5, direction = "horizontal", gap = unit("1", "cm"))


Heatmap(as.matrix(plot.data), col = c("Diploid"="rosybrown","Complex"="palevioletred", "Not Complex"="lightpink",
                                      "del7"="steelblue3","del7q"="steelblue3","del5q"="steelblue3", "Trisomy8"="steelblue3", "del17_17p"="steelblue3","No"="lightgray",
                                      "1" = "blue", "0"="lightgray", 
                                      "Favorable"="darkseagreen3", "Intermediate"="khaki", "Adverse"="paleturquoise", 
                                      "Secondary"= "wheat3", "Denovo"="seashell1"), 
        row_title=" ", row_split=c(rep("A", 5), rep("B",26), rep("C",2)),border=FALSE, width=unit(5, "in"),
        height=unit(6,"in"), rect_gp = gpar(col = "black", lwd = 1), show_heatmap_legend = FALSE)
draw(pd, just=c("center","bottom"), y=unit(4, "cm"))


############################################
# Fig2C dotplot
#################################
features=c("CD34","CD33","KIT","CFD","AZU1","CD3E","CD8A","CD19","CD79A","NCAM1","NKG7","CD68","CD14","ITGAX","HBA1","HBB","GATA1","CD1C","PLD4")
feat_AML=c("CD34","KIT","CD68")
feat_prog=c("CD33","CFD","AZU1")
feat_T=c("CD3E","CD8A");feat_B=c("CD19","CD79A");feat_NK=c("NCAM1","NKG7")
feat_mono=c("CD14","ITGAX","S100A12")
feat_ery=c("HBA1","HBB");feat_dc=c("CD1C","PLD4")
dot_plot_x <- function(seurat_object, features, lim=NULL, group="class9", feats=NULL, min=-1.5, max=1.5){
  dot_plot <- DotPlot(object = seurat_object, features=features, group.by=group, col="RdYlBu", col.min = min, col.max = max, scale.max = 50) + 
    theme(legend.position="none",axis.ticks.x = element_blank(), axis.text.x = element_text(size=10, angle=60, hjust=1), 
          axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
    labs(x=feats, y=NULL) +  scale_y_discrete(limits=lim) + coord_flip()
  return(dot_plot)
}

dot_plot_x(allcell, rev(features), group = "class2")  ##allcell is the seurat object





############################################
# Fig2D barplot
#################################
library(scales)
library(RColorBrewer)
ColAssign <- function(Var,palettes="Classic 20"){
  require(ggthemes);require(RColorBrewer)
  pal <- tableau_color_pal(palette = palettes,direction = 1,type="regular")
  if (length(Var) > 20) {
    palOut <- colorRampPalette(pal(20))(length(Var))
    names(palOut) <- Var
  } else if (length(Var) == 20) {
    palOut <- pal(20)
    names(palOut) <- Var
  } else if (length(Var) < 20) {
    palOut <- pal(20)
    palOut <- setdiff(palOut,c("#7F7F7F","#C7C7C7"))# remove grey colors
    #palOut <- sample(palOut)
    palOut <- c(palOut,c("#7F7F7F","#C7C7C7"))
    palOut <- palOut[1:length(Var)]
    names(palOut) <- Var
  }
  return(palOut)
}

show_col(ColAssign(letters[1:10]))
meta=allcell@meta.data[,c("orig.ident","cyto4","class2")]
fraction=meta %>% group_by(cyto4,class2)  %>% summarise(count=n()) 
fraction=fraction %>% group_by(cyto4) %>% mutate(Percent = round(count/sum(count),4))
fraction=as.data.frame(fraction)
fraction$orig.ident=factor(fraction$orig.ident,levels = c('PT20A',"PT23A","PT25A","PT12A","PT19A","PT30A","PT32A","PT9A", "PT13A","PT17A","PT22A","PT26A","PT14A","PT15A","PT21A","PT27A","PT28A","PT10A","PT16A","PT29A"))
fraction$cyto4=factor(fraction$cyto4,levels = c("Diploid-Nonmono","Diploid-mono","Del5.5q","Del7.7q","Double del"))
fraction$class2=factor(fraction$class2,levels = c("AML","Progenitor","T","B","NK","Monocytes","Erythroids","DC","Other"))
ggplot(fraction, aes(fill=class2, y=Percent, x=cyto4)) + theme_classic() + xlab("") +
  geom_bar(position="stack",stat="identity") + labs(fill='Cell type') +
  scale_fill_manual(values = as.character(ColAssign(letters[1:10])))+
  theme(plot.title = element_text(hjust = 0.5),axis.title.y = element_text(size = 18),
        axis.text.y =element_text(size=16), axis.text.x =element_text(angle = 45,hjust = 1,size=15))


############################################
# Fig2G Correlation
#################################
metadata=as.data.frame(metadata[metadata$Bofei_AML_cellAnalysis=="Include in Analysis",])
metadata=metadata[metadata$Timepoint=="A",]
fraction=merge(fraction,metadata[,c("SampleID","% Blast by IHC","% Blast by FC","% Blast by aspirate smear")],by.x="orig.ident",by.y='SampleID')
fraction$Percent=fraction$Percent * 100
ggscatter(fraction, x="Percent", y="% Blast by FC", add="reg.line",conf.int = FALSE,shape=1,
          cor.coef = TRUE,cor.coef.size = 4,cor.method = "pearson",
          ylab = "%  by Flow Cytometry", xlab = "% by scRNA-seq",size=2) + 
  theme(legend.position = "right", axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))

