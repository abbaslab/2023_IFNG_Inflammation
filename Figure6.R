library(pROC)
library(caret)
library(e1071)
library(glmnet)
library(lars)
set.seed(7)


#############################################################
###fig 6A, B
#############################################################
drug.auc=read.csv("Research/AML/info_files/public_AMLdata/Beat_AML/BEAT_drug_AUC.csv")
drug.auc=drug.auc[,c("Patient","Venetoclax")]
drug.auc=merge(drug.auc,bulk.meta[,c("HALLMARK.IFNG.Parsimony","Sample")],by.x="Patient",by.y="Sample") 
drug.auc=drug.auc[!is.na(drug.auc$Venetoclax),]
ggscatter(as.data.frame(drug.auc), x = "HALLMARK.IFNG.Parsimony", y = "Venetoclax", add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",xlab = "Parsimonious IFNG score", ylab = "Venetoclax resistance",size=2,
          cor.coef.size=7) + #,cor.coef.coord=c(4000,300)
  theme(legend.position = "none", axis.title = element_text(size = 26),
        axis.text.x = element_text(size = 23),axis.text.y = element_text(size = 23),
        axis.ticks = element_blank()) 



#############################################################
###fig 6F,G Parsimonious score
############################################################
dat=merge(dat,bulk.meta[,c("Sample","OS_days","Cohort")])
dat=dat[dat$Cohort %in% c("BEAT1","BEAT2"),];dat=dat[,-203]
rownames(dat)=dat$Sample; dat=dat[,-1]
dat=dat[!is.na(dat$OS_days),]
dat$OS_days=ifelse(dat$OS_days>median(dat$OS_days),0,1)
y=as.matrix(dat[,201])
x=as.matrix(dat[,-201])
cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial")
#plot(cv.lasso)
fit.beat=glmnet(x,y,family="binomial",alpha=1,lambda=cv.lasso$lambda.min)
coef(fit.beat)
row.names(coef(fit.beat))[which(coef(fit.beat)!=0)][-1]


#############################################################
###apply lasso logistic regression(glmnet) on OS_days (below/above) to reduce genes
dat=merge(dat,bulk.meta[,c("Sample","OS_days")])
rownames(dat)=dat$Sample; dat=dat[,-1]
dat=dat[!is.na(dat$OS_days),]
dat$OS_days=ifelse(dat$OS_days>median(dat$OS_days),0,1)
y=as.matrix(dat[,201])
x=as.matrix(dat[,-201])
cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial")
plot(cv.lasso)
fit=glmnet(x,y,family="binomial",alpha=1,lambda=cv.lasso$lambda.min)
coef(fit)
row.names(coef(fit))[which(coef(fit)!=0)][-1]

#####gsva score 
gsva_score=gsva(as.matrix(bulk.tpm), list(hallmark.parsi.beat=row.names(coef(fit.beat))[which(coef(fit.beat)!=0)][-1]),method="ssgsea",abs.ranking=F)
gsva_score=as.data.frame(t(gsva_score))
gsva_score$Sample=rownames(gsva_score)
bulk.meta=merge(bulk.meta,gsva_score,by="Sample")



#####survival curve by new IFNG score
survival=bulk.meta[,c("Sample","OS_days","Dead_Alive","HALLMARK.IFNG.Parsimony","hallmark.parsi.beat","Cohort","hla1")]
survival[(!is.na(survival$OS_days) & survival$OS_days> 5*365 ),"Dead_Alive"]="Alive"
survival[(!is.na(survival$OS_days) & survival$OS_days> 5*365 ),"OS_days"]= 5*365
#survival=survival[survival$Cohort=="BEAT2",]
fit=survfit(Surv(OS_days,Dead_Alive == "Dead") ~ HALLMARK.IFNG.Parsimony>quantile(HALLMARK.IFNG.Parsimony,0.5), survival)
p=ggsurvplot(fit,data=survival,pval=T,palette = c("blue","red"),xlab="Survival Time (Years)",
             xscale=365,xlim=c(0,5*365),break.x.by = 1*365,surv.median.line = "hv",
             legend.title="Parsimonious IFNG", risk.table = TRUE,#title="BEAT2, >0.5",
             legend.labs = c("Below Median", "Above Median"),pval.coord = c(1000,0.8),pval.size =10) 

p$plot+theme(axis.text.x=element_text(size=28),axis.title.y = element_text(size = 27),
             axis.text.y=element_text(size = 28),axis.title.x = element_text(size = 27),
             legend.text=element_text(size = 25),legend.title=element_text(size = 26)) 



####HR plot for parsimonous score
survival[survival$Dead_Alive=="Alive",'Dead_Alive']=0
survival[survival$Dead_Alive=="Dead",'Dead_Alive']=1
survival$Dead_Alive=as.numeric(survival$Dead_Alive)
smoothCoxph(survival$OS_days,survival$Dead_Alive,survival$HALLMARK.IFNG.Parsimony,xlab="IFNG Parsimony",main="IFNG Parsimony") 
abline(v = median(survival$HALLMARK.IFNG.Parsimony), col="grey", lwd=3, lty=2)
gene.cox <- coxph(Surv(OS_days, Dead_Alive==1) ~ HALLMARK.IFNG.Parsimony, survival)
summary(gene.cox)

##########################################
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
###apply lasso cox regression(glmnet) on OS_days (continuous) to reduce proteins 
dat=merge(dat,bulk.meta[,c("Sample","OS_days","Dead_Alive")])
rownames(dat)=dat$Sample; dat=dat[,-1]
dat=dat[!is.na(dat$OS_days),]
dat$Dead_Alive=ifelse(dat$Dead_Alive=='Alive',0,1)
y=as.matrix(dat[,c(201,202)])
x=as.matrix(dat[,-c(201,202)])
fit2=coxph(as.formula(paste("Surv(OS_days, Dead_Alive==1)", paste(colnames(x), collapse=" + "), sep=" ~ ")), dat)
cv.lasso2 <- cv.glmnet(x, y, alpha = 1, family = "cox")
fit2=glmnet(x,y,family="cox",alpha=1,lambda=cv.lasso2$lambda.min)
coef(fit2)
row.names(coef(fit2))[which(coef(fit2)!=0)][-1]


















