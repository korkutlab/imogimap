#This code reproduces assessment of combinatorial interactions for T-cell dysfunction signature genes in Uterine Corpus Endometrial Carcinoma (UCEC).

# BAsed on ImogiMap package. Developed by Behnaz Bozorgui for the Korkut lab at MDA.

##Install from github-----------------------------

library(remotes)
install_github('korkutlab/imogimap')
library(imogimap)

#Input--------------------------------------------

#Set directory
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)

#UCEC T cell dysfunction genes
df_TCDF_genes<- read.csv("TCell_dys_genes_ucec.csv",header = T)



#Calculate and evaluate synergy scores------------


#Specificity and robustness evaluation of scores may take some time.
#To skip these evaluations set specificity and sensitivity(Robustness) to FALSE.
my_scores <- im_syng_tcga( onco_gene=df_TCDF_genes$Gene,
                           icp_gene = icp_gene_list,
                           cohort="ucec",
                           select_iap = "IFNGscore",
                           specificity = F,
                           sensitivity = F,
                           N_iteration_specificity = 1000,
                           N_iteration_sensitivity = 1000)

my_scores <- my_scores[-which(is.na(my_scores$Synergy_score)),]


#Compare-----------------------------------------

#To compare your results with ours
df_scores <- read.csv("TCell_dys_scores.csv")
df_scores <- df_scores[-which(is.na(df_scores$Synergy_score)),]


library(ggplot2)
library(circlize)
#Plot2A---------------------------------------------
df <- df_scores
df$r <- -log(df$sensitivity_R)

while (!is.null(dev.list()))  dev.off()
ggplot(df,aes(Synergy_score,r))+
  geom_abline(slope =0,intercept = 1,color="red",size=2)+
  geom_point(cex=1.5)+  theme_bw()+
  xlim(c(-0.8,0.8))+
  #ylim(c(-3,3))+
  labs(title= "", x= "Synergy scores", y= "-log(R) Robustness",tag="1")+
  guides(x.sec = "axis", y.sec = "axis") +
  theme( axis.text.x = element_text(size=30,angle = 90),
         axis.text.y = element_text(size=30,angle = 90),
         axis.title = element_text(size=30),
         plot.margin=unit(c(1,1,1.5,1.2),"cm"),
         plot.tag.position = c(.05, 1.15),
         plot.tag = element_text(size=60,colour="red",hjust = -0.3,vjust = 9),
         axis.line.x.top  = element_line( size = 1),
         axis.line.x.bottom = element_line( size = 1),
         axis.line.y.left  = element_line( size = 1),
         axis.line.y.right  = element_line( size = 1),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())

#Plot2B---------------------------------------------
df <- df_scores
df$r <- -log(df$sensitivity_R)
df <- df[df$r>1,]
n <- nrow(df)
df$Q <- p.adjust(df$wilcox_pvalue,method="BH",n =n)

while (!is.null(dev.list()))  dev.off()
par(mar = c(5,7,10,20))
ggplot(df,aes(Synergy_score,-log((Q))))+
  geom_abline(slope =0,intercept = -log(0.1),color="red",size=2)+
  geom_point(size=1.5)+  theme_bw()+xlim(c(-0.8,0.8))+
  labs(title= "", x= "Synergy scores", y= " -log(Q) Significance\n",tag = "-log(0.1)")+
  guides(x.sec = "axis", y.sec = "axis") +
  theme( axis.text.x = element_text(size=30,angle = 90),
         axis.text.y = element_text(size=30,angle = 90),
         axis.title = element_text(size=30),
         axis.line.x.top  = element_line( size = 1),
         axis.line.x.bottom = element_line( size = 1),
         axis.line.y.left  = element_line( size = 1),
         axis.line.y.right  = element_line( size = 1),
         plot.margin=unit(c(1,1,1.5,1.2),"cm"),
         plot.tag.position = c(.14, .4),
         plot.tag = element_text(size=30,colour="red",hjust = -0.1,vjust = 1.5),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())

#Plot2C---------------------------------------------
df <- df_scores
df$r <- -log(df$sensitivity_R)
df <- df[df$r>1,]
n <- nrow(df)
df$Q <- p.adjust(df$wilcox_pvalue,method="BH",n =n)
df <- df[df$Q<0.1,]
n<- nrow(df)
df$Q <- p.adjust(df$specificity_pvalue,method="BH",n =n)

while (!is.null(dev.list()))  dev.off()
par(mar = c(5,7,10,20))
ggplot(df,aes(Synergy_score,-log((Q))))+
  geom_abline(slope =0,intercept = -log(0.1),color="red",size=2)+
  geom_point(size=2)+  theme_bw()+xlim(c(-0.8,0.8))+
  guides(x.sec = "axis", y.sec = "axis") +
  labs(title= "", x= "Synergy scores", y= "-log(Q) Specificity",tag = "-log(0.1)")+
  theme( axis.text.x = element_text(size=30,angle = 90),
         axis.text.y = element_text(size=30,angle = 90),
         axis.title = element_text(size=30),
         axis.line.x.top  = element_line( size = 1),
         axis.line.x.bottom = element_line( size = 1),
         axis.line.y.left  = element_line( size = 1),
         axis.line.y.right  = element_line( size = 1),
         plot.margin=unit(c(1,1,1.5,1.2),"cm"),
         plot.tag.position = c(.1, .65),
         plot.tag = element_text(size=30,colour="red",hjust = -0.1,vjust = 1.),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())


#plot2D---------------------------------------------

df <- df_scores
df$r <- -log(df$sensitivity_R)
df <- df[df$r>1,]
n <- nrow(df)
df$Q <- p.adjust(df$wilcox_pvalue,method="BH",n =n)
df <- df[df$Q<0.1,]
n<- nrow(df)
df$Q <- p.adjust(df$specificity_pvalue,method="BH",n =n)
df <- df[df$Q<0.1,]

#change seed to explore different layouts
while (!is.null(dev.list()))  dev.off()
im_netplot(df =df, Immune_phenotype  ="IFNGscore",cutoff = 0.,seed=3)


#Plot2E---------------------------------------------
while (!is.null(dev.list()))  dev.off()
im_boxplot_tcga(onco_gene = "SERPINB9",
                icp_gene = "CTLA4",
                cohort = "ucec",
                Immune_Feature ="IFNGscore",
                logtrans = F)
dev.off()

#Plot2F---------------------------------------------
while (!is.null(dev.list()))  dev.off()
im_boxplot_tcga(onco_gene = "CD70",
                icp_gene = "CD86",
                cohort = "ucec",
                Immune_Feature ="IFNGscore",
                logtrans = F)


#Plot2G---------------------------------------------
library(survival)
library(survminer)

#Read TCGA UCEC survival values
#Survival values are taken from Liu J., et al. An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics. Cell 173,400 (2018).
clin <- read.csv("TCGA_PFIsurvival_UCEC.csv",header = T)

#Read TCGA UCEC expression
df <-curatedTCGAData::curatedTCGAData( diseaseCode = "ucec",
                                       assays = c("RNASeq2GeneNorm"),
                                       dry.run = F)@ExperimentList@listData[[1]]
data_expression <- df@assays$data@listData[[1]]
colnames(data_expression)<-  substr(colnames(data_expression), 1, 15)

#Select genes
df_selected <- t(data_expression[rownames(data_expression ) %in% c("CD86","CD70"),])

#Log2p1 transform
df_selected <- scale(log2(df_selected+1),center = T,scale = T)

#Take mean of multiple samples
df_selected <- as.data.frame(df_selected)
df_selected$bcr_patient_barcode <- substr(rownames(df_selected), 1, 12)
df_selected <- aggregate(list(df_selected$CD70,df_selected$CD86),
                         list(df_selected$bcr_patient_barcode), FUN=mean)
rownames(df_selected)<- df_selected$Group.1
df_selected$Group.1<-NULL
colnames(df_selected)<- c("CD70","CD86")
df_selected <- as.matrix(df_selected)

#Get quantile ranks for each sample
df_select_qr <- get_quantile_rank(df_selected)
df_select_qr$bcr_patient_barcode <- df_select_qr$Tumor_Sample_ID

#Select samples with  expression values in Q1 or Q3 quantile.
df_select_qr <- cbind(df_select_qr,as.integer((df_select_qr[,2]*2+df_select_qr[,3])/3))
df_select_qr <- df_select_qr[df_select_qr[,5]%in% c(1,4),]
colnames(df_select_qr)[5]<-"status"

#Merge data
df_surv <- merge(df_select_qr ,clin,by="bcr_patient_barcode")
df_surv$pcevent <- df_surv$PFI
df_surv$pctime <- df_surv$PFI.time
df_surv <- df_surv[!is.na(df_surv$pctime),]
df_surv <- df_surv[,c("Tumor_Sample_ID","pcevent","pctime","status")]

#Survival analysis
cox1 <-  coxph( Surv(pctime,pcevent) ~ status , data = df_surv )
hr <- exp(cox1$coefficients)
fit1 <- survfit( Surv(pctime,pcevent) ~ status , data = df_surv ,type="kaplan-meier")

pval <- surv_pvalue(
  fit1,
  data = df_surv,
  method = "Log-rank"
)$pval

#plot
while (!is.null(dev.list()))  dev.off()
ggsurv <- ggsurvplot(fit1, data = df_surv, risk.table=F,  size = 2,
                     palette = c("blue","red"),
                     xlab="Time (Days)",  legend.title="",
                     legend.labs = c("Low CD86 Low CD70","High CD86 High CD70"),
                     legend=c(0.7,0.4),font.legend=30,
                     #pval=T,pval.method=T,
                     font.x=30,font.y=30,font.tickslab=30)
ggsurv$plot <- ggsurv$plot +
  ggplot2::annotate(
    "text",
    x = 0, y = 0,
    vjust = 0.1, hjust = 0,
    label = paste0("Log-Rank Pvalue", pval  ,"\nHazard Ratio  ,",hr),
    size = 10
  )

ggsurv


#Supplementary figure----------------------------------------

df_all<-list()

#Read TCGA mRNA data for all disease cohorts
for(cohortID in 1:length(TCGA_disease_list)){
  disease <- TCGA_disease_list[cohortID]

  message("\nReading TCGA ",toupper(disease)," data\n")
  df <-curatedTCGAData::curatedTCGAData( diseaseCode =disease,
                                         assays = c("RNASeq2GeneNorm"), dry.run = F)@ExperimentList@listData[[1]]

  df <- df@assays$data@listData[[1]]
  colnames(df)<-  substr(colnames(df), 1, 15)
  df_all[[cohortID]] <- df
}

df <- lapply(df_all, as.data.frame)
df <- lapply(df, function(x){x$gene <- rownames(x); x})
df <-Reduce(merge,df)
rownames(df)<-df$gene
df$gene<-NULL


#Calculate IAP values
df1 <- get_emt_score(df)
df2 <- get_angio_score(df)
df3 <- get_ifng_score(df)
df4 <-TCGA_Leukocyte_fraction
df5 <-TCGA_TMB[,c("Tumor_Sample_ID","TMB_Non.silent_per_Mb")]
df6 <- TCGA_TMB[,c("Tumor_Sample_ID","TMB_Silent_per_Mb")]

data_feature <- get_features(df)

d1<- merge(df1, data_feature[,c("Tumor_Sample_ID" ,"EMTscore")],by= "Tumor_Sample_ID")
d2<- merge(df2, data_feature[,c("Tumor_Sample_ID" ,"AGscore")],by= "Tumor_Sample_ID")
d3<- merge(df3, data_feature[,c("Tumor_Sample_ID" ,"IFNGscore")],by= "Tumor_Sample_ID")
d4<- merge(df4, data_feature[,c("Tumor_Sample_ID" ,"Leukocyte_fraction")],by= "Tumor_Sample_ID")
d5<- merge(df5, data_feature[,c("Tumor_Sample_ID" ,"TMB_Non.silent_per_Mb")],by= "Tumor_Sample_ID")
d6<- merge(df6, data_feature[,c("Tumor_Sample_ID" ,"TMB_Silent_per_Mb")],by= "Tumor_Sample_ID")

#plot
while (!is.null(dev.list()))  dev.off()
plot(d1$EMTscore.x,d1$EMTscore.y,
     xlab="EMT score",ylab="Transformed\n EMT score")
plot(d2$AGscore.x,d2$AGscore.y,
     xlab="Vascularization",ylab="Transformed\n Vascularization")
plot(d3$IFNGscore.x,d3$IFNGscore.y,
     xlab="IFNG score",ylab="Transformed\n IFNG score")
plot(d4$Leukocyte_fraction.x,d4$Leukocyte_fraction.y,
     xlab="Leukocyte fraction",ylab="Transformed\n Leukocyte fraction")
plot(d5$TMB_Non.silent_per_Mb.x,d5$TMB_Non.silent_per_Mb.y,
     xlab="Non.silent TMB",ylab="Transformed\n Non.silent TMB")
plot(d6$TMB_Silent_per_Mb.x,d6$TMB_Silent_per_Mb.y,
     xlab="Silent TMB",ylab="Transformed\n Silent TMB")

hist(df1$EMTscore,breaks=20,xlab="EMT score",main = "")
hist(df2$AGscore,breaks=20,xlab="Vascularization",main = "")
hist(df3$IFNGscore,breaks=20,xlab="IFNG score",main = "")
hist(df4$Leukocyte_fraction,breaks=20,xlab="Leukocyte fraction",main = "")
hist(df5$TMB_Non.silent_per_Mb,breaks=20,xlab="Non.silent TMB",main = "")
hist(df6$TMB_Silent_per_Mb,breaks=20,xlab="Silent TMB",main = "")

hist(data_feature$EMTscore,breaks=20,xlab="Transformed EMT score",main = "")
hist(data_feature$AGscore,breaks=20,xlab="Transformed Vascularization",main = "")
hist(data_feature$IFNGscore,breaks=20,xlab="Transformed IFNG score",main = "")
hist(data_feature$Leukocyte_fraction,breaks=20,xlab="Transformed Leukocyte fraction",main = "")
hist(data_feature$TMB_Non.silent_per_Mb,breaks=20,xlab="Transformed Non.silent TMB",main = "")
hist(data_feature$TMB_Silent_per_Mb,breaks=20,xlab="Transformed Silent TMB",main = "")

