#############  R Code for Score tests for scale effects, with application to genomic analysis
#############  David Soave, dsoave@wlu.ca, Date: June 2020

source('../Functions/funcs1.R')

library(SummarizedExperiment)
library(DESeq2)
library(sva)
library(parallel)
library(xtable)

#################################################################################
### Download data from TCGA
### This may be skipped - the data file brcaExp.rda is already available within this folder

#BiocManager::install(c("TCGAbiolinks"))
library(TCGAbiolinks)

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")

GDCdownload(query)
exp.brca<-GDCprepare(query = query, save = TRUE, save.filename = "../data/brcaExp.rda")
#################################################################################
citation('TCGAbiolinks')



load(file = "../data/brcaExp.rda")

# Gender selection for female (Remove missing gender)
dim(data)
data_f=data[,!is.na(data$gender)]
data_f=data_f[,data_f$gender=="female"]
table(data_f$gender)
table(table(colData(data_f)$patient))

# Remove missing Tumour Stage (NA)
table(colData(data_f)$paper_pathologic_stage)
data_f=data_f[,!is.na(data_f$paper_pathologic_stage)]
data_f=data_f[,!(data_f$paper_pathologic_stage=="NA")]
dim(data_f)

table(data_f$gender)
table(table(colData(data_f)$patient))

# Remove duplicated patients (keep data for tissue that was formalin-fixed paraffin-embedded (FFPE))
duplicated<-names(which(table(data_f$patient)>1))
data_f2<-data_f[,-which(pmatch(data_f$patient,duplicated,dup=TRUE)>0&!data_f$is_ffpe)]
dim(data_f2)
table(table(colData(data_f2)$patient))
table(colData(data_f2)$paper_pathologic_stage)

table(colData(data_f2)$paper_BRCA_Subtype_PAM50)
table(colData(data_f2)$paper_BRCA_Subtype_PAM50,colData(data_f2)$paper_pathologic_stage)


#################
### New stage variable Stage IV==1 Stage I,II,II==0
table(data_f2$paper_pathologic_stage)
data_f2=data_f2[,grep("Stage_I",data_f2$paper_pathologic_stage)]  ## remove all but stages I-IV
data_f2$stage=0
data_f2[,grep("Stage_IV",data_f2$paper_pathologic_stage)]$stage=1
table(data_f2$paper_pathologic_stage,data_f2$stage)
#################


# Remove genes where more than 5% of samples have zero count
dim(data_f2)
countzeros<-apply(assay(data_f2),1,function(x){sum(x==0)});summary(countzeros)
cut_percent<-dim(data_f2)[2]*0.05;cut_percent
data_f2<-data_f2[(countzeros<cut_percent),]
dim(data_f2)


# Normalize expression values with respect to library size and then use log transformation
dds <- DESeqDataSet(se = data_f2,design = ~ 1,ignoreRank = T)
dds<-estimateSizeFactors(dds)
dds<-counts(dds,normalized=TRUE)
exprs.vsd<-dds
exprs.vsd<-log(exprs.vsd+1)


# Create two (2) surrogate variables
X<-cbind(data_f2$stage)
X<-as.matrix(X)
Z0<-NULL
mod = model.matrix(~X)
mod0 = model.matrix(~1)

svobj.2= sva(exprs.vsd,mod,mod0=NULL,vfilter=2000,n.sv=2)
Z<-cbind(Z0,svobj.2$sv)
head(Z);dim(Z)

#saveRDS(exprs.vsd, "../Results/exprs_vsd.rds")
#saveRDS(X, "../Results/cov_X.rds")
#saveRDS(Z, "../Results/cov_Z.rds")
exprs.vsd<-readRDS("../Results/exprs_vsd.rds")
X<-readRDS("../Results/cov_X.rds")
Z<-readRDS("../Results/cov_Z.rds")



# Function to conduct score tests and estimate parameters/coefficients of interest
testf<-function(y,x,z=NULL){
  y<-as.numeric(y)
  x<-as.numeric(as.factor(x))
  z=z
  W=x;Z=cbind(x,z)
  
  #Scale coefficient (gamma_x)
  rq1<-rq(y~x+z,tau=.5)
  d1 <- abs(resid(rq1))
  lms<-lm(d1~x)
  coef.s<-summary(lms)$coef[2,1]
  #Scale test p-values
  Scale.Score.perm<-Scale_Test(y=y,X=NULL,W=W,Z=Z)
  p.s<-Scale.Score.perm$VPL_p

  lm0<-lm(y~z)
  lm1<-lm(y~x+z)
  #Location coefficient (beta_x)
  coef.l<-summary(lm1)$coef[2,1]
  #Location test p-values
  p.l<-anova(lm0,lm1,test="Rao")[2,5]

  return(c(coef.l,p.l,coef.s,p.s))
}  

#geneidx=which(rownames(exprs.vsd)=="ENSG00000247317") #LY6E-DT
#geneidx=which(rownames(exprs.vsd)=="ENSG00000237854") #LINC00674
#testf(y=exprs.vsd[geneidx,],x=X,z=Z)
#rq1<-rq(exprs.vsd[geneidx,]~X+Z,tau=.5);d1 <- abs(resid(rq1));summary(lm(d1~X)) 
#summary(lm(exprs.vsd[geneidx,]~X+Z));anova(lm(exprs.vsd[geneidx,]~Z),lm(exprs.vsd[geneidx,]~X+Z),test="Rao")[2,5]


# Run the analysis on expression values for each gene
dim1<-dim(exprs.vsd)[1]
num.cores<-20
df1<-data.frame(do.call(rbind,mclapply(seq(1:dim1), function(x) {testf(y=(exprs.vsd[x,]),x=X,z=Z)},mc.cores=num.cores)))
names(df1)<-c("cl","pl","cs","ps")
row.names(df1)<-row.names(exprs.vsd)
saveRDS(df1, "../Results/Stage4vAll.rds")




################################################################
################################################################
### Produce Application Results for Figure 2 and Web Table 13
### User may begin here with file "../Results/Stage4vAll.rds" already available


df1<-readRDS("../Results/Stage4vAll.rds");head(df1)

############# Generate FDR q-values, select significant genes based on FDR<0.1
qthresh<-0.1
df1$pl_q<-p.adjust(df1$pl,method="fdr")
df1$ps_q<-p.adjust(df1$ps,method="fdr")

candidate_genes.L<-row.names(df1[df1$pl_q<qthresh,])
length(candidate_genes.L)
candidate_genes.S<-row.names(df1[df1$ps_q<qthresh,])
length(candidate_genes.S)



##### Figure 2(b) Plot of scale test results: -log10(p) versus gamma_x using Laplacian-model pemutation approximation test V_[P(L)] 
X11(width=8.5, height=7,type="cairo")
par(mar=c(5, 4+2, 4-3, 2) + 0.1)   #Needed for this plot

plot((df1$cs),-log10(df1$ps),pch=NA,xlab="",ylab="",cex.lab=2,cex.axis=2,ylim=c(0,6),xlim=c(-0.95,0.95))
thresh1<-max(df1[which(df1$ps_q<qthresh),]$ps);thresh1
abline(h=-log10(thresh1))
points((df1$cs),-log10(df1$ps),col=gray.colors(n=1,.7))
points((df1[rownames(df1)%in%candidate_genes.S,]$cs),-log10(df1[rownames(df1)%in%candidate_genes.S,]$ps),col=1,pch=19)
mtext(expression(paste(-log[10](p[S]))), side=2, line=3,cex=2)
mtext(expression(paste("Coefficient ",hat(gamma)[X])), side=1, line=4,cex=2)
dev.copy2pdf(device = x11,file='../Results/F2b_VPL.pdf', out.type = "pdf")


##### Figure 2(a) Plot of location test results: -log10(p) versus beta_x
X11(width=8.5, height=7,type="cairo")
par(mar=c(5, 4+2, 4-3, 2) + 0.1)   #Needed for this plot

plot(df1$cl,-log10(df1$pl),pch=NA,xlab="",ylab="",cex.lab=2,cex.axis=2,ylim=c(0,6))
points(df1$cl,-log10(df1$pl),col=gray.colors(n=1,.7))
mtext(expression(paste(-log[10](p[L]))), side=2, line=3,cex=2)
mtext(expression(paste("Coefficient ",hat(beta)[X])), side=1, line=4,cex=2)
dev.copy2pdf(device = x11,file='../Results/F2a.pdf', out.type = "pdf")


### Results Web Table 13
### Note: These results must be transfered manually to the Supplementary latex document
df1<-df1[order(df1$ps),]
head(df1)
mat1<-df1[row.names(df1)%in%candidate_genes.S,c(2,1,4,3)]
mat1<-mat1[order(mat1$ps),]
mat1
library(xtable)
print(xtable(mat1,digits=c(0,-1,2,-1,2),include.rownames=FALSE))



#####################################################################
#####################################################################
#####################################################################
######### Model-Based Results (Application Web Figure 4 and Web Table 14)

exprs.vsd<-readRDS("../Results/exprs_vsd.rds")
X<-readRDS("../Results/cov_X.rds")
Z<-readRDS("../Results/cov_Z.rds")

X11(width=8, height=7,type="cairo")
par(mar=c(5, 4+1, 4-1, 2) + 0.1) 


#library(gamlss)
geneidx1=which(rownames(exprs.vsd)=="ENSG00000247317") #LY6E-DT
yy=exprs.vsd[geneidx1,]

fitNO<-gamlss(yy~X+Z,sigma.formula = ~X,family=NO)
summary(fitNO)
Figurei="S43a"
qqnorm(resid(fitNO),main=expression(paste(italic("LY6E-DT")," / Gaussian Fit")));abline(0,1,col=2)
dev.copy2pdf(device = x11,file=paste("../Results/F",Figurei,".pdf",sep=""), out.type = "pdf")

fitEGB2<-gamlss(yy~X+Z,sigma.formula = ~X,family=EGB2)
summary(fitEGB2)
qresid1<-pEGB2(yy, mu = fitted(fitEGB2,"mu"), sigma = fitted(fitEGB2,"sigma"), nu = fitted(fitEGB2,"nu"), tau = fitted(fitEGB2,"tau") )
Figurei="S43b"
qqnorm(qnorm(qresid1),main=expression(paste(italic("LY6E-DT")," / EGB2 Fit")));abline(0,1,col=2)
dev.copy2pdf(device = x11,file=paste("../Results/F",Figurei,".pdf",sep=""), out.type = "pdf")


geneidx2=which(rownames(exprs.vsd)=="ENSG00000237854") #LINC00674
yy=exprs.vsd[geneidx2,]

fitNO<-gamlss(yy~X+Z,sigma.formula = ~X,family=NO)
summary(fitNO)
Figurei="S43c"
qqnorm(resid(fitNO),main=expression(paste(italic("LINC00674")," / Gaussian Fit")));abline(0,1,col=2)
dev.copy2pdf(device = x11,file=paste("../Results/F",Figurei,".pdf",sep=""), out.type = "pdf")

fitEGB2<-gamlss(yy~X+Z,sigma.formula = ~X,family=EGB2)
summary(fitEGB2)
qresid1<-pEGB2(yy, mu = fitted(fitEGB2,"mu"), sigma = fitted(fitEGB2,"sigma"), nu = fitted(fitEGB2,"nu"), tau = fitted(fitEGB2,"tau") )
Figurei="S43d"
qqnorm(qnorm(qresid1),main=expression(paste(italic("LINC00674")," / EGB2 Fit")));abline(0,1,col=2)
dev.copy2pdf(device = x11,file=paste("../Results/F",Figurei,".pdf",sep=""), out.type = "pdf")


######### Figures (Distributions) Web Figure5
X11(width=8, height=7,type="cairo")
par(mar=c(5, 4+1, 4-1, 2) + 0.1) 

Figurei="S5a"
geneidx1=which(rownames(exprs.vsd)=="ENSG00000247317") #LY6E-DT
yy=exprs.vsd[geneidx1,]
boxplot(yy~X,axes=F,xlab='Group',ylab='Expression',cex.lab=1.5,lwd=1.5)
axis(1, 1:2, c("Stages I-III","Stage IV"),las=1,cex.axis=1.5)
axis(2,cex.axis=1.5)
box()
title(expression(paste(italic("LY6E-DT"))),cex.main=2)
dev.copy2pdf(device = x11,file=paste("../Results/F",Figurei,".pdf",sep=""), out.type = "pdf")

Figurei="S5b"
geneidx2=which(rownames(exprs.vsd)=="ENSG00000237854") #LINC00674
yy=exprs.vsd[geneidx2,]
boxplot(yy~X,axes=F,xlab='Group',ylab='Expression',cex.lab=1.5,lwd=1.5)
axis(1, 1:2, c("Stages I-III","Stage IV"),las=1,cex.axis=1.5)
axis(2,cex.axis=1.5)
box()
title(expression(paste(italic("LINC00674"))),cex.main=2)
dev.copy2pdf(device = x11,file=paste("../Results/F",Figurei,".pdf",sep=""), out.type = "pdf")


