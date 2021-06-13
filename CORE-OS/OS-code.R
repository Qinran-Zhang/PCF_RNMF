########## sample #######
library(survival);library(survminer);library(rms);library(timeROC);library(polspline)
library(data.table);library(MASS)
setwd("C:/Users/MSI-NB/Desktop/code/CORE-OS")
load("CORE-OS.Rdata")

source("RNMF.R")
RNMF_data=list(train_data1,train_data2,train_data3)
K_select=7;lambda_select=10^(-2)
set.seed(2191)
t1=Sys.time()
result=RNMF(RNMF_data,K=K_select,lambda=lambda_select,maxiter=10000,speak=T)
result=RNMF(RNMF_data,K=K_select,lambda=lambda_select,maxiter=10000,speak=T)
Sys.time()-t1
result$iter;result$loss
clin<-a[,c(1,3,4,rep(1,ncol(result$W)))];dim(clin)
rownames(clin)=clin[,1]
clin=clin[co_sample,];dim(clin)
clin=clin[c(sample_train),]
for(i in 1:ncol(result$W)){clin[,3+i]=result$W[,i]}
clin<-na.omit(clin);dim(clin)
colnames(clin)<-c("sampleID","status","time",paste("module",1:ncol(result$W),sep=""))
clin=clin[,-1]
ph_df <- del_nonPH(gene_list = colnames(clin)[c(3:ncol(clin))], survival_info_df = clin)
res.cox <- coxph(Surv(time, status) ~., data = clin[,c("status","time",as.character(ph_df[,1]))])
n=nrow(clin)
res.cox=stepAIC(res.cox,k = log(n))
summary(res.cox)
test.ph <- cox.zph(res.cox)
test.ph
#ROC
pre.result=predict(res.cox,newdata=clin)
ROC.time1=timeROC(clin$time, clin$status, pre.result, other_markers = NULL, cause=1,
	    weighting = "marginal",
		times=c(365,365*2,365*3,365*4,365*5),
		 ROC = TRUE, iid = F) #
plot(ROC.time1,time=365*1,title=F)     
title(main=paste("AUC at 1-year =",as.numeric(round(ROC.time1$AUC[1],3)),seq=""),cex=2,line=1)
#KM-Plot
pre.result= predict(res.cox, newdata=clin)
med=median(pre.result)
class=rep("low-risk",length(pre.result))
class[which(pre.result>med)]="high-risk"
clin_km<- clin[,c(1:ncol(clin),1)];dim(clin_km)
clin_km[,ncol(clin_km)]=class
colnames(clin_km)[ncol(clin_km)]<-"predict"
clin_km$time=clin_km$time/30
library("survival")
library("survminer")
fit <- survfit(Surv(time, status)~ predict, data = clin_km)
ggsurvplot(fit,data=clin_km,legend.title="Group",xlab="Months",
       pval = TRUE, #conf6.int = TRUE,
	 #pval = round(p2.val,2),
	 palette = c("#E41A1C","#4DAF4A","#377EB8"),
       #risk.table = TRUE, # Add risk table
       #risk.table.col = "strata", # Change risk table color by groups
       linetype = "strata", # Change line type by groups
       #ggtheme = theme_bw(), # Change ggplot2 theme
        )
surv_diff <- survdiff(Surv(time, status) ~ predict, data = clin_km)
p.val <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1) ##p-value from survdiff function
p3=ggsurvplot(fit,data=clin_km,legend.title="Group",xlab="Months",
       #pval = FALSE, #conf.int = TRUE,
	 #pval = round(p2.val,2),
	 palette = c("#E41A1C","#4DAF4A","#377EB8"),
       #risk.table = TRUE, # Add risk table
       #risk.table.col = "strata", # Change risk table color by groups
       linetype = "strata", # Change line type by groups
       #ggtheme = theme_bw(), # Change ggplot2 theme
        )
objCoxph <- coxph(Surv(time, status) ~ predict, data = clin_km)
p3$plot = p3$plot + ggplot2::annotate("text",x = 25, y = 0.15,size = 4,
                             label = paste0("HR = ", round(1/summary(objCoxph)$conf.int[1],2),
	"(",round(1/summary(objCoxph)$conf.int[4],2),"-",round(1/summary(objCoxph)$conf.int[3],2),")")) +
 	 ggplot2::annotate("text",x = 25, y = 0.25,size = 4,
      	              label = paste("p =",round(p.val,5)))
p3

ginvH=ginv(result$H[[1]])
W_new=test_data1%*%ginvH
dim(W_new)

ginvH=ginv(result$H[[2]])
W_new=test_data2%*%ginvH
dim(W_new)

ginvH=ginv(result$H[[3]])
W_new=test_data3%*%ginvH
dim(W_new)

clin_new<-a[,c(1,3,4,rep(1,ncol(W_new)))];dim(clin_new)
rownames(clin_new)=clin_new[,1]
clin_new=clin_new[co_sample,];dim(clin_new)
clin_new=clin_new[sample_test,]
for(i in 1:ncol(W_new)){clin_new[,3+i]=W_new[,i]}
clin_new<-na.omit(clin_new);dim(clin_new)
colnames(clin_new)<-c("sampleID","status","time",paste("module",1:ncol(W_new),sep=""))
clin_new=clin_new[,-1]
dd=datadist(clin)
options(datadist="dd")
concordance(cph(Surv(clin$time, clin$status) ~ predict(res.cox, newdata=clin)), newdata=clin)
concordance(cph(Surv(clin_new$time, clin_new$status) ~ predict(res.cox, newdata=clin_new)), newdata=clin_new) 
#KM-Plot
pre.result= predict(res.cox, newdata=clin_new)
clin_km<- clin_new[,c(1:ncol(clin_new),1)];dim(clin_km)
 med_zz=median(pre.result)
class=rep("low-risk",length(pre.result))
class[which(pre.result>med)]="high-risk"
clin_km[,ncol(clin_km)]=class
colnames(clin_km)[ncol(clin_km)]<-"predict"
clin_km$time=clin_km$time/30
fit <- survfit(Surv(time, status)~ predict, data = clin_km)
ggsurvplot(fit,data=clin_km,legend.title="Group",xlab="Months",
       pval = TRUE, #conf.int = TRUE,
	 #pval = round(p2.val,2),
	 palette = c("#E41A1C","#4DAF4A","#377EB8"),
       #risk.table = TRUE, # Add risk table
       #risk.table.col = "strata", # Change risk table color by groups
       linetype = "strata", # Change line type by groups
       #ggtheme = theme_bw(), # Change ggplot2 theme
        )
surv_diff <- survdiff(Surv(time, status) ~ predict, data = clin_km)
p.val <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1) ##p-value from survdiff function
p3=ggsurvplot(fit,data=clin_km,legend.title="Group",xlab="Months",
       #pval = FALSE, #conf.int = TRUE,
	 #pval = round(p2.val,2),
	 palette = c("#E41A1C","#4DAF4A","#377EB8"),
       #risk.table = TRUE, # Add risk table
       #risk.table.col = "strata", # Change risk table color by groups
       linetype = "strata", # Change line type by groups
       #ggtheme = theme_bw(), # Change ggplot2 theme
        )
objCoxph <- coxph(Surv(time, status) ~ predict, data = clin_km)
p3$plot = p3$plot + ggplot2::annotate("text",x = 25, y = 0.15,size = 4,
                             label = paste0("HR = ", round(1/summary(objCoxph)$conf.int[1],2),
	"(",round(1/summary(objCoxph)$conf.int[4],2),"-",round(1/summary(objCoxph)$conf.int[3],2),")")) +
 	 ggplot2::annotate("text",x = 25, y = 0.25,size = 4,
      	              label = paste("p =",round(p.val,4)))
p3
##########################  3-omics predict    ############################
# Gene
ginvH=ginv(result$H[[1]])
W_new=train_data1%*%ginvH
dim(W_new)
clin_new<-a[,c(1,3,4,rep(1,ncol(W_new)))];dim(clin_new)
rownames(clin_new)=clin_new[,1]
clin_new=clin_new[co_sample,];dim(clin_new)
clin_new=clin_new[c(sample_train),]  #
for(i in 1:ncol(W_new)){clin_new[,3+i]=W_new[,i]}
clin_new<-na.omit(clin_new);dim(clin_new)
colnames(clin_new)<-c("sampleID","status","time",paste("module",1:ncol(W_new),sep=""))
clin_new=clin_new[,-1]
dd=datadist(clin)
options(datadist="dd")
c1=concordance(cph(Surv(clin_new$time, clin_new$status) ~ predict(res.cox, newdata=clin_new)), newdata=clin_new)$concordance 
# Exon
ginvH=ginv(result$H[[2]])
W_new=train_data2%*%ginvH
dim(W_new)
clin_new<-a[,c(1,3,4,rep(1,ncol(W_new)))];dim(clin_new)
rownames(clin_new)=clin_new[,1]
clin_new=clin_new[co_sample,];dim(clin_new)
clin_new=clin_new[c(sample_train),]  #
for(i in 1:ncol(W_new)){clin_new[,3+i]=W_new[,i]}
clin_new<-na.omit(clin_new);dim(clin_new)
colnames(clin_new)<-c("sampleID","status","time",paste("module",1:ncol(W_new),sep=""))
clin_new=clin_new[,-1]
c2=concordance(cph(Surv(clin_new$time, clin_new$status) ~ predict(res.cox, newdata=clin_new)), newdata=clin_new)$concordance
# DM
ginvH=ginv(result$H[[3]])
W_new=train_data3%*%ginvH
dim(W_new)
clin_new<-a[,c(1,3,4,rep(1,ncol(W_new)))];dim(clin_new)
rownames(clin_new)=clin_new[,1]
clin_new=clin_new[co_sample,];dim(clin_new)
clin_new=clin_new[c(sample_train),]  #
for(i in 1:ncol(W_new)){clin_new[,3+i]=W_new[,i]}
clin_new<-na.omit(clin_new);dim(clin_new)
colnames(clin_new)<-c("sampleID","status","time",paste("module",1:ncol(W_new),sep=""))
clin_new=clin_new[,-1]
c3=concordance(cph(Surv(clin_new$time, clin_new$status) ~ predict(res.cox, newdata=clin_new)), newdata=clin_new)$concordance
k1=kappa(result$H[[1]],exact=T,norm="2")
k2=kappa(result$H[[2]],exact=T,norm="2")
k3=kappa(result$H[[3]],exact=T,norm="2")

alpha=c(c1,c2,c3)/sum(c(c1,c2,c3))
ginvH=ginv(result$H[[1]])
W_new1=test_data1%*%ginvH
ginvH2=ginv(result$H[[2]])
W_new2=test_data2%*%ginvH2
ginvH3=ginv(result$H[[3]])
W_new3=test_data3%*%ginvH3
W_new=(W_new1*alpha[1]+W_new2*alpha[2]+W_new3*alpha[3])

clin_new<-a[,c(1,3,4,rep(1,ncol(W_new)))];dim(clin_new)
rownames(clin_new)=clin_new[,1]
clin_new=clin_new[co_sample,];dim(clin_new)
clin_new=clin_new[sample_test,]
for(i in 1:ncol(W_new)){clin_new[,3+i]=W_new[,i]}
clin_new<-na.omit(clin_new);dim(clin_new)
colnames(clin_new)<-c("sampleID","status","time",paste("module",1:ncol(W_new),sep=""))
clin_new=clin_new[,-1]
dd=datadist(clin)
options(datadist="dd")
concordance(cph(Surv(clin$time, clin$status) ~ predict(res.cox, newdata=clin)), newdata=clin)
concordance(cph(Surv(clin_new$time, clin_new$status) ~ predict(res.cox, newdata=clin_new)), newdata=clin_new) 
# KM-Plot
pre.result= predict(res.cox, newdata=clin_new)
clin_km<- clin_new[,c(1:ncol(clin_new),1)];dim(clin_km)
med_zz=median(pre.result)
class=rep("low-risk",length(pre.result))
class[which(pre.result>med)]="high-risk"
clin_km[,ncol(clin_km)]=class
colnames(clin_km)[ncol(clin_km)]<-"predict"
clin_km$time=clin_km$time/30
fit <- survfit(Surv(time, status)~ predict, data = clin_km)
ggsurvplot(fit,data=clin_km,legend.title="Group",xlab="Months",
       pval = TRUE, #conf.int = TRUE,
	 #pval = round(p2.val,2),
	 palette = c("#E41A1C","#4DAF4A","#377EB8"),
       #risk.table = TRUE, # Add risk table
       #risk.table.col = "strata", # Change risk table color by groups
       linetype = "strata", # Change line type by groups
       #ggtheme = theme_bw(), # Change ggplot2 theme
        )

############  Combined with clinical  ###########################################
clin<-a[,c(1,3,4,rep(1,ncol(result$W)))];dim(clin)
rownames(clin)=clin[,1]
clin=clin[co_sample,];dim(clin)
clin=clin[c(sample_train),]
for(i in 1:ncol(result$W)){clin[,3+i]=result$W[,i]}
clin<-na.omit(clin);dim(clin)
colnames(clin)<-c("sampleID","status","time",paste("module",1:ncol(result$W),sep=""))
clin=clin[,-1]
ph_df <- del_nonPH(gene_list = colnames(clin)[c(3:ncol(clin))], survival_info_df = clin)
res.cox <- coxph(Surv(time, status) ~., data = clin[,c("status","time",as.character(ph_df[,1]))])
n=nrow(clin)
res.cox=stepAIC(res.cox,k = log(n))
summary(res.cox)
test.ph <- cox.zph(res.cox)
test.ph

pre.result= predict(res.cox, newdata=clin)
### Supplementary clinical information ###
clin_roc<-a[,c(1,3,4)];dim(clin_roc)
rownames(clin_roc)=clin_roc[,1]
clin_roc=clin_roc[co_sample,]
clin_roc=clin_roc[c(sample_train),]
dim(clin_roc)
clin_roc<-na.omit(clin_roc);dim(clin_roc)
colnames(clin_roc)<-c("sampleID","status","time")
clin_roc=clin_roc[,-1]
clin_roc2=cbind(clin_roc,pre.result)
a1<-read.csv("COADREAD_clinicalMatrix.csv",head=T,stringsAsFactors = FALSE)
clin_roc<-a1[,c(1,80,22)];dim(clin_roc)
rownames(clin_roc)<-clin_roc[,1]
clin_roc=clin_roc[co_sample,]
clin_roc=clin_roc[c(sample_train),]
clin_roc=cbind(clin_roc,clin_roc2)
clin_roc=na.omit(clin_roc)
clin_roc=clin_roc[,-1];dim(clin_roc)
clin_roc[which(clin_roc[,1]=="Stage I"),1]<-"Stage I-II"
clin_roc[which(clin_roc[,1]=="Stage II"),1]<-"Stage I-II"
clin_roc[which(clin_roc[,1]=="Stage III"),1]<-"Stage III-IV"
clin_roc[which(clin_roc[,1]=="Stage IV"),1]<-"Stage III-IV"
clin_roc[which(clin_roc[,1]=="Stage IA"),1]<-"Stage I-II"
clin_roc[which(clin_roc[,1]=="Stage IIA"),1]<-"Stage I-II"
clin_roc[which(clin_roc[,1]=="Stage IIB"),1]<-"Stage I-II"
clin_roc[which(clin_roc[,1]=="Stage IIC"),1]<-"Stage I-II"
clin_roc[which(clin_roc[,1]=="Stage IIIA"),1]<-"Stage III-IV"
clin_roc[which(clin_roc[,1]=="Stage IIIB"),1]<-"Stage III-IV"
clin_roc[which(clin_roc[,1]=="Stage IIIC"),1]<-"Stage III-IV"
clin_roc[which(clin_roc[,1]=="Stage IVA"),1]<-"Stage III-IV"
clin_roc[which(clin_roc[,1]=="Stage IVB"),1]<-"Stage III-IV"
ind_noTNM=which(clin_roc[,1]=="[Discrepancy]"|clin_roc[,1]=="")
clin_roc=clin_roc[-ind_noTNM,]
table(clin_roc[,1])
ROC.time2=timeROC(clin_roc$time, clin_roc$status, as.numeric(factor(clin_roc$pathologic_stage)), other_markers = NULL, cause=1,
	    weighting = "marginal",
		times=c(365*3,365*5),
		 ROC = TRUE, iid = F)

ROC.time_1=timeROC(clin_roc$time, clin_roc$status, pre.result[-ind_noTNM], other_markers = NULL, cause=1,
	    weighting = "marginal", times=c(365,365*2,365*3,365*4,365*5), ROC = TRUE, iid = TRUE)	
ROC.time_2=timeROC(clin_roc$time, clin_roc$status, as.numeric(factor(clin_roc$pathologic_stage)), other_markers = NULL, cause=1,
	    weighting = "marginal",
		times=c(365,365*2,365*3,365*4,365*5),ROC = TRUE, iid = TRUE)
clin_cl_train=clin_roc
res.cox_cl <- coxph(Surv(time, status) ~., data = clin_cl_train)
summary(res.cox_cl)
library(MASS)
res.cox_cl=stepAIC(res.cox_cl)
summary(res.cox_cl)
test.ph <- cox.zph(res.cox_cl)
test.ph

cl.predist=predict(res.cox_cl,data=clin_cl_train)
ROC.time_cl=timeROC(clin_cl_train$time, clin_cl_train$status,cl.predist , other_markers = NULL, cause=1,
	    weighting = "marginal",
		times=c(365,365*2,365*3,365*4,365*5),
		 ROC = TRUE, iid = T)
nyears=1
plot(ROC.time_cl,time=365*nyears,col="blue",title=F) 
plot(ROC.time_1,time=365*nyears,add=TRUE,col="red")
plot(ROC.time_2,time=365*nyears,add=TRUE,col="green") 
lebel_2=paste(c("AUC of combine all  = ","AUC of score = ","AUC of stage = "),
		c(round(ROC.time_cl$AUC[nyears],2),round(ROC.time_1$AUC[nyears],2),round(ROC.time_2$AUC[nyears],2)),sep=" ")
legend("bottomright",lebel_2,col=c("blue","red","green"),lty=1,lwd=2)
title(main="ROC at 1-year",cex=1,line=0.8)
#KM-Plot
pre.result= predict(res.cox_cl, newdata=clin_cl_train)
med_cl=median(pre.result)
class=rep("low-risk",length(pre.result))
class[which(pre.result>med_cl)]="high-risk"
clin_km<- clin_cl_train[,c(1:ncol(clin_cl_train),1)];dim(clin_km)
clin_km[,ncol(clin_km)]=class
colnames(clin_km)[ncol(clin_km)]<-"predict"
clin_km$time=clin_km$time/30
fit <- survfit(Surv(time, status)~ predict, data = clin_km)
surv_diff <- survdiff(Surv(time, status) ~ predict, data = clin_km)
p.val <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1) ##p-value from survdiff function
ggsurvplot(fit,data=clin_km,legend.title="Group",xlab="Months",
       pval = TRUE, #conf.int = TRUE,
	 #pval = round(p2.val,2),
	 palette = c("#E41A1C","#4DAF4A","#377EB8"),
       #risk.table = TRUE, # Add risk table
       #risk.table.col = "strata", # Change risk table color by groups
       linetype = "strata", # Change line type by groups
       #ggtheme = theme_bw(), # Change ggplot2 theme
        )
surv_diff <- survdiff(Surv(time, status) ~ predict, data = clin_km)
p.val <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1) ##p-value from survdiff function
p3=ggsurvplot(fit,data=clin_km,legend.title="Group",xlab="Months",
       #pval = FALSE, #conf.int = TRUE,
	 #pval = round(p2.val,2),
	 palette = c("#E41A1C","#4DAF4A","#377EB8"),
       #risk.table = TRUE, # Add risk table
       #risk.table.col = "strata", # Change risk table color by groups
       linetype = "strata", # Change line type by groups
       #ggtheme = theme_bw(), # Change ggplot2 theme
        )
objCoxph <- coxph(Surv(time, status) ~ predict, data = clin_km)
p3$plot = p3$plot + ggplot2::annotate("text",x = 25, y = 0.15,size = 4,
                             label = paste0("HR = ", round(1/summary(objCoxph)$conf.int[1],2),
	"(",round(1/summary(objCoxph)$conf.int[4],2),"-",round(1/summary(objCoxph)$conf.int[3],2),")")) +
 	 ggplot2::annotate("text",x = 25, y = 0.25,size = 4,
      	              label = paste0("p < 0.0001"))
p3

# test
ginvH=ginv(result$H[[1]])
W_new1=test_data1%*%ginvH
ginvH2=ginv(result$H[[2]])
W_new2=test_data2%*%ginvH2
ginvH3=ginv(result$H[[3]])
W_new3=test_data3%*%ginvH3
W_new=(W_new1*alpha[1]+W_new2*alpha[2]+W_new3*alpha[3])
clin_new<-a[,c(1,3,4,rep(1,ncol(W_new)))];dim(clin_new)
rownames(clin_new)=clin_new[,1]
clin_new=clin_new[co_sample,];dim(clin_new)
clin_new=clin_new[sample_test,]
for(i in 1:ncol(W_new)){clin_new[,3+i]=W_new[,i]}
clin_new<-na.omit(clin_new);dim(clin_new)
colnames(clin_new)<-c("sampleID","status","time",paste("module",1:ncol(W_new),sep=""))
clin_new=clin_new[,-1]
pre.result= predict(res.cox, newdata=clin_new)
### Supplementary clinical information ###
clin_roc<-a[,c(1,3,4)];dim(clin_roc)
rownames(clin_roc)=clin_roc[,1]
clin_roc=clin_roc[co_sample,]
clin_roc=clin_roc[sample_test,]
dim(clin_roc)
clin_roc<-na.omit(clin_roc);dim(clin_roc)
colnames(clin_roc)<-c("sampleID","status","time")
clin_roc=clin_roc[,-1]
clin_roc2=cbind(clin_roc,pre.result)
clin_roc<-a1[,c(1,80,22,84)];dim(clin_roc)
rownames(clin_roc)<-clin_roc[,1]
clin_roc=clin_roc[co_sample,]
clin_roc=clin_roc[sample_test,]
clin_roc=cbind(clin_roc,clin_roc2)
clin_roc=na.omit(clin_roc)
clin_roc=clin_roc[,-1];dim(clin_roc)
clin_roc[which(clin_roc[,1]=="Stage I"),1]<-"Stage I-II"
clin_roc[which(clin_roc[,1]=="Stage II"),1]<-"Stage I-II"
clin_roc[which(clin_roc[,1]=="Stage III"),1]<-"Stage III-IV"
clin_roc[which(clin_roc[,1]=="Stage IV"),1]<-"Stage III-IV"
clin_roc[which(clin_roc[,1]=="Stage IA"),1]<-"Stage I-II"
clin_roc[which(clin_roc[,1]=="Stage IIA"),1]<-"Stage I-II"
clin_roc[which(clin_roc[,1]=="Stage IIB"),1]<-"Stage I-II"
clin_roc[which(clin_roc[,1]=="Stage IIC"),1]<-"Stage I-II"
clin_roc[which(clin_roc[,1]=="Stage IIIA"),1]<-"Stage III-IV"
clin_roc[which(clin_roc[,1]=="Stage IIIB"),1]<-"Stage III-IV"
clin_roc[which(clin_roc[,1]=="Stage IIIC"),1]<-"Stage III-IV"
clin_roc[which(clin_roc[,1]=="Stage IVA"),1]<-"Stage III-IV"
clin_roc[which(clin_roc[,1]=="Stage IVB"),1]<-"Stage III-IV"

ind_noTNM=which(clin_roc[,1]=="[Discrepancy]"|clin_roc[,1]=="")
clin_roc=clin_roc[-ind_noTNM,]
table(clin_roc[,1])


ROC.time_1=timeROC(clin_roc$time, clin_roc$status, pre.result[-ind_noTNM], other_markers = NULL, cause=1,
	    weighting = "marginal", times=c(365,365*2,365*3,365*4,365*5), ROC = TRUE, iid = TRUE)	
ROC.time_2=timeROC(clin_roc$time, clin_roc$status, as.numeric(factor(clin_roc$pathologic_stage)), other_markers = NULL, cause=1,
	    weighting = "marginal",
		times=c(365,365*2,365*3,365*4,365*5),ROC = TRUE, iid = TRUE)
clin_cl_test=clin_roc
cl.predist=predict(res.cox_cl,newdata=clin_cl_test)
ROC.time_cl=timeROC(clin_cl_test$time, clin_cl_test$status,cl.predist , other_markers = NULL, cause=1,
	    weighting = "marginal",
		times=c(365,365*2,365*3,365*4,365*5),
		 ROC = TRUE, iid = T)
nyears=1
plot(ROC.time_cl,time=365*nyears,col="blue",title=F) 
plot(ROC.time_1,time=365*nyears,add=TRUE,col="red")
plot(ROC.time_2,time=365*nyears,add=TRUE,col="green") 
title(main="ROC at 1-year",cex=2,line=0.8)
lebel_2=paste(c("AUC of combine all = ","AUC of score = ","AUC of stage = "),
		c(round(ROC.time_cl$AUC[nyears],2),round(ROC.time_1$AUC[nyears],2),round(ROC.time_2$AUC[nyears],2)),sep=" ")
legend("bottomright",lebel_2,col=c("blue","red","green"),lty=1,lwd=2)
#KM-Plot
pre.result= predict(res.cox_cl, newdata=clin_cl_test)
class=rep("low-risk",length(pre.result))
class[which(pre.result>med_cl)]="high-risk"
clin_km<- clin_cl_test[,c(1:ncol(clin_cl_test),1)];dim(clin_km)
clin_km[,ncol(clin_km)]=class
colnames(clin_km)[ncol(clin_km)]<-"predict"
clin_km$time=clin_km$time/30
fit <- survfit(Surv(time, status)~ predict, data = clin_km)
surv_diff <- survdiff(Surv(time, status) ~ predict, data = clin_km)
p.val <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1) ##p-value from survdiff function
ggsurvplot(fit,data=clin_km,legend.title="Group",xlab="Months",
       pval = TRUE, #conf.int = TRUE,
	 #pval = round(p2.val,2),
	 palette = c("#E41A1C","#4DAF4A","#377EB8"),
       #risk.table = TRUE, # Add risk table
       #risk.table.col = "strata", # Change risk table color by groups
       linetype = "strata", # Change line type by groups
       #ggtheme = theme_bw(), # Change ggplot2 theme
        )
dd=datadist(clin_cl_train)
options(datadist="dd")
concordance(cph(Surv(clin_cl_train$time, clin_cl_train$status) ~ predict(res.cox_cl, newdata=clin_cl_train)), newdata=clin_cl_train)
concordance(cph(Surv(clin_cl_test$time, clin_cl_test$status) ~ predict(res.cox_cl, newdata=clin_cl_test)), newdata=clin_cl_test) 
###################### Compare LASSO
I=3  # I=  1,2,3
ls_data=RNMF_data3
library(glmnet)
set.seed(2200)
	clin<-a[,c(1,3,4)];dim(clin)
	rownames(clin)=clin[,1]
	clin=clin[co_sample,1:3];dim(clin)
	clin=clin[sample_train,] # sample.index
	data.out=Surv(clin$OS.time+1,clin$OS)
	if(I==1) {data=RNMF_data1[sample_train,]} 
	if(I==2) {data=RNMF_data2[sample_train,]} 
	if(I==3) {data=RNMF_data3[sample_train,]} 
	fit=cv.glmnet(data,data.out,family="cox",trace.it=1,nfolds=3)
	plot(fit)
	fit.coef.lambda.lse=coef(fit,s=fit$lambda.min)
	fit.1se.out=fit.coef.lambda.lse[which(fit.coef.lambda.lse!=0)]
	fit.1se.out=round(fit.1se.out,4)
	fit.1se.out2=matrix(fit.1se.out,length(fit.1se.out),1)
	rownames(fit.1se.out2)=rownames(fit.coef.lambda.lse)[which(fit.coef.lambda.lse!=0)]
	colnames(fit.1se.out2)=c("coef")
	lasso_id=rownames(fit.1se.out2)
	lasso_id


clin<-a[,c(1,3,4)];dim(clin)
rownames(clin)=clin[,1]
clin=clin[sample_train,];dim(clin)
clin=cbind(clin,ls_data[sample_train,lasso_id])
clin<-na.omit(clin);dim(clin)
colnames(clin)[1:3]<-c("sampleID","status","time")
clin=clin[,-1]
dd=datadist(clin)
options(datadist="dd")
pre_tcga=as.numeric(as.matrix(ls_data[sample_train,])%*% coef(fit,s=fit$lambda.min))
concordance(cph(Surv(clin$time, clin$status) ~ pre_tcga), newdata=clin)
pre_new=as.numeric(as.matrix(ls_data[sample_test,])%*% coef(fit,s=fit$lambda.min))
concordance(cph(Surv(clin_new$time, clin_new$status) ~ pre_new), newdata=clin_new)



