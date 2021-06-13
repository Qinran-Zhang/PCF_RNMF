########## sample #######
library(survival);library(survminer);library(rms);library(timeROC);library(polspline)
library(data.table);library(MASS)
setwd("C:/Users/MSI-NB/Desktop/code/CORE-DFI")
load("CORE-DFI.Rdata")
source("RNMF.R")
RNMF_data=list(train_data1,train_data2,train_data3)
K_select=9;lambda_select=10^(-8)

set.seed(204754)
t1=Sys.time()
result=RNMF(RNMF_data,K=K_select,lambda=lambda_select,maxiter=10000,speak=T)
Sys.time()-t1
clin<-a[,c(1,2,3,rep(1,ncol(result$W)))];dim(clin)
rownames(clin)=clin[,1]
clin=clin[co_sample,];dim(clin)
for(i in 1:ncol(result$W)){clin[,3+i]=result$W[,i]}
clin<-na.omit(clin);dim(clin)
colnames(clin)<-c("sampleID","status","time",paste("module",1:ncol(result$W),sep=""))
clin=clin[,-1]
ph_df <- del_nonPH(gene_list = colnames(clin)[c(3:ncol(clin))], survival_info_df = clin)
res.cox <- coxph(Surv(time, status) ~., data = clin[,c("status","time",as.character(ph_df [,1]))])
n=nrow(clin)
res.cox=stepAIC(res.cox,k = log(n))
summary(res.cox)
test.ph <- cox.zph(res.cox)
test.ph
#timeROC
pre.result=predict(res.cox,newdata=clin)
ROC.time1=timeROC(clin$time, clin$status, pre.result, other_markers = NULL, cause=1,
	    weighting = "marginal",
		times=c(365*1,365*2,365*3,365*4,365*5),
		 ROC = TRUE, iid = F) 
plot(ROC.time1,time=365*3,col="red",title=F)     
title(main=paste("AUC of 3-years =",as.numeric(round(ROC.time1$AUC[3],3)),seq=""),cex=2)
#KM-Plot
pre.result= predict(res.cox, newdata=clin)
med=median(pre.result)
class=rep("low-risk",length(pre.result))
class[which(pre.result>med)]="high-risk"
clin_km<- clin[,c(1:ncol(clin),1)];dim(clin_km)
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
p3$plot = p3$plot + ggplot2::annotate("text",x = 30, y = 0.15,size = 4,
                             label = paste0("HR = ", round(1/summary(objCoxph)$conf.int[1],2),
	"(",round(1/summary(objCoxph)$conf.int[4],2),"-",round(1/summary(objCoxph)$conf.int[3],2),")")) +
 	 ggplot2::annotate("text",x = 30, y = 0.25,size = 4,
      	              label = paste("p =",round(p.val,5)))
p3

ZZ_data=ZZ_data[,colnames(RNMF_data[[3]])];dim(ZZ_data)
ginvH=ginv(result$H[[3]])
zz_data=as.matrix(ZZ_data)
zz_data=t((t(ZZ_data[,colnames(RNMF_data[[3]])])-std_min3)/(std_max3-std_min3)*1000/F3)
W_new=zz_data%*%ginvH
dim(W_new)
clin_new<-b[,c(1,21,22,rep(1,ncol(W_new)))];dim(clin_new)
rownames(clin_new)=clin_new[,1]
clin_new=clin_new[rownames(W_new),]
dim(clin_new)
for(i in 1:ncol(W_new)){
clin_new[,3+i]=W_new[,i]}
clin_new<-na.omit(clin_new);dim(clin_new)
colnames(clin_new)<-c("sampleID","status","time",paste("module",1:ncol(W_new),sep=""))
clin_new=clin_new[,-1]
clin_new$time=clin_new$time*30
dd=datadist(clin)
options(datadist="dd")
concordance(cph(Surv(clin$time, clin$status) ~ predict(res.cox, newdata=clin)), newdata=clin)
concordance(cph(Surv(clin_new$time, clin_new$status) ~ predict(res.cox, newdata=clin_new)), newdata=clin_new)
#ROC
library(timeROC)
pre.result= predict(res.cox, newdata=clin_new)
ROC.time=timeROC(clin_new$time,clin_new$status, pre.result, other_markers = NULL, cause=1,
	    weighting = "marginal",
		times=c(365*1,365*2,365*3,365*4,365*5),
		 ROC = TRUE, iid = F) ## 
plot(ROC.time,time=365*3,title=F)     
title(main=paste("AUC of 3-years =",as.numeric(round(ROC.time$AUC[3],3)),seq=""),cex=2)
#KM-Plot
pre.result= predict(res.cox, newdata=clin_new)
clin_km<- clin_new[,c(1:ncol(clin_new),1)];dim(clin_km)
 med_zz=median(pre.result)
class=rep("low-risk",length(pre.result))
class[which(pre.result>med_zz)]="high-risk"
clin_km[,ncol(clin_km)]=class
colnames(clin_km)[ncol(clin_km)]<-"predict"
clin_km$time=clin_km$time/30
fit <- survfit(Surv(time, status)~ predict, data = clin_km)
surv_diff <- survdiff(Surv(time, status) ~ predict, data = clin_km)
p.val <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1) ##p-value from survdiff function
ggsurvplot(fit,data=clin_km,legend.title="Group",xlab="Months",
       pval = TRUE, #conf.int = TRUE,
	 #pval = round(p.val,5),
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

n=3
fev3 <- cph(Surv(clin$time, clin$status) ~ predict(res.cox, newdata=clin), x=T, y=T, surv=T, data=clin, time.inc=365*n)
calev3 <- calibrate(fev3, cmethod="KM", method="boot", u=365*n, m=42)
plot(calev3,subtitles=F) #,xlim=c(0.5,1),ylim=c(0.5,1)
fev3 <- cph(Surv(clin_new$time, clin_new$status) ~ predict(res.cox, newdata=clin_new), x=T, y=T, surv=T, data=clin_new, time.inc=365*n)
calev3 <- calibrate(fev3, cmethod="KM", method="boot", u=365*n, m=38)
plot(calev3,subtitles=F) #,xlim=c(0.8,1),ylim=c(0.8,1)
######################Further screening##############
summary(res.cox)
module=findmodule(result$W,result$H,tred=c(1,3))
idn=c(4,9);ID3=c();for(idni in idn){ID3=c(ID3,module[[4]][[idni]])};ID3=unique(ID3)
library(glmnet)
set.seed(153363) #153363
	clin<-a[,c(1,2,3)];dim(clin)
	rownames(clin)=clin[,1]
	clin=clin[co_sample,1:3];dim(clin)
	data.out=Surv(clin$DFI.time,clin$DFI)
	 {data=RNMF_data3[,ID3];k=1} #0.9
	fit=cv.glmnet(data,data.out,family="cox",trace.it=1,nfolds=3)#,foldid=foldid
	plot(fit)
	fit.coef.lambda.lse=coef(fit,s=fit$lambda.1se)
	#if(I==2) fit.coef.lambda.lse=coef(fit,s=fit$lambda.1se)
	fit.1se.out=fit.coef.lambda.lse[which(fit.coef.lambda.lse!=0)]
	fit.1se.out=round(fit.1se.out,4)
	fit.1se.out2=matrix(fit.1se.out,length(fit.1se.out),1)
	rownames(fit.1se.out2)=rownames(fit.coef.lambda.lse)[which(fit.coef.lambda.lse!=0)]
	colnames(fit.1se.out2)=c("coef")
	lasso_id=rownames(fit.1se.out2)
	lasso_id

clin<-a[,c(1,2,3)];dim(clin)
rownames(clin)=clin[,1]
clin=clin[co_sample,];dim(clin)
clin=cbind(clin,RNMF_data3[co_sample,lasso_id])
clin<-na.omit(clin);dim(clin)
colnames(clin)[1:3]<-c("sampleID","status","time")
clin=clin[,-1]
res.cox <- coxph(Surv(time, status) ~., data = clin)
test.ph <- cox.zph(res.cox)
test.ph
ph_df  <- del_nonPH(gene_list = colnames(clin)[c(3:ncol(clin))], survival_info_df = clin)
n=nrow(clin)
res.cox=stepAIC(res.cox,k = log(n))
summary(res.cox)
test.ph <- cox.zph(res.cox)
test.ph
ZZ_data=ZZ_data[,colnames(RNMF_data[[3]])];dim(ZZ_data)
ginvH=ginv(result$H[[3]])
zz_data=as.matrix(ZZ_data)
clin_new<-b[,c(1,21,22)];dim(clin_new)
rownames(clin_new)=clin_new[,1]
clin_new=clin_new[rownames(W_new),]
dim(clin_new)
clin_new=cbind(clin_new,zz_data)
clin_new<-na.omit(clin_new);dim(clin_new)
colnames(clin_new)[1:3]<-c("sampleID","status","time")
clin_new=clin_new[,-1]
clin_new$time=clin_new$time*30
dd=datadist(clin)
options(datadist="dd")
concordance(cph(Surv(clin$time, clin$status) ~ predict(res.cox, newdata=clin)), newdata=clin)
concordance(cph(Surv(clin_new$time, clin_new$status) ~ predict(res.cox, newdata=clin_new)), newdata=clin_new)

#########Combined with clinical####################
pre.result11=RNMF_data3[,c("cg25053798","cg01117833")]
pre.result22=ZZ_data[,c("cg25053798","cg01117833")]
pre.result=pre.result11
### Supplementary clinical information ###
clin_roc<-a[,c(1,2,3)];dim(clin_roc)
rownames(clin_roc)=clin_roc[,1]
clin_roc=clin_roc[co_sample,]
dim(clin_roc)
clin_roc<-na.omit(clin_roc);dim(clin_roc)
colnames(clin_roc)<-c("sampleID","status","time")
clin_roc=clin_roc[,-1]
clin_roc2=cbind(clin_roc,pre.result)
a1<-read.csv("COADREAD_clinicalMatrix.csv",head=T,stringsAsFactors = FALSE)
clin_roc<-a1[,c(1,78,79,83,106,80,22,42)];dim(clin_roc)# ,22,84,85,86,99,105
rownames(clin_roc)<-clin_roc[,1]
clin_roc=clin_roc[co_sample,]
clin_roc=cbind(clin_roc,clin_roc2)
clin_roc=na.omit(clin_roc)
clin_roc=clin_roc[,-1];dim(clin_roc)
clin_roc[which(clin_roc[,1]=="N0"),1]<-"N1"
clin_roc[which(clin_roc[,1]=="N1"),1]<-"N1"
clin_roc[which(clin_roc[,1]=="N1a"),1]<-"N1"
clin_roc[which(clin_roc[,1]=="N1b"),1]<-"N1"
clin_roc[which(clin_roc[,1]=="N1c"),1]<-"N1"
clin_roc[which(clin_roc[,1]=="N2"),1]<-"N2"
clin_roc[which(clin_roc[,1]=="N2a"),1]<-"N2"
clin_roc[which(clin_roc[,1]=="N2b"),1]<-"N2"
ind_noTNM=which(clin_roc[,1]=="NX")
clin_roc=clin_roc[-ind_noTNM,]
table(clin_roc[,1])
clin_roc[which(clin_roc[,2]=="T4a"),2]<-"T4"
clin_roc[which(clin_roc[,2]=="T4b"),2]<-"T4"
clin_roc[which(clin_roc[,2]=="T1"),2]<-"T1-2"
clin_roc[which(clin_roc[,2]=="T2"),2]<-"T1-2"
clin_roc[which(clin_roc[,2]=="T3"),2]<-"T3-4"
clin_roc[which(clin_roc[,2]=="T4"),2]<-"T3-4"
table(clin_roc[,2])
clin_roc[which(clin_roc[,5]=="Stage I"),5]<-"I"
clin_roc[which(clin_roc[,5]=="Stage II"),5]<-"II"
clin_roc[which(clin_roc[,5]=="Stage IIA"),5]<-"II"
clin_roc[which(clin_roc[,5]=="Stage IIB"),5]<-"II"
clin_roc[which(clin_roc[,5]=="Stage IIC"),5]<-"II"
clin_roc[which(clin_roc[,5]=="Stage III"),5]<-"III"
clin_roc[which(clin_roc[,5]=="Stage IIIA"),5]<-"III"
clin_roc[which(clin_roc[,5]=="Stage IIIB"),5]<-"III"
clin_roc[which(clin_roc[,5]=="Stage IIIC"),5]<-"III"
table(clin_roc[,5])
ind_noTNM=which(clin_roc[,5]==""|clin_roc[,5]=="[Discrepancy]")
clin_roc=clin_roc[-ind_noTNM,]
colnames(clin_roc)[1:7]=c("N_stage","T_stage","PI","VI","stage","Age","gender")
clin_cl_train=clin_roc
res.cox_cl <- coxph(Surv(time, status) ~cg25053798+cg01117833+N_stage+PI+VI+stage+T_stage,
	 data = clin_cl_train)
res.cox_cl=stepAIC(res.cox_cl) 
summary(res.cox_cl)
test.ph <- cox.zph(res.cox_cl)
test.ph
# test
pre.result=pre.result22
### Supplementary clinical information ###
clin_roc<-b[,c(1,10,11,13,14,15,21,22,16,17,10)];dim(clin_roc)
rownames(clin_roc)=clin_roc[,1]
clin_roc=clin_roc[rownames(W_new),]
colnames(clin_roc)<-c("sampleID","Gender","Age","stage","N_stage","T_stage","status","time","VI","PI","gender")
clin_roc=clin_roc[,-1]
clin_roc=cbind(clin_roc,pre.result)
clin_roc[which(clin_roc[,4]=="0"),4]<-"N1"
clin_roc[which(clin_roc[,4]=="1b"),4]<-"N1"
clin_roc[which(clin_roc[,4]=="1a"),4]<-"N1"
clin_roc[which(clin_roc[,4]=="2a"),4]<-"N2"
clin_roc[which(clin_roc[,4]=="2b"),4]<-"N2"
table(clin_roc[,4])
clin_roc[which(clin_roc[,5]=="1"),5]<-"T1"
clin_roc[which(clin_roc[,5]=="2"),5]<-"T2"
clin_roc[which(clin_roc[,5]=="3"),5]<-"T3"
clin_roc[which(clin_roc[,5]=="4a"),5]<-"T4"
clin_roc[which(clin_roc[,5]=="4b"),5]<-"T4"
clin_roc[which(clin_roc[,5]=="T1"),5]<-"T1-2"
clin_roc[which(clin_roc[,5]=="T2"),5]<-"T1-2"
clin_roc[which(clin_roc[,5]=="T3"),5]<-"T3-4"
clin_roc[which(clin_roc[,5]=="T4"),5]<-"T3-4"
table(clin_roc[,5])
table(clin_roc[,8])
table(clin_roc[,9])
clin_roc[which(clin_roc[,3]=="IIB"),3]<-"II"
clin_roc[which(clin_roc[,3]=="IIA"),3]<-"II"
clin_roc[which(clin_roc[,3]=="IIC"),3]<-"II"
clin_roc[which(clin_roc[,3]=="IIIA"),3]<-"III"
clin_roc[which(clin_roc[,3]=="IIIB"),3]<-"III"
clin_roc[which(clin_roc[,3]=="IIIC"),3]<-"III"
ind_noTNM=which(clin_roc[,3]=="IVB")
clin_roc=clin_roc[-ind_noTNM,]
clin_roc[,8]=factor(ifelse(clin_roc[,8]==2,"YES","NO"))
clin_roc[,9]=factor(ifelse(clin_roc[,9]==2,"YES","NO"))
clin_cl_test=clin_roc
dd=datadist(clin_cl_train)
options(datadist="dd")
concordance(cph(Surv(clin_cl_train$time, clin_cl_train$status) ~ predict(res.cox_cl, newdata=clin_cl_train)), newdata=clin_cl_train)
concordance(cph(Surv(clin_cl_test$time, clin_cl_test$status) ~ predict(res.cox_cl, newdata=clin_cl_test)), newdata=clin_cl_test) 
clin_cl_test$time=clin_cl_test$time*30
pre.result_plot=predict(res.cox_cl, newdata=clin_cl_test)
ROC.time=timeROC(clin_cl_test$time, clin_cl_test$status, pre.result_plot, other_markers = NULL, cause=1,
	    weighting = "marginal",
		times=c(365*1,365*2,365*3,365*4,365*5),
		 ROC = TRUE, iid = F) 
plot(ROC.time,time=365*3,title=F)     
title(main=paste("AUC of 3-years =",as.numeric(round(ROC.time$AUC[3],3)),seq=""),cex=2)
dd=datadist(clin_cl_train)
f=cph(Surv(clin_cl_train$time, clin_cl_train$status) ~cg25053798+cg01117833+N_stage, data = clin_cl_train,
	x=TRUE,y=TRUE,surv=TRUE)

	survival=Survival(f)
	survival1=function(x) survival(365*1,x)
	survival2=function(x) survival(365*3,x)
	survival3=function(x) survival(365*5,x)
	nom=nomogram(f,fun=list(survival1,survival2,survival3),lp=FALSE,
		fun.at=c(0,seq(0.1,1,by=0.2),1),
		funlabel=c("1 year survival","3 years survival","5 years survival")
		)
plot(nom)
nyears=3
fev3 <- cph(Surv(clin_cl_test$time, clin_cl_test$status) ~ predict(res.cox_cl, newdata=clin_cl_test), x=T, y=T, surv=T, data=clin_cl_test, time.inc=365*nyears)
calev3 <- calibrate(fev3, cmethod="KM", method="boot", u=365*nyears, m=38,surv=TRUE, time.inc=u)
plot(calev3,subtitles=F) 
fev3 <- cph(Surv(clin_cl_train$time, clin_cl_train$status) ~ predict(res.cox_cl, newdata=clin_cl_train), x=T, y=T, surv=T, data=clin_cl_train, time.inc=365*nyears)
calev3 <- calibrate(fev3, cmethod="KM", method="boot", u=365*nyears, m=39,surv=TRUE, time.inc=u)
plot(calev3,subtitles=F) 
# KM-Plot
pre.result= predict(res.cox_cl, newdata=clin_cl_train)
clin_km<- clin_cl_train;dim(clin_km)
 med=median(pre.result)
class=rep("low-risk",length(pre.result))
class[which(pre.result>med)]="high-risk"
clin_km[,ncol(clin_km)]=class
clin_km$time=clin_km$time/30
colnames(clin_km)[ncol(clin_km)]<-"predict"
fit <- survfit(Surv(time, status)~ predict, data = clin_km)
summary(fit)$table
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
p3$plot = p3$plot + ggplot2::annotate("text",x = 30, y = 0.15,size = 4,
                             label = paste0("HR = ", round(1/summary(objCoxph)$conf.int[1],2),
	"(",round(1/summary(objCoxph)$conf.int[4],2),"-",round(1/summary(objCoxph)$conf.int[3],2),")")) +
 	 ggplot2::annotate("text",x = 30, y = 0.25,size = 4,
      	              label = paste("p =",round(p.val,4)))
p3


pre.result= predict(res.cox_cl, newdata=clin_cl_test)
clin_km<- clin_cl_test[,c(1:ncol(clin_cl_test),1)];dim(clin_km)
 med_zz=median(pre.result)
class=rep("low-risk",length(pre.result))
class[which(pre.result>med_zz)]="high-risk"
clin_km[,ncol(clin_km)]=class
colnames(clin_km)[ncol(clin_km)]<-"predict"
clin_km$time=clin_km$time/30
fit <- survfit(Surv(time, status)~ predict, data = clin_km)
summary(fit)$table
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
p3$plot = p3$plot + ggplot2::annotate("text",x = 30, y = 0.15,size = 4,
                             label = paste0("HR = ", round(1/summary(objCoxph)$conf.int[1],2),
	"(",round(1/summary(objCoxph)$conf.int[4],2),"-",round(1/summary(objCoxph)$conf.int[3],2),")")) +
 	 ggplot2::annotate("text",x = 30, y = 0.25,size = 4,
      	              label = paste("p =",round(p.val,5)))
p3
pre.result=predict(res.cox_cl, newdata=clin_cl_train)
ROC.time1=timeROC(clin_cl_train$time, clin_cl_train$status, pre.result, other_markers = NULL, cause=1,
	    weighting = "marginal",
		times=c(365*1,365*2,365*3,365*4,365*5),
		 ROC = TRUE, iid = F) 
plot(ROC.time1,time=365*3,col="red",title=F)     
title(main=paste("AUC of 3-years =",as.numeric(round(ROC.time1$AUC[3],3)),seq=""),cex=2)
pre.result= predict(res.cox_cl, newdata=clin_cl_test)
#clin_cl_test$time=clin_cl_test$time*30
ROC.time=timeROC(clin_cl_test$time,clin_cl_test$status, pre.result, other_markers = NULL, cause=1,
	    weighting = "marginal",
		times=c(365*1,365*2,365*3,365*4,365*5),
		 ROC = TRUE, iid = F) 
plot(ROC.time,time=365*3,title=F)     
title(main=paste("AUC of 3-years =",as.numeric(round(ROC.time$AUC[3],3)),seq=""),cex=2)
#######ROC
colnames(clin_cl_train)
res.cox0 <- coxph(Surv(time, status) ~., data =clin_cl_train )
summary(res.cox0)
dd=datadist(clin_cl_train)
options(datadist="dd")
concordance(cph(Surv(clin_cl_train$time, clin_cl_train$status) ~ predict(res.cox0, newdata=clin_cl_train)), newdata=clin_cl_train )
concordance(cph(Surv(clin_cl_test$time, clin_cl_test$status) ~ predict(res.cox0, newdata=clin_cl_test)), newdata=clin_cl_test )
res.cox1 <- coxph(Surv(time, status) ~., data = clin_cl_train[,c(1,8:9,10,11)] )
summary(res.cox1)
dd=datadist(clin_cl_train)
options(datadist="dd")
concordance(cph(Surv(clin_cl_train$time, clin_cl_train$status) ~ predict(res.cox1, newdata=clin_cl_train)), newdata=clin_cl_train)
concordance(cph(Surv(clin_cl_test$time, clin_cl_test$status) ~ predict(res.cox1, newdata=clin_cl_test)), newdata=clin_cl_test)
res.cox_cg <- coxph(Surv(time, status) ~., data = clin_cl_train[,c(8:11)] )
summary(res.cox_cg)
dd=datadist(clin_cl_train)
options(datadist="dd")
concordance(cph(Surv(clin_cl_train$time, clin_cl_train$status) ~ predict(res.cox_cg, newdata=clin_cl_train)), newdata=clin_cl_train)
concordance(cph(Surv(clin_cl_test$time, clin_cl_test$status) ~ predict(res.cox_cg, newdata=clin_cl_test)), newdata=clin_cl_test)
res.cox_stage <- coxph(Surv(time, status) ~., data = clin_cl_train[,c(5,8:9)] )
summary(res.cox_stage)
dd=datadist(clin_cl_train)
options(datadist="dd")
concordance(cph(Surv(clin_cl_train$time, clin_cl_train$status) ~ predict(res.cox_stage, newdata=clin_cl_train)), newdata=clin_cl_train)
concordance(cph(Surv(clin_cl_test$time, clin_cl_test$status) ~ predict(res.cox_stage, newdata=clin_cl_test)), newdata=clin_cl_test)
res.cox_clin <- coxph(Surv(time, status) ~., data = clin_cl_train[,c(1,8:9)] )
summary(res.cox_clin)
dd=datadist(clin_cl_train)
options(datadist="dd")
concordance(cph(Surv(clin_cl_train$time, clin_cl_train$status) ~ predict(res.cox_clin, newdata=clin_cl_train)), newdata=clin_cl_train)
concordance(cph(Surv(clin_cl_test$time, clin_cl_test$status) ~ predict(res.cox_clin, newdata=clin_cl_test)), newdata=clin_cl_test)
library("RColorBrewer")
#display.brewer.all()
display.brewer.pal(n = 12, name = 'Paired')
mycolor<-brewer.pal(12,"Paired")
par(mfrow=c(1,2))
pre.result_cg=predict(res.cox_cg, newdata=clin_cl_train)
pre.result_combine=predict(res.cox1, newdata=clin_cl_train)
pre.result_clin=predict(res.cox_clin, newdata=clin_cl_train)
pre.result_stage=predict(res.cox_stage, newdata=clin_cl_train)
ROC.time_1=timeROC(clin_cl_train$time, clin_cl_train$status, pre.result_cg, other_markers = NULL, cause=1,
	    weighting = "marginal", times=c(365,365*2,365*3,365*4,365*5), ROC = TRUE, iid = TRUE)	
ROC.time_2=timeROC(clin_cl_train$time, clin_cl_train$status, pre.result_combine, other_markers = NULL, cause=1,
	    weighting = "marginal",times=c(365,365*2,365*3,365*4,365*5),ROC = TRUE, iid = TRUE)
ROC.time_3=timeROC(clin_cl_train$time, clin_cl_train$status, pre.result_clin, other_markers = NULL, cause=1,
	    weighting = "marginal", times=c(365,365*2,365*3,365*4,365*5), ROC = TRUE, iid = TRUE)	
ROC.time_4=timeROC(clin_cl_train$time, clin_cl_train$status, pre.result_stage, other_markers = NULL, cause=1,
	    weighting = "marginal", times=c(365,365*2,365*3,365*4,365*5), ROC = TRUE, iid = TRUE)	
nyears=3
plot(ROC.time_1,time=365*nyears,col=mycolor[2],lwd=2,title=F)
plot(ROC.time_2,time=365*nyears,add=TRUE,col=mycolor[6],lwd=2) 
plot(ROC.time_3,time=365*nyears,add=TRUE,col=mycolor[4],lwd=2) 
plot(ROC.time_4,time=365*nyears,add=TRUE,col=mycolor[8],lwd=2) 
lebel_2=paste(c("AUC of combine=","AUC of 2 makers=","AUC of N_stage=",
			"AUC of stage="),
		c(round(ROC.time_2$AUC[nyears],2),round(ROC.time_1$AUC[nyears],2),round(ROC.time_3$AUC[nyears],2),
		round(ROC.time_4$AUC[nyears],2)
			),sep=" ")
legend("bottomright",lebel_2,col=mycolor[c(6,2,4,8)],lty=1,lwd=2)
compare(ROC.time_1,ROC.time_2) #compute p-values of comparison tests
### test set
pre.result_cg=predict(res.cox_cg, newdata=clin_cl_test)
pre.result_combine=predict(res.cox1, newdata=clin_cl_test)
pre.result_clin=predict(res.cox_clin, newdata=clin_cl_test)
pre.result_stage=predict(res.cox_stage, newdata=clin_cl_test)
ROC.time_1=timeROC(clin_cl_test$time, clin_cl_test$status, pre.result_cg, other_markers = NULL, cause=1,
	    weighting = "marginal", times=c(365,365*2,365*3,365*4,365*5), ROC = TRUE, iid = TRUE)	
ROC.time_2=timeROC(clin_cl_test$time, clin_cl_test$status, pre.result_combine, other_markers = NULL, cause=1,
	    weighting = "marginal",times=c(365,365*2,365*3,365*4,365*5),ROC = TRUE, iid = TRUE)
ROC.time_3=timeROC(clin_cl_test$time, clin_cl_test$status, pre.result_clin, other_markers = NULL, cause=1,
	    weighting = "marginal", times=c(365,365*2,365*3,365*4,365*5), ROC = TRUE, iid = TRUE)	
ROC.time_4=timeROC(clin_cl_test$time, clin_cl_test$status, pre.result_stage, other_markers = NULL, cause=1,
	    weighting = "marginal", times=c(365,365*2,365*3,365*4,365*5), ROC = TRUE, iid = TRUE)	
plot(ROC.time_1,time=365*nyears,col=mycolor[2],lwd=2,title=F)
plot(ROC.time_2,time=365*nyears,add=TRUE,col=mycolor[6],lwd=2) 
plot(ROC.time_3,time=365*nyears,add=TRUE,col=mycolor[4],lwd=2) 
plot(ROC.time_4,time=365*nyears,add=TRUE,col=mycolor[8],lwd=2) 
lebel_2=paste(c("AUC of combine=","AUC of 2 makers=","AUC of N_stage=",
			"AUC of stage="),
		c(round(ROC.time_2$AUC[nyears],2),round(ROC.time_1$AUC[nyears],2),round(ROC.time_3$AUC[nyears],2),
		round(ROC.time_4$AUC[nyears],2)
			),sep=" ")
legend("bottomright",lebel_2,col=mycolor[c(6,2,4,8,3)],lty=1,lwd=2)
dev.off()

#################Compared with LASSO########################
library(glmnet)
set.seed(105)
	clin<-a[,c(1,2,3)];dim(clin)
	rownames(clin)=clin[,1]
	clin=clin[co_sample,1:3];dim(clin)
	data.out=Surv(clin$DFI.time,clin$DFI)
	 {data=RNMF_data3[,];k=1} #0.9
	fit=cv.glmnet(data,data.out,family="cox",trace.it=1,nfolds=3)#,foldid=foldid
	plot(fit)
	fit.coef.lambda.lse=coef(fit,s=fit$lambda.min)
	#if(I==2) fit.coef.lambda.lse=coef(fit,s=fit$lambda.1se)
	fit.1se.out=fit.coef.lambda.lse[which(fit.coef.lambda.lse!=0)]
	fit.1se.out=round(fit.1se.out,4)
	fit.1se.out2=matrix(fit.1se.out,length(fit.1se.out),1)
	rownames(fit.1se.out2)=rownames(fit.coef.lambda.lse)[which(fit.coef.lambda.lse!=0)]
	colnames(fit.1se.out2)=c("coef")
	lasso_id=rownames(fit.1se.out2)
	lasso_id

clin<-a[,c(1,2,3)];dim(clin)
rownames(clin)=clin[,1]
clin=clin[co_sample,];dim(clin)
clin=cbind(clin,RNMF_data3[co_sample,lasso_id])
clin<-na.omit(clin);dim(clin)
colnames(clin)[1:3]<-c("sampleID","status","time")
clin=clin[,-1]
ZZ_data=ZZ_data[,colnames(RNMF_data[[3]])];dim(ZZ_data)
zz_data=as.matrix(ZZ_data)
clin_new<-b[,c(1,21,22)];dim(clin_new)
rownames(clin_new)=clin_new[,1]
clin_new=clin_new[rownames(ZZ_data),]
dim(clin_new)
clin_new=cbind(clin_new)
clin_new<-na.omit(clin_new);dim(clin_new)
colnames(clin_new)[1:3]<-c("sampleID","status","time")
clin_new=clin_new[,-1]
clin_new$time=clin_new$time*30
dd=datadist(clin)
options(datadist="dd")
pre_tcga=as.numeric(as.matrix(RNMF_data3[,])%*% coef(fit,s=fit$lambda.min))
concordance(cph(Surv(clin$time, clin$status) ~ pre_tcga), newdata=clin)
pre_new=as.numeric(as.matrix(zz_data[,])%*% coef(fit,s=fit$lambda.min))
concordance(cph(Surv(clin_new$time, clin_new$status) ~ pre_new), newdata=clin_new)



