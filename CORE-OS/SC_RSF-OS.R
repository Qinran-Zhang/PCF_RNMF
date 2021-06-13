
###########################################################################
rm(list=ls())
library(data.table)
library("survival")
#library("survminer")
library("data.table")
co_sample=read.csv("co_sample.csv",head=T,stringsAsFactors = FALSE)[,2]
co_sample<-gsub(".", "-", co_sample, fixed = TRUE)
group=substr(co_sample,14,15)
co_sample=co_sample[which(group=="01")]
a<-read.csv("COADREAD_survival.csv",head=T,stringsAsFactors = FALSE)
clin=a[which(a[,1] %in% co_sample),c(1,3,4)]
clin=na.omit(clin)
rownames(clin)=clin[,1]
dim(clin)
co_sample=clin[,1]
co_sample=co_sample[order(co_sample)]
clin=clin[co_sample,]
threshold.p=0.05

library(data.table)
eps=2.2204e-16
##########################################################
 sample.index=read.csv("sample_train.csv",stringsAsFactors=F)[,2]
# sample.index=read.csv("sample_train2.csv",stringsAsFactors=F)[,2]

#sample_train1=read.csv("sample_train1.csv",stringsAsFactors=F)[,2]
#sample_train2=read.csv("sample_train2.csv",stringsAsFactors=F)[,2]
#sample.index=c(sample_train1,sample_train2)
# sample.index=co_sample 
clin=clin[sample.index,];dim(clin)
mySurv =Surv(clin$OS.time, clin$OS)


# exon
ind0=function(x){ return(sum(x==0)>=0.9*length(x) ) }   ###########

file2<-fread("HiSeqV2_exon",header=TRUE, 
	sep = "\t", stringsAsFactors = F, na.strings = "", data.table = T)
# sep="\t"   sep=","
gene=data.frame(file2);rm(file2);gc()
data2<-gene[,-1]
rownames(data2)<-gene[,1];# rm(gene);gc()
## sample
colnames(data2)<-gsub(".", "-", colnames(data2), fixed = TRUE)
dim(data2)
data2=data2[,sample.index]
## geneID
ind=apply(data2,1,ind0)
data2<-data2[!ind,]
dim(data2)
geneID=rownames(data2)
Var=apply(data2,1,var)
if(length(which(Var==0))!=0) data2=data2[-which(Var==0),]
library(randomForestSRC)
	rf.data=cbind(as.data.frame(t(data2)),time=clin$OS.time,status=clin$OS)
	pbc.obj <- rfsrc(Surv(time, status) ~ .,rf.data,importance=F,
		block.size=1,ntree=2000,seed=(-c(1:6000)) )  #

	print(pbc.obj)				 #  
	# head(pbc.obj$forest$seed)     

	plot(pbc.obj)
	#if(I==1) ID1=max.subtree(pbc.obj)$topvars
	#if(I==2) ID2=max.subtree(pbc.obj)$topvars
	#if(I==3) ID3=max.subtree(pbc.obj)$topvars
	var_select=var.select(object = pbc.obj,method="vh" )
	ID=var_select$topvars;var_freq=var_select$varselect
	#imp=pbc.obj$importance;imp=imp[order(imp,decreasing=T)]
	#ID=names(imp)[which(imp>0)]
length(ID)
data2=t(data2[ID,])
write.csv(data2,"RNMF_data2-OS.csv")


### Methylation 24w
NAto0=function(x){ x[which(x=="NA")]<-0 ; return(x) }
NAind=function(x){ return(sum(x=="NA")>=0.9*length(x) ) }  ####################
file3<-fread("HumanMethylation450",header=TRUE, 
	sep = "\t", stringsAsFactors = F, na.strings = "", data.table = T)
gene=data.frame(file3);rm(file3);gc()
data3<-gene[,-1]
name3=gene[,1]
rownames(data3)<-gene[,1];rm(gene);gc()
dim(data3)
ID450_850=read.csv("850K&450K_cgID.csv",header=T)[,2]
data3=data3[which(rownames(data3) %in% ID450_850),]
dim(data3)
## sample
colnames(data3)<-gsub(".", "-", colnames(data3), fixed = TRUE)
data3=data3[,sample.index]
colname3=colnames(data3)
## cgID
 NAind=function(x){ return(sum(x=="NA")>0 ) }
ind=apply(data3,1,NAind)
data3<-data3[!ind,]
dim(data3)
rowname3=rownames(data3)

data3=t(data3)
#data3<-apply(data3,2,NAto0)
#sum(is.na(data3))
data3<-apply(data3,2,as.numeric)
data3=t(data3)
rownames(data3)<-rowname3
colnames(data3)<-colname3
Var=apply(data3,1,var)
if(length(which(Var==0))!=0) data3=data3[-which(Var==0),]

	rf.data=cbind(as.data.frame(t(data3)),time=clin$OS.time,status=clin$OS)
	pbc.obj <- rfsrc(Surv(time, status) ~ .,rf.data,importance=F,
		block.size=1,ntree=2000,seed=(-c(1:6000)) )  #

	print(pbc.obj)				 #  
	# head(pbc.obj$forest$seed)     

	plot(pbc.obj)
	#if(I==1) ID1=max.subtree(pbc.obj)$topvars
	#if(I==2) ID2=max.subtree(pbc.obj)$topvars
	#if(I==3) ID3=max.subtree(pbc.obj)$topvars
	var_select=var.select(object = pbc.obj,method="vh" )
	ID=var_select$topvars;var_freq=var_select$varselect
	#imp=pbc.obj$importance;imp=imp[order(imp,decreasing=T)]
	#ID=names(imp)[which(imp>0)]
length(ID)
data3=t(data3[ID,])
write.csv(data3,"RNMF_data3-OS.csv")
 # quit()

###################################



################
library(data.table)
RNMF_data3=fread("RNMF_data3-OS.csv",header=TRUE, 
	sep = ",", stringsAsFactors = F, na.strings = "", data.table = T)
RNMF_data3=data.frame(RNMF_data3)
rownames(RNMF_data3)=RNMF_data3[,1]
cgID=colnames(RNMF_data3)
library(data.table)
file3<-fread("betas_dt.csv",header=TRUE, 
	sep = ",", stringsAsFactors = F, na.strings = "", data.table = T)
gene=data.frame(file3);rm(file3);gc()

data=gene[,which(colnames(gene) %in% cgID)]
rownames(data)<-gene[,1]

write.csv(data,"RNMF_data4-OS.csv")


###  RNA-seq
file1<-fread("HiSeqV2",header=TRUE, 
	sep = "\t", stringsAsFactors = F, na.strings = "", data.table = T)
# sep="\t"   sep=","
gene=data.frame(file1);rm(file1);gc()
data1<-gene[,-1]
rownames(data1)<-gene[,1];# rm(gene);gc()
colnames(data1)<-gsub(".", "-", colnames(data1), fixed = TRUE)
data1=data1[,sample.index]
geneID=rownames(data1)
dim(data1)
Var=apply(data1,1,var)
if(length(which(Var==0))!=0) data1=data1[-which(Var==0),]
ind0=function(x){ return(sum(x==0)>=0.9*length(x) ) } #
ind=apply(data1,1,ind0) #
data1=data1[!ind,] #
library(randomForestSRC)
	rf.data=cbind(as.data.frame(t(data1)),time=clin$OS.time,status=clin$OS)
	pbc.obj <- rfsrc(Surv(time, status) ~ .,rf.data,importance=F,
		block.size=1,ntree=2000,seed=(-c(1:6000)) )  #

	print(pbc.obj)				 #  
	# head(pbc.obj$forest$seed)     

	plot(pbc.obj)
	#if(I==1) ID1=max.subtree(pbc.obj)$topvars
	#if(I==2) ID2=max.subtree(pbc.obj)$topvars
	#if(I==3) ID3=max.subtree(pbc.obj)$topvars
	var_select=var.select(object = pbc.obj,method="vh" )
	ID=var_select$topvars;var_freq=var_select$varselect
	#imp=pbc.obj$importance;imp=imp[order(imp,decreasing=T)]
	#ID=names(imp)[which(imp>0)]
length(ID)
data1=t(data1[ID,])
write.csv(data1,"RNMF_data1-OS.csv")



##########################################################################
rm(list=ls())
library(data.table)
library("survival")
#library("survminer")
library("data.table")
co_sample=read.csv("co_sample.csv",head=T,stringsAsFactors = FALSE)[,2]
co_sample<-gsub(".", "-", co_sample, fixed = TRUE)
group=substr(co_sample,14,15)
co_sample=co_sample[which(group=="01")]
a<-read.csv("COADREAD_survival.csv",head=T,stringsAsFactors = FALSE)
clin=a[which(a[,1] %in% co_sample),c(1,3,4)]
clin=na.omit(clin)
rownames(clin)=clin[,1]
dim(clin)
co_sample=clin[,1]
co_sample=co_sample[order(co_sample)]
clin=clin[co_sample,]
threshold.p=0.05
eps=2.2204e-16

### RNA_seq    # 1.5,0.1
RNMF_data1=read.csv("RNMF_data1-OS.csv",head=T,stringsAsFactors = FALSE,check.names=F)
#colnames(RNMF_data1)<-gsub(".", "-", colnames(RNMF_data1), fixed = TRUE)
ID1=colnames(RNMF_data1)[-1]

dim(RNMF_data1)
### miRNA 2k   # 2, 0.1
RNMF_data2=fread("RNMF_data2-OS.csv",sep = ",",head=T,stringsAsFactors = FALSE,check.names=F)
RNMF_data2=data.frame(RNMF_data2,check.names=F)
#colnames(RNMF_data2)<-gsub(".", "-", colnames(RNMF_data2), fixed = TRUE)
ID2=colnames(RNMF_data2)[-1]

### Methylation # 1.5, 0.1
RNMF_data3=fread("RNMF_data3-OS.csv",sep = ",",head=T,stringsAsFactors = FALSE,check.names=F)
RNMF_data3=data.frame(RNMF_data3,check.names=F)
#colnames(RNMF_data3)<-gsub(".", "-", colnames(RNMF_data3), fixed = TRUE)
ID3=colnames(RNMF_data3)[-1]

# exon

file2<-fread("HiSeqV2_exon",header=TRUE, 
	sep = "\t", stringsAsFactors = F, na.strings = "", data.table = T)
# sep="\t"   sep=","
gene=data.frame(file2);rm(file2);gc()
data2<-gene[,-1]
rownames(data2)<-gene[,1];# rm(gene);gc()
## sample
colnames(data2)<-gsub(".", "-", colnames(data2), fixed = TRUE)
dim(data2)
length(ID2)
data2=t(data2[ID2,])
write.csv(data2,"RNMF_data2-OS_end.csv")


file1<-fread("HiSeqV2",header=TRUE, 
	sep = "\t", stringsAsFactors = F, na.strings = "", data.table = T)
# sep="\t"   sep=","
gene=data.frame(file1);rm(file1);gc()
data1<-gene[,-1]
rownames(data1)<-gene[,1];# rm(gene);gc()
colnames(data1)<-gsub(".", "-", colnames(data1), fixed = TRUE)
length(ID1)
data1=t(data1[ID1,])
write.csv(data1,"RNMF_data1-OS_end.csv")


### Methylation 24w

NAto0=function(x){ x[which(x=="NA")]<-0 ; return(x) }
NAind=function(x){ return(sum(x=="NA")>=0.9*length(x) ) }  ####################
file3<-fread("HumanMethylation450",header=TRUE, 
	sep = "\t", stringsAsFactors = F, na.strings = "", data.table = T)
gene=data.frame(file3);rm(file3);gc()
data3<-gene[,-1]
name3=gene[,1]
rownames(data3)<-gene[,1];rm(gene);gc()
dim(data3)
ID450_850=read.csv("850K&450K_cgID.csv",header=T)[,2]
data3=data3[which(rownames(data3) %in% ID450_850),]
dim(data3)
## sample
colnames(data3)<-gsub(".", "-", colnames(data3), fixed = TRUE)
data3=data3[,]
colname3=colnames(data3)
## cgID
data3<-data3[ID3,]
dim(data3)
rowname3=rownames(data3)

data3=t(data3)
data3<-apply(data3,2,as.numeric)
data3=t(data3)
rownames(data3)<-rowname3
colnames(data3)<-colname3
Var=apply(data3,1,var)
if(length(which(Var==0))!=0) data3=data3[-which(Var==0),]
length(ID3)
data3=t(data3[ID3,])
write.csv(data3,"RNMF_data3-OS_end.csv")
 # quit()

###################################

file3<-fread("betas_dt.csv",header=TRUE, 
	sep = ",", stringsAsFactors = F, na.strings = "", data.table = T)
gene=data.frame(file3);rm(file3);gc()

data=gene[,which(colnames(gene) %in% ID3)]
rownames(data)<-gene[,1]

write.csv(data,"RNMF_data4-OS_end.csv")


