eps=2.2204e-256
doublecol=function(X){
	m=dim(X)[1]
	n=dim(X)[2]
	Y=matrix(0,nrow=m,ncol=2*n)
	Y[,seq(1,2*n,2)]=X*(X>0)
	Y[,seq(2,2*n,2)]=-X*(X<0)
	return(Y)
}
doublecol_inv=function(Y,para=NULL){
	#    UNTITLED Summary of this function goes here
	#    Detailed explanation goes here
	#    [m,n]=size(Y);
	if(is.null(para)) para=0
	if(para==0){ X=Y[,seq(1,ncol(Y),2)]-Y[,seq(2,ncol(Y),2)]}
		else { Y1=Y[,seq(1,ncol(Y),2)]
			 Y2=Y[,seq(2,ncol(Y),2)]
			 flag=abs(Y1)>abs(Y2)
			 Y1[!flag]=0
			 Y2[flag]=0
			 X=Y1-Y2	
			}	
	return(X)
}	
nmf_euclidean_dist=function(data,W,H,lambda){
	I=length(data)
	err=0
	for(i in 1:I){
		err = err+norm(data[[i]]-W%*%H[[i]],"F")^2  # 计算F-范数
		}
	#err = err+lambda*norm(H[[1]],"F")^2+lambda*norm(H[[2]],"F")^2+lambda*norm(H[[3]],"F")^2
	err = err+lambda*norm(W,"F")^2
	return(err)
}

RNMF=function(data,K,lambda=NULL,maxiter=NULL,speak=NULL){
	if(!is.list(data))	data=list(data)
	if(is.null(lambda))     lambda=0
	if(is.null(maxiter))    maxiter=3000
	print_iter=maxiter/5	# print_iter is what?
	if(is.null(speak))      speak=TRUE

	I=length(data)
	Dist=0;nall=0
	for(Ii in 1:I){	Dist=Dist+norm(data[[Ii]],"F")^2;nall=nall+ncol(data[[Ii]]) }

	m=dim(data[[1]])[1]
	n=dim(data[[1]])[2]
	H=list()
	for(i in 1:I){  n[i] = dim(data[[i]])[2] 
	H= c(H,list( matrix(max(data[[i]])*runif((K*n[i])),nrow=K,ncol=n[i]) ))# 生成1*3元包数组/列表
	}	## 乘最大值 max(X)
	W=matrix(max(data[[i]])*runif(m*K),nrow=m,ncol=K)   # 生成m行K列随机数矩阵
	# str(H)
	
	R=m
	Omega=matrix(runif(R*m),nrow=R,ncol=m)
	Y=list();Q=list();XR=list();HR=list()
	for(i in 1:I){
		Y[[i]] = Omega %*% data[[i]]
		qrresult = qr( t(Y[[i]]))  #进行qr分解 #,LAPACK = T
		Q[[i]] = qr.Q(qrresult)
		XR[[i]]= doublecol(data[[i]] %*% Q[[i]])
		HR[[i]]= matrix( runif(K*dim(XR[[i]])[2]),nrow=K )
	}
	t1=Sys.time()
	for(iter in 1:maxiter){
		# 更新各个H
		for(i in 1:I){
			# matlab中 ./   .*    *    /
			#    R  中  /    *   %*%  %/%
			HR[[i]] = HR[[i]]*(t(W)%*%XR[[i]])/ ((t(W)%*%W)%*%HR[[i]]+eps)
			}
		# 更新W
		HXt=matrix(0,nrow=K,ncol=m)
		HHt=matrix(0,nrow=K,ncol=K)
		for(i in 1:I){
			HXt = HXt+HR[[i]] %*% t(XR[[i]])
			HHt = HHt+HR[[i]] %*% t(HR[[i]])
			}
		W=W*t(HXt)/(W%*%HHt+eps+lambda*W)
		if(iter%%print_iter==0){
			for(i in 1:I){
				H[[i]] = doublecol_inv(HR[[i]])%*% t(Q[[i]])
				H[[i]][which(H[[i]]<0)] = 0
				#HR[[i]]= doublecol(H[[i]]%*%Q[[i]])
				}
			if(speak){
				eucl_dist = nmf_euclidean_dist(data,W,H,lambda)
				print(list(Iter=iter,eucl_dist=eucl_dist))				
				}
			}
#
		if(iter%%100==0){
			for(i in 1:I){
				H[[i]] = doublecol_inv(HR[[i]])%*% t(Q[[i]])
				H[[i]][which(H[[i]]<0)] = 0
				#HR[[i]]= doublecol(H[[i]]%*%Q[[i]])
				}
			if(iter!=100){
				eucl_dist = nmf_euclidean_dist(data,W,H,lambda)
				d_eucl=abs(eucl_dist-old_eucl)
				old_eucl=eucl_dist
				print(old_eucl)
				#if(d_eucl<10^(-5)*Dist)	 break;	
				if(d_eucl< 10^(-3)*nall*nrow(data[[1]]))	 {print(d_eucl);break;}	
				#if(d_eucl< 10^(-5))	 {print(d_eucl);break;}
				}
			else{old_eucl = nmf_euclidean_dist(data,W,H,lambda)}
			}
		}
	time = Sys.time()-t1
	Wmax = apply(W,2,max)
	W = W %*% diag(1/(Wmax))
	for(i in 1:I){
		H[[i]] = diag(Wmax) %*% H[[i]]
		}
#	WF=c();WFHF=c();HF=list();length(HF)=I
#	for(i in 1:ncol(W)){
#		WF[i]=norm(t(W[,i]),"F")
#		WFHF[i]=WF[i]
#		for(j in 1:I){
#			HF[[j]][i]=norm(t(H[[j]][i,]),"F")
#			WFHF[i]=WFHF[i]*HF[[j]][i]
#		}
#	}
	eucl_dist=nmf_euclidean_dist(data,W,H,lambda)  # error
#	index=order(WFHF,decreasing = TRUE)
#	for(indexi in 1:length(H)){ H[[indexi]]=H[[indexi]][index,]}
	return(list(time=time,iter=iter,
			W=W,H=H,loss=eucl_dist))

}


findmodule=function(W,H,tred){
	n = dim(W)[1]
	K = dim(W)[2]
	I=length(H)
	module=list();length(module)=I+1
	meanW=colMeans(W)
	sdW=apply(W,2,sd)
	meanH=list();sdH=c()
	for(i in 1:I){
		meanH[[i]]=rowMeans(H[[i]])
		sdH[[i]]=apply(H[[i]],1,sd)
		}
	for(i in 1:K){
		module[[1]]=c(module[[1]], list(which( W[,i]>(meanW[i]+tred[1]*sdW[i]) )))
		for(j in 1:I){
			module[[j+1]]=c(module[[j+1]], 
				# list( which(H[[j]][i,] > meanH[[j]][i]+tred[2]*sdH[[j]][i] )) )
				list( which(H[[j]][i,] > meanH[[j]][i]+tred[2]*sdH[[j]][i] )) )
			}
		}
	return(module)
}


 