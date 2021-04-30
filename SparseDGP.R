
install.packages("MASS")


library(MASS)

#The Data generating process


## n is the number of observations 
## p is the dimension of the random design matrix X
## theta is the equi-correlation parameter between the columns of the random design matrix X.
## s is the sparsity percentage. 
## mu_x is multivariate mean of the random design matrix. If TRUE, it is set to 0.
## rho is the AR(1) parameter. Unless specified it is set to zero.
## equicor is for the case of equicorrelated errors

SparseDGP<-function(n,p,s,theta=0,mu_x=rep(0,times=p),rho=0,equicor=FALSE,EV=FALSE,BIV=FALSE,GARCH=FALSE){


	#Generate random sparse coefficeints, where parameter s is the sparsity percentage. Finally, add a random intercept.
	ss<-round((1-s)*p)
	coefs<-rnorm(p+1,mean=0,sd=1)
	toReplace<-sample(p,ss)
	sCoefs<-replace(coefs,toReplace,0)
	


	#Generate the random design matrix. If mu_x==TRUE, then the means are zero.
	if(mu_x!=rep(0,times=p)){
		mu_x<-rep(0,times=p)	
	}
	
	#An equi-correlated random design matrix. If theta=0, the columns of matrix are drawn with i.i.d N(0,1) entries.
	sigma_x<-(1-theta)*diag(p)+theta*matrix(data=1,ncol=p,nrow=p)
	x<-mvrnorm(n,mu=mu_x,Sigma=sigma_x)

	#Add a column of ones as the first column of the matrix X.
	X<-cbind(rep(1,times=n),x)


	#If equicor==FALSE, the errors follow an AR(1) process with parameter rho. When rho=0, the errors are i.i.d N(0,1). 
	#Otherwise equicorrelated errors are generated.
	
	if(equicor==TRUE){

		Omega<-(1-rho)*diag(n)+rho*matrix(data=1,ncol=n,nrow=n)
		e<-mvrnorm(n=1,mu=rep(0,times=n),Sigma=Omega)


	}else if(EV==TRUE){


		Sqnc<-exp(0.5*(1:n))
		sqnc<-Sqnc^2
		Omega<-diag(sqnc)
		e<-mvrnorm(n=1,mu=rep(0,times=n),Sigma=Omega)


	}else if(BIV==TRUE){


		Sigma<-diag(n)
		t<-round(n/2)
		Sigma[t,t]<-1000
		Omega<-Sigma
		e<-mvrnorm(n=1,mu=rep(0,times=n),Sigma=Omega)



	}else if(GARCH==TRUE){

		
		e<-rep(0,times=n)
		e[1]<-rnorm(1)
		Sigma<-rep(0,times=n)
		Sigma[1]<-1

		for(t in 2:n){

			Sigma[t]<-0.00037+0.0888*e[t-1]^2+0.9024*Sigma[t-1]
			e[t]<-rnorm(1,mean=0,sd=sqrt(Sigma[t]))

		}

		Omega<-diag(Sigma)

	}else{


		sqnc<-rho^seq(0,n,by=1)
		Sigma<-toeplitz(sqnc[1:n])
		Omega<-(1/(1-rho^2))*Sigma
		e<-mvrnorm(n=1,mu=rep(0,times=n),Sigma=Omega)

	}



	#Finally the sparse linear model is generated as follow

	y<-X%*%sCoefs+e

	X_NI<-X[,2:p]

	data_list<-list(y,X_NI,sCoefs,e,Omega)
	names(data_list)<-c("Y_Gen","X_Gen","sparseCoefs","eps","Omega")
	dataSparse<-return(data_list)

}