
install.packages("MASS")


library(MASS)

##The Data generating process

# n is the number of observations 
# p is the dimension of the random design matrix X
# s is the sparsity percentage. 
# theta is the equi-correlation parameter for X ensemble. It is set to 0 unless specified otherwise.
# mu_x is multivariate mean of the X ensemble. If TRUE, it is set to a 0 vector.
# noC specifies the presence (or the lack of) an intercept in the model. When TRUE the DGP has no intercept.
# EquiProb follow Wainwright by setting beta={-0.5,+0.5} with equal probability. When TRUE the sparse beta vector assumes only the values of -0.5, +0.5 and 0. 
# rho is the AR(1) parameter. Unless specified it is set to zero.
# epStd is the standard deviation of the error term in the case of independent errors. Set to 1 unless specified otherwise.
# equicor is the case of equicorrelated errors. Set to FALSE unless specified otherwise.	
# EV allows errors to exhibit exponential variance. Set to FALSE unless specified otherwise.
# BIV allows error to exhibit break in their variances across observations. Set to FALSE unless specified otherwise.
# GARCH allows the errors to follow a stationary GARCH(1,1) process. Set to FALSE unless specified otherwise.


SparseDGP <- function(n,p,s,theta=0,mu_x=rep(0,times=p),noC=TRUE,equiProb=TRUE,rho=0,epStd=1,equicor=FALSE,EV=FALSE,BIV=FALSE,GARCH=FALSE){


	# Generate a random vector of coefficeints.


	if(equiProb==TRUE && noC==FALSE){

		w <- rbinom(n=p+1,size=1,prob=0.5)
		coefs <- w*0.5-(1-w)*(0.5)

	}else if(equiProb==FALSE && noC==FALSE){

		coefs <- rnorm(p+1,mean=0,sd=1)

	}else if(equiProb==TRUE && noC==TRUE){

		w <- rbinom(n=p,size=1,prob=0.5)
		coefs <- w*0.5-(1-w)*(05)

	}else if(equiProb==FALSE && noC==TRUE){

		coefs <- rnorm(p,mean=0,sd=1)

	}
	

	# Randomly replace (1-s)% of the coefficients with zero

	ss <- round((1-s)*p)

	toReplace <- sample(p,ss)
	sCoefs <- replace(coefs,toReplace,0)
	


	# Generate the random design matrix. If mu_x is a zero vetor unless specified otherwise.
	# A non-zero theta allows for an equi-correlated random design matrix. 
	# If theta=0, the columns of matrix are drawn with i.i.d N(0,1) entries.
	
	sigma_x <- (1-theta)*diag(p)+theta*matrix(data=1,ncol=p,nrow=p)
	x <- mvrnorm(n,mu=mu_x,Sigma=sigma_x)

	# If noC=FALSE, add a column of ones as the first column of matrix X.

	if(noC==FALSE){

		X <- cbind(rep(1,times=n),x)

	}else{

		X <- x

	}

	# If equicor==TRUE, equicorrelated errors are generated. 
	# If rho is non-zero, an AR(1) process with parameter rho is generated. 
	# When rho=0, the errors are i.i.d N(0,1). 
	# We may also generate different forms of heteroskedasticity.
	
	if(equicor==TRUE){

		Omega <- (1-rho)*diag(n)+rho*matrix(data=1,ncol=n,nrow=n)
		e <- mvrnorm(n=1,mu=rep(0,times=n),Sigma=Omega)


	}else if(EV==TRUE){


		Sqnc <- exp(0.5*(1:n))
		sqnc <- Sqnc^2
		Omega <- diag(sqnc)
		e <- mvrnorm(n=1,mu=rep(0,times=n),Sigma=Omega)


	}else if(BIV==TRUE){


		Sigma <- diag(n)
		t <- round(n/2)
		Sigma[t,t] <- 1000
		Omega <- Sigma
		e <- mvrnorm(n=1,mu=rep(0,times=n),Sigma=Omega)



	}else if(GARCH==TRUE){

		
		e <- rep(0,times=n)
		e[1] <- rnorm(1)
		Sigma <- rep(0,times=n)
		Sigma[1] <- 1

		for(t in 2:n){

			Sigma[t] <- 0.00037+0.0888*e[t-1]^2+0.9024*Sigma[t-1]
			e[t] <- rnorm(1,mean=0,sd=sqrt(Sigma[t]))

		}

		Omega <- diag(Sigma)

	}else if(rho==0 && epStd==1){


		sqnc <- rho^seq(0,n,by=1)
		Sigma <- toeplitz(sqnc[1:n])
		Omega <- (1/(1-rho^2))*Sigma
		e <- mvrnorm(n=1,mu=rep(0,times=n),Sigma=Omega)

	}else if(rho==0 && epStd!=1){

		Omega <- diag(epStd^2,n,n)
		e <- mvrnorm(n=1,mu=rep(0,times=n),Sigma=Omega)

	}

	# Finally a sparse linear model is generated as follows:

	y <- X%*%sCoefs+e

	# In presence of an intercept, - i.e. noC=FALSE, the first column of X constituing of ones is removed for estimation. 

	if(noC==FALSE){

		X_NI <- X[,-1]

	}else{

		X_NI <- X
	}

	# The outcomes are stored in dataSparse. For instance, to obtain the generated Y vector, write dataSparse$Y_Gen.

	data_list <-list(y,X_NI,sCoefs,e,Omega)
	names(data_list) <- c("Y_Gen","X_Gen","sparseCoefs","eps","Omega")
	dataSparse <<- data_list
}




SparseDGP(n=20,p=40,s=0.1,epStd=0.5)  