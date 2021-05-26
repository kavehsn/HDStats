install.packages("lars")
install.packages("glmnet")
install.packages("ADMM")
install.packages("tidyverse")

library(LowRankQP)
library(tidyverse)
library(latex2exp)

## Replication of Fig 2 of Wainwright's Sharp Thresholds for HD and Noisy Sparsity Recovery 

# 'n_min' is the minimum and initial number of observations for the simulation.
# 'n_max' is the maximum and final number of observations for the simulation.
# 'p_rng' specifies the set of dimension sizes to consider. This argument can be scalar, or a vector of dimensions.
# 'sss' is the sparsity percentage.
# 'step' is the increment size for sequencing from 'n_min' to 'n_max'. e.g., 'n_min = 5', 'n_max = 20', 'step = 10', -> n={5,15}.	
# 'iter' is the number of iterations for the simulations. Set to '1000' unless specified otherwise. 
# 'errStd' is the "noise level", corresponding to the standard deviation of the error terms. Set to '1' unless specifised otherwise.
# 'Intcpt' is the intercept for the DGP. Set to 'FALSE' unless specified otherwise.




BPSuccess.QP<-function(n_min,n_max,p_rng,sss,step,iter=1000,errStd=1,intcpt=FALSE,exact=FALSE,GLS=FALSE,...){


	m <-length(p_rng)

	B <- iter

	# 'matRow' is the maximum number of samples considered.

	matRow <- length(seq(n_min,n_max,by=step))

	# 'n_sizes' lists the samples from 'n_min' to 'n_max' sequenced by increments of size 'step'.

	n_sizes <- seq(n_min,n_max,by=step)


	# counter is a 'matRow x m' matrix of zeros. Each element of the matrix corresponds to the number of sign/support recovery successes in 'iter' simulations.

	counter <- matrix(data=0,nrow=matRow,ncol=m)


	# Loop through different dimensions.


	for(t in 1:m){


		# Loop through different sample sizes for each dimension 't'.


		for(j in 1:matRow){


			# Frequent messages that show the progress of the code, while it is running.


			if(j%%5==0){

				cat("Obs", j, "column", t,"\n")

			}


			# For each dimension and observations size simulate 'iter' times 


			set.seed(100)


			for(l in 1:B){

				d <- p_rng[t]
				n_sim <- n_sizes[j]

				if(intcpt==FALSE){

					# For each 'iter' generate the sparse model.

					DGP_list <- SparseDGP(n=n_sim,p=d,s=sss,epStd=errStd,...)

				}else{

					# For each 'iter' generate the sparse model.

					DGP_list <- SparseDGP(n=n_sim,p=d,s=sss,epStd=errStd,noC=FALSE,...)

				}


				# Extract the variables

				Y <- DGP_list$Y_Gen
				X <- DGP_list$X_Gen
				Beta <- DGP_list$sparseCoefs
				Sigma <- DGP_list$Omega
				


				if(GLS==FALSE){



					constraint <- sum(abs(Beta))

				#Estimate the Lasso model using 'glmnet'. For Lasso argument 'lambda' is set to '1'.
					Betahat <- lasso.qr(X=X,y=Y,constraint)


				}else{


					Rotation(Y=Y,X=X,Sigma=Sigma)
					Betahat <- lasso.qr(X=xStar,y=yStar,constraint)

				}





				if(intcpt==TRUE){

					# If intercept is present in the DGP, remove the first value of Beta. 

					Beta_noC <- Beta[-1]


				}else{


					Beta_noC <- Beta

				}

				# Remove the first value of 'Betahat' corresponding to the intercept. 

				Betahat_noC <- Betahat[-1]

				# Compare the signs of 'Betahat' and 'Beta'.

				if(exact==TRUE){

					if(all(sign(Betahat_noC)==sign(Beta_noC))){

					# If the the equality holds, add one to the counter of at row 'j' and column 't'.

						counter[j,t] <- counter[j,t]+1

					}

				}else{

					if(isSubset(Betahat_noC,Beta_noC)==TRUE){

					# If the the subset condition holds, add one to the counter of at row 'j' and column 't'.

						counter[j,t] <- counter[j,t]+1

					}				

				}

				

			}
		}


	}

		# The success probablity is calculated as 'counter / B'. 

	successRate <- counter/B
	counterTbl <<- counter
	probSuccess.QP <<- successRate 
}



isSubset <- function(a,b){
  # looks at support (subset) recovery
  # is supp(a) subset supp(b)
	if(sum(abs(a))>0){

		s1 <- which(abs(a)>0)
		s2 <- which(abs(b)>0)
		return(all(s1 %in% s2))

	}else{

		return(FALSE)
	}
}


Rotation <- function(Y,X,Sigma){

	#Cholesky factorization of the covariance matrix
	U <- chol(Sigma)
	L <- t(U)

	#Transform the variables
	yStar <- solve(L)%*%Y
	xStar <- solve(L)%*%X


	#Create a list of the variables
	data_list <- list(yStar, xStar)
	names(data_list) <- c("yStar", "xStar")
	rotatedVars <<- data_list

}


## Example

BPSuccess(n_min=350,n_max=350,p_rng=c(128,256,512),errStd=0.5,iter=20,sss=0.1,step=10,exact=TRUE) 
BPSuccess(n_min=5,n_max=100,p_rng=c(32,64,128),sss=0.1,step=5,iter=1000,GLS=FALSE,exact=TRUE, rho=0.95)  


x_axis<- seq(5,100,by=5)



png(filename="/home/kaveh/Desktop/NoGLSExactRho0.95.png")

mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 2, 0, 0)) 

#plot(x_axis, probSuccess[,1], type = "b", frame = FALSE, pch = 19, 
	#col = "red", xlab = "Number of observations", ylab = TeX(r'($P\[supp(\hat{\beta})\subseteq supp(\beta^*)\]$,  $\rho =0.95$)'))


plot(x_axis, probSuccess[,1], type = "b", frame = FALSE, pch = 19, 
	col = "red", xlab = "Number of observations", ylab = TeX(r'($P\[S_{\pm}(\hat{\beta})=S_{\pm}(\hat{\beta^*})\]$,  $\rho =0.5$)'))



lines(x_axis, probSuccess[,2], pch = 18, col = "blue", type = "b", lty = 2)

lines(x_axis, probSuccess[,3], pch = 18, col = "green", type = "b", lty = 2)


legend("topleft", legend=c(TeX("$d=32$"), TeX("$d=64$"),TeX("$d=128$")),
	col=c("red", "blue","green"), lty = 1:2, cex=0.8)


dev.off()