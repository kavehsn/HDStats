install.packages("lars")
install.packages("glmnet")
install.packages("ADMM")
install.packages("tidyverse")


library(lars)
library(glmnet)
library(ADMM)
library(tidyverse)

## Replication of Fig 2 of Wainwright's Sharp Thresholds for HD and Noisy Sparsity Recovery 

# 'n_min' is the minimum and initial number of observations for the simulation.
# 'n_max' is the maximum and final number of observations for the simulation.
# 'p_rng' specifies the set of dimension sizes to consider. This argument can be scalar, or a vector of dimensions.
# 'sss' is the sparsity percentage.
# 'step' is the increment size for sequencing from 'n_min' to 'n_max'. e.g., 'n_min = 5', 'n_max = 20', 'step = 10', -> n={5,15}.	
# 'iter' is the number of iterations for the simulations. Set to '1000' unless specified otherwise. 
# 'errStd' is the "noise level", corresponding to the standard deviation of the error terms. Set to '1' unless specifised otherwise.
# 'Intcpt' is the intercept for the DGP. Set to 'FALSE' unless specified otherwise.




BPSuccess<-function(n_min,n_max,p_rng,sss,step,iter=1000,errStd=1,intcpt=FALSE,exact=FALSE,GLS=FALSE,lambda_star=FALSE,...){


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

				# Number of non-zero elements of the random coefficient vector.					

				k_p <- round(sss*d)


				#### FOR MY OWN INFORMATION
				#theta <- runif(n=1,min=0.1,max=2.3)
				#theta <- 1
				#n_lambda <- 2*theta*log(d-k_p)
				####
				
				# Following Wainwright's paper, calculate the scaled 'lambda' for each dimension 'd', sample size 'n' and 'errStd'. 


				if(lambda_star==TRUE){
					

					if(GLS==FALSE){


						lambda_martin <- sqrt(((2*(errStd^2)*log(k_p)*log(d-k_p))/n_sim))

				#Estimate the Lasso model using 'glmnet'. For Lasso argument 'lambda' is set to '1'.
						lasso_model <- glmnet(X, Y, alpha = 1, lambda = lambda_martin, standardize = FALSE)
						Betahat <- coef(lasso_model)


					}else{


						#***THE TRUE LAMBDA MUST CHANGE ONCE THERE IS SUFFICIENT THEORY***

						lambda_martin <- sqrt(((2*(errStd^2)*log(k_p)*log(d-k_p))/n_sim))

						Rotation(Y=Y,X=X,Sigma=Sigma)

						lasso_model <- glmnet(X=xStar, Y=yStar, alpha = 1, lambda = lambda_martin, standardize = FALSE)
						Betahat <- coef(lasso_model)

					}


				}else{

					if(GLS==FALSE){


						lasso_model <- lars(x=X,y=Y,type="lasso",use.Gram=FALSE) 
						Betahat <- coef(lasso_model)



					}else{


						objRot <- Rotation(Y=Y,X=X,Sigma=Sigma)

						xRot <- objRot$xStar
						yRot <- objRot$yStar

						lasso_model <- lars(x=xRot,y=yRot,type="lasso",use.Gram=FALSE) 
						Betahat <- coef(lasso_model)
						
					}

				}

				
				if(lambda_star==TRUE){

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

				}else{


					dimBetahat <- dim(Betahat)
					nRowBetahat <- dimBetahat[1]

					if(intcpt==TRUE){

					# If intercept is present in the DGP, remove the first value of Beta. 

						Beta_noC <- Beta[-1]


					}else{


						Beta_noC <- Beta

					}

					if(exact==TRUE){


						for(jj in 1:nRowBetahat){


							if(all(sign(Betahat[jj,])==sign(Beta_noC))){

								counter[j,t] <- counter[j,t]+1

								break

							}


						}


					}else{

						for(kk in 1:nRowBetahat){

							if(isSubset(Betahat[kk,],Beta_noC)==TRUE){

					# If the the subset condition holds, add one to the counter of at row 'j' and column 't'.

								counter[j,t] <- counter[j,t]+1

								break

							}

						}


					}


				}
			}


		}

		# The success probablity is calculated as 'counter / B'. 

		successRate <- counter/B
		counterTbl <<- counter
		probSuccess <<- successRate 
	}

}


isSubset <- function(a,b){
  # looks at support (subset) recovery
  # is supp(a) subset supp(b)
	s1 <- which(abs(a)>0)
	s2 <- which(abs(b)>0)
	return(all(s1 %in% s2))
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

BPSuccess(n_min=350,n_max=350,p_rng=c(128,256,512),errStd=0.5,iter=20,sss=0.1,step=10) 
BPSuccess(n_min=5,n_max=350,p_rng=c(128,256,512),sss=0.1,step=10,iter=100,GLS=TRUE,exact=TRUE, rho=0.5)  



x_axis<- seq(5,350,by=10)



plot(x_axis, probSuccess[,1], type = "b", frame = FALSE, pch = 19, 
	col = "red", xlab = "Number of observations", ylab = "Probablity of success")


lines(x_axis, probSuccess[,2], pch = 18, col = "blue", type = "b", lty = 2)

lines(x_axis, probSuccess[,3], pch = 18, col = "green", type = "b", lty = 2)


legend("topleft", legend=c("d=128", "d=256","d=512"),
	col=c("red", "blue","green"), lty = 1:2, cex=0.8)

png(filename="/home/kaveh/Desktop/WainwrightSim.png")