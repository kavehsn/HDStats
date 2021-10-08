

install.packages("lars")
install.packages("tidyverse")


library(glmnet)
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




LassoEstError <- function(n_min, n_max, p_rng, sss, step, iter=1000, errStd=1, intcpt=FALSE, exact=FALSE, GLS=FALSE, lambda_star=FALSE,...){


	m <- length(p_rng)

	B <- iter

	# 'matRow' is the maximum number of samples considered.

	matRow <- length(seq(n_min, n_max, by=step))

	# 'n_sizes' lists the samples from 'n_min' to 'n_max' sequenced by increments of size 'step'.

	n_sizes <- seq(n_min, n_max, by=step)


	# counter is a 'matRow x m' matrix of zeros. Each element of the matrix corresponds to the number of sign/support recovery successes in 'iter' simulations.

	counter <- matrix(data=0, nrow=matRow, ncol=m*3)


	# Loop through different dimensions.


	FTT <- 1

	for(t in 1:m){


		# Loop through different sample sizes for each dimension 't'.


		for(j in 1:matRow){


			# Frequent messages that show the progress of the code, while it is running.


			if(j%%5==0){

				cat("Obs", j, "column", t,"\n")

			}


			# For each dimension and observations size simulate 'iter' times 


			l2Error <- matrix(data=0, nrow=B, ncol=1)


			for(l in 1:B){

				d <- p_rng[t]
				n_sim <- n_sizes[j]

				if(intcpt==FALSE){

					# For each 'iter' generate the sparse model.

					DGP_list <- SparseDGP(n=n_sim, p=d, s=sss, epStd=errStd, ...)

				}else{

					# For each 'iter' generate the sparse model.

					DGP_list <- SparseDGP(n=n_sim, p=d, s=sss, epStd=errStd, noC=FALSE,...)

				}


				# Extract the variables

				Y <- DGP_list$Y_Gen
				X <- DGP_list$X_Gen
				Beta <- DGP_list$sparseCoefs

				# Number of non-zero elements of the random coefficient vector.					

				k_p <- round(sss*d)


				#### FOR MY OWN INFORMATION
				#theta <- runif(n=1,min=0.1,max=2.3)
				#theta <- 1
				#n_lambda <- 2*theta*log(d-k_p)
				####
				
				# Following Wainwright's paper, calculate the scaled 'lambda' for each dimension 'd', sample size 'n' and 'errStd'. 




				if(GLS==FALSE){


					cvfit <- glmnet::cv.glmnet(x=X, y=Y) 
					Betahat <- coef(cvfit, s = "lambda.1se")



				}else{


					cvfit1 <- glmnet::cv.glmnet(x=X, y=Y) 

					Betahat1 <- coef(cvfit1, s = "lambda.1se")

					X_intercept <- matrix(data=c(rep(1, times = n), X), nrow = n, ncol = d+1) 

					errHat <- as.numeric(Y - X_intercept%*%Betahat1)

					RhoHat <- ar.ols(errHat, aic = FALSE, order.max = 1, na.action = na.fail, demean = TRUE)

					rhoHat <- as.numeric(RhoHat[2])

					objRot <- Rotation(Y=Y, X=X, parma=rhoHat)

					xRot <- objRot$xStar
					yRot <- objRot$yStar

					cvfit2 <- glmnet::cv.glmnet(x=xRot, y=yRot) 
					Betahat2 <- coef(cvfit2)

				}

				#Calculate the L_2 error
				if(GLS==FALSE){

					l2Error[l] <- norm(Beta-Betahat[2:length(Betahat)], type = "2")

				}else{


					l2Error[l] <- norm(Beta-Betahat2[2:length(Betahat2)], type = "2")

				}



			}

			

			counter[j, FTT] <- min(l2Error)
			counter[j, FTT+1] <- mean(l2Error)
			counter[j, FTT+2] <- max(l2Error)


		}

		FTT <- FTT + 3
	}


	l2ErrorTbl <<- counter


}


LassoEstError(n_min=50, n_max=100, p_rng=c(50, 100), sss=0.1, step=50, iter=50, errStd=1, intcpt=FALSE, exact=FALSE, GLS=FALSE, lambda_star=FALSE)



Rotation <- function(Y, X, parma){

	sqnc <- parma^seq(0, n, by=1)
	Sigma <- toeplitz(sqnc[1:n])
	Omega <- (1/(1-parma^2))*Sigma

	#Cholesky factorization of the covariance matrix
	U <- chol(Omega)
	L <- t(U)

	#Transform the variables
	yStar <- solve(L)%*%Y
	xStar <- solve(L)%*%X


	#Create a list of the variables
	data_list <- list(yStar, xStar)
	names(data_list) <- c("yStar", "xStar")
	rotatedVars <<- data_list

}


