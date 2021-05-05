install.packages("lars")
install.packages("glmnet")


library(lars)
library(glmnet)

## Replication of Fig 2 of Wainwright's Sharp Thresholds for HD and Noisy Sparsity Recovery 

# 'n_min' is the minimum and initial number of observations for the simulation.
# 'n_max' is the maximum and final number of observations for the simulation.
# 'p_rng' specifies the set of dimension sizes to consider. This argument can be scalar, or a vector of dimensions.
# 'sss' is the sparsity percentage.
# 'step' is the increment size for sequencing from 'n_min' to 'n_max'. e.g., 'n_min = 5', 'n_max = 20', 'step = 10', -> n={5,15}.	
# 'iter' is the number of iterations for the simulations. Set to '1000' unless specified otherwise. 
# 'errStd' is the "noise level", corresponding to the standard deviation of the error terms. Set to '1' unless specifised otherwise.
# 'Intcpt' is the intercept for the DGP. Set to 'FALSE' unless specified otherwise.




BPSuccess<-function(n_min,n_max,p_rng,sss,step,iter=1000,errStd=1,intcpt=FALSE,...){


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

				lambda_martin <- sqrt(((2*(errStd^2)*log(k_p)*log(d-k_p))/n_sim))

				#Estimate the Lasso model using 'glmnet'. For Lasso argument 'lambda' is set to '1'.


				lasso_model <- glmnet(X, Y, alpha = 1, lambda = lambda_martin, standardize = TRUE)
				Betahat <- coef(lasso_model)
				
				if(intcpt==TRUE){

					# If intercept is present in the DGP, remove the first value of Beta. 

					Beta_noC <- Beta[-1]


				}else{


					Beta_noC <- Beta

				}

				# Remove the first value of 'Betahat' corresponding to the intercept. 

				Betahat_noC <- Betahat[-1] 

				# Compare the signs of 'Betahat' and 'Beta'.

				if(all(sign(Beta_noC)==sign(Betahat_noC))){


					# If the the equality holds, add one to the counter of at row 'j' and column 't'.

					counter[j,t] <- counter[j,t]+1

				}
			}


		}

		# The success probablity is calculated as 'counter / B'. 

		successRate <- counter/B
		counterTbl <<-counter
		probSuccess <<- successRate 
	}

}


## Example

BPSuccess(n_min=5,n_max=350,p_rng=c(128,256,512),errStd=0.5,iter=250,sss=0.1,step=10) 