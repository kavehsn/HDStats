install.packages("lars")
install.packages("glmnet")


library(lars)
library(glmnet)

#Replication of Figure 7.3 of Martin Wainwright's High-Dimensional Statistics



BPSuccess<-function(n_min,n_max,p_rng,sss,step,iter=1000,...){


	m<-length(p_rng)

	B<-iter

	DGP_list <- vector(mode = "list", length = m)

	matRow<-length(seq(n_min,n_max,by=step))

	counter<-matrix(data=0,nrow=matRow,ncol=m)


	for(t in 1:m){


		for(j in seq(n_min,n_max,by=step)){


			for(l in 1:B){

				d<-p_rng[t]
				DGP_list[[t]]<-SparseDGP(n=j,p=d,s=sss,...)
				dataLasso<-DGP_list[[t]][c('Y_Gen','X_Gen','sparseCoefs','Omega')]


				Y<-dataLasso$Y_Gen
				X<-dataLasso$X_Gen
				Beta<-dataLasso$sparseCoefs
				Sigma<-dataLasso$Omega

				k_p<-round(s*d)
				lambda_martin<-sqrt((2*log(k_p)*log(d-k_p)/n))

				lasso_model <- glmnet(X, Y, alpha = 1, lambda = lambda_martin, standardize = TRUE)
				Betahat<-coef(lasso_model)


				if(sign(Beta)==sign(Betahat)){


					counter[j,t]<-counter[j,t]+1

				}
			}
		}


	}



}

BPSuccess(n_min=10,n_max=30,p_rng=c(40,50,60),s=0.1,step=10) 