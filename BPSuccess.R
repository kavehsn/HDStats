install.packages("lars")


library(lars)

#Replication of Figure 7.3 of Martin Wainwright's High-Dimensional Statistics



BPSuccess<-function(n_max,p_rng,iter=1000,...){

	n_rng<-1:n_max
	m<-length(p_rng)

	DGP_list <- vector(mode = "list", length = m)

	for(t in 1:m){

		d<-p_rng[t]
		DGP_list[[t]]<-SparseDGP(n=n_max,p=d,...)
		dataLasso<-DGP_list[[t]][c('Y_Gen','X_Gen','sparseCoefs')]


		Y<-dataLasso$Y_Gen
		X<-dataLasso$X_Gen
		Beta<-dataLasso$sparseCoefs


		Beta_hat<-lars(x=X,y=Y,type="lasso")

	}

}