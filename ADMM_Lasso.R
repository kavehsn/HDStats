library(Matrix)
library(matrixcalc)



ADMM_Lasso <- function(A,b,lambda=1,rho=1,alpha=1,MAX_ITER=1000,ABSTOL=1e-4, RELTOL=1e-2){



	dim_A <- dim(A)

	r_dim <- dim_A[1]
	c_dim <- dim_A[2]

	z <- rep(0,times=c_dim)
	u <- rep(0,times=c_dim)
	x <- rep(0,times=c_dim)

	Atb <- t(A)%*%b


	objval <- rep(0,times=MAX_ITER)
	r_norm <- rep(0,times=MAX_ITER)
	s_norm <- rep(0,times=MAX_ITER)
	eps_pri <- rep(0,times=MAX_ITER)
	eps_dual <- rep(0,times=MAX_ITER)


	factor(A,rho)

	for(k in 1:MAX_ITER){

		q <- Atb + rho*(z-u)

	#x-update
		if(r_dim>=c_dim){

			x <- solve(U)%*%(solve(L)%*%q)

		}else{

			x <- (q/rho) - ((t(A)%*%(solve(U)%*%((solve(L)%*%(A%*%q))))/(rho^2)))

		}


	#z-update

		zold <- z
		x_hat <- alpha*x+(1-alpha)*zold
		z <- shrinkage(x_hat+u, lambda/rho)
		z <- matrix(data=z,nrow=length(z),ncol=1)


	#u-update

		u <- u+(x_hat-z)


	#Diagnostics, reporting, termination checks
		objval[k] <- objective(A,b,lambda,x,z)

		r_norm[k] <- norm(x-z, type="F")
		s_norm[k] <- norm(-rho*(z-zold), type="F")

		eps_pri[k] <- sqrt(c_dim)*ABSTOL + RELTOL*max(norm(x,type="F"),norm(-z,type="F"))
		eps_dual[k] <- sqrt(c_dim)*ABSTOL + RELTOL*norm(rho*u,type="F")


		if((r_norm[k] < eps_pri[k]) && (s_norm[k] < eps_dual[k])){

			objval <- head(objval, k)
			r_norm <- head(r_norm, k)
			s_norm <- head(s_norm, k)
			eps_pri <- head(eps_pri, k)
			eps_dual <- head(eps_dual, k)

			break

		}
	}


	History <<- list(objval, r_norm, s_norm, eps_pri, eps_dual)
	names(History) <<- c("objval","r_norm", "s_norm", "eps_pri", "eps_dual")

	z_opt <<- z
	x_opt <<- x

}


#The shrinkage parameter

shrinkage <- function(x,kappa){

	z <- pmax(0,x-kappa)-pmax(0,-x-kappa)
}



	#Cholesky decomposition of (A'A+rho*I)

factor <- function(A,rho){

	dim_A <- dim(A)

	r_dim <- dim_A[1]
	c_dim <- dim_A[2]

	if(r_dim>=c_dim){

		U <- chol(t(A)%*%A + rho*diag(c_dim))

	}else{

		U <- chol(diag(r_dim)+(1/rho)*(A%*%t(A)))

	}


	L <<- Matrix(t(U),sparse=TRUE)
	U <<- Matrix(U,sparse=TRUE)

}

	#Opbjective function for the ADMM_LASSO optimization problem

objective <- function(A,b,lambda,x,z){

	dim_A <- dim(A)

	r_dim <- dim_A[1]
	c_dim <- dim_A[2]

	p <- (1/(2*r_dim))*sum(((A%*%x)-b)^2)+(lambda*norm(z,type="1")) 

}

