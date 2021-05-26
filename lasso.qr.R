install.packages("LowRankQP")

library(LowRankQP)


#### Note: This function has been copied from https://www.erikdrysdale.com/l1_solution/

lasso.qr <- function(X,y,t) {

  y <- as.vector(y)	
  p <- dim(X)[2]
  #Get negative and psoitive side of X
  Vn <- X %*% cbind (cbind (diag(p), -diag(p)), 0)

  #Compute YX
  yX <- apply(sweep (X, MARGIN=1, -y, '*'), MARGIN=2, FUN=sum)
  Zn <- c (2*yX, -2*yX, 0)

  #Get upper bound for coefficients
  mod <- lm.fit(X, y)
  bOls <- mod$coefficients %>% replace_na(max(mod$coefficients, na.rm=TRUE))

  #Create matrices and vectors
  u <- c(abs(bOls), abs(bOls), sum(abs(bOls)))
  A <- matrix (c(rep (1, 2*p), 1), nrow=1)
  b <- c(min(t, sum(abs(bOls),na.rm=TRUE)))

  # Solve
  soln <- LowRankQP::LowRankQP(Vmat=sqrt(2)*t(Vn), dvec=Zn, Amat=A, bvec=b, uvec=u, method="LU")
  return(round(soln$alpha[1:p] - soln$alpha[(p+1):(2*p)], digits=5))
}