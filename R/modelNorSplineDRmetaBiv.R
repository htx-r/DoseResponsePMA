#******* jags model of spline dose-response model with binomial likelihood for OR

modelNorSplineDRmetaBiv <- function(){
  b[1] <-0
  for (i in 1:ns) { ## for each study

    # multivariate normal likelihood of the response Y = logOR or logRR for the *reference* and *non-refernce* dose in a study i
    Y[i,1:(nd[i]-1)]  ~ dmnorm(mean[i,1:(nd[i]-1)], prec[(b[i]+1):(b[i]+nd[i]-1),1:(nd[i]-1)])

    # the restricted cubic spline dose-response model
    mean[i,1:(nd[i]-1)] <-  beta[i,1]*((X1[i, 2:(nd[i])]-X1[i, 1]))+ beta[i,2]*((X2[i, 2:(nd[i])]-X2[i,1])) #+beta3[i]*(X3[i, 2:(nd[i]+1)]-X3[i,1])

    # window to change the index in precision matrix
    b[i+1] <- b[i]+ nd[i]-1
  }


  # distribution of random effects
  for(i in 1:ns) {
    beta[i,1:2]~dmnorm(beta.pooled[1:2],inv.det*(tau.sq*idmat - rho*idmati))
    u[i]~dnorm(0,0.001)
  }

  # prior distribution for heterogenity
  tau~ dnorm(0,1)%_%T(0,)
  tau.sq <- tau^2
  inv.det <- 1/(tau.sq^2 + rho^2)
  #tau~dwish(R,df)
  rho ~ dunif(-1,1)

  # prior distribution for both regression coeff beta1 and beta2

  beta1.pooled <- beta.pooled[1]
  beta2.pooled <- beta.pooled[2]

  beta.pooled[1] ~ dnorm(0,0.001)
  beta.pooled[2] ~ dnorm(0,0.001)

}

#
















