#******* jags model of spline dose-response model with binomial likelihood for OR

modelBinSplineDRmetaORBiv <- function(){
  for (i in 1:ns) { ## for each study
    # binomial likelihood of number of events in the *refernce* dose level in a study i
    r[i,1] ~ dbinom(p[i,1],n[i,1])

    # logit parametrization of probabilities at each *refernce* dose level: by that exp(beta)= OR
    logit(p[i,1]) <- u[i]

    for (j in 2:(nd[i])) { ## for each dose
      # binomial likelihood of number of events for the *non-refernce* dose in a study i
      r[i,j] ~ dbinom(p[i,j],n[i,j])

      # logit parametrization of probabilities at each *non-refernce* dose level: by that exp(beta)= OR
      logit(p[i,j]) <- u[i] + delta[i,j]
      delta[i,j] <-   beta[i,1]*(X1[i,j]-X1[i,1]) + beta[i,2]*(X2[i,j]-X2[i,1])
    }

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
  beta.pooled[1] ~ dnorm(0,0.001)
  beta.pooled[2] ~ dnorm(0,0.001)

}

#
















