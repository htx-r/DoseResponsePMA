#******* Linear dose-response jags model using binomial likelihood

modelBinLinearDRmetaRR <- function(){

  for (i in 1:ns) { ## for each study
    r[i,1] ~ dbinom(p[i,1],n[i,1])
    log(p[i,1])<- u[i]
    for (j in 2:(nd[i])) { ## for each dose
      # Likelihood
      r[i,j] ~ dbinom(p[i,j],n[i,j])
      log(p[i,j])<- u[i] + delta[i,j]
      delta[i,j] <-  beta[i]*(X[i,(j)]-X[i,1])

    }
  }
  # Random effect

  for(i in 1:ns) {
    beta[i]~dnorm(beta.pooled,prec.beta)
    u[i]~dnorm(0,0.1)%_%T(,0)
  }

  # Priors
  prec.beta<-1/variance
  variance<-tau*tau
  tau~ dnorm(0,4)%_%T(0,)
  beta.pooled ~ dnorm(0,0.1)
}
