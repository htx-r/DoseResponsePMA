#******* Linear dose-response jags model using binomial likelihood

modelBinSplineDRmetaOR <- function(){

  for (i in 1:ns) { ## for each study
    r[i,1] ~ dbinom(p[i,1],n[i,1])
    logit(p[i,1])<- u[i]
    for (j in 2:(nd[i])) { ## for each dose
      # Likelihood
      r[i,j] ~ dbinom(p[i,j],n[i,j])
      logit(p[i,j])<- u[i] + delta[i,j]
      delta[i,j] <-  beta1[i]*(X1[i,(j)]-X1[i,1])+beta2[i]*(X2[i,(j)]-X2[i,1])

    }
  }
  # Random effect

  for(i in 1:ns) {
    # beta[i] <- xi*eta[i]
    # eta[i]~dnorm(beta.pooled,tau.eta)
    beta1[i] ~dnorm(beta1.pooled,prec.tau)
    beta2[i] ~dnorm(beta2.pooled,prec.tau)

    u[i]~dnorm(0,0.001)
  }

  # Priors

  # xi~dnorm(0,0.1)
  # tau.eta~dgamma(5,5)
  # tau <- abs(xi)/sqrt(tau.eta)
  # log(tau) <- log.tau
  # log.tau~ dunif(-20,20)
  prec.tau<-1/variance
  variance<-tau*tau
  tau~dnorm(0,400)%_%T(0,)
  beta1.pooled ~ dnorm(0,0.001)
  beta2.pooled ~ dnorm(0,0.001)

}

