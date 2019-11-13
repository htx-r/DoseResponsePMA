#******* Linear dose-response jags model using binomial likelihood

modelBinSplineDRmetaOR <- function(){

  for (i in 1:ns) { ## for each study
    logit(p[i,1])<- u[i]
    r[i,1] ~ dbinom(p[i,1],n[i,1])

    for (j in 2:(nd[i])) { ## for each dose
      # Likelihood
      r[i,j] ~ dbinom(p[i,j],n[i,j])
      # logit(p[i,1])<- u[i] - delta[i,2]/3 - delta[i,3]/3
      # logit(p[i,2])<- u[i] + 2*delta[i,2]/3 - delta[i,3]/3
      # logit(p[i,3])<- u[i] - delta[i,2]/3 +2*delta[i,3]/3
      #
      # delta[i,2] <-  beta1[i]*(X1[i,2]-X1[i,1]) + beta2[i]*(X2[i,2]-X2[i,1])
      # delta[i,3] <-  beta1[i]*(X1[i,3]-X1[i,1]) + beta2[i]*(X2[i,3]-X2[i,1])

      logit(p[i,j])<- u[i] + delta[i,j]
      delta[i,j] <-   beta1[i]*(X1[i,j]-X1[i,1]) +beta2[i]*(X2[i,j]-X2[i,1]) #

    }
  }
  # Random effect

  for(i in 1:ns) {
    beta1[i]~dnorm(beta1.pooled,prec.beta)
    beta2[i]~dnorm(beta2.pooled,prec.beta)
    u[i]~dnorm(0,0.001)

  }

  # Priors
  beta1.pooled ~ dnorm(0,0.001)
  beta2.pooled ~ dnorm(0,0.001)

  prec.beta<-1/variance
  variance<-tau*tau
  tau~ dnorm(0,1)%_%T(0,)

}

