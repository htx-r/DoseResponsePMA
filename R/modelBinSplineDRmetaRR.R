#******* Linear dose-response jags model using binomial likelihood

modelBinSplineDRmetaRR <- function(){

  for (i in 1:ns) { ## for each study
    log(p[i,1])<- u[i]
    r[i,1] ~ dbinom(p[i,1],n[i,1])

    for (j in 2:(nd[i])) { ## for each dose
      # Likelihood
      r[i,j] ~ dbinom(p[i,j],n[i,j])
      log(p[i,j])<- u[i] + delta[i,j]
      delta[i,j] <-  beta1[i]*(X1[i,j]-X1[i,1]) + beta2[i]*(X2[i,j]-X2[i,1])#+beta3[i]*(X3[i,j]-X3[i,1])

    }
  }
  # Random effect

  for(i in 1:ns) {
    beta1[i]~dnorm(beta1.pooled,prec.beta)
    beta2[i]~dnorm(beta2.pooled,prec.beta)
    #beta3[i]~dnorm(beta3.pooled,prec.beta)

    u[i]~dnorm(0,0.0001)%_%T(,0)

  }

  # Priors
  prec.beta<-1/variance
  variance<-tau*tau
  beta1.pooled ~ dnorm(0,0.01)
  beta2.pooled ~ dnorm(0,0.01)
  #beta3.pooled ~ dnorm(0,0.01)

  tau~ dnorm(0,4)%_%T(0,)
  # for (i in 1:new.n) {
  #   newbeta1[i]~dnorm(beta1.pooled,prec.beta)
  #   newbeta2[i]~dnorm(beta2.pooled,prec.beta)
  #   newbeta3[i]~dnorm(beta3.pooled,prec.beta)
  #
  #   newY[i] <-newbeta1[i]*new.dose1[i]+newbeta2[i]*new.dose2[i]+newbeta3[i]*new.dose3[i]
  #   newOR[i]<-exp(newY[i])
  # }

}

#
