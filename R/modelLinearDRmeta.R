#******* Linear dose-response jags model with normal likelihood for logRR=Y

modelLinearDRmeta <- function(){
  b[1] <-0
  for (i in 1:ns) { ## for each study

    # Likelihood

    Y[i,1:(nd[i])]  ~ dmnorm(mean[i,1:(nd[i])], prec[(b[i]+1):(b[i]+nd[i]),1:(nd[i])])

    mean[i,1:(nd[i])] <-  beta[i]*(X[i, 2:(nd[i]+1)]-X[i, 1])

    b[i+1] <- b[i]+ nd[i]
  }

  # Random effect

  for(i in 1:ns) {
    beta[i]~dnorm(beta.pooled,prec.tau)

  }

  # Priors
  prec.tau<-1/variance
  variance<-tau*tau
  beta.pooled ~ dnorm(0,0.1)

  tau~ dnorm(0,10)%_%T(0,)

  ## Predictions

   for (i in 1:new.n) {
  newbeta[i]~dnorm(beta.pooled,prec.tau)
 newY[i] <-newbeta[i]*new.dose[i]
 newRR[i]<-exp(newY[i])
   }
}

