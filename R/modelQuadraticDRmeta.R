#******* Quadratic dose-response jags model with normal likelihood for logRR=Y

modelQuadraticDRmeta <- function(){
  b[1] <-0
  for (i in 1:ns) { ## for each study

    # Likelihood

    Y[i,1:(nd[i])]  ~ dmnorm(mean[i,1:(nd[i])], prec[(b[i]+1):(b[i]+nd[i]),1:(nd[i])])

    mean[i,1:(nd[i])] <-  beta1[i]*(X[i, 1:(nd[i])]-Xref[i])+ beta2[i]*(X[i, 1:(nd[i])]^2-Xref[i]^2)

    b[i+1] <- b[i]+ nd[i]
  }

  # Random effect

  for(i in 1:ns) {
    beta1[i]~dnorm(beta1.pooled,prec.tau)
    beta2[i]~dnorm(beta2.pooled,prec.tau)

  }

  # Priors
  prec.tau<-1/variance
  variance<-tau*tau
  beta1.pooled ~ dnorm(0,0.1)
  beta2.pooled ~ dnorm(0,0.1)

  tau~ dnorm(0,0.01)%_%T(0,)

}
