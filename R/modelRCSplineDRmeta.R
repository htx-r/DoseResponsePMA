#******* 3. Dose-Response Restricted Cubic Splines jags model

modelRCSplineDRmeta <- function(){
  b[1] <-0
  for (i in 1:ns) { ## for each study

    # Likelihood

    Y[i,1:(nd[i])]  ~ dmnorm(mean[i,1:(nd[i])], prec[(b[i]+1):(b[i]+nd[i]),1:(nd[i])])

    mean[i,1:(nd[i])] <-  beta1[i]*(X1[i, 1:(nd[i])]-X1ref[i]) + beta2[i]*(X2[i, 1:(nd[i])]-X2ref[i])#+beta3[i]*(X3[i, 1:(nd[i])]-X3ref[i])


    b[i+1] <- b[i]+ nd[i]
  }

  # Random effect

  for(i in 1:ns) {
    beta1[i]~dnorm(beta1.pooled,prec.tau)
    beta2[i]~dnorm(beta2.pooled,prec.tau)
    #beta3[i]~dnorm(beta3.pooled,prec.tau)
  }

  # Priors
  prec.tau<-1/variance
  variance<-tau*tau
  beta1.pooled ~ dnorm(0,0.1)
  beta2.pooled ~ dnorm(0,1)
  #beta3.pooled ~ dnorm(0,0.1)

  tau~ dnorm(0,10)%_%T(0,)

}
