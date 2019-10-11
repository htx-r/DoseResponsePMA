#******* Restricted cubic splines dose-response  jags model with normal likelihood for the logRR

modelNorSplineDRmeta <- function(){
  b[1] <-0
  for (i in 1:ns) { ## for each study

    # Likelihood

    Y[i,1:(nd[i])]  ~ dmnorm(mean[i,1:(nd[i])], prec[(b[i]+1):(b[i]+nd[i]),1:(nd[i])])

    mean[i,1:(nd[i])] <-  beta1[i]*(X1[i, 2:(nd[i]+1)]-X1[i, 1])+ beta2[i]*(X2[i, 2:(nd[i]+1)]-X2[i,1]) #+beta3[i]*(X3[i, 2:(nd[i]+1)]-X3[i,1])


    b[i+1] <- b[i]+ nd[i]
  }

  # Random effect

  for(i in 1:ns) {
    beta1[i]~dnorm(beta1.pooled,prec.beta)
    beta2[i]~dnorm(beta2.pooled,prec.beta)
  }

  #Priors
  prec.beta<-1/variance
  variance<-tau*tau
  tau~ dnorm(0,0.1)%_%T(0,)

  # prec.beta2<-1/variance2
  # variance2<-tau2*tau2
  # tau2~ dnorm(0,20)%_%T(0,)


  beta1.pooled ~ dnorm(0,0.1)
  beta2.pooled ~ dnorm(0,0.1)

  ## Predictions

  # for (i in 1:new.n) {
  #   newbeta1[i]~dnorm(beta1.pooled,prec.beta)
  #   newbeta2[i]~dnorm(beta2.pooled,prec.beta)
  #   newbeta3[i]~dnorm(beta3.pooled,prec.beta)
  #
  #   newY[i] <-newbeta1[i]*new.dose1[i]+newbeta2[i]*new.dose2[i]+newbeta3[i]*new.dose3[i]
  #   newOR[i]<-exp(newY[i])
  # }
}
