#******* jags model of spline dose-response model with normal likelihood for OR and RR

modelNorSplineDRmeta <- function(){
  b[1] <-0
  for (i in 1:ns) { ## for each study

    # multivariate normal likelihood of the response Y = logOR or logRR for the *reference* and *non-refernce* dose in a study i
    Y[i,1:(nd[i]-1)]  ~ dmnorm(mean[i,1:(nd[i]-1)], prec[(b[i]+1):(b[i]+nd[i]-1),1:(nd[i]-1)])

    # the restricted cubic spline dose-response model
    mean[i,1:(nd[i]-1)] <-  beta1[i]*((X1[i, 2:(nd[i])]-X1[i, 1]))+ beta2[i]*((X2[i, 2:(nd[i])]-X2[i,1])) #+beta3[i]*(X3[i, 2:(nd[i]+1)]-X3[i,1])

    # window to change the index in precision matrix
    b[i+1] <- b[i]+ nd[i]-1
  }

  # distribution of random effects
  for(i in 1:ns) {
    beta1[i]~dnorm(beta1.pooled,prec.beta)
    beta2[i]~dnorm(beta2.pooled,prec.beta)
  }

  # prior distribution for heterogenity
  prec.beta<-1/variance
  variance<-tau*tau
  tau~ dnorm(0,1)%_%T(0,)

  # prior distribution for both regression coeff beta1 and beta2
  beta1.pooled ~ dnorm(0,10)
  beta2.pooled ~ dnorm(0,10)
}
