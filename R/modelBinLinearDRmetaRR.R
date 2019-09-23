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
    u[i]~dnorm(0,0.01)%_%T(,0)
  }

  # Priors
  prec.beta<-1/variance
  variance<-tau*tau
  tau~ dnorm(0,10)%_%T(0,)
  beta.pooled ~ dnorm(0,0.01)
}

# modelBinomialLinearDRmeta <- function(){
#
#   for (i in 1:ns) { ## for each study
#     log(p[i,1])<- u[i]
#     for (j in 2:(nd[i])) { ## for each dose
#       # Likelihood
#       r[i,j] ~ dbinom(p[i,j],n[i,j])
#       log(p[i,j])<- u[i] + delta[i,j]
#       delta[i,j] <-  beta[i]*(X[i,(j-1)]-Xref[i])
#
#     }
#   }
#   # Random effect
#
#   for(i in 1:ns) {
#     beta[i]~dnorm(beta.pooled,prec.beta)
#     u[i]~dnorm(0,0.1)%_%T(,0)
#   }
#
#   # Priors
#   prec.beta<-1/variance
#   variance<-tau*tau
#   tau~ dnorm(0,0.01)%_%T(0,)
#   beta.pooled ~ dnorm(0,1e-6)
#
#   for (i in 1:ns) { ## for each study
#     ones[i,1] <- 1
#     for (j in 2:(nd[i])) { ## for each dose
#       ones[i,j]  ~ dbern(C1[i,j])
#       C1[i,j] <- step(1-p[i,j])
#
#     }
#   }
#
# }
