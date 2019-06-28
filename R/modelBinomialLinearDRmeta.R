#******* Linear dose-response jags model using binomial likelihood

modelBinomialLinearDRmeta <- function(){

for (i in 1:ns) { ## for each study
  log(p[i,1])<- u[i]
  for (j in 2:(nd[i])) { ## for each dose
    # Likelihood
r[i,j] ~ dbinom(p[i,j],n[i,j])
log(p[i,j])<- u[i] + delta[i,j]
delta[i,j] <-  beta[i]*(X[i,(j-1)]-Xref[i])

  }
}
  # Random effect

  for(i in 1:ns) {
    beta[i]~dnorm(beta.pooled,prec.tau)
 u[i]~dnorm(0,0.1)%_%T(,0)
  }

  # Priorsx
  prec.tau<-1/variance
  variance<-tau*tau
  beta.pooled ~ dnorm(0,0.1)

  tau~ dnorm(0,0.01)%_%T(0,)

}



# I tried to fit poission distribution
# model2LinearDRmeta <- function(){
#
#   for (i in 1:ns) { ## for each study
#     log(p[i,1])<- u[i]
#     for (j in 2:(nd[i])) { ## for each dose
#       # Likelihood
#       r[i,j] ~ dpois(p[i,j]*n[i,j])
#       log(p[i,j])<- u[i] + delta[i,j]
#       delta[i,j] <-  beta[i]*(X[i,(j-1)]-Xref[i])
#
#     }
#   }
#   # Random effect
#
#   for(i in 1:ns) {
#     beta[i]~dnorm(beta.pooled,prec.tau)
#     u[i]~dnorm(0,0.1)
#   }
#
#   # Priors
#   prec.tau<-1/variance
#   variance<-tau*tau
#   beta.pooled ~ dnorm(0,0.1)
#
#   tau~ dnorm(0,0.01)%_%T(0,)
#
# }
#
