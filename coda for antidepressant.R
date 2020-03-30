load('antidepORspline')
library(coda)
library(R2jags)
 sim.beta1 <-  doseresORsplineBin$BUGSoutput$sims.array[,,'beta1.pooled']
 sim.beta2 <-  doseresORsplineBin$BUGSoutput$sims.array[,,'beta2.pooled']

m <-as.mcmc(doseresORsplineBin)

# for beta1.pooled , beta2.pooled and tau, the following are some coda plots and measures of them

# traceplots
coda::traceplot(m[,'beta1.pooled',drop(FALSE)])
coda::traceplot(m[,'beta2.pooled',drop(FALSE)])
coda::traceplot(m[,'tau',drop(FALSE)])

# effective size
effectiveSize(m[,'beta1.pooled',drop(FALSE)])
effectiveSize(m[,'beta2.pooled',drop(FALSE)])
effectiveSize(m[,'tau',drop(FALSE)])

# Raftery diagnosis measure
raftery.diag(m[,'beta1.pooled',drop(FALSE)])
raftery.diag(m[,'beta1.pooled',drop(FALSE)])


# geweke plot
geweke.plot(m[,'beta1.pooled',drop(FALSE)])
geweke.plot(m[,'beta2.pooled',drop(FALSE)])
geweke.plot(m[,'tau',drop(FALSE)])

# gelman plot
gelman.plot(m[,'beta1.pooled',drop(FALSE)],ylim=c(1,3.8))
gelman.plot(m[,'beta2.pooled',drop(FALSE)])
gelman.plot(m[,'tau',drop(FALSE)],ylim=c(1,3.8))

# gelman reduction factor R_hat
gelman.diag(m[,'beta1.pooled',drop(FALSE)],autoburnin = TRUE)
gelman.diag(m[,'beta2.pooled',drop(FALSE)],autoburnin = TRUE)
gelman.diag(m[,'tau',drop(FALSE)],autoburnin = TRUE)

#
heidel.diag(m[,'beta1.pooled',drop(FALSE)])
heidel.diag(m[,'beta2.pooled',drop(FALSE)])
heidel.diag(m[,'tau',drop(FALSE)])


data(line)
x1 <- line[[1]]                    #Select first chain
x2 <- line[,1, drop=FALSE]         #Select first var from all chains
varnames(x2) == varnames(line)[1]  #TRUE
