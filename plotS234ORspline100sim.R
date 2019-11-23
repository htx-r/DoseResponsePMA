library(MASS)
load('S2ORspline')
load('S3ORspline')
load('S4ORspline')
rval234 <- rbind.data.frame(S2ORspline$res1,S3ORspline$res1,S4ORspline$res1)
## compute our X: doses
d <- 0:10
knots<-unlist(round(quantile(d,c(0.25,0.5,0.75))))
new.dose<-rcspline.eval(d,knots,inclx = T)
new.dose1 <- c(new.dose[,1])
new.dose2 <- c(new.dose[,2])
n.d <- length(d)

## compute our Y = OR; in DR model: OR = exp(beta1*new.dose1+beta2*new.dose2)

# beta1
beta1 <- cbind(rval234[c('true.beta1')],bayesB.beta1=rval234[c('true.beta1')]+rval234[c('BayesB1bias')],
                 bayesN.beta1= rval234[c('true.beta1')]+rval234[c('BayesN1bias')],
                 freq.beta1=rval234[c('true.beta1')]+rval234[c('Freq1bias')])
colnames(beta1) <- c('true.beta1','bayesB.beta1','bayesN.beta1','freq.beta1')


# beta2
beta2 <- cbind(rval234[c('true.beta2')],bayesB.beta1=rval234[c('true.beta2')]+rval234[c('BayesB2bias')],
                 bayesN.beta1= rval234[c('true.beta2')]+rval234[c('BayesN2bias')],
                 freq.beta1=rval234[c('true.beta2')]+rval234[c('Freq2bias')])
colnames(beta2) <- c('true.beta2','bayesB.beta2','bayesN.beta2','freq.beta2')

# finally, y
y <- array(NA,dim=c(n.d,4,nrow(beta1)))
for (k in 1:nrow(beta1)) {
  for (j in 1:4) {
    y[1:n.d,j,k] <- exp(beta1[k,j]*new.dose1+beta2[k,j]*new.dose2)
  }
}



####################################################################################
####################################################################################
##########  Plot the 4 DR curves: true, bayesBin, BayesNor, Freq ###################
####################################################################################
####################################################################################

# unified plots' settings
ylim=c(1,3.5)
ylab='OR'
xlab='dose'
col=1:4
lty=1
beta1[1:4,1]
beta2[1:4,1]
#
par(mfrow=c(2,2),las=1)

matplot(new.dose1,cbind(y[,,1]),type='l',main=expression(paste( beta, "1 = 0.04, ", beta, "2= 0"))
        ,lwd=3,col=col,lty=lty,ylab=ylab,xlab=xlab,ylim=ylim)
matplot(new.dose1,cbind(y[,,2]),type='l',main=expression(paste( beta, "1 = 0.1, ", beta, "2= 0.03"))
        ,lwd=3,col=col,lty=lty,ylab=ylab,xlab=xlab,ylim=ylim)
matplot(new.dose1,cbind(y[,,3]),type='l',main=expression(paste( beta, "1 = 0.2, ", beta, "2= -0.2"))
        ,lwd=3,col=col,lty=lty,ylab=ylab,xlab=xlab,ylim=ylim)

legend('topleft',legend=c('true','binomial Bayes', 'normal Bayes','Freq' ),
       col=1:4,horiz = F,lty=1, bty='n',xjust = 0.1,cex = 0.5,lwd=3,xpd = 'NA')


##########  ##########  ##########  ##########
##  Rhat distribution for binomial and normal
##########  ##########  ##########  ##########
par(mfrow=c(2,2))
truehist(t(S2ORspline$res2[c('RhatB1'),]))
truehist(t(S2ORspline$res2[c('RhatB2'),]))
truehist(t(S2ORspline$res2[c('RhatN1'),]))
truehist(t(S2ORspline$res2[c('RhatN2'),]))

par(mfrow=c(2,2))
truehist(t(S3ORspline$res2[c('RhatB1'),]))
truehist(t(S3ORspline$res2[c('RhatB2'),]))
truehist(t(S3ORspline$res2[c('RhatN1'),]))
truehist(t(S3ORspline$res2[c('RhatN2'),]))

par(mfrow=c(2,2))
truehist(t(S4ORspline$res2[c('RhatB1'),]))
truehist(t(S4ORspline$res2[c('RhatB2'),]))
truehist(t(S4ORspline$res2[c('RhatN1'),]))
truehist(t(S4ORspline$res2[c('RhatN2'),]))



