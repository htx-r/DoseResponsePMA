## In this file I plot the 4 spline curves with beta1 and beta2
    # 1. true 2. binomial Bayes 3. normal Bayes 4. one-stage

rval<-read.csv("~/Google Drive/DoseResponseNMA/DoseResponseNMA/2019-11-24resORspline40sim1000.csv")
rval <- rval[1:8,]
round(rval[,c('true.beta1',   'BayesB1bias',   'BayesN1bias',   'Freq1bias','BayesB1mse',   'BayesN1mse',     'Freq1mse',
         'true.beta2',   'BayesB2bias',   'BayesN2bias',   'Freq2bias','BayesB2mse',   'BayesN2mse',     'Freq2mse')],4)
## compute our X: doses
d <- 0:10
knots<-unlist(round(quantile(d,c(0.25,0.5,0.75))))
new.dose<-rcspline.eval(d,knots,inclx = T)
new.dose1 <- c(new.dose[,1])
new.dose2 <- c(new.dose[,2])
n.d <- length(d)

## compute our Y = OR; in DR model: OR = exp(beta1*new.dose1+beta2*new.dose2)

# beta1
mybeta1 <- cbind(resORspline40sim1000[c('true.beta1')],bayesB.beta1=resORspline40sim1000[c('true.beta1')]+resORspline40sim1000[c('BayesB1bias')],
                 bayesN.beta1= resORspline40sim1000[c('true.beta1')]+resORspline40sim1000[c('BayesN1bias')],
                 freq.beta1=resORspline40sim1000[c('true.beta1')]+resORspline40sim1000[c('Freq1bias')])
colnames(mybeta1) <- c('true.beta1','bayesB.beta1','bayesN.beta1','freq.beta1')
beta1 <- mybeta1[-c(1,5),]

# beta2
mybeta2 <- cbind(resORspline40sim1000[c('true.beta2')],bayesB.beta1=resORspline40sim1000[c('true.beta2')]+resORspline40sim1000[c('BayesB2bias')],
                 bayesN.beta1= resORspline40sim1000[c('true.beta2')]+resORspline40sim1000[c('BayesN2bias')],
                 freq.beta1=resORspline40sim1000[c('true.beta2')]+resORspline40sim1000[c('Freq2bias')])
colnames(mybeta2) <- c('true.beta2','bayesB.beta2','bayesN.beta2','freq.beta2')
beta2 <- mybeta2[-c(1,5),]

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
ylim=c(0.5,3.5)
ylab='OR'
xlab='dose'
col=1:4
lty=1
beta1[1:4,1]
beta2[1:4,1]
#
par(mfcol=c(3,2),las=1)

matplot(new.dose1,cbind(y[,,1]),type='l'#,main=expression(paste( beta, "1 = 0.04, ", beta, "2= 0"))
        ,lwd=3,col=col,lty=c(1,5,7,10),ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4,ylim=ylim)
matplot(new.dose1,cbind(y[,,2]),type='l'#,main=expression(paste( beta, "1 = 0.1, ", beta, "2= 0.03"))
        ,lwd=3,col=col,lty=1:4,ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4,ylim=ylim)
matplot(new.dose1,cbind(y[,,3]),type='l'#,main=expression(paste( beta, "1 = 0.2, ", beta, "2= -0.2"))
        ,lwd=3,col=col,lty=lty,ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4,ylim=ylim)

# legend(-2,9,legend=c('true','binomial Bayes', 'normal Bayes','Freq' ),
#        col=1:4,horiz = F,lty=1, bty='n',xjust = 0.1,cex = 0.5,lwd=3,xpd = 'NA')

matplot(new.dose1,cbind(y[,,4]),type='l'#,main=expression(paste( beta, "1 = 0.04, ", beta, "2= 0"))
        ,lwd=3,col=col,lty=lty,ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4,ylim=ylim)
matplot(new.dose1,cbind(y[,,5]),type='l'#,main=expression(paste( beta, "1 = 0.1, ", beta, "2= 0.03"))
        ,lwd=3,col=col,lty=lty,ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4,ylim=ylim)
matplot(new.dose1,cbind(y[,,6]),type='l'#,main=expression(paste( beta, "1 = 0.2, ", beta, "2= -0.2"))
        ,lwd=3,col=col,lty=lty,ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4,ylim=ylim)


# legend(-2,9,legend=c('true','binomial Bayes', 'normal Bayes','Freq' ),
#        col=1:4,horiz = F,lty=1, bty='n',xjust = 0.1,cex = 0.5,lwd=3,xpd = 'NA')

## investigate Rhat distribution for binomial and normal.
load('resORspline40sim1000ALL')
dim(resORspline40sim1000ALL)
rownames(resORspline40sim1000ALL)
t(S4ORspline$res2[c('RhatB1','RhatB2','RhatN1','RhatN2'),])
truehist(t(S4ORspline$res2[c('RhatB1'),]))
truehist(t(S4ORspline$res2[c('RhatB2'),]))
truehist(t(S4ORspline$res2[c('RhatN1'),]))
truehist(t(S4ORspline$res2[c('RhatN2'),]))
RhatB2 <- t(S4ORspline$res2[c('RhatB2'),])
mean(S4ORspline$res2['BayesB2',][S4ORspline$res2['RhatB2',]<1.2])
mean(S4ORspline$res2['BayesB2',])
mean(S4ORspline$res2['BayesN2',][S4ORspline$res2['RhatB2',]<1.2])
sum(S4ORspline$res2['RhatN2',]<1.2)
mean(S1ORspline$res2['BayesB2',][S1ORspline$res2['RhatB2',]<1.2])
mean(S1ORspline$res2['BayesB2',])
