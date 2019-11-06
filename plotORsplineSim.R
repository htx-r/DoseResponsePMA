## I need to set common 1. legend 2. title and ylim
# I may use the matplot

load('resORspline40sim1000')

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
beta1 <- mybeta1[-c(1,6),]

# beta2
mybeta2 <- cbind(resORspline40sim1000[c('true.beta2')],bayesB.beta1=resORspline40sim1000[c('true.beta2')]+resORspline40sim1000[c('BayesB2bias')],
               bayesN.beta1= resORspline40sim1000[c('true.beta2')]+resORspline40sim1000[c('BayesN2bias')],
               freq.beta1=resORspline40sim1000[c('true.beta2')]+resORspline40sim1000[c('Freq2bias')])
colnames(mybeta2) <- c('true.beta2','bayesB.beta2','bayesN.beta2','freq.beta2')
beta2 <- mybeta2[-c(1,6),]

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
ylim=c(0.5,6.5)
ylab='OR'
xlab='dose'
col=1:4
lty=1
beta1[1:4,1]
beta2[1:4,1]
#
par(mfrow=c(2,2),las=1)

matplot(new.dose1,cbind(y[,,1]),type='l',main=expression(paste( beta, "1 = 0.04, ", beta, "2= 0"))
        ,lwd=3,col=col,lty=lty,ylab=ylab,xlab=xlab)
matplot(new.dose1,cbind(y[,,2]),type='l',main=expression(paste( beta, "1 = 0.1, ", beta, "2= 0.03"))
        ,lwd=3,col=col,lty=lty,ylab=ylab,xlab=xlab)
matplot(new.dose1,cbind(y[,,3]),type='l',main=expression(paste( beta, "1 = 0.2, ", beta, "2= -0.2"))
        ,lwd=3,col=col,lty=lty,ylab=ylab,xlab=xlab)
matplot(new.dose1,cbind(y[,,4]),type='l',main=expression(paste( beta, "1 = 0.2, ", beta, "2= -0.3"))
        ,lwd=3,col=col,lty=lty,ylab=ylab,xlab=xlab)

legend(-2,9,legend=c('true','binomial Bayes', 'normal Bayes','Freq' ),
       col=1:4,horiz = F,lty=1, bty='n',xjust = 0.1,cex = 0.5,lwd=3,xpd = 'NA')

par(mfrow=c(2,2))
matplot(new.dose1,cbind(y[,,5]),type='l',main=expression(paste( beta, "1 = 0.04, ", beta, "2= 0"))
        ,lwd=3,col=col,lty=lty,ylab=ylab,xlab=xlab)
matplot(new.dose1,cbind(y[,,6]),type='l',main=expression(paste( beta, "1 = 0.1, ", beta, "2= 0.03"))
        ,lwd=3,col=col,lty=lty,ylab=ylab,xlab=xlab)
matplot(new.dose1,cbind(y[,,7]),type='l',main=expression(paste( beta, "1 = 0.2, ", beta, "2= -0.2"))
        ,lwd=3,col=col,lty=lty,ylab=ylab,xlab=xlab)
matplot(new.dose1,cbind(y[,,8]),type='l',main=expression(paste( beta, "1 = 0.2, ", beta, "2= -0.3"))
        ,lwd=3,col=col,lty=lty,ylab=ylab,xlab=xlab)

legend(-2,9,legend=c('true','binomial Bayes', 'normal Bayes','Freq' ),
       col=1:4,horiz = F,lty=1, bty='n',xjust = 0.1,cex = 0.5,lwd=3,xpd = 'NA')


####################################################################################
####################################################################################
##########  Plot the true DR curve with its confidence interval ####################
####################################################################################
####################################################################################

beta1<- c(0.04,0.1,0.2,0.2)
beta2 <- c(0,0.03,-0.2,-0.3 )
tau <- c(0.001,0.01)

## larger tau
par(mfrow=c(2,2))

# S2
plot(new.dose1,exp(beta1[1]*new.dose1+beta2[1]*new.dose2),col=1,type='l',ylim = c(0.5,4)
     ,las=1,ylab='OR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) #  true
lines(new.dose1,exp((beta1[1]+2*tau[1])*new.dose1+(beta2[1]+2*tau[1])*new.dose2),lty=2)
lines(new.dose1,exp((beta1[1]-2*tau[1])*new.dose1+(beta2[1]-2*tau[1])*new.dose2),lty=2)

lines(new.dose1,exp((beta1[1]+2*tau[2])*new.dose1+(beta2[1]+2*tau[2])*new.dose2),lty=2,col=2)
lines(new.dose1,exp((beta1[1]-2*tau[2])*new.dose1+(beta2[1]-2*tau[2])*new.dose2),lty=2,col=2)

# S3
plot(new.dose1,exp(beta1[2]*new.dose1+beta2[2]*new.dose2),col=1,type='l',ylim = c(0.5,6)
     ,las=1,ylab='OR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) #  true
lines(new.dose1,exp((beta1[2]+2*tau[1])*new.dose1+(beta2[2]+2*tau[1])*new.dose2),lty=2)
lines(new.dose1,exp((beta1[2]-2*tau[1])*new.dose1+(beta2[2]-2*tau[1])*new.dose2),lty=2)

lines(new.dose1,exp((beta1[2]+2*tau[2])*new.dose1+(beta2[2]+2*tau[2])*new.dose2),lty=2,col=2)
lines(new.dose1,exp((beta1[2]-2*tau[2])*new.dose1+(beta2[2]-2*tau[2])*new.dose2),lty=2,col=2)

# S4
plot(new.dose1,exp(beta1[3]*new.dose1+beta2[3]*new.dose2),col=1,type='l',ylim = c(0.5,4)
     ,las=1,ylab='OR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) #  true
lines(new.dose1,exp((beta1[3]+2*tau[1])*new.dose1+(beta2[3]+2*tau[1])*new.dose2),lty=2)
lines(new.dose1,exp((beta1[3]-2*tau[1])*new.dose1+(beta2[3]-2*tau[1])*new.dose2),lty=2)

lines(new.dose1,exp((beta1[3]+2*tau[2])*new.dose1+(beta2[3]+2*tau[2])*new.dose2),lty=2,col=2)
lines(new.dose1,exp((beta1[3]-2*tau[2])*new.dose1+(beta2[3]-2*tau[2])*new.dose2),lty=2,col=2)




## simulate the antidepressant scenario for OR: beta1=0.25

load('antidepSimulation')
# new dose for predictions and plot
library(rms)
knots = c(10,20,50)
new.dose <- seq(0,80,1)
new.dose1 <- c(rcs(new.dose,knots)[,1])
new.dose2 <- c(rcs(new.dose,knots)[,2])

par(mfrow=c(1,2))
plot(new.dose1,exp(antidepSim20$true.beta1*new.dose1+antidepSim20$true.beta2*new.dose2),col=1,type='l',ylim = c(0,8)
     ,las=1,ylab='OR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) # true
lines(new.dose1,exp((antidepSim20$true.beta1+antidepSim20$BayesB1bias)*new.dose1+(antidepSim20$true.beta2+antidepSim20$BayesB2bias)*new.dose2)
      ,col=2,lwd=3) # bayes binomial
lines(new.dose1,exp((antidepSim20$true.beta1+antidepSim20$BayesN1bias)*new.dose1+(antidepSim20$true.beta2+antidepSim20$BayesN2bias)*new.dose2)
      ,col=3,lwd=3) # bayes normal
lines(new.dose1,exp((antidepSim20$true.beta1+antidepSim20$Freq1bias)*new.dose1+(antidepSim20$true.beta2+antidepSim20$Freq2bias)*new.dose2)
      ,col=4,lwd=3) # freq



plot(new.dose1,exp(antidepSim40$true.beta1*new.dose1+antidepSim40$true.beta2*new.dose2),col=1,type='l',ylim = c(0,2)
     ,las=1,ylab='OR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) # true
lines(new.dose1,exp((antidepSim40$true.beta1+antidepSim40$BayesB1bias)*new.dose1+(antidepSim40$true.beta2+antidepSim40$BayesB2bias)*new.dose2)
      ,col=2,lwd=3) # bayes binomial
lines(new.dose1,exp((antidepSim40$true.beta1+antidepSim40$BayesN1bias)*new.dose1+(antidepSim40$true.beta2+antidepSim40$BayesN2bias)*new.dose2)
      ,col=3,lwd=3) # bayes normal
lines(new.dose1,exp((antidepSim40$true.beta1+antidepSim40$Freq1bias)*new.dose1+(antidepSim40$true.beta2+antidepSim40$Freq2bias)*new.dose2)
      ,col=4,lwd=3) # freq

 legend('topright',legend=c('True','binomial Bayes', 'normal Bayes','Freq' ),
        col=1:4,horiz = F,lty=1, bty='n',xjust = 0,cex = 0.7,lwd=3)

