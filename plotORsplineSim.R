
## I need to set common 1. legend 2. title and ylim
# I may use the matplot
load('resORspline')

## the doses
d <- 0:10
knots<-unlist(round(quantile(d,c(0.25,0.5,0.75))))
new.dose<-rcspline.eval(d,knots,inclx = T)
new.dose1 <- c(new.dose[,1])
new.dose2 <- c(new.dose[,2])
beta1 <- c(0.04,0.1,0.2,0.2)
beta2 <- c(0,0.03,-0.2,-0.3)
tau <- c(0.01,0.001)
## larger tau
par(mfrow=c(2,2))

plot(new.dose1,exp(beta1[1]*new.dose1+beta2[1]*new.dose2),col=1,type='l',ylim = c(0.5,4)
     ,las=1,ylab='OR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) #  true
lines(new.dose1,exp((beta1[1]+2*tau[1])*new.dose1+(beta2[1]+2*tau[1])*new.dose2),lty=2)
lines(new.dose1,exp((beta1[1]-2*tau[1])*new.dose1+(beta2[1]-2*tau[1])*new.dose2),lty=2)

lines(new.dose1,exp((beta1[1]+2*tau[2])*new.dose1+(beta2[1]+2*tau[2])*new.dose2),lty=2,col=2)
lines(new.dose1,exp((beta1[1]-2*tau[2])*new.dose1+(beta2[1]-2*tau[2])*new.dose2),lty=2,col=2)

plot(new.dose1,exp(beta1[2]*new.dose1+beta2[2]*new.dose2),col=1,type='l',ylim = c(0.5,6)
     ,las=1,ylab='OR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) #  true
lines(new.dose1,exp((beta1[2]+2*tau[1])*new.dose1+(beta2[2]+2*tau[1])*new.dose2),lty=2)
lines(new.dose1,exp((beta1[2]-2*tau[1])*new.dose1+(beta2[2]-2*tau[1])*new.dose2),lty=2)

lines(new.dose1,exp((beta1[2]+2*tau[2])*new.dose1+(beta2[2]+2*tau[2])*new.dose2),lty=2,col=2)
lines(new.dose1,exp((beta1[2]-2*tau[2])*new.dose1+(beta2[2]-2*tau[2])*new.dose2),lty=2,col=2)

plot(new.dose1,exp(beta1[3]*new.dose1+beta2[3]*new.dose2),col=1,type='l',ylim = c(0.5,4)
     ,las=1,ylab='OR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) #  true
lines(new.dose1,exp((beta1[3]+2*tau[1])*new.dose1+(beta2[3]+2*tau[1])*new.dose2),lty=2)
lines(new.dose1,exp((beta1[3]-2*tau[1])*new.dose1+(beta2[3]-2*tau[1])*new.dose2),lty=2)

lines(new.dose1,exp((beta1[3]+2*tau[2])*new.dose1+(beta2[3]+2*tau[2])*new.dose2),lty=2,col=2)
lines(new.dose1,exp((beta1[3]-2*tau[2])*new.dose1+(beta2[3]-2*tau[2])*new.dose2),lty=2,col=2)


## S2
plot(new.dose1,exp(S2ORspline$true.beta1*new.dose1+S2ORspline$true.beta2*new.dose2),col=1,type='l',ylim = c(0.5,8)
     ,las=1,ylab='OR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) #  true
lines(new.dose1,exp((S2ORspline$true.beta1+S2ORspline$BayesB1bias)*new.dose1+(S2ORspline$true.beta2+S2ORspline$BayesB2bias)*new.dose2)
      ,col=2,lwd=3) # bayes binomial
lines(new.dose1,exp((S2ORspline$true.beta1+S2ORspline$BayesN1bias)*new.dose1+(S2ORspline$true.beta2+S2ORspline$BayesN2bias)*new.dose2)
      ,col=3,lwd=3) # bayes normal
lines(new.dose1,exp((S2ORspline$true.beta1+S2ORspline$Freq1bias)*new.dose1+(S2ORspline$true.beta2+S2ORspline$Freq2bias)*new.dose2)
      ,col=4,lwd=3) # freq
# legend(0,1.5,legend=c('True','binomial Bayes', 'normal Bayes','Freq' ),
#        col=1:4,horiz = T,lty=1, bty='n',xjust = 0,cex = 0.7,lwd=3)


# S3
plot(new.dose1,exp(S3ORspline$true.beta1*new.dose1+S3ORspline$true.beta2*new.dose2),col=1,type='l',ylim = c(0.5,8)
     ,las=1,ylab='OR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) # true
lines(new.dose1,exp((S3ORspline$true.beta1+S3ORspline$BayesB1bias)*new.dose1+(S3ORspline$true.beta2+S3ORspline$BayesB2bias)*new.dose2)
      ,col=2,lwd=3) # bayes binomial
lines(new.dose1,exp((S3ORspline$true.beta1+S3ORspline$BayesN1bias)*new.dose1+(S3ORspline$true.beta2+S3ORspline$BayesN2bias)*new.dose2)
      ,col=3,lwd=3) # bayes normal
lines(new.dose1,exp((S3ORspline$true.beta1+S3ORspline$Freq1bias)*new.dose1+(S3ORspline$true.beta2+S3ORspline$Freq2bias)*new.dose2)
      ,col=4,lwd=3) # freq
# legend(0,0.7,legend=c('True','binomial Bayes', 'normal Bayes','Freq' ),
#        col=1:4,horiz = T,lty=1, bty='n',xjust = 0,cex = 0.7,lwd=3)

# S4
plot(new.dose1,exp(S4ORspline$true.beta1*new.dose1+S4ORspline$true.beta2*new.dose2),col=1,type='l',ylim = c(0.5,8)
     ,las=1,ylab='OR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) # true
lines(new.dose1,exp((S4ORspline$true.beta1+S4ORspline$BayesB1bias)*new.dose1+(S4ORspline$true.beta2+S4ORspline$BayesB2bias)*new.dose2)
      ,col=2,lwd=3) # bayes binomial
lines(new.dose1,exp((S4ORspline$true.beta1+S4ORspline$BayesN1bias)*new.dose1+(S4ORspline$true.beta2+S4ORspline$BayesN2bias)*new.dose2)
      ,col=3,lwd=3) # bayes normal
lines(new.dose1,exp((S4ORspline$true.beta1+S4ORspline$Freq1bias)*new.dose1+(S4ORspline$true.beta2+S4ORspline$Freq2bias)*new.dose2)
      ,col=4,lwd=3) # freq
# legend(0,,legend=c('True','binomial Bayes', 'normal Bayes','Freq' ),
#        col=1:4,horiz = T,lty=1, bty='n',xjust = 0,cex = 0.7,lwd=3)


# S5
plot(new.dose1,exp(S5ORspline$true.beta1*new.dose1+S5ORspline$true.beta2*new.dose2),col=1,type='l',ylim = c(0.5,8)
     ,las=1,ylab='OR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) # true
lines(new.dose1,exp((S5ORspline$true.beta1+S5ORspline$BayesB1bias)*new.dose1+(S5ORspline$true.beta2+S5ORspline$BayesB2bias)*new.dose2)
      ,col=2,lwd=3) # bayes binomial
lines(new.dose1,exp((S5ORspline$true.beta1+S5ORspline$BayesN1bias)*new.dose1+(S5ORspline$true.beta2+S5ORspline$BayesN2bias)*new.dose2)
      ,col=3,lwd=3) # bayes normal
lines(new.dose1,exp((S5ORspline$true.beta1+S5ORspline$Freq1bias)*new.dose1+(S5ORspline$true.beta2+S5ORspline$Freq2bias)*new.dose2)
      ,col=4,lwd=3) # freq


# S6
plot(new.dose1,exp(S6ORspline$true.beta1*new.dose1+S6ORspline$true.beta2*new.dose2),col=1,type='l',ylim = c(0.5,8)
     ,las=1,ylab='OR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) # true
lines(new.dose1,exp((S6ORspline$true.beta1+S6ORspline$BayesB1bias)*new.dose1+(S6ORspline$true.beta2+S6ORspline$BayesB2bias)*new.dose2)
      ,col=2,lwd=3) # bayes binomial
lines(new.dose1,exp((S6ORspline$true.beta1+S6ORspline$BayesN1bias)*new.dose1+(S6ORspline$true.beta2+S6ORspline$BayesN2bias)*new.dose2)
      ,col=3,lwd=3) # bayes normal
lines(new.dose1,exp((S6ORspline$true.beta1+S6ORspline$Freq1bias)*new.dose1+(S6ORspline$true.beta2+S6ORspline$Freq2bias)*new.dose2)
      ,col=4,lwd=3) # freq

# S7
plot(new.dose1,exp(S7ORspline$true.beta1*new.dose1+S7ORspline$true.beta2*new.dose2),col=1,type='l',ylim = c(0.5,8)
     ,las=1,ylab='OR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) # true
lines(new.dose1,exp((S7ORspline$true.beta1+S7ORspline$BayesB1bias)*new.dose1+(S7ORspline$true.beta2+S7ORspline$BayesB2bias)*new.dose2)
      ,col=2,lwd=3) # bayes binomial
lines(new.dose1,exp((S7ORspline$true.beta1+S7ORspline$BayesN1bias)*new.dose1+(S7ORspline$true.beta2+S7ORspline$BayesN2bias)*new.dose2)
      ,col=3,lwd=3) # bayes normal
lines(new.dose1,exp((S7ORspline$true.beta1+S7ORspline$Freq1bias)*new.dose1+(S7ORspline$true.beta2+S7ORspline$Freq2bias)*new.dose2)
      ,col=4,lwd=3) # freq




## larger tau


par(mfrow=c(2,3))
## S9
plot(new.dose1,exp(S9ORspline$true.beta1*new.dose1+S9ORspline$true.beta2*new.dose2),col=1,type='l',ylim = c(0.5,2)
     ,las=1,ylab='OR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) #  true
lines(new.dose1,exp((S9ORspline$true.beta1+S9ORspline$BayesB1bias)*new.dose1+(S9ORspline$true.beta2+S9ORspline$BayesB2bias)*new.dose2)
      ,col=2,lwd=3) # bayes binomial
lines(new.dose1,exp((S9ORspline$true.beta1+S9ORspline$BayesN1bias)*new.dose1+(S9ORspline$true.beta2+S9ORspline$BayesN2bias)*new.dose2)
      ,col=3,lwd=3) # bayes normal
lines(new.dose1,exp((S9ORspline$true.beta1+S9ORspline$Freq1bias)*new.dose1+(S9ORspline$true.beta2+S9ORspline$Freq2bias)*new.dose2)
      ,col=4,lwd=3) # freq


# S10
plot(new.dose1,exp(S10ORspline$true.beta1*new.dose1+S10ORspline$true.beta2*new.dose2),col=1,type='l',ylim = c(0.5,2)
     ,las=1,ylab='OR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) # true
lines(new.dose1,exp((S10ORspline$true.beta1+S10ORspline$BayesB1bias)*new.dose1+(S10ORspline$true.beta2+S10ORspline$BayesB2bias)*new.dose2)
      ,col=2,lwd=3) # bayes binomial
lines(new.dose1,exp((S10ORspline$true.beta1+S10ORspline$BayesN1bias)*new.dose1+(S10ORspline$true.beta2+S10ORspline$BayesN2bias)*new.dose2)
      ,col=3,lwd=3) # bayes normal
lines(new.dose1,exp((S10ORspline$true.beta1+S10ORspline$Freq1bias)*new.dose1+(S10ORspline$true.beta2+S10ORspline$Freq2bias)*new.dose2)
      ,col=4,lwd=3) # freq


# S11
plot(new.dose1,exp(S11ORspline$true.beta1*new.dose1+S11ORspline$true.beta2*new.dose2),col=1,type='l',ylim = c(0.5,12)
     ,las=1,ylab='OR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) # true
lines(new.dose1,exp((S11ORspline$true.beta1+S11ORspline$BayesB1bias)*new.dose1+(S11ORspline$true.beta2+S11ORspline$BayesB2bias)*new.dose2)
      ,col=2,lwd=3) # bayes binomial
lines(new.dose1,exp((S11ORspline$true.beta1+S11ORspline$BayesN1bias)*new.dose1+(S11ORspline$true.beta2+S11ORspline$BayesN2bias)*new.dose2)
      ,col=3,lwd=3) # bayes normal
lines(new.dose1,exp((S11ORspline$true.beta1+S11ORspline$Freq1bias)*new.dose1+(S11ORspline$true.beta2+S11ORspline$Freq2bias)*new.dose2)
      ,col=4,lwd=3) # freq
# legend(0,,legend=c('True','binomial Bayes', 'normal Bayes','Freq' ),
#        col=1:4,horiz = T,lty=1, bty='n',xjust = 0,cex = 0.7,lwd=3)


# S12
plot(new.dose1,exp(S12ORspline$true.beta1*new.dose1+S12ORspline$true.beta2*new.dose2),col=1,type='l',ylim = c(0.5,12)
     ,las=1,ylab='OR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) # true
lines(new.dose1,exp((S12ORspline$true.beta1+S12ORspline$BayesB1bias)*new.dose1+(S12ORspline$true.beta2+S12ORspline$BayesB2bias)*new.dose2)
      ,col=2,lwd=3) # bayes binomial
lines(new.dose1,exp((S12ORspline$true.beta1+S12ORspline$BayesN1bias)*new.dose1+(S12ORspline$true.beta2+S12ORspline$BayesN2bias)*new.dose2)
      ,col=3,lwd=3) # bayes normal
lines(new.dose1,exp((S12ORspline$true.beta1+S12ORspline$Freq1bias)*new.dose1+(S12ORspline$true.beta2+S12ORspline$Freq2bias)*new.dose2)
      ,col=4,lwd=3) # freq


# S13
plot(new.dose1,exp(S13ORspline$true.beta1*new.dose1+S13ORspline$true.beta2*new.dose2),col=1,type='l',ylim = c(0.5,12)
     ,las=1,ylab='OR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) # true
lines(new.dose1,exp((S13ORspline$true.beta1+S13ORspline$BayesB1bias)*new.dose1+(S13ORspline$true.beta2+S13ORspline$BayesB2bias)*new.dose2)
      ,col=2,lwd=3) # bayes binomial
lines(new.dose1,exp((S13ORspline$true.beta1+S13ORspline$BayesN1bias)*new.dose1+(S13ORspline$true.beta2+S13ORspline$BayesN2bias)*new.dose2)
      ,col=3,lwd=3) # bayes normal
lines(new.dose1,exp((S13ORspline$true.beta1+S13ORspline$Freq1bias)*new.dose1+(S13ORspline$true.beta2+S13ORspline$Freq2bias)*new.dose2)
      ,col=4,lwd=3) # freq

# S14
plot(new.dose1,exp(S14ORspline$true.beta1*new.dose1+S14ORspline$true.beta2*new.dose2),col=1,type='l',ylim = c(0.5,12)
     ,las=1,ylab='OR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) # true
lines(new.dose1,exp((S14ORspline$true.beta1+S14ORspline$BayesB1bias)*new.dose1+(S14ORspline$true.beta2+S14ORspline$BayesB2bias)*new.dose2)
      ,col=2,lwd=3) # bayes binomial
lines(new.dose1,exp((S14ORspline$true.beta1+S14ORspline$BayesN1bias)*new.dose1+(S14ORspline$true.beta2+S14ORspline$BayesN2bias)*new.dose2)
      ,col=3,lwd=3) # bayes normal
lines(new.dose1,exp((S14ORspline$true.beta1+S14ORspline$Freq1bias)*new.dose1+(S14ORspline$true.beta2+S14ORspline$Freq2bias)*new.dose2)
      ,col=4,lwd=3) # freq


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

