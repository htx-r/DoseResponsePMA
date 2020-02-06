## In this file I plot the 4 spline curves with beta1 and beta2 esimated in one of the three approaches
# 1. true 2. binomial Bayes 3. normal Bayes 4. one-stage

# load the results as csv file
rvalwo48 <- read.csv("~/Google Drive/DoseResponseNMA/DoseResponsePMA/results used in the published article/RR spline with 40 studies 1000 sim/2020-01-09resRRspline40sim100048.csv")
rval48 <- read.csv("~/Google Drive/DoseResponseNMA/DoseResponsePMA/results used in the published article/RR spline with 40 studies 1000 sim/2020-02-05resRRsplineONLY4840sim1000.csv")
rval <- rbind.data.frame(rvalwo48[1:3,],rval48[1,],rvalwo48[4:6,],rval48[2,])
## compute our X: doses
d <- 0:10
knots<-unlist(round(quantile(d,c(0.25,0.5,0.75))))
new.dose<-rcspline.eval(d,knots,inclx = T)
new.dose1 <- c(new.dose[,1])
new.dose2 <- c(new.dose[,2])
n.d <- length(d)

## compute our Y = OR; in DR model: OR = exp(beta1*new.dose1+beta2*new.dose2)

# our beta1
mybeta1 <- cbind(rval[c('true.beta1')],bayesB.beta1=rval[c('true.beta1')]+rval[c('BayesB1bias')],
                 bayesN.beta1= rval[c('true.beta1')]+rval[c('BayesN1bias')],
                 freq.beta1=rval[c('true.beta1')]+rval[c('Freq1bias')])
colnames(mybeta1) <- c('true.beta1','bayesB.beta1','bayesN.beta1','freq.beta1')
beta1 <- mybeta1[-c(1,5),]

# our beta2
mybeta2 <- cbind(rval[c('true.beta2')],bayesB.beta1=rval[c('true.beta2')]+rval[c('BayesB2bias')],
                 bayesN.beta1= rval[c('true.beta2')]+rval[c('BayesN2bias')],
                 freq.beta1=rval[c('true.beta2')]+rval[c('Freq2bias')])
colnames(mybeta2) <- c('true.beta2','bayesB.beta2','bayesN.beta2','freq.beta2')
beta2 <- mybeta2[-c(1,5),]

# finally, y as an arry with dim: doses x 4 ways (true,bayesB,BayesN,freq) +  different scenarios (8)
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

# unified plot settings
ylim=c(0.5,3.5)
ylab=''
xlab=''
col=c('black','green','darkorchid1','blue')
lty=1
beta1[1:4,1]
beta2[1:4,1]

#
par(mar=c(3.2,3.2,3.2,3.2),las=1)

##
matplot(new.dose1,cbind(y[,,1]),type='l'#,main=expression(paste( beta, "1 = 0.04, ", beta, "2= 0"))
        ,lwd=3,col=col,lty=1,ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4,ylim=ylim)
matplot(new.dose1[9:11],cbind(y[9:11,,1]),type='l'#,main=expression(paste( beta, "1 = 0.04, ", beta, "2= 0"))
        ,lwd=3,col=col,lty=1,ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4)


matplot(new.dose1,cbind(y[,,2]),type='l'#,main=expression(paste( beta, "1 = 0.1, ", beta, "2= 0.03"))
        ,lwd=3,col=col,lty=1,ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4,ylim=ylim)
matplot(new.dose1[9:11],cbind(y[9:11,,2]),type='l'#,main=expression(paste( beta, "1 = 0.04, ", beta, "2= 0"))
        ,lwd=3,col=col,lty=1,ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4)


matplot(new.dose1,cbind(y[,,3]),type='l'#,main=expression(paste( beta, "1 = 0.2, ", beta, "2= -0.2"))
        ,lwd=3,col=col,lty=lty,ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4,ylim=ylim)
matplot(new.dose1[9:11],cbind(y[9:11,,3]),type='l'#,main=expression(paste( beta, "1 = 0.04, ", beta, "2= 0"))
        ,lwd=3,col=col,lty=1,ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4)
# legend(-2,9,legend=c('true','binomial Bayes', 'normal Bayes','Freq' ),
#        col=1:4,horiz = F,lty=1, bty='n',xjust = 0.1,cex = 0.5,lwd=3,xpd = 'NA')

matplot(new.dose1,cbind(y[,,4]),type='l'#,main=expression(paste( beta, "1 = 0.04, ", beta, "2= 0"))
        ,lwd=3,col=col,lty=lty,ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4,ylim=ylim)
matplot(new.dose1,cbind(y[,,5]),type='l'#,main=expression(paste( beta, "1 = 0.1, ", beta, "2= 0.03"))
        ,lwd=3,col=col,lty=lty,ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4,ylim=ylim)
matplot(new.dose1,cbind(y[,,6]),type='l'#,main=expression(paste( beta, "1 = 0.2, ", beta, "2= -0.2"))
        ,lwd=3,col=col,lty=lty,ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4,ylim=ylim)



############################################################################
# Histogram of both betas in normal and binomial models over all scenarios
####
rval <- read.csv("~/Google Drive/DoseResponseNMA/DoseResponsePMA/2019-12-22resRRspline40sim1000ALL.csv",header = TRUE)
rval1 <- rval[1:25,-1]
rval2 <- rval[26:50,-1]
rval3 <- rval[51:75,-1]
rval4 <- rval[76:100,-1]
rval5 <- rval[101:125,-1]
rval6 <- rval[126:150,-1]
rval7 <- rval[151:175,-1]
rval8 <- rval[176:200,-1]

#rval1$X == rval2$X == rval3$X == rval4$X == rval5$X == rval6$X == rval7$X == rval8$X
#rownames(rval1) <- rownames(rval2) <- rownames(rval3) <- rownames(rval4) <- rownames(rval5) <- rownames(rval6) <- rownames(rval7) <- rownames(rval8) <- (rval$X[1:25])

par(mfcol=c(1,1),las=1,cex.axis=1.3)
beta1.true<- c(0,0.04,0.1,0.2)
beta2.true <- c(0,0,0.03,-0.2 )
col2 <- 'limegreen'
lwd=3

######### Binomial
# beta1
col='orchid3'

truehist(unlist(rval1[1,]),col=col,xlab='')
abline(v=beta1.true[1],lwd=lwd,col=col2)

truehist(unlist(rval2[1,]),col=col,xlab = '')
abline(v=beta1.true[2],lwd=lwd,col=col2)

truehist(unlist(rval3[1,]),col=col,xlab = '')
abline(v=beta1.true[3],lwd=lwd,col=col2)

truehist(unlist(rval4[1,]),col=col,xlab = '')
abline(v=beta1.true[4],lwd=lwd,col=col2)

truehist(unlist(rval5[1,]),col=col,xlab = '')
abline(v=beta1.true[1],lwd=lwd,col=col2)

truehist(unlist(rval6[1,]),col=col,xlab = '')
abline(v=beta1.true[2],lwd=lwd,col=col2)

truehist(unlist(rval7[1,]),col=col,xlab = '')
abline(v=beta1.true[3],lwd=lwd,col=col2)

truehist(unlist(rval8[1,]),col=col,xlab = '')
abline(v=beta1.true[4],lwd=lwd,col=col2)


## beta2
col='orchid4'

truehist(unlist(rval1[12,]),col=col,xlab='')
abline(v=beta2.true[1],lwd=lwd,col=col2)

truehist(unlist(rval2[12,]),col=col,xlab = '')
abline(v=beta2.true[2],lwd=lwd,col=col2)

truehist(unlist(rval3[12,]),col=col,xlab = '')
abline(v=beta2.true[3],lwd=lwd,col=col2)

truehist(unlist(rval4[12,]),col=col,xlab = '')
abline(v=beta2.true[4],lwd=lwd,col=col2)

truehist(unlist(rval5[12,]),col=col,xlab = '')
abline(v=beta2.true[1],lwd=lwd,col=col2)

truehist(unlist(rval6[12,]),col=col,xlab = '')
abline(v=beta2.true[2],lwd=lwd,col=col2)

truehist(unlist(rval7[12,]),col=col,xlab = '')
abline(v=beta2.true[3],lwd=lwd,col=col2)

truehist(unlist(rval8[12,]),col=col,xlab = '')
abline(v=beta2.true[4],lwd=lwd,col=col2)






######### Normal
# beta1
par(mfcol=c(1,1),las=1,cex.axis=1.3)


col='orchid3'

truehist(unlist(rval1[2,]),col=col,xlab='')
abline(v=beta1.true[1],lwd=lwd,col=col2)

truehist(unlist(rval2[2,]),col=col,xlab = '')
abline(v=beta1.true[2],lwd=lwd,col=col2)

truehist(unlist(rval3[2,]),col=col,xlab = '')
abline(v=beta1.true[3],lwd=lwd,col=col2)

truehist(unlist(rval4[2,]),col=col,xlab = '')
abline(v=beta1.true[4],lwd=lwd,col=col2)

truehist(unlist(rval5[2,]),col=col,xlab = '')
abline(v=beta1.true[1],lwd=lwd,col=col2)

truehist(unlist(rval6[2,]),col=col,xlab = '')
abline(v=beta1.true[2],lwd=lwd,col=col2)

truehist(unlist(rval7[2,]),col=col,xlab = '')
abline(v=beta1.true[3],lwd=lwd,col=col2)

truehist(unlist(rval8[2,]),col=col,xlab = '',ylim=c(0,25))
abline(v=beta1.true[4],lwd=lwd,col=col2)


## beta2
col='orchid4'

truehist(unlist(rval1[13,]),col=col,xlab='')
abline(v=beta2.true[1],lwd=lwd,col=col2)

truehist(unlist(rval2[13,]),col=col,xlab = '')
abline(v=beta2.true[2],lwd=lwd,col=col2)

truehist(unlist(rval3[13,]),col=col,xlab = '')
abline(v=beta2.true[3],lwd=lwd,col=col2)

truehist(unlist(rval4[13,]),col=col,xlab = '')
abline(v=beta2.true[4],lwd=lwd,col=col2)

truehist(unlist(rval5[13,]),col=col,xlab = '')
abline(v=beta2.true[1],lwd=lwd,col=col2)

truehist(unlist(rval6[13,]),col=col,xlab = '',ylim=c(0,10))
abline(v=beta2.true[2],lwd=lwd,col=col2)

truehist(unlist(rval7[13,]),col=col,xlab = '',ylim=c(0,10))
abline(v=beta2.true[3],lwd=lwd,col=col2)

truehist(unlist(rval8[13,]),col=col,xlab = '')
abline(v=beta2.true[4],lwd=lwd,col=col2)


#


