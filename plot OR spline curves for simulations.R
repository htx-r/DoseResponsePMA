## In this file I plot the 4 spline curves with beta1 and beta2 esimated in one of the three approaches
    # 1. true 2. binomial Bayes 3. normal Bayes 4. one-stage
library(rms)
# load the results as csv file
rval <- read.csv("~/Google Drive/DoseResponseNMA/DoseResponsePMA/2019-12-11resORspline40sim1000.csv")
# rval <- rval[1:8,]
# round(rval[,c('true.beta1',   'BayesB1bias',   'BayesN1bias',   'Freq1bias','BayesB1mse',   'BayesN1mse',     'Freq1mse',
#          'true.beta2',   'BayesB2bias',   'BayesN2bias',   'Freq2bias','BayesB2mse',   'BayesN2mse',     'Freq2mse')],4)

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

# lets try to use ggplot .....
library(ggplot2)
library(viridis)
library(gridExtra)
# S2
df1 <- data.frame(new.dose1=new.dose1,y1 =cbind(y[,1,1]),y2=cbind(y[,2,1]),
                      y3=cbind(y[,3,1]),y4=cbind(y[,4,1]))
df1 <- data.frame(new.dose1=new.dose1,y1 =c(y[,1,1],y[,2,1],
                  y[,3,1],y[,4,1]),method=rep(c('true values','binomial Bayesian','normal Bayesian', 'one-stage (freq)'),each=11))

theme_set(
  theme_minimal() +
    theme(legend.position = "right",axis.text.x = element_text(size=12),axis.text.y = element_text(size=12))
)
ggplot(data = df1,aes(x=new.dose1,y=y1,color=method)) +
  # geom_line(aes(y =y1),color='black' , size=1.3) +
  # geom_line(aes(y = y2),color='lightblue2' , size=1.3)+
  # geom_line(aes(y = y3),color='steelblue' , size=1.3)+
  # geom_line(aes(y = y4),color='darkred' , size=1.3)+
  geom_line(size=1.3)+
  xlab('')+
  ylab('')+
  ylim(0.5,3.5)+
  scale_color_manual(values=c('lightblue2','steelblue','darkred','black'))+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        legend.position = c(0.72,0.90), legend.key.size = unit(3, "cm"), legend.text = element_text(size=12),
        legend.key.height  = unit(0.5,"cm"),legend.title = element_blank(), legend.background = element_blank(),
       axis.text.x = element_text(face='bold',size=14),axis.text.y = element_text(face='bold',size=14))

ggplot(data = df1[9:11,],aes(x=new.dose1)) +
  geom_line(aes(y =y1),color='black' , size=1.3) +
  geom_line(aes(y = y2),color='lightblue2' , size=1.3)+
  geom_line(aes(y = y3),color='steelblue' , size=1.3)+
  geom_line(aes(y = y4),color='darkred' , size=1.3)+
  xlab('')+
  ylab('')+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        legend.position = "right",axis.text.x = element_text(face='bold',size=14),axis.text.y =element_text(face='bold',size=14))

# S3
df2 <- data.frame(new.dose1=new.dose1,y1 =cbind(y[,1,2]),y2=cbind(y[,2,2]),
                  y3=cbind(y[,3,2]),y4=cbind(y[,4,2]))


ggplot(data = df2,aes(x=new.dose1)) +
  geom_line(aes(y =y1),color='black' , size=1.3) +
  geom_line(aes(y = y2),color='lightblue2' , size=1.3)+
  geom_line(aes(y = y3),color='steelblue' , size=1.3)+
  geom_line(aes(y = y4),color='darkred' , size=1.3)+
  xlab('')+
  ylab('')+
  ylim(0.5,3.5)+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        legend.position = "right",axis.text.x = element_text(face='bold',size=14),axis.text.y = element_text(face='bold',size=14))

ggplot(data = df2[9:11,],aes(x=new.dose1)) +
  geom_line(aes(y =y1),color='black' , size=1.3) +
  geom_line(aes(y = y2),color='lightblue2' , size=1.3)+
  geom_line(aes(y = y3),color='steelblue' , size=1.3)+
  geom_line(aes(y = y4),color='darkred' , size=1.3)+
  xlab('')+
  ylab('')+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        legend.position = "right",axis.text.x = element_text(face='bold',size=14),axis.text.y =element_text(face='bold',size=14))


# S4
df3 <- data.frame(new.dose1=new.dose1,y1 =cbind(y[,1,3]),y2=cbind(y[,2,3]),
                  y3=cbind(y[,3,3]),y4=cbind(y[,4,3]))


ggplot(data = df3,aes(x=new.dose1)) +
  geom_line(aes(y =y1),color='black' , size=1.3) +
  geom_line(aes(y = y2),color='lightblue2' , size=1.3)+
  geom_line(aes(y = y3),color='steelblue' , size=1.3)+
  geom_line(aes(y = y4),color='darkred' , size=1.3)+
  xlab('')+
  ylab('')+
  ylim(0.5,3.5)+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        legend.position = "right",axis.text.x = element_text(face='bold',size=14),axis.text.y =element_text(face='bold',size=14))

ggplot(data = df3[9:11,],aes(x=new.dose1)) +
  geom_line(aes(y =y1),color='black' , size=1.3) +
  geom_line(aes(y = y2),color='lightblue2' , size=1.3)+
  geom_line(aes(y = y3),color='steelblue' , size=1.3)+
  geom_line(aes(y = y4),color='darkred' , size=1.3)+
  xlab('')+
  ylab('')+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        legend.position = "right",axis.text.x = element_text(face='bold',size=14),axis.text.y =element_text(face='bold',size=14))


#theme(panel.background = element_rect(fill = 'snow')) # beige lavenderblush
#"#00AFBB", "#E7B800", "#FC4E07"
# legend(-2,9,legend=c('true','binomial Bayes', 'normal Bayes','Freq' ),
#        col=1:4,horiz = F,lty=1, bty='n',xjust = 0.1,cex = 0.5,lwd=3,xpd = 'NA')

## investigate Rhat distribution for binomial and normal.
load('resORspline40sim1000ALL')
dim(resORspline40sim1000ALL)
rownames(resORspline40sim1000ALL)

t(S4ORspline[c('RhatB1','RhatB2','RhatN1','RhatN2'),])
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



####################################################################################
####################################################################################
##########  Plot the true DR curve with its confidence interval ####################
####################################################################################
####################################################################################

beta1.true<- c(0.04,0.1,0.2)
beta2.true <- c(0,0.03,-0.2 )
tau <- c(0.001,0.01)

## larger tau
par(mfrow=c(1,1),las=1)

# S2
plot(new.dose1,exp(beta1.true[1]*new.dose1+beta2.true[1]*new.dose2),col=1,type='l',ylim = c(0.5,5)
     ,las=1,ylab='OR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) #  true
lines(new.dose1,exp((beta1.true[1]+2*tau[1])*new.dose1+(beta2.true[1]+2*tau[1])*new.dose2),lty=2,col=2,lwd=3)
lines(new.dose1,exp((beta1.true[1]-2*tau[1])*new.dose1+(beta2.true[1]-2*tau[1])*new.dose2),lty=2,col=2,lwd=3)

lines(new.dose1,exp((beta1.true[1]+2*tau[2])*new.dose1+(beta2.true[1]+2*tau[2])*new.dose2),lty=2,col=3,lwd=3)
lines(new.dose1,exp((beta1.true[1]-2*tau[2])*new.dose1+(beta2.true[1]-2*tau[2])*new.dose2),lty=2,col=3,lwd=3)

# S3
plot(new.dose1,exp(beta1.true[2]*new.dose1+beta2.true[2]*new.dose2),col=1,type='l',ylim = c(0.5,5)
     ,las=1,ylab='OR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) #  true
lines(new.dose1,exp((beta1.true[2]+2*tau[1])*new.dose1+(beta2.true[2]+2*tau[1])*new.dose2),lty=2,col=2,lwd=3)
lines(new.dose1,exp((beta1.true[2]-2*tau[1])*new.dose1+(beta2.true[2]-2*tau[1])*new.dose2),lty=2,col=2,lwd=3)

lines(new.dose1,exp((beta1.true[2]+2*tau[2])*new.dose1+(beta2.true[2]+2*tau[2])*new.dose2),lty=2,col=3,lwd=3)
lines(new.dose1,exp((beta1.true[2]-2*tau[2])*new.dose1+(beta2.true[2]-2*tau[2])*new.dose2),lty=2,col=3,lwd=3)

# S4
plot(new.dose1,exp(beta1.true[3]*new.dose1+beta2.true[3]*new.dose2),col=1,type='l',ylim = c(0.5,5)
     ,las=1,ylab='OR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) #  true
lines(new.dose1,exp((beta1.true[3]+2*tau[1])*new.dose1+(beta2.true[3]+2*tau[1])*new.dose2),lty=2,col=2,lwd=3)
lines(new.dose1,exp((beta1.true[3]-2*tau[1])*new.dose1+(beta2.true[3]-2*tau[1])*new.dose2),lty=2,col=2,lwd=3)

lines(new.dose1,exp((beta1.true[3]+2*tau[2])*new.dose1+(beta2.true[3]+2*tau[2])*new.dose2),lty=2,col=3,lwd=3)
lines(new.dose1,exp((beta1.true[3]-2*tau[2])*new.dose1+(beta2.true[3]-2*tau[2])*new.dose2),lty=2,col=3,lwd=3)

############################################################################
# Histogram of both betas in normal and binomial models over all scenarios
####
rval <- read.csv("~/Google Drive/DoseResponseNMA/DoseResponsePMA/2019-12-13resORspline40sim1000ALL.csv",header = TRUE)
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
















# dim(rval)
# rvalbeta1 <- unlist(rval[1,-1])
# rvalsd1 <- unlist(rval[17,-1])
#
# rvalbeta2 <- unlist(rval[12,-1])
#
# truehist(rvalbeta1,prob = TRUE,col='yellow')
# range(rvalsd1)
# truehist(rvalsd1,prob = TRUE,h=0.002)
#
# x=seq(-0.05,0.05,l=1000)
# y=dnorm(x,mean = mean(rvalbeta1),sd=mean(rvalsd1))
# lines(x,y)
# rval[rval$X=='BayesB1',][1,]
# rval[1,]
#
# myds <- data.frame(beta1=rvalbeta1,beta2=rvalbeta2)
# mydss <- myds[order(myds$beta1),]
#
# x <- rvalbeta1[order(rvalbeta1)]
# y <- rvalbeta2[order(rvalbeta2)]
# z <- as.matrix(expand.grid(x,y))
# contour(x,y,z)




 # #
 # par(mar=c(3.2,3.2,3.2,3.2),las=1)
 #
 # ##
 # matplot(new.dose1,cbind(y[,,1]),type='l'#,main=expression(paste( beta, "1 = 0.04, ", beta, "2= 0"))
 #         ,lwd=3,col=col,lty=1,ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4,ylim=ylim)
 # matplot(new.dose1[9:11],cbind(y[9:11,,1]),type='l'#,main=expression(paste( beta, "1 = 0.04, ", beta, "2= 0"))
 #         ,lwd=3,col=col,lty=1,ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4)
 #
 #
 # matplot(new.dose1,cbind(y[,,2]),type='l'#,main=expression(paste( beta, "1 = 0.1, ", beta, "2= 0.03"))
 #         ,lwd=3,col=col,lty=1,ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4,ylim=ylim)
 # matplot(new.dose1[9:11],cbind(y[9:11,,2]),type='l'#,main=expression(paste( beta, "1 = 0.04, ", beta, "2= 0"))
 #         ,lwd=3,col=col,lty=1,ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4)
 #
 #
 # matplot(new.dose1,cbind(y[,,3]),type='l'#,main=expression(paste( beta, "1 = 0.2, ", beta, "2= -0.2"))
 #         ,lwd=3,col=col,lty=lty,ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4,ylim=ylim)
 # matplot(new.dose1[9:11],cbind(y[9:11,,3]),type='l'#,main=expression(paste( beta, "1 = 0.04, ", beta, "2= 0"))
 #         ,lwd=3,col=col,lty=1,ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4)
 # # legend(-2,9,legend=c('true','binomial Bayes', 'normal Bayes','Freq' ),
 # #        col=1:4,horiz = F,lty=1, bty='n',xjust = 0.1,cex = 0.5,lwd=3,xpd = 'NA')
 #
 # matplot(new.dose1,cbind(y[,,4]),type='l'#,main=expression(paste( beta, "1 = 0.04, ", beta, "2= 0"))
 #         ,lwd=3,col=col,lty=lty,ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4,ylim=ylim)
 # matplot(new.dose1,cbind(y[,,5]),type='l'#,main=expression(paste( beta, "1 = 0.1, ", beta, "2= 0.03"))
 #         ,lwd=3,col=col,lty=lty,ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4,ylim=ylim)
 # matplot(new.dose1,cbind(y[,,6]),type='l'#,main=expression(paste( beta, "1 = 0.2, ", beta, "2= -0.2"))
 #         ,lwd=3,col=col,lty=lty,ylab=ylab,xlab=xlab,cex.axis=1.4,cex.lab=1.4,ylim=ylim)
 #



