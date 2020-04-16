###########################################################
# plots in the appendix: Figure1 , Figure2a and Figure2b


# Figure 1: RR vs dose for simulations S2-S4 with zoom-in plots
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

# plot x vs y = exp(beta1*x1+beta2*x2)
df11 <- data.frame(new.dose1=new.dose1,y1 =c(y[,1,1],y[,2,1],
                                             y[,3,1],y[,4,1]),method=rep(c('true values','binomial Bayesian','normal Bayesian', 'one-stage (freq)'),each=11))
df1 <- data.frame(new.dose1=new.dose1,y1 =cbind(y[,1,1]),y2=cbind(y[,2,1]),
                  y3=cbind(y[,3,1]),y4=cbind(y[,4,1]))
theme_set(
  theme_minimal() +
    theme(legend.position = "right",axis.text.x = element_text(size=12),axis.text.y = element_text(size=12))
)

# S2
ggplot(data = df11,aes(x=new.dose1,y=y1,color=method)) +
  geom_line(size=1.3)+
  xlab('')+
  ylab('')+
  ylim(0.5,3.5)+
  scale_color_manual(values=c('lightblue2','steelblue','darkred','black'))+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        legend.position = c(0.72,0.90), legend.key.size = unit(3, "cm"), legend.text = element_text(size=12),
        legend.key.height  = unit(0.5,"cm"),legend.title = element_blank(), legend.background = element_blank(),
        axis.text.x = element_text(face='bold',size=14),axis.text.y = element_text(face='bold',size=14))

# zoom-in
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

# zoom-in
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

# zoom-in

ggplot(data = df3[9:11,],aes(x=new.dose1)) +
  geom_line(aes(y =y1),color='black' , size=1.3) +
  geom_line(aes(y = y2),color='lightblue2' , size=1.3)+
  geom_line(aes(y = y3),color='steelblue' , size=1.3)+
  geom_line(aes(y = y4),color='darkred' , size=1.3)+
  xlab('')+
  ylab('')+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        legend.position = "right",axis.text.x = element_text(face='bold',size=14),axis.text.y =element_text(face='bold',size=14))

# Figure 2: true dosres curves for all scenarios

beta1.true<- c(0.04,0.1,0.2)
beta2.true <- c(0,0.03,-0.2 )
tau <- c(0.001,0.01)

d <- 0:10
knots<-unlist(round(quantile(d,c(0.25,0.5,0.75))))
new.dose<-rcspline.eval(d,knots,inclx = T)
new.dose1 <- c(new.dose[,1])
new.dose2 <- c(new.dose[,2])
n.d <- length(d)
## larger tau
par(mfrow=c(1,1),las=1)

# S2
f1 <- function(beta1.true,beta2.true) exp(beta1.true*new.dose1+beta2.true*new.dose2)
f2 <- function(beta1.true,beta2.true,tau)exp((beta1.true+2*tau)*new.dose1+(beta2.true+2*tau)*new.dose2)
f3 <- function(beta1.true,beta2.true,tau)exp((beta1.true-2*tau)*new.dose1+(beta2.true-2*tau)*new.dose2)

df1 <- data.frame(new.dose1=new.dose,y1 =f1(beta1.true[1],beta2.true[1]),y2=f2(beta1.true[1],beta2.true[1],tau[1]),
                  y3=f3(beta1.true[1],beta2.true[1],tau[1]),y4=f2(beta1.true[1],beta2.true[1],tau[2]),y5=f3(beta1.true[1],beta2.true[1],tau[2]))

theme_set(
  theme_minimal() +
    theme(legend.position = "right",axis.text.x = element_text(size=12),axis.text.y = element_text(size=12))
)
ggplot(data = df1,aes(x=new.dose1))+
  geom_line(aes(y =y1),color='grey20' , size=1.15) +
  geom_line(aes(y = y2),color='coral3' , size=1.2,linetype = "dashed")+
  geom_line(aes(y = y3),color='coral3' , size=1.15,linetype = "dashed")+
  geom_line(aes(y = y4),color='steelblue' , size=1.2,linetype = "dashed")+
  geom_line(aes(y = y5),color='steelblue' , size=1.15,linetype = "dashed")+
  geom_hline(yintercept =  1,color='grey',size=1.5)+
  xlab('')+
  ylab('')+
  #ggtitle(mytitle)+
  ylim(1,5)+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        axis.text.x = element_text(face='bold',size=14),axis.text.y = element_text(face='bold',size=14),
        axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),
        plot.title = element_text(size = 16, face = "bold"))

# S3
df2 <- data.frame(new.dose1=new.dose,y1 =f1(beta1.true[2],beta2.true[2]),y2=f2(beta1.true[2],beta2.true[2],tau[1]),
                  y3=f3(beta1.true[2],beta2.true[2],tau[1]),y4=f2(beta1.true[2],beta2.true[2],tau[2]),y5=f3(beta1.true[2],beta2.true[2],tau[2]))

theme_set(
  theme_minimal() +
    theme(legend.position = "right",axis.text.x = element_text(size=12),axis.text.y = element_text(size=12))
)
ggplot(data = df2,aes(x=new.dose1))+
  geom_line(aes(y =y1),color='grey20' , size=1.15) +
  geom_line(aes(y = y2),color='coral3' , size=1.2,linetype = "dashed")+
  geom_line(aes(y = y3),color='coral3' , size=1.15,linetype = "dashed")+
  geom_line(aes(y = y4),color='steelblue' , size=1.2,linetype = "dashed")+
  geom_line(aes(y = y5),color='steelblue' , size=1.15,linetype = "dashed")+
  geom_hline(yintercept =  1,color='grey',size=1.5)+
  xlab('')+
  ylab('')+
  #ggtitle(mytitle)+
  ylim(1,5)+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        axis.text.x = element_text(face='bold',size=14),axis.text.y = element_text(face='bold',size=14),
        axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),
        plot.title = element_text(size = 16, face = "bold"))

# S4
df3 <- data.frame(new.dose1=new.dose,y1 =f1(beta1.true[3],beta2.true[3]),y2=f2(beta1.true[3],beta2.true[3],tau[1]),
                  y3=f3(beta1.true[3],beta2.true[3],tau[1]),y4=f2(beta1.true[3],beta2.true[3],tau[2]),y5=f3(beta1.true[3],beta2.true[3],tau[2]))

theme_set(
  theme_minimal() +
    theme(legend.position = "right",axis.text.x = element_text(size=12),axis.text.y = element_text(size=12))
)
ggplot(data = df3,aes(x=new.dose1))+
  geom_line(aes(y =y1),color='grey20' , size=1.15) +
  geom_line(aes(y = y2),color='coral3' , size=1.2,linetype = "dashed")+
  geom_line(aes(y = y3),color='coral3' , size=1.15,linetype = "dashed")+
  geom_line(aes(y = y4),color='steelblue' , size=1.2,linetype = "dashed")+
  geom_line(aes(y = y5),color='steelblue' , size=1.15,linetype = "dashed")+
  geom_hline(yintercept =  1,color='grey',size=1.5)+
  xlab('')+
  ylab('')+
  #ggtitle(mytitle)+
  ylim(1,5)+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        axis.text.x = element_text(face='bold',size=14),axis.text.y = element_text(face='bold',size=14),
        axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),
        plot.title = element_text(size = 16, face = "bold"))



# Figure3: 60 dose-response curves with the averaged curve
load('antidepORsplineFINAL')

# dose
knots = c(10,20,50)
new.dose <- seq(0,80,1)
new.dose1 <- c(rcs(new.dose,knots)[,1])
new.dose2 <- c(rcs(new.dose,knots)[,2])
n.d <- length(new.dose)

# beta1.pooled and beta2.pooled
beta1bOR <- doseresORsplineBin$BUGSoutput$mean$beta1.pooled
beta2bOR <- doseresORsplineBin$BUGSoutput$mean$beta2.pooled

# y
beta1 <- doseresORsplineBin$BUGSoutput$mean$beta1
beta2 <- doseresORsplineBin$BUGSoutput$mean$beta2
y <-matrix(NA,60,n.d)
for (j in 1:60) {
  y[j,1:n.d] <- exp(beta1[j]*new.dose1+beta2[j]*new.dose2)
}

# plot ..
colNo <- c(seq(1,100,2),seq(4,100,4)[1:10])
cc <- paste0('gray',colNo[order(colNo)])  # 60 colors
matplot(new.dose,t(y),col=cc,type='l',lty=1,lwd=2,xlab='Dose',ylab = 'Odds Ratio',
        cex.axis=1.3,las=1,axes=FALSE) # 60 DR curves
axis(1,at=seq(0,80,20))
axis(2,at=seq(0.5,2.5,0.5))
lines(new.dose,exp(beta1bOR*new.dose1+beta2bOR*new.dose2),col='orchid3',lwd=4) # vertical line at the true value
grid(col='gray')
box(col='gray50')


# Histograms for beta1, beta2 and tau:  OR simulations
library(MASS)
rval <- read.csv("~/Google Drive/DoseResponseNMA/DoseResponsePMA/2020-04-13resORspline40sim1000ALL.csv")
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
tau.true <- c(0.001,0.01)
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

## tau
col='orchid1'

truehist(unlist(rval1[8,]),col=col,xlab='')
abline(v=tau.true[1],lwd=lwd,col=col2)

truehist(unlist(rval2[8,]),col=col,xlab = '')
abline(v=tau.true[1],lwd=lwd,col=col2)

truehist(unlist(rval3[8,]),col=col,xlab = '')
abline(v=tau.true[1],lwd=lwd,col=col2)

truehist(unlist(rval4[8,]),col=col,xlab = '')
abline(v=tau.true[1],lwd=lwd,col=col2)

truehist(unlist(rval5[8,]),col=col,xlab = '')
abline(v=tau.true[2],lwd=lwd,col=col2)

truehist(unlist(rval6[8,]),col=col,xlab = '')
abline(v=tau.true[2],lwd=lwd,col=col2)

truehist(unlist(rval7[8,]),col=col,xlab = '')
abline(v=tau.true[2],lwd=lwd,col=col2)

truehist(unlist(rval8[8,]),col=col,xlab = '')
abline(v=tau.true[2],lwd=lwd,col=col2)



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


## tau
col='orchid1'

truehist(unlist(rval1[7,]),col=col,xlab='')
abline(v=tau.true[1],lwd=lwd,col=col2)

truehist(unlist(rval2[7,]),col=col,xlab = '')
abline(v=tau.true[1],lwd=lwd,col=col2)

truehist(unlist(rval3[7,]),col=col,xlab = '')
abline(v=tau.true[1],lwd=lwd,col=col2)

truehist(unlist(rval4[7,]),col=col,xlab = '')
abline(v=tau.true[1],lwd=lwd,col=col2)

truehist(unlist(rval5[7,]),col=col,xlab = '')
abline(v=tau.true[2],lwd=lwd,col=col2)

truehist(unlist(rval6[7,]),col=col,xlab = '')
abline(v=tau.true[2],lwd=lwd,col=col2)

truehist(unlist(rval7[7,]),col=col,xlab = '')
abline(v=tau.true[2],lwd=lwd,col=col2)

truehist(unlist(rval8[7,]),col=col,xlab = '')
abline(v=tau.true[2],lwd=lwd,col=col2)


###### freq(one-stage)
# beta1
par(mfcol=c(1,1),las=1,cex.axis=1.3)


col='orchid3'

truehist(unlist(rval1[3,]),col=col,xlab='')
abline(v=beta1.true[1],lwd=lwd,col=col2)

truehist(unlist(rval2[3,]),col=col,xlab = '')
abline(v=beta1.true[2],lwd=lwd,col=col2)

truehist(unlist(rval3[3,]),col=col,xlab = '')
abline(v=beta1.true[3],lwd=lwd,col=col2)

truehist(unlist(rval4[3,]),col=col,xlab = '')
abline(v=beta1.true[4],lwd=lwd,col=col2)

truehist(unlist(rval5[3,]),col=col,xlab = '')
abline(v=beta1.true[1],lwd=lwd,col=col2)

truehist(unlist(rval6[3,]),col=col,xlab = '')
abline(v=beta1.true[2],lwd=lwd,col=col2)

truehist(unlist(rval7[3,]),col=col,xlab = '')
abline(v=beta1.true[3],lwd=lwd,col=col2)

truehist(unlist(rval8[3,]),col=col,xlab = '')
abline(v=beta1.true[4],lwd=lwd,col=col2)


## beta2
col='orchid4'

truehist(unlist(rval1[14,]),col=col,xlab='')
abline(v=beta2.true[1],lwd=lwd,col=col2)

truehist(unlist(rval2[14,]),col=col,xlab = '')
abline(v=beta2.true[2],lwd=lwd,col=col2)

truehist(unlist(rval3[14,]),col=col,xlab = '')
abline(v=beta2.true[3],lwd=lwd,col=col2)

truehist(unlist(rval4[14,]),col=col,xlab = '')
abline(v=beta2.true[4],lwd=lwd,col=col2)

truehist(unlist(rval5[14,]),col=col,xlab = '')
abline(v=beta2.true[1],lwd=lwd,col=col2)

truehist(unlist(rval6[14,]),col=col,xlab = '')
abline(v=beta2.true[2],lwd=lwd,col=col2)

truehist(unlist(rval7[14,]),col=col,xlab = '')
abline(v=beta2.true[3],lwd=lwd,col=col2)

truehist(unlist(rval8[14,]),col=col,xlab = '')
abline(v=beta2.true[4],lwd=lwd,col=col2)

#

