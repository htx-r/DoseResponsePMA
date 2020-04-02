# load libraries
library(rms) # for rcs()
library(MASS) # for truehist()
library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponsePMA",force=TRUE)
library(DoseResponseNMA)
library(meta)
source('Functions needed for dosres MA antidep.R')

########################################
#     load data an prepare

# load and exclude single arm studies
mydata <-  read.csv('~/Google Drive/DoseResponseNMA/DoseResponsePMA/DOSEmainanalysis.csv')
antidep=mydata[mydata$exc==F,]
sum(antidep$No_randomised)
#
antidep$studyid <- as.numeric(as.factor(antidep$Study_No))
antidep$nonResponders <- antidep$No_randomised- antidep$Responders

# apply the function above to all studies
logORmat <- sapply(unique(antidep$studyid),function(i) createORreference.fun(antidep$Responders[antidep$studyid==i],antidep$No_randomised[antidep$studyid==i]),simplify = FALSE)
logORmat <- do.call(rbind,logORmat)
antidep$logOR <- c(logORmat[,1])
antidep$selogOR <- c(logORmat[,2])

# restricted cubic spline transformation doses
knots = c(10,20,50)
antidep$dose1 <- as.matrix(rcs(antidep$hayasaka_ddd,knots))[,1]
antidep$dose2 <- as.matrix(rcs(antidep$hayasaka_ddd,knots))[,2]

# transform data into jags format
jagsdataORspline<- makejagsDRmeta(studyid=studyid,logOR,dose1=dose1,dose2=dose2,cases=Responders,noncases=nonResponders,se=selogOR,type=type,data=antidep,splines=T)
#


########## ##### ##### #########################
# 1. Bivariate normal prior for beta1 and beta2 without residual heterogeneity

# binomial
# jagsdataORspline$idmat <- diag(1,2)
# jagsdataORspline$idmati <- matrix(c(0,1,1,0), nrow = 2,ncol = 2)

doseresORsplineBinBiv <- jags.parallel(data = jagsdataORspline,inits=NULL,parameters.to.save = c('beta.pooled','tau','Z','p.drug','p.drug3020','p.drug4030','beta','rho'),model.file = modelBinSplineDRmetaORBiv,
                                       n.chains=3,n.iter = 100000,n.burnin = 30000,DIC=F,n.thin = 5)
cor_beta1_beta2 <- (doseresORsplineBinBiv$BUGSoutput$mean$rho)/doseresORsplineBinBiv$BUGSoutput$mean$tau^2
doseresORsplineBinBiv$BUGSoutput$mean$beta.pooled

# load the results
load('doseresORsplineBinBiv')

# normal
doseresORsplineNorBiv <- jags.parallel(data = jagsdataORspline,inits=NULL,parameters.to.save = c('beta.pooled','beta.pooled','tau','rho'),model.file = modelNorSplineDRmetaBiv,
                                       n.chains=3,n.iter = 100000,n.burnin = 20000,DIC=F,n.thin = 3)
load('doseresORsplineNorBiv')
# freq
doseresORsplineFreq <- dosresmeta(formula=logOR~rcs(hayasaka_ddd,knots), proc="1stage",id=Study_No, type=type,cases=Responders,n=No_randomised,se=selogOR,data=antidep,method = 'reml')

## use ggplot
beta1fOR <- coef(doseresORsplineFreq)[1]
beta2fOR <- coef(doseresORsplineFreq)[2]

beta1nOR <- doseresORsplineNorBiv$BUGSoutput$mean$beta.pooled[1]
beta2nOR <- doseresORsplineNorBiv$BUGSoutput$mean$beta.pooled[2]

beta1bOR <- doseresORsplineBinBiv$BUGSoutput$mean$beta.pooled[1]
beta2bOR <- doseresORsplineBinBiv$BUGSoutput$mean$beta.pooled[2]

# new dose range to plot the results of the three curves
new.dose <- seq(0,80,1)
new.dose1 <- c(rcs(new.dose,knots)[,1])
new.dose2 <- c(rcs(new.dose,knots)[,2])
df2 <- data.frame(new.dose1=new.dose1,y1 =c(exp(beta1fOR*new.dose1+beta2fOR*new.dose2)
                                            ,exp(beta1nOR*new.dose1+beta2nOR*new.dose2),
                                            exp(beta1bOR*new.dose1+beta2bOR*new.dose2)),
                  method=rep(c('binomial Bayesian','normal Bayesian', 'one-stage (freq)'),each=81))

ggplot(data = df2,aes(x=new.dose1,y=y1,color=method)) +
  geom_line(size=1.3)+ # c('lightblue2','steelblue','darkred')+
  scale_color_manual(values=c('lightblue2','steelblue','darkred'))+
  # geom_line(aes(y = y3),color='lightblue2' , size=1.3)+
  # geom_line(aes(y = y2),color='steelblue' , size=1.3)+
  # geom_line(aes(y = y1),color='darkred' , size=1.3)+
  xlab('')+
  ylab('')+
  ylim(0.5,2)+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        legend.position = c(0.72,0.95), legend.key.size = unit(3, "cm"), legend.text = element_text(size=12),
        legend.key.height  = unit(0.5,"cm"),legend.title = element_blank(), legend.background = element_blank(),
        axis.text.x = element_text(face='bold',size=14),
        axis.text.y = element_text(face='bold',size=14))

p1+scale_colour_manual(name="legend",labels=c('binomial Bayesian','normal Bayesian', 'one-stage (freq)'), values=c('lightblue2', "steelblue",'darkred'))

#** Figure
########## ##### ##### #########################
# 2. dosres MA model with residual heterogeneity

# compute sigma matrix per study
nd    <- as.numeric(table(antidep$studyid)) ## number of all doses (with zero dose)
max.nd <- max(nd)                        ## maximum number of doses
ns <- length(unique(antidep$studyid))       ## number of studies
tncomp <- sum(as.numeric(table(antidep$studyid))-1) ## total number of non-zero comparisons
sigma <- sapply(1:ns, sigmaMat)

# and residual het. matrix per study
resMat <- sapply(1:ns, function(i) diag(1,nrow = nd[i]-1)+(1-diag(1,nrow = nd[i]-1))*0.5)

#
sigmamat <- matrix(NA,tncomp,max.nd-1)
s <- matrix(NA, ns,max.nd-1)
b <- no.d <-vector()
index <- 1:tncomp
for (i in 1:ns) {
  b[1] <- 0
  no.d[i] <- as.numeric(table(antidep$studyid)[i])-1
  sigmamat[(b[i]+1):(b[i]+no.d[i]),1:(no.d[i])] <- sigma[[i]]
  s[i,1:no.d[i]] <- index[(b[i]+1):(b[i]+no.d[i])]
  b[i+1] <- b[i]+ no.d[i]
  sigmamat
}

resmat <- matrix(NA,tncomp,max.nd-1)
s <- matrix(NA, ns,max.nd-1)
b <- no.d <-vector()
index <- 1:tncomp
for (i in 1:ns) {
  b[1] <- 0
  no.d[i] <- as.numeric(table(antidep$studyid)[i])-1
  resmat[(b[i]+1):(b[i]+no.d[i]),1:(no.d[i])] <- resMat[[i]]
  s[i,1:no.d[i]] <- index[(b[i]+1):(b[i]+no.d[i])]
  b[i+1] <- b[i]+ no.d[i]
  resmat
}

jagsdataORspline$sigmamat <- sigmamat
jagsdataORspline$resmat <- resmat


doseresORsplineBin2 <- jags.parallel(data = jagsdataORspline,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau','tau.res'),model.file = modelBinSplineDRmetaORwithRes,
                                    n.chains=2,n.iter = 100,n.burnin = 20,DIC=F,n.thin = 1)

doseresORsplineBin2$BUGSoutput$mean$tau




