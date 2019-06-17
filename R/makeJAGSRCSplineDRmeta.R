makeJAGSRCSplineDRmeta <- function(mydata){
library(rms)
knots <- c(10,20,30)
t <-rcs(mydata$hayasaka_ddd,knots=knots)

# Construct from each column a vector and add it to the dataset
mydata1 <- mydata
mydata1$X1 <- as.vector(t[,1])
mydata1$X2 <- as.vector(t[,2])
mydata1$X3 <- as.vector(t[,3])#### maybe skip this????/@@@@@@@@@@@@@@@@@@@
max.nd <- max(as.numeric(table(mydata1$Study_No)))
ns <- length(unique(mydata$Study_No))        ## number of studies
study_id <- unique(mydata$Study_No)

Ymat <- matrix(NA,ns,max.nd)
for (i in 1:ns) {
  Ymat[i,1:as.numeric(table(mydata$Study_No)[i])] <- mydata$logRR[mydata$Study_No == study_id[i]]
}

pr <- sapply(unique(mydata$Study_No), function(i) invVarcov(cases=mydata$Responders[mydata$Study_No==i],
                                                            controls=mydata$No_randomised[mydata$Study_No==i]-mydata$Responders[mydata$Study_No==i]))
tncomp <- sum(as.numeric(table(mydata$Study_No))-1) ## total number of comparisons


precmat <- matrix(NA,tncomp,max.nd-1)
b <- vector()
nd <- vector()
for (i in 1:ns) {
  b[1] <- 0
  nd[i] <- as.numeric(table(mydata$Study_No)[i])-1
  precmat[(b[i]+1):(b[i]+nd[i]),1:(nd[i])] <- pr[[i]]
  b[i+1] <- b[i]+ nd[i]
  precmat
}

## construct matrix with
X1mat <- matrix(NA,ns,max.nd)
for (i in 1:ns) {
  X1mat[i,1:as.numeric(table(mydata1$Study_No)[i])] <- mydata1$X1[mydata1$Study_No == study_id[i]]
}

X2mat <- matrix(NA,ns,max.nd)
for (i in 1:ns) {
  X2mat[i,1:as.numeric(table(mydata1$Study_No)[i])] <- mydata1$X2[mydata1$Study_No == study_id[i]]
}

X3mat <- matrix(NA,ns,max.nd)
for (i in 1:ns) {
  X3mat[i,1:as.numeric(table(mydata1$Study_No)[i])] <- mydata1$X3[mydata1$Study_No == study_id[i]]
}

# Data
dat_rcs   <- list(Y=Ymat,X1=X1mat[,-1],X1ref=X1mat[,1],X2=X2mat[,-1],X2ref=X2mat[,1],X3=X3mat[,-1],nd=nd,ns=ns,prec=precmat) # X3ref=X3mat[,1],
}
