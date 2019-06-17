## Function to create JAGS data so one can use it easily in the model.


makeJAGSLinearQuadraticDRmeta <- function(mydata){

# data$studyid <-  eval(substitute(studyid), data)
# data$logrr <- eval(substitute(logrr), data)
# data$dose <- eval(substitute(dose), data)
# study_id <- data$studyid


study_id <- unique(mydata$Study_No)

# 9 studies are excluded

#
nd    <- as.numeric(table(mydata$Study_No)) ## number of all doses (with zero dose)
max.nd <- max(nd)                             ## maximum number of doses
ns <- length(unique(mydata$Study_No))        ## number of studies


## 1.Matrix for the effects 'log RR' where each row refers to study and the columns refers to doses.
Ymat <- matrix(NA,ns,max.nd)
for (i in 1:ns) {
  Ymat[i,1:as.numeric(table(mydata$Study_No)[i])] <- mydata$logRR[mydata$Study_No == study_id[i]]
}

## 2.Matrix for the doses 'hayasaka_ddd' where each row refers to a study and the columns refers to doses.

Xmat <- matrix(NA,ns,max.nd)
for (i in 1:ns) {
  Xmat[i,1:as.numeric(table(mydata$Study_No)[i])] <- mydata$hayasaka_ddd[mydata$Study_No == study_id[i]]
}

## 3. Find the variance covariance matrix; then use a multivariate normal distributon


#pr is a list of the within-study variance-covariance matrices between LogRR@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

pr <- sapply(unique(mydata$Study_No), function(i) invVarcov(cases=mydata$Responders[mydata$Study_No==i],
                                                          controls=mydata$No_randomised[mydata$Study_No==i]-mydata$Responders[mydata$Study_No==i]))
tncomp <- sum(as.numeric(table(mydata$Study_No))-1) ## total number of comparisons


precmat <- matrix(NA,tncomp,max.nd-1)
b <- vector()
for (i in 1:ns) {
  b[1] <- 0
  nd[i] <- as.numeric(table(mydata$Study_No)[i])-1
  precmat[(b[i]+1):(b[i]+nd[i]),1:(nd[i])] <- pr[[i]]
  b[i+1] <- b[i]+ nd[i]
  precmat
}


# reference dosages: either 0 or minimum dose
Xref <- ifelse(unique(mydata$Study_No) %in% mydata$Study_No[mydata$hayasaka_ddd == 0] ,0,min(mydata$hayasaka_ddd[which(!mydata$Study_No %in% mydata$Study_No[mydata$hayasaka_ddd ==0])]))

# JAGSdata
dat    <- list(Y=Ymat[,-1],X=Xmat[,-1],Xref=Xref,nd=nd,ns=ns,prec=precmat)
}
