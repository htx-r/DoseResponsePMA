# Libraries
library(rms)
library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force=TRUE)
library(DoseResponseNMA)

#########################################
# Data
#########################################

antiDep <-  read.csv('~/Desktop/TasnimPhD/DoseResponseNMA/DoseResponseNMA/DOSEmainanalysis.csv')
#antiDep <-  read.csv("~/Google Drive/_mydrive/HTx/HTx-R/DoseResponseNMA/DOSEmainanalysis.csv")
NAstudyid <- antiDep$Study_No[is.na(antiDep$logRR)]
antiDep1 <-antiDep[!antiDep$Study_No %in% NAstudyid,]

# add noncases to the data and studyid
antiDep1$nonResponders <- antiDep1$No_randomised- antiDep1$Responders
antiDep1$studyid <- as.numeric(as.factor(antiDep1$Study_No))

max.nd <- max(as.numeric(table(antiDep1$Study_No)))
ns <- length(unique(antiDep1$Study_No))
nd <- as.numeric(table(antiDep1$Study_No))-1

# Compute the OR and selogOR
ORfun <- function(p0,p1) p1*(1-p0)/(p0*(1-p1))
selogORfun <- function(c0,n0,c1,n1) sqrt(1/c0+1/n0+1/c1+1/n1)
OR <- sapply(unique(antiDep1$studyid), function(i) c(1,ORfun(p0=antiDep1$Responders[antiDep1$studyid==i&antiDep1$hayasaka_ddd==0]/antiDep1$No_randomised[antiDep1$studyid==i&antiDep1$hayasaka_ddd==0],p1=antiDep1$Responders[antiDep1$studyid==i&antiDep1$hayasaka_ddd!=0]/antiDep1$No_randomised[antiDep1$studyid==i&antiDep1$hayasaka_ddd!=0])),simplify = F)
selogOR <- sapply(unique(antiDep1$studyid), function(i) c(NA,selogORfun(c0=antiDep1$Responders[antiDep1$studyid==i&antiDep1$hayasaka_ddd==0],n0=antiDep1$No_randomised[antiDep1$studyid==i&antiDep1$hayasaka_ddd==0],c1=antiDep1$Responders[antiDep1$studyid==i&antiDep1$hayasaka_ddd!=0],n1=antiDep1$No_randomised[antiDep1$studyid==i&antiDep1$hayasaka_ddd!=0])),simplify = F)

# compute OR and SE for study '56' by hand, since it does not have 0 dose, I will do it later more generic
p0 <- 47/111
p1 <- 90/110
OR56 <- ORfun(p0,p1)
OR[[56]] <- c(1, OR56)
selogOR[[56]] <- c(NA, selogORfun(47,11,90,110))

antiDep1$logOR <- log(unlist(OR))
antiDep1$selogOR <- unlist(selogOR)

# dose rcs transformation
dose  <- as.matrix(rcs(antiDep1$hayasaka_ddd,knots=c(10,20,30,50)))
rcspline.eval(antiDep1$hayasaka_ddd,knots=c(10,20,30,50),inclx = T)
antiDep1$dose1 <- dose[,1]
antiDep1$dose2 <- dose[,2]
antiDep1$dose3 <-dose[,3]
study_id <- unique(antiDep1$studyid)


###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Linear: 1. odds ratio (OR) 2. risk ratio (RR)
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrep <- 3
