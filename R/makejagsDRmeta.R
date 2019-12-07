
makejagsDRmeta <- function(studyid, y,dose1,dose2,cases,noncases,se,type,data,splines=F){
  # This function represent the dataset format so it can be input in all dose-response meta-analysis JAGS model
  # different settings are covered in this function: linear, spline, OR or RR

  # The arguments: all of them are vectors that are assumed to be elements in a dataframe ''data''
  # studyid: study id or number
  # y: a dose-specific outomce as log RR = log p_nonzero/p_zero or log OR
  # dose1: the doses in all studies
  # dose2: the spline transformed doses in all studies
  # cases: the dose-specific number of events
  # noncases: the dose-specific number of non events
  # se: standard error of y: log RR or log OR
  # type: (characterstic) either 'cc' when y= logOR or 'ci' when y= logRR
  # data: a dataframe that contains all the above arguments
  # splines: logical (T/F) to indicate whether we want a jags data for splines, or linear

  library(dosresmeta)

  # evaluate the arguments in the data
  data$studyid <-  eval(substitute(studyid), data)
  data$y <- eval(substitute(y), data)
  data$v <- eval(substitute(se), data)^2
  data$type <-eval(substitute(type), data)
  data$dose1 <- eval(substitute(dose1), data)
  data$dose2 <- eval(substitute(dose2), data)
  data$cases <- eval(substitute(cases), data)
  data$noncases <- eval(substitute(noncases), data)
  data$n <- data$cases + data$noncases


  # additional needed quantities
  study_id <- unique(data$studyid)         ## a vector of the study id
  nd    <- as.numeric(table(data$studyid)) ## number of all doses (with zero dose)
  max.nd <- max(nd)                        ## maximum number of doses
  ns <- length(unique(data$studyid))       ## number of studies
  tncomp <- sum(as.numeric(table(data$studyid))-1) ## total number of non-zero comparisons

  ######################################################################  ######################################################################
  # convert the information from vectors format into matrices with study-specific values in columns and dose levels in rows
  ######################################################################  ######################################################################

  ## Matrix for the number of cases/events 'r' where each row refers to study and the columns refers to the dose levels.
  rmat <- matrix(NA,ns,max.nd)
  for (i in 1:ns) {
    rmat[i,1:as.numeric(table(data$studyid)[i])] <- data$cases[data$studyid == study_id[i]]
  }

  ## Matrix for the sample size 'n' where each row refers to study and the columns refers to the dose levels.

  nmat <- matrix(NA,ns,max.nd)
  for (i in 1:ns) {
    nmat[i,1:as.numeric(table(data$studyid)[i])] <- data$n[data$studyid == study_id[i]]
  }

  ## Matrix for the effects log RR or log OR where each row refers to study and the columns refers to the dose levels.
  Ymat <- matrix(NA,ns,max.nd)
  for (i in 1:ns) {
    Ymat[i,1:as.numeric(table(data$studyid)[i])] <- data$y[data$studyid == study_id[i]]
  }

  ## Matrix for the doses  where each row refers to a study and the columns refers to the dose levels.
  Xmat <- matrix(NA,ns,max.nd)
  for (i in 1:ns) {
    Xmat[i,1:as.numeric(table(data$studyid)[i])] <- data$dose1[data$studyid == study_id[i]]
  }

  ## Find the inverse of the variance-covariance matrix of the dose-specific effects y within each study;
    # so within each study we get a matrix (nd-1)x(nd-1) then we bind all the matrices into one big matrix so we can input it in jags model

  Slist <- sapply(unique(data$studyid), function(i) covar.logrr( cases = cases, n = cases+noncases, y=y,v=v,type = type,data = data[data$studyid==i,])
                  ,simplify = F)

  precmat <- matrix(NA,tncomp,max.nd-1)
  s <- matrix(NA, ns,max.nd-1)
  b <- no.d <-vector()
  index <- 1:tncomp
  for (i in 1:ns) {
    b[1] <- 0
    no.d[i] <- as.numeric(table(data$studyid)[i])-1
    precmat[(b[i]+1):(b[i]+no.d[i]),1:(no.d[i])] <- Slist[[i]]
   s[i,1:no.d[i]] <- index[(b[i]+1):(b[i]+no.d[i])]
    b[i+1] <- b[i]+ no.d[i]
    precmat
  }



  ######################################################################
  #  The final JAGSdataset format
  ######################################################################

  if (splines==TRUE) { # for splines
    ## add the spline dose transformation to the list
    X2mat <- matrix(NA,ns,max.nd)
    for (i in 1:ns) {
      X2mat[i,1:as.numeric(table(data$studyid)[i])] <- data$dose2[data$studyid == study_id[i]]
    }
    JAGSdata <- list(Y=Ymat[,-1],r=rmat,n=nmat,X1=Xmat,X2=X2mat,nd=nd,ns=ns,prec=precmat,s=s) # X3=X3mat[,-1], X3ref=X3mat[,1],

  }else { # for linear
    JAGSdata<- list(Y=Ymat[,-1],r=rmat,n=nmat,X=Xmat,nd=nd,ns=ns,prec=precmat,s=s)
  }

  return(JAGSdata)

}

# End






























######################################################################
##%% 2. For Restricted Cubic Splines model
######################################################################
# if(splines){
#
#   # Construct from each column a vector and add it to the dataset
#   data$X1 <- data$dose1#as.vector(t[,1])
#   data$X2 <- data$dose2#as.vector(t[,2])
#
#   #
#   X1mat <- matrix(NA,ns,max.nd)
#   for (i in 1:ns) {
#     X1mat[i,1:as.numeric(table(data$studyid)[i])] <- data$X1[data$studyid == study_id[i]]
#   }
#
#   #
#   X2mat <- matrix(NA,ns,max.nd)
#   for (i in 1:ns) {
#     X2mat[i,1:as.numeric(table(data$studyid)[i])] <- data$X2[data$studyid == study_id[i]]
#   }
#
# }


# b <- 1
# ndose <- nd-1
# index <- 1:tncom
# maxndose <- max.nd-1
# s <- matrix(NA, ns,maxndose)
# for (i in 1:ns) {
#   s[i,1:ndose[i]] <- index[(b[i]):(b[i]+ndose[i]-1)]
#   b[i+1] <- b[i]+ ndose[i]
# }
