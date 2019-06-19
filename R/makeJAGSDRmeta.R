# This function reconstruct the dataset so it fits the JAGS model.
 # Depends on the model you want to fit (Linear, Quadratic or Cubic spline) the JAGS dataset differ for
 # Cubic spline compared to Linear and Quadratic.

makeJAGSDRmeta <- function(studyid, logrr,dose,cases,controls,data,Splines=F,knots){
  library(rms) ## contains rcs function

  #
  data$studyid <-  eval(substitute(studyid), data)
  data$logrr <- eval(substitute(logrr), data)
  data$dose <- eval(substitute(dose), data)
  data$cases <- eval(substitute(cases), data)
  data$controls <- eval(substitute(controls), data)


  #
  study_id <- unique(data$studyid)         ## a vector of the study id
  nd    <- as.numeric(table(data$studyid)) ## number of all doses (with zero dose)
  max.nd <- max(nd)                        ## maximum number of doses
  ns <- length(unique(data$studyid))       ## number of studies

######################################################################
##%% 1. For Linear and Quadratic models
######################################################################

  ## Matrix for the effects 'log RR' where each row refers to study and the columns refers to doses.
  Ymat <- matrix(NA,ns,max.nd)
  for (i in 1:ns) {
        Ymat[i,1:as.numeric(table(data$studyid)[i])] <- data$logrr[data$studyid == study_id[i]]
      }

  ## Matrix for the doses  where each row refers to a study and the columns refers to doses.
  Xmat <- matrix(NA,ns,max.nd)
  for (i in 1:ns) {
        Xmat[i,1:as.numeric(table(data$studyid)[i])] <- data$dose[data$studyid == study_id[i]]
      }

  ## Find the inverse of the variance covariance matrix for the doses within each study

  pr <- sapply(unique(data$studyid), function(i) invVarcov(cases=data$cases[data$studyid==i],
                                                           controls=data$controls[data$studyid==i]))

  tncomp <- sum(as.numeric(table(data$studyid))-1) ## total number of non-zero comparisons


  precmat <- matrix(NA,tncomp,max.nd-1)
  b <- vector()

  for (i in 1:ns) {
        b[1] <- 0
        nd[i] <- as.numeric(table(data$studyid)[i])-1
        precmat[(b[i]+1):(b[i]+nd[i]),1:(nd[i])] <- pr[[i]]
        b[i+1] <- b[i]+ nd[i]
        precmat
      }

######################################################################
##%% 2. For Restricted Cubic Splines model
######################################################################
  #knots <- c(10,20,30)
  if(Splines){
      t <-rcs(data$dose,knots=knots)

      # Construct from each column a vector and add it to the dataset
      data$X1 <- as.vector(t[,1])
      data$X2 <- as.vector(t[,2])
      #data$X3 <- as.vector(t[,3])#### maybe skip this????/@@@@@@@@@@@@@@@@@@@

      X1mat <- matrix(NA,ns,max.nd)
      for (i in 1:ns) {
        X1mat[i,1:as.numeric(table(data$studyid)[i])] <- data$X1[data$studyid == study_id[i]]
      }

      X2mat <- matrix(NA,ns,max.nd)
      for (i in 1:ns) {
        X2mat[i,1:as.numeric(table(data$studyid)[i])] <- data$X2[data$studyid == study_id[i]]
      }

     # X3mat <- matrix(NA,ns,max.nd)
    #  for (i in 1:ns) {
    #    X3mat[i,1:as.numeric(table(data$studyid)[i])] <- data$X3[data$studyid == study_id[i]]
    #  }

  }
  ######################################################################
  # 3. The final JAGSdataset
  ######################################################################
  if (Splines) {
    JAGSdata <- list(Y=Ymat,X1=X1mat[,-1],X1ref=X1mat[,1],X2=X2mat[,-1],X2ref=X2mat[,1],nd=nd,ns=ns,prec=precmat) # X3=X3mat[,-1], X3ref=X3mat[,1],

  }else {
    JAGSdata<- list(Y=Ymat[,-1],X=Xmat[,-1],Xref=Xmat[,1],nd=nd,ns=ns,prec=precmat)
  }

  return(JAGSdata)

}
