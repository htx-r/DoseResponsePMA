
makeJAGSDRmeta <- function(studyid, logrr,dose1,dose2,cases,noncases,data,Splines=F,new.dose.range=NULL){
  # This function reconstruct the dataset so it fits the JAGS model based on normal likelihood of the logRR.
  # Depends on the model you want to fit (Linear, Quadratic or Cubic spline) the JAGS dataset differs for
  # Cubic spline compared to Linear and Quadratic.

  # The arguments: all of them assumed to be within a dataframe ''data''
  # studyid: study id or number
  # logrr: a vector of the log RR = log p_nonzero/p_zero for each dose within each study
  # dose: a vector of all doses in all studies
  # cases:a vector of the number of people that have the outcome
  # noncases: a vector of the number of people that don't have the outcome
  # data: a dataframe that contains all the arguments
  # Splines: logical (T/F) to indicate whether we want a jags data for splines, or not (linear or quadratic)
  # knots: if Splines=T then we need to specify the position of knots that represented in the spline function

  library(rms) ## contains rcs function

  #
  data$studyid <-  eval(substitute(studyid), data)
  data$logrr <- eval(substitute(logrr), data)
  data$dose1 <- eval(substitute(dose1), data)
  data$dose2 <- eval(substitute(dose2), data)
  data$cases <- eval(substitute(cases), data)
  data$noncases <- eval(substitute(noncases), data)


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
        Xmat[i,1:as.numeric(table(data$studyid)[i])] <- data$dose1[data$studyid == study_id[i]]
      }

  ## Find the inverse of the variance covariance matrix for the doses within each study

  #pr <- sapply(unique(data$studyid), function(i) logRRprecmatix(cases=data$cases[data$studyid==i&data$dose1!=0],
   #                                                        casesRef =data$cases[data$studyid==i&data$dose1==0]),simplify = F)

  tncomp <- sum(as.numeric(table(data$studyid))-1) ## total number of non-zero comparisons


  precmat <- matrix(NA,tncomp,max.nd-1)
  # b <- vector()
  #
  # for (i in 1:ns) {
  #       b[1] <- 0
  #       nd[i] <- as.numeric(table(data$studyid)[i])-1
  #       precmat[(b[i]+1):(b[i]+nd[i]),1:(nd[i])] <- pr[[i]]
  #       b[i+1] <- b[i]+ nd[i]
  #       precmat
  #     }

######################################################################
##%% 2. For Restricted Cubic Splines model
######################################################################
  #knots <- c(10,20,30)
  if(Splines){
      #t <-rcs(data$dose,knots=knots)

      # Construct from each column a vector and add it to the dataset
      data$X1 <- data$dose1#as.vector(t[,1])
      data$X2 <- data$dose2#as.vector(t[,2])
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
  ## for the predictions
  if(is.null(new.dose.range)){
  new.dose <- min(data$dose1):max(data$dose1)
  }else{
    new.dose <- new.dose.range[1]:new.dose.range[2]
  }

  if (Splines) {
    JAGSdata <- list(Y=Ymat[,-1],X1=X1mat,X2=X2mat,nd=nd,ns=ns,prec=precmat,new.dose=new.dose,new.n=length(new.dose)) # X3=X3mat[,-1], X3ref=X3mat[,1],

  }else {
    JAGSdata<- list(Y=Ymat[,-1],X=Xmat,nd=nd,ns=ns,prec=precmat,new.dose=new.dose,new.n=length(new.dose))
  }

  return(JAGSdata)

}
