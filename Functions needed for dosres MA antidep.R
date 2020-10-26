# function to compute the relative odds ratio as: odds of non-reference dose / odds of refernce dose using metabin
createORreference.fun=function(r,n)
{

  logOR=c(0)
  selogOR=c(NA)

  for(i in 2:c(length(n)))
  {
    calculate=metabin(r[i],n[i],r[1],n[1],sm="OR")
    logOR=c(logOR,calculate$TE)
    selogOR=c(selogOR,calculate$seTE)

  }
  return(cbind(logOR=logOR,selogOR=selogOR))
}

# create dataset per drug
data_per_drug<-function(data, name_of_drug){
  data2<-data
  data2$count<-1
  data_drug<-data2%>%filter(data2$Drug==name_of_drug)%>%select(Study_No)
  data2 <-data2 %>% filter(Study_No%in%data_drug$Study_No)%>%filter(Drug=='placebo'||Drug==name_of_drug)
  # data_count<-data2 %>% group_by(Study_No, Study_No) %>% summarise(arms=sum(count)) %>% filter(arms>2)
  # data2<-data2 %>% filter(Study_No%in%data_count$Study_No)
  data2
}

sample_per_drug <- function(data,name_of_drug){
  dd <- data_per_drug(antidep,name_of_drug)
  sample_per_drug <- sum(dd %>% filter(dd$Drug==name_of_drug)%>%select(No_randomised))
  return(sample_per_drug)
}












sigmaMat <- function(i){
  m <- t(combn(1:nd[i],2))
  x1mat <- sapply(1:ns, function(i) antidep$dose1[antidep$studyid==i])
  x2mat <- sapply(1:ns, function(i) antidep$dose2[antidep$studyid==i])
  A <- function(j,k) (x1mat[[i]][j]-x1mat[[i]][k])^2 + (x2mat[[i]][j]-x2mat[[i]][k])^2
  Amat1 <- diag( mapply(A, j=m[(1:nd[i]-1),1], k=m[(1:nd[i]-1),2]),nd[i]-1)
  for (k in 1:(nd[i]-2)) {
    Amat1[row(Amat1)>col(Amat1) & col(Amat1)==k] <- Amat1[k,k]
  }

  if(nd[i]>2){
  Amat2 <- diag(0,nd[i]-1)
  for (k in (nd[i]-2):2) {
    Amat2[row(Amat2)>col(Amat2) & row(Amat2)==k] <- Amat1[k,k]
  }

  Amat3 <- diag(0,nd[i]-1)
  Amat3[row(Amat3)>col(Amat3)] <- mapply(A, j=m[nd[i]:nrow(m),1], k=m[nd[i]:nrow(m),2])

  Amat0 <- Amat1+Amat2-Amat3
  Amat0[row(Amat0)<col(Amat0)] <- t(Amat0)[row(Amat0)<col(Amat0)]
  Amat <- Amat0
  }else{
  Amat <- Amat1
  }
  return(Amat)
}


tauMatF <- function(i){
  m <- t(combn(1:nd[i],2))
  x1mat <- sapply(1:ns, function(i) antidep$dose1[antidep$studyid==i])
  x2mat <- sapply(1:ns, function(i) antidep$dose2[antidep$studyid==i])
  A <- function(j,k) (x1mat[[i]][j]-x1mat[[i]][k])^2 + (x2mat[[i]][j]-x2mat[[i]][k])^2
  Amat1 <- diag( mapply(A, j=m[(1:nd[i]-1),1], k=m[(1:nd[i]-1),2]),nd[i]-1)
  for (k in 1:(nd[i]-2)) {
    Amat1[row(Amat1)>col(Amat1) & col(Amat1)==k] <- Amat1[k,k]
  }


 if(nd[i]>2){
    Amat2 <- diag(0,nd[i]-1)
    for (k in (nd[i]-2):2) {
      Amat2[row(Amat2)>col(Amat2) & row(Amat2)==k] <- Amat1[k,k]
    }

    Amat3 <- diag(0,nd[i]-1)
    Amat3[row(Amat3)>col(Amat3)] <- mapply(A, j=m[nd[i]:nrow(m),1], k=m[nd[i]:nrow(m),2])

    Amat0 <- Amat1+Amat2-Amat3
    Amat0[row(Amat0)<col(Amat0)] <- t(Amat0)[row(Amat0)<col(Amat0)]
    Amat0[row(Amat0)==col(Amat0)] <- 2*Amat0[row(Amat0)==col(Amat0)]
    Amat0[row(Amat0)!=col(Amat0)] <- 0.5*Amat0[row(Amat0)!=col(Amat0)]

    Amat <- Amat0
  }else{
    Amat <- Amat1
  }
  return(Amat)
}


rhoMatF <- function(i){
  m <- t(combn(1:nd[i],2))
  x1mat <- sapply(1:ns, function(i) antidep$dose1[antidep$studyid==i])
  x2mat <- sapply(1:ns, function(i) antidep$dose2[antidep$studyid==i])
  A <- function(j,k) (x1mat[[i]][j]-x1mat[[i]][k])*(x2mat[[i]][j]-x2mat[[i]][k])
  Amat1 <- diag( mapply(A, j=m[(1:nd[i]-1),1], k=m[(1:nd[i]-1),2]),nd[i]-1)
  for (k in 1:(nd[i]-2)) {
    Amat1[row(Amat1)>col(Amat1) & col(Amat1)==k] <- Amat1[k,k]
  }

  if(nd[i]>2){
    Amat2 <- diag(0,nd[i]-1)
    for (k in (nd[i]-2):2) {
      Amat2[row(Amat2)>col(Amat2) & row(Amat2)==k] <- Amat1[k,k]
    }

    Amat3 <- diag(0,nd[i]-1)
    Amat3[row(Amat3)>col(Amat3)] <- mapply(A, j=m[nd[i]:nrow(m),1], k=m[nd[i]:nrow(m),2])

    Amat0 <- Amat1+Amat2-Amat3
    Amat0[row(Amat0)<col(Amat0)] <- t(Amat0)[row(Amat0)<col(Amat0)]
    Amat0[row(Amat0)==col(Amat0)] <- 2*Amat0[row(Amat0)==col(Amat0)]
    Amat0[row(Amat0)!=col(Amat0)] <- 0.5*Amat0[row(Amat0)!=col(Amat0)]
    Amat <- Amat0
  }else{
    Amat <- Amat1
  }
  return(Amat)
}


