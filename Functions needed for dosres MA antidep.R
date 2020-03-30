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


i <- 1
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





