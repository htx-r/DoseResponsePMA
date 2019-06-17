# This function computes the inverse of variance covariance matrix (precision matrix)
 # for the different doses within each study based on multinomial distribution for counts
 # and using delta method for large sample inference (see Riley article and chapter page 11)

invVarcov <- function(cases,controls){
  c_mat <-cbind(cases,controls)
  c_mat_inv <-1/c_mat
  ncomp <- nrow(c_mat)-1
  off_diag <- c_mat_inv[1, 1] + c_mat_inv[1, 2]
  C1 <- matrix(off_diag,nrow=ncomp,ncol=ncomp)
  diag(C1) <- c_mat_inv[1:ncomp, 1]+ c_mat_inv[1:ncomp, 2] + off_diag
  return(solve(C1))
}

#















#precF(cases=antiDep$Responders[antiDep$Study_No==1],
#    controls=antiDep$No_randomised[antiDep$Study_No==1]-antiDep2$Responders[antiDep$Study_No==1])

# varcov <- function(cases,controls){
#   c_mat <-cbind(cases,controls)
#   c_mat_inv <-1/c_mat
#   ncomp <- nrow(c_mat)-1
#   off_diag <- c_mat_inv[1, 1] + c_mat_inv[1, 2]
#   C1 <- matrix(off_diag,nrow=ncomp,ncol=ncomp)
#   diag(C1) <- c_mat_inv[1:ncomp, 1]+ c_mat_inv[1:ncomp, 2] + off_diag
#   return(C1)
# }

#varcov(cases=antiDep$Responders[antiDep$Study_No==1],
 #      controls=antiDep$No_randomised[antiDep$Study_No==1]-antiDep2$Responders[antiDep$Study_No==1])
#mm
