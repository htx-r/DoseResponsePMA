# This function computes the inverse of variance covariance matrix (precision matrix)
# for the different doses within each study based on multinomial distribution for counts
# and using delta method for large sample inference (see Orsini article and on his chapter, see page 11)

# The arguments:
# cases:a vector of the number of people that have the outcome
# controls: a vector of the number of people that don't have the outcome

logRRprecmatix <- function(cases,casesRef){
  # c_mat <-cbind(cases,controls)
  # c_mat_inv <-1/c_mat
  # ncomp <- nrow(c_mat)-1
  # off_diag <- c_mat_inv[1, 1] + c_mat_inv[1, 2]
  # C1 <- matrix(off_diag,nrow=ncomp,ncol=ncomp)
  # diag(C1) <- c_mat_inv[1:ncomp, 1]+ c_mat_inv[1:ncomp, 2] + off_diag

  precmatrix <- matrix(1/casesRef, nrow = length(cases),ncol = length(cases))
  diag(precmatrix) <- 1/casesRef+ 1/cases
  return(precmatrix)
}
#















