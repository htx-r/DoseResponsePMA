# This function computes the inverse of variance covariance matrix (precision matrix)
# for the different doses within each study based on multinomial distribution for counts
# and using delta method for large sample inference (see Orsini article and on his chapter, see page 11)

# The arguments:
# cases:a vector of the number of people that have the outcome
# controls: a vector of the number of people that don't have the outcome

logRRprecmatix <- function(cases,casesRef){
  precmatrix <- matrix(1/casesRef , nrow = length(cases),ncol = length(cases))
  diag(precmatrix) <- 1/casesRef+ 1/cases
  return(solve(precmatrix))
}
#















