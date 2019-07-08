#MultiArgSimulateDRlinear <- function(nrep=3,beta.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Linear
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta.pooled <- c(0,0.05,0.1)
tau <- c(0.0001,0.1)
ns <- 20
doserange <- c(0,10)
samplesize<- 200
nrep=3
res <- mapply(MultiRunSimulateDRlinear, beta.pooled=beta.pooled,tau=tau,MoreArgs=list(ns=ns,doserange=doserange,samplesize=samplesize,nrep=nrep),SIMPLIFY =F)
ret.obj <- do.call(rbind,res)
rownames(ret.obj) <- paste0('Scenario ',1:nrow(ret.obj))
ret.obj
#}
## Run all scenarios

