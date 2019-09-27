## Explore

library(Hmisc)
rcspline.eval (seq(0,1, .01),
                   knots =seq(.05 ,.95 ,length =3), inclx =T,fractied = 2)
library(rms)
rcs(seq(0,1, .01),
    knots =seq(.05 ,.95 ,length =3), inclx =T)


antiDep <-  read.csv('~/Desktop/TasnimPhD/DoseResponseNMA/DoseResponseNMA/DOSEmainanalysis.csv')
#antiDep <-  read.csv("~/Google Drive/_mydrive/HTx/HTx-R/DoseResponseNMA/DOSEmainanalysis.csv")
NAstudyid <- antiDep$Study_No[is.na(antiDep$logRR)]
antiDep1 <-antiDep[!antiDep$Study_No %in% NAstudyid,]

# add noncases to the data and studyid
antiDep1$nonResponders <- antiDep1$No_randomised- antiDep1$Responders
antiDep1$studyid <- as.numeric(as.factor(antiDep1$Study_No))






options(knots=4, poly.degree=2)
# To get the old behavior of rcspline.eval knot placement (which didnt' handle
# clumping at the lowest or highest value of the predictor very well):
# options(fractied = 1.0)   # see rcspline.eval for details
country <- factor(country.codes)
blood.pressure <- cbind(sbp=systolic.bp, dbp=diastolic.bp)
fit <- lrm(Y ~ sqrt(x1)*rcs(x2) + rcs(x3,c(5,10,15)) +
             lsp(x4,c(10,20)) + country + blood.pressure + poly(age,2))

