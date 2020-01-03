library(knitr)
rval<-read.csv("~/Google Drive/DoseResponseNMA/DoseResponseNMA/2019-12-11resORspline40sim1000.csv")
rval <- rval[1:8,]
kable(round(rval,4))
rmarkdown::render("analysis.R", "pdf_document")
