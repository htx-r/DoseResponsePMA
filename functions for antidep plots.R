data_per_drug<-function(data, name_of_drug){
  # data <- antidep
  # name_of_drug <- 'escitalopram'
   data2<-data
  data2$count<-1
  data_drug<-data2%>%filter(data2$Drug==name_of_drug)#%>%select(Study_No)
  data2 <-data2 %>% filter(Study_No%in%data_drug$Study_No)%>%filter(Drug=='placebo'||Drug==name_of_drug)
  # data_count<-data2 %>% group_by(Study_No, Study_No) %>% summarise(arms=sum(count)) %>% filter(arms>2)
  # data2<-data2 %>% filter(Study_No%in%data_count$Study_No)
  data2
}
dd <-data_per_drug(antidep,'escitalopram')
name_of_drug <- 'escitalopram'
sample_per_drug <- function(data,name_of_drug){
  dd <- data_per_drug(data,name_of_drug)
  sample_per_drug <- sum(dd %>% filter(dd$Drug==name_of_drug)%>%dplyr::select(No_randomised))
return(sample_per_drug)
  }

# Dose: cubic spline transformation
run3modelsDRmeta <- function(drug_name){
  #drug_name <- 'escitalopram'
  dataDrug <- data_per_drug(antidep,drug_name)
  knots = quantile(dataDrug$Dose_delivered_mean,probs=c(0.10,0.50,0.90))
  dataDrug$dose1 <- as.matrix(rcs(dataDrug$Dose_delivered_mean,knots))[,1]
  dataDrug$dose2 <- as.matrix(rcs(dataDrug$Dose_delivered_mean,knots))[,2]

  # data_per_drug(antidep,'citalopram')
  # data_per_drug(antidep,'escitalopram')
  # data_per_drug(antidep,'sertraline')

  #

  #
  dataDrug$studyid <- as.numeric(as.factor(dataDrug$Study_No))
  dataDrug$nonResponders <- dataDrug$No_randomised- dataDrug$Responders

  # Response:  odds ratio
  logORmat <- sapply(unique(dataDrug$studyid),function(i) createORreference.fun(dataDrug$Responders[dataDrug$studyid==i],dataDrug$No_randomised[dataDrug$studyid==i]),simplify = FALSE)
  logORmat <- do.call(rbind,logORmat)
  dataDrug$logOR <- c(logORmat[,1])
  dataDrug$selogOR <- c(logORmat[,2])

  antidep$studyid <- as.numeric(as.factor(antidep$Study_No))
  study_id <- unique(antidep$studyid)         ## a vector of the study id
  max.nd <- max(as.numeric(table(antidep$studyid)))
  ns <- length(unique(antidep$studyid))
  rr <- matrix(NA,ns,max.nd)
  for (i in 1:ns) {
    rr[i,1:as.numeric(table(antidep$studyid)[i])] <- antidep$Responders[antidep$studyid == study_id[i]]
  }

  ## Matrix for the sample size 'n' where each row refers to study and the columns refers to the dose levels.

  nn <- matrix(NA,ns,max.nd)
  for (i in 1:ns) {
    nn[i,1:as.numeric(table(antidep$studyid)[i])] <- antidep$No_randomised[antidep$studyid == study_id[i]]
  }

  # transform data to jags format
  jagsdataORspline<- makejagsDRmeta(studyid=studyid,logOR,dose1=dose1,dose2=dose2,cases=Responders,noncases=nonResponders,se=selogOR,type=type,data=dataDrug,splines=T)

  # additional arguments into jagsdata to compute the absolute response for the placebo and drug arms
  jagsdataORspline$np <- ns#length(unique(dataDrug$Study_No))#58 # sum(jagsdataORspline$X1[,1]==0)
  jagsdataORspline$nn <- nn
  jagsdataORspline$rr <- rr
  jagsdataORspline$new.dose <-1:max(dataDrug$Dose_delivered_mean)
  jagsdataORspline$f.new.dose <- rcspline.eval(jagsdataORspline$new.dose,knots,inclx = T)[,2]
  jagsdataORspline$nd.new <- length(jagsdataORspline$new.dose)


  ## Frequentist: one-stage model using dosresmeta
  doseresORsplineFreq <- dosresmeta(formula=logOR~dose1+dose2, proc="1stage",id=Study_No, type=type,cases=Responders,n=No_randomised,se=selogOR,data=dataDrug,method = 'reml')

  # Bayes with normal likelihood
  doseresORsplineNor <- jags.parallel(data = jagsdataORspline,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelNorSplineDRmeta,
                                      n.chains=3,n.iter = 100000,n.burnin = 40000,DIC=F,n.thin = 1)
  # Bayes with binomial likelihood
  doseresORsplineBin <- jags.parallel(data = jagsdataORspline,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau','Z','p.drug','p.drug3020','p.drug4030'),model.file = modelBinSplineDRmetaORdrug,
                                      n.chains=3,n.iter = 1000,n.burnin = 400,DIC=F,n.thin = 1)

  return(list(doseresORsplineBin=doseresORsplineBin,doseresORsplineNor=doseresORsplineNor, doseresORsplineFreq=doseresORsplineFreq))
}

plotdata.fun <- function(dataDrug,drug_name,doseresORsplineBin){
  # new dose range to plot the results of the three curves
  # dataDrug <- data_per_drug(antidep,'citalopram')
  # drug_name <- 'citalopram'
  # doseresORsplineBin <- doseresORsplineBin
  df <- data.frame(obsDose =dataDrug$Dose_delivered_mean,drug=drug_name)
  new.dose <- 1:max(dataDrug$Dose_delivered_mean)
  # knots = quantile(dataDrug$Dose_delivered_mean,probs=c(0.10,0.50,0.90))
  # new.dose1 <- c(rcs(new.dose,knots)[,1])
  # new.dose2 <- c(rcs(new.dose,knots)[,2])

  # Figure 2a:  absoulte responses vs  dose and add placebo response effect
  #plotdata <- as.data.frame(doseresORsplineBin[["BUGSoutput"]][["summary"]])

  p.placebo<- exp(doseresORsplineBin$BUGSoutput$mean$Z)/(1+exp(doseresORsplineBin$BUGSoutput$mean$Z))
  p.drug <- doseresORsplineBin$BUGSoutput$mean$p.drug
  l.ci <-  doseresORsplineBin$BUGSoutput$summary[paste0('p.drug[',1:max(new.dose),']'),'2.5%']
  u.ci <- doseresORsplineBin$BUGSoutput$summary[paste0('p.drug[',1:max(new.dose),']'),'97.5%']
  lp.ci <-  exp(doseresORsplineBin$BUGSoutput$summary['Z','2.5%'])/(1+exp(doseresORsplineBin$BUGSoutput$summary['Z','2.5%']))
  up.ci <- exp(doseresORsplineBin$BUGSoutput$summary['Z','97.5%'])/(1+exp(doseresORsplineBin$BUGSoutput$summary['Z','97.5%']))
plotdata <- data.frame(p.drug=p.drug,l.ci=l.ci,u.ci=u.ci ,drug=drug_name,p.placebo=p.placebo,lp.ci=lp.ci,up.ci=up.ci,dose=new.dose)
plotdata = list(plotdata=plotdata,df=df)
return(plotdata)
}

myggplot <- function(plotdata){
  df <- plotdata$df
  plotdata <- plotdata$plotdata

  # add observed doses to plotdata
m <- plotdata%>%as.data.frame()%>%group_split(drug)
z <-sapply(1:5, function(i){
  k <- as.data.frame(m[[i]])
  k$obsDose <- rep_len(df[df$drug==unique(as.data.frame(m[[i]])$drug),]$obsDose,nrow(as.data.frame(m[[i]])))
  k
  },simplify = F)
#df %>%dplyr::select()
#unique(rep(df$obsDose['drug'==],each=as.numeric(table(plotdata[,'drug']))[1]))
plotdata <- do.call(rbind,z)

  theme_set(
    theme_minimal()
  )
  g <- ggplot2::ggplot(plotdata, ggplot2::aes(x=dose)) +
    #ggplot2::geom_line(ggplot2::aes(linetype="Posterior Median"))+
    # ggplot2::geom_line(ggplot2::aes(y=`2.5%`, linetype="95% CrI")) +
    # ggplot2::geom_line(ggplot2::aes(y=`97.5%`, linetype="95% CrI"))+
    geom_line(aes(y=p.drug,color='treatment'))+
    geom_line(aes(y=p.placebo,color='placebo'))+
    scale_color_manual(values=c('steelblue','coral4'),name="")+
    guides(color = guide_legend(override.aes = list(size = 1.5)))+
    ggplot2::geom_smooth(ggplot2::aes(x=dose,y=p.drug,ymin=l.ci,ymax=u.ci),color='coral4',fill='coral1',
                         data=plotdata, stat="identity")+
    geom_smooth(aes(x=dose, y=p.placebo, ymin=lp.ci,
                    ymax=up.ci),color='steelblue',fill='steelblue2',
                data=plotdata, stat="identity")+geom_rug(data=plotdata,mapping=aes(x=obsDose),inherit.aes = F,color='red',length = unit(0.05, "npc"))
  #coord_cartesian(ylim = c(0, 0.8))

  #ggplot2::geom_line(ggplot2::aes(y=`97.5%`, linetype="95% CrI"))+
  #geom_line(aes(y=p.placebo,color='coral'))


  g <- g + ggplot2::facet_wrap(~drug,scales = 'free_x')

  g <- g + ggplot2::labs(y="Predicted absolute response", x="Dose")

  g <- g + ggplot2::scale_linetype_manual(name="")

  g <- g+theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
          legend.position = "right",
          legend.key.size = unit(2, "cm"), legend.text = element_text(size=16), # c(0.72,0.90)
          legend.key.height  = unit(1,"cm"),legend.title = element_blank(),
          axis.text.x = element_text(face='bold',size=14),axis.text.y = element_text(face='bold',size=14),
          axis.title.x=element_text(size=16,face = "bold"),axis.title.y=element_text(size=16,face = "bold"),
          # plot.title = element_text(size = 16, face = "bold"),
          strip.background =element_rect(fill="snow3"),
          strip.text.x = element_text(size = 16))

  g
  }

plotdata.fun2 <- function(dataDrug,drug_name,rval){
  df <- data.frame(obsDose =dataDrug$Dose_delivered_mean)

  #
  doseresORsplineFreq <- rval$doseresORsplineFreq
  doseresORsplineNor <- rval$doseresORsplineNor
  doseresORsplineBin <- rval$doseresORsplineBin
  # beta1 and beta2
  beta1fOR <- coef(doseresORsplineFreq)[1]
  beta2fOR <- coef(doseresORsplineFreq)[2]

  beta1nOR <- doseresORsplineNor$BUGSoutput$mean$beta1.pooled
  beta2nOR <- doseresORsplineNor$BUGSoutput$mean$beta2.pooled

  beta1bOR <- doseresORsplineBin$BUGSoutput$mean$beta1.pooled
  beta2bOR <- doseresORsplineBin$BUGSoutput$mean$beta2.pooled

  # new dose range to plot the results of the three curves
  new.dose <- 1:max(dataDrug$Dose_delivered_mean)
  knots = quantile(dataDrug$Dose_delivered_mean,probs=c(0.10,0.50,0.90))
  new.dose1 <- c(rcs(new.dose,knots)[,1])
  new.dose2 <- c(rcs(new.dose,knots)[,2])

  f.new.dose1 <- c(exp(beta1fOR*new.dose1+beta2fOR*new.dose2),exp(beta1nOR*new.dose1+beta2nOR*new.dose2),
                  exp(beta1bOR*new.dose1+beta2bOR*new.dose2))
  method <- rep(c('one-stage (freq)','normal Bayesian','binomial Bayesian' ),each=length(new.dose))
  dose <- rep(new.dose,times=3)
  plotdata <- data.frame(drug=drug_name,dose=dose,y=f.new.dose1,method=method)
  plotdata = list(plotdata=plotdata,df=df)
  return(plotdata)
}


myggplot2 <- function(plotdata){
  df <- plotdata$df
  plotdata <- plotdata$plotdata

  g <- ggplot(data = plotdata,aes(x=dose,y=y,color=method)) +
    geom_line(size=1.3)+ # c('lightblue2','steelblue','darkred')+
    scale_color_manual(values=c('lightblue2','steelblue','darkred'))
    #ylim(0.5,2)+
    # theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
    #       legend.position = c(0.72,0.95), legend.key.size = unit(3, "cm"), legend.text = element_text(size=12),
    #       legend.key.height  = unit(0.5,"cm"),legend.title = element_blank(), legend.background = element_blank(),
    #       axis.text.x = element_text(face='bold',size=14),
    #       axis.text.y = element_text(face='bold',size=14))

  g <- g + ggplot2::facet_wrap(~drug,scales = 'free')

  g <- g + ggplot2::labs(y="odds ratio (OR)", x="Dose")

  g <- g + ggplot2::scale_linetype_manual(name="")

  g <- g+theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
          axis.text.x = element_text(face='bold',size=14),axis.text.y = element_text(face='bold',size=14),
          axis.title.x=element_text(size=16,face = "bold"),axis.title.y=element_text(size=16,face = "bold"),
          # plot.title = element_text(size = 16, face = "bold"),
          strip.background =element_rect(fill="snow3"),
          strip.text.x = element_text(size = 16),
          legend.text = element_text(size=12))

  g#+geom_rug(data=df,mapping=aes(x=obsDose),inherit.aes = F,color='red',length = unit(0.05, "npc"))
}
#







# Figure 2b
# df2 <-  plotdata2(dataDrug = data_per_drug(antidep,drug_name[1]),drug_name=drug_name[1],res[['fluoxetine']])
# names(df2)
# ggplot(data = df2,aes(x=dose,y=y,color=method)) +
#   geom_line(size=1.3)+ # c('lightblue2','steelblue','darkred')+
#   scale_color_manual(values=c('lightblue2','steelblue','darkred'))+
#   xlab('')+
#   ylab('')+
#   #ylim(0.5,2)+
#   theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
#         legend.position = c(0.72,0.95), legend.key.size = unit(3, "cm"), legend.text = element_text(size=12),
#         legend.key.height  = unit(0.5,"cm"),legend.title = element_blank(), legend.background = element_blank(),
#         axis.text.x = element_text(face='bold',size=14),
#         axis.text.y = element_text(face='bold',size=14))

