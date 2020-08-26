plotdata <- function(dataDrug,doseresORsplineBin){
  # new dose range to plot the results of the three curves
  new.dose <- 0:max(dataDrug$Dose_delivered_mean)
  knots = quantile(dataDrug$Dose_delivered_mean,probs=c(0.10,0.50,0.90))
  new.dose1 <- c(rcs(new.dose,knots)[,1])
  new.dose2 <- c(rcs(new.dose,knots)[,2])

  # Figure 2a:  absoulte responses vs  dose and add placebo response effect

  p.placebo<- exp(doseresORsplineBin$BUGSoutput$mean$Z)/(1+exp(doseresORsplineBin$BUGSoutput$mean$Z))
  p.drug <- doseresORsplineBin$BUGSoutput$mean$p.drug
  l.ci <-  doseresORsplineBin$BUGSoutput$summary[paste0('p.drug[',1:max(new.dose),']'),'2.5%']
  u.ci <- doseresORsplineBin$BUGSoutput$summary[paste0('p.drug[',1:max(new.dose),']'),'97.5%']
  lp.ci <-  exp(doseresORsplineBin$BUGSoutput$summary['Z','2.5%'])/(1+exp(doseresORsplineBin$BUGSoutput$summary['Z','2.5%']))
  up.ci <- exp(doseresORsplineBin$BUGSoutput$summary['Z','97.5%'])/(1+exp(doseresORsplineBin$BUGSoutput$summary['Z','97.5%']))
plotdata <- data.frame(p.drug=p.drug, drug=drug_name,p.placebo=p.placebo,lp.ci=lp.ci,up.ci=up.ci)
plotdata
}

# beta1 and beta2
beta1fOR <- coef(doseresORsplineFreq)[1]
beta2fOR <- coef(doseresORsplineFreq)[2]

beta1nOR <- doseresORsplineNor$BUGSoutput$mean$beta1.pooled
beta2nOR <- doseresORsplineNor$BUGSoutput$mean$beta2.pooled

beta1bOR <- doseresORsplineBin$BUGSoutput$mean$beta1.pooled
beta2bOR <- doseresORsplineBin$BUGSoutput$mean$beta2.pooled

# new dose range to plot the results of the three curves
new.dose <- 0:max(dataDrug$Dose_delivered_mean)
knots = quantile(dataDrug$Dose_delivered_mean,probs=c(0.10,0.50,0.90))
new.dose1 <- c(rcs(new.dose,knots)[,1])
new.dose2 <- c(rcs(new.dose,knots)[,2])

# Figure 2a:  absoulte responses vs  dose and add placebo response effect

p.placebo<- exp(doseresORsplineBin$BUGSoutput$mean$Z)/(1+exp(doseresORsplineBin$BUGSoutput$mean$Z))
p.drug <- doseresORsplineBin$BUGSoutput$mean$p.drug
l.ci <-  doseresORsplineBin$BUGSoutput$summary[paste0('p.drug[',1:max(new.dose),']'),'2.5%']
u.ci <- doseresORsplineBin$BUGSoutput$summary[paste0('p.drug[',1:max(new.dose),']'),'97.5%']
lp.ci <-  exp(doseresORsplineBin$BUGSoutput$summary['Z','2.5%'])/(1+exp(doseresORsplineBin$BUGSoutput$summary['Z','2.5%']))
up.ci <- exp(doseresORsplineBin$BUGSoutput$summary['Z','97.5%'])/(1+exp(doseresORsplineBin$BUGSoutput$summary['Z','97.5%']))



df1 <- data.frame(new.dose=new.dose[-1],y1 =p.drug,y2=l.ci,
                  y3=u.ci,yp1=p.placebo,yp2=lp.ci,yp3=up.ci)

ggplot(data = df1,aes(x=new.dose)) +
  geom_line(aes(y=y1,color='treatment'))+
  geom_line(aes(y=yp1,color='placebo'))+
  scale_color_manual(values=c('steelblue','coral4'))+
  geom_smooth(aes(x=new.dose, y=y1, ymin=y2,
                  ymax=y3),color='coral4',fill='coral1',
              data=df1, stat="identity")+
  geom_smooth(aes(x=new.dose, y=yp1, ymin=yp2,
                  ymax=yp3),color='steelblue',fill='steelblue2',
              data=df1, stat="identity")+
  xlab('')+
  ylab('')+
 # ylim(0.2,0.6)+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        legend.position = c(0.72,0.95), legend.key.size = unit(3, "cm"), legend.text = element_text(size=12),
        legend.key.height  = unit(0.5,"cm"),legend.title = element_blank(), legend.background = element_blank(),
        axis.text.x = element_text(face='bold',size=14),
        axis.text.y = element_text(face='bold',size=14))

# Figure 2b
df2 <- data.frame(new.dose1=new.dose1,y1 =c(exp(beta1fOR*new.dose1+beta2fOR*new.dose2)
                                            ,exp(beta1nOR*new.dose1+beta2nOR*new.dose2),
                                            exp(beta1bOR*new.dose1+beta2bOR*new.dose2)),
                  method=rep(c('binomial Bayesian','normal Bayesian', 'one-stage (freq)'),each=length(new.dose)))

ggplot(data = df2,aes(x=new.dose1,y=y1,color=method)) +
  geom_line(size=1.3)+ # c('lightblue2','steelblue','darkred')+
  scale_color_manual(values=c('lightblue2','steelblue','darkred'))+
  xlab('')+
  ylab('')+
  #ylim(0.5,2)+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        legend.position = c(0.72,0.95), legend.key.size = unit(3, "cm"), legend.text = element_text(size=12),
        legend.key.height  = unit(0.5,"cm"),legend.title = element_blank(), legend.background = element_blank(),
        axis.text.x = element_text(face='bold',size=14),
        axis.text.y = element_text(face='bold',size=14))

