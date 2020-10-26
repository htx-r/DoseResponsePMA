library("gridExtra")

# load
load('S8ORspline2.6')
load('S5ORspline2')
load('S3ORspline2')

ORdose_singleplot <- function(rval,dmax=10, knots=NULL,ymin=0.5,ymax=3.5){

  # x = dose
  d <- seq(0,dmax,l=100)

  if(is.null(knots)){
    knots <- unlist(round(quantile(d,c(0.25,0.5,0.75))))
  } else{
    knots <- knots
  }

  new.dose<-rcspline.eval(d,knots,inclx = T)
  new.dose1 <- c(new.dose[,1])
  new.dose2 <- c(new.dose[,2])
  n.d <- length(d)

  # beta1
  beta1 <- cbind(bayesB.beta1=rval[c('BayesB1bias')],
                 bayesN.beta1= rval[c('BayesN1bias')],
                 freq.beta1=rval[c('Freq1bias')])
  beta1 <- unlist(cbind(rval[c('true.beta1')],beta1+c(rval[c('true.beta1')])))

  # beta2
  beta2 <- cbind(bayesB.beta1=rval[c('BayesB2bias')],
                 bayesN.beta1= rval[c('BayesN2bias')],
                 freq.beta1=rval[c('Freq2bias')])
  beta2 <-unlist(cbind(rval[c('true.beta2')],beta2+c(rval[c('true.beta2')])))


  #  y = OR
  y <- matrix(NA,n.d,length(beta1))

  for (j in 1:length(beta1)) {
    y[,j] <- exp(beta1[j]*new.dose1+beta2[j]*new.dose2)
  }


  # plot x = dose vs y= OR = exp(beta1*x1+beta2*x2)
  df <- data.frame(dose=new.dose1,y=c(y[,1],y[,2],
                                      y[,3],y[,4]),method=rep(c('true curve','binomial Bayesian','normal Bayesian', 'one-stage (freq)'),each=100))
  theme_set(
    theme_minimal() +
      theme(legend.position = "top",axis.text.x = element_text(size=12),axis.text.y = element_text(size=12))
  )

  # S1
  ggplot(data = df,aes(x=dose,y=y,color=method)) +
    geom_line(size=1.3)+
    xlab('')+
    ylab('')+
    ylim(ymin,ymax)+
    scale_color_manual(values=c('lightblue2','steelblue','darkred','black'))+
    theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
          #legend.position = c(-100,100),
          legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=16), # c(0.72,0.90)
          legend.key.height  = unit(0.5,"cm"),legend.title = element_blank(),
          axis.text.x =element_blank(),axis.text.y =element_text(face='bold',size=14)) #element_text(face='bold',size=14)
}
#load('S8ORspline2.8')
# 1. Create the 4 plots
#+++++++++++++++++++++++
p1 <- ORdose_singleplot(rval=S8ORspline2.5$res1,ymax = 3.5)
p2 <- ORdose_singleplot(rval=S8ORspline2.6$res1,ymax = 3.5)
p3 <- ORdose_singleplot(rval=S8ORspline2.7$res1,ymax = 3.5)
p4 <- ORdose_singleplot(rval=S8ORspline2.8$res1,ymax = 3.5)

# 2. Save the legend
#+++++++++++++++++++++++
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

legend <- g_legend(p1)

# 4. Produce the multiplot
#+++++++++++++++++++++++
grid.arrange(legend,arrangeGrob(p1+ theme(legend.position="none"), p2+ theme(legend.position="none"),p3+ theme(legend.position="none"),p4+ theme(legend.position="none"),
            nrow=2), nrow = 2, heights = c(0.5,7))

#


# S5: quadratic
load('S16ORspline2')
rval <- colMeans(t(S16ORspline2))
#rval <- S2ORspline2$res1
#  our x : doses
d <- seq(0,10,l=1000)

#knots<-unlist(round(quantile(d,c(0.10,0.5,0.9))))

  dose<-rchisq(1000,2)
 knots <- quantile(dose,c(0.1,0.5,0.9))

knots <- c(0,1,3) #c(0,1,3)c(0,0.5,3)
new.dose<-rcspline.eval(d,knots,inclx = T)
new.dose1 <- c(new.dose[,1])
new.dose2 <- c(new.dose[,2])
n.d <- length(d)

# our beta1
beta1 <- cbind(1,bayesB.beta1=rval[c('BayesB1')],
               bayesN.beta1= rval[c('BayesN1')],
               freq.beta1=rval[c('Freq1')])
colnames(beta1) <- c('true.beta1','bayesB.beta1','bayesN.beta1','freq.beta1')
#beta1 <- mybeta1[-c(1,5),]

# our beta2
beta2 <- cbind(1,bayesB.beta2=rval[c('BayesB2')],
               bayesN.beta2= rval[c('BayesN2')],
               freq.beta2=rval[c('Freq2')])
colnames(beta2) <- c('true.beta2','bayesB.beta2','bayesN.beta2','freq.beta2')
#beta2 <- mybeta2[-c(1,5),]

y <- array(NA,dim=c(n.d,4,nrow(beta1)))
for (k in 1:nrow(beta1)) {
  #y[1:n.d,1,k] <- exp(new.dose1/sqrt(1+(new.dose1)^2)) # beta2[k,1]*ifelse(new.dose1==0,0,log(new.dose1))

  y[1:n.d,1,k] <- log(new.dose1+1)+1
  for (j in 2:4) {
    y[1:n.d,j,k] <- exp(beta1[k,j]*new.dose1+beta2[k,j]*new.dose2)
    #y[1:n.d,j,k] <- beta1[k,j]*new.dose1+beta2[k,j]*new.dose2
  }
}
var.logOR <- new.dose1^2*rval['sdBin1']^2+new.dose2^2*rval['sdBin2']^2
low.ci <- exp(log(y[,2,1]) - (1.96*sqrt(var.logOR)))
upp.ci <-exp(log( y[,2,1]) + (1.96*sqrt(var.logOR)))

df11 <- data.frame(new.dose1=new.dose1,y=c(y[,1,1],y[,2,1],y[,3,1],y[,4,1]), # y=y[,1,1],y1 =y[,2,1],
                   low.ci=low.ci,upp.ci=upp.ci,
                   method=rep(c('true curve','binomial Bayesian','normal Bayesian', 'one-stage (freq)'),each=1000))
# df1 <- data.frame(new.dose1=new.dose1,y1 =cbind(y[,1,1]),y2=cbind(y[,2,1]),
#                   y3=cbind(y[,3,1]),y4=cbind(y[,4,1]))
theme_set(
  theme_minimal() +
    theme(legend.position = "top",axis.text.x = element_text(size=12),axis.text.y = element_text(size=12))
)

# S1
# ggplot(data = df11,aes(x=new.dose1)) +
#   geom_line(size=1.3,aes(y=y))+
#   geom_smooth(ggplot2::aes(x=new.dose1,y=y1,ymin=low.ci,ymax=upp.ci),color='coral4',fill='coral1',
#              data=df11, stat="identity")+
# p1 <-ggplot(data = df11,aes(x=new.dose1,y=y,color=method)) +
#   geom_line(size=1.3)+
# xlab('')+
#   ylab('')+
#   # xlim(2.5,10)+
#   ylim(1,4.2)+
#   scale_color_manual(values=c('lightblue2','steelblue','darkred','black'))+
#   theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
#         #legend.position = c(-100,100),
#         legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=16), # c(0.72,0.90)
#         legend.key.height  = unit(0.5,"cm"),legend.title = element_blank(),
#         axis.text.x =element_blank(),axis.text.y =element_text(face='bold',size=14)) #element_text(face='bold',size=14)
# p2 <-ggplot(data = df11,aes(x=new.dose1,y=y,color=method)) +
#   geom_line(size=1.3)+
#   xlab('')+
#   ylab('')+
#   # xlim(2.5,10)+
#   ylim(1,4.2)+
#   scale_color_manual(values=c('lightblue2','steelblue','darkred','black'))+
#   theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
#         #legend.position = c(-100,100),
#         legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=16), # c(0.72,0.90)
#         legend.key.height  = unit(0.5,"cm"),legend.title = element_blank(),
#         axis.text.x =element_text(face='bold',size=14),axis.text.y =element_text(face='bold',size=14)) #element_text(face='bold',size=14)
# p3 <-ggplot(data = df11,aes(x=new.dose1,y=y,color=method)) +
#   geom_line(size=1.3)+
#   xlab('')+
#   ylab('')+
#   # xlim(2.5,10)+
#   ylim(1,4.2)+
#   scale_color_manual(values=c('lightblue2','steelblue','darkred','black'))+
#   theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
#         #legend.position = c(-100,100),
#         legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=16), # c(0.72,0.90)
#         legend.key.height  = unit(0.5,"cm"),legend.title = element_blank(),
#         axis.text.x =element_blank(),axis.text.y =element_blank()) #element_text(face='bold',size=14)
#
# p4 <-ggplot(data = df11,aes(x=new.dose1,y=y,color=method)) +
#   geom_line(size=1.3)+
#   xlab('')+
#   ylab('')+
#   # xlim(2.5,10)+
#   ylim(1,4.2)+
#   scale_color_manual(values=c('lightblue2','steelblue','darkred','black'))+
#   theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
#         #legend.position = c(-100,100),
#         legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=16), # c(0.72,0.90)
#         legend.key.height  = unit(0.5,"cm"),legend.title = element_blank(),
#         axis.text.x =element_text(face='bold',size=14),axis.text.y =element_blank()) #element_text(face='bold',size=14)

# legend <- g_legend(p4)
#
# # 4. Produce the multiplot
# #+++++++++++++++++++++++
# grid.arrange(legend,arrangeGrob(p1+ theme(legend.position="none"), p3+ theme(legend.position="none"),p2+ theme(legend.position="none"),p4+ theme(legend.position="none"),
#                                 nrow=2), nrow = 2, heights = c(0.5,7))
#

# non RCS chi2
pp1 <-ggplot(data = df11,aes(x=new.dose1,y=y,color=method)) +
  geom_line(size=1.3)+
  xlab('')+
  ylab('')+
  # xlim(2.5,10)+
  ylim(1,4.2)+
  scale_color_manual(values=c('lightblue2','steelblue','darkred','black'))+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        #legend.position = c(-100,100),
        legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=16), # c(0.72,0.90)
        legend.key.height  = unit(0.5,"cm"),legend.title = element_blank(),
        axis.text.x =element_blank(),axis.text.y =element_text(face='bold',size=14)) #element_text(face='bold',size=14)
pp2 <-ggplot(data = df11,aes(x=new.dose1,y=y,color=method)) +
  geom_line(size=1.3)+
  xlab('')+
  ylab('')+
  # xlim(2.5,10)+
  ylim(1,4.2)+
  scale_color_manual(values=c('lightblue2','steelblue','darkred','black'))+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        #legend.position = c(-100,100),
        legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=16), # c(0.72,0.90)
        legend.key.height  = unit(0.5,"cm"),legend.title = element_blank(),
        axis.text.x =element_text(face='bold',size=14),axis.text.y =element_text(face='bold',size=14)) #element_text(face='bold',size=14)
pp3 <-ggplot(data = df11,aes(x=new.dose1,y=y,color=method)) +
  geom_line(size=1.3)+
  xlab('')+
  ylab('')+
  # xlim(2.5,10)+
  ylim(1,4.2)+
  scale_color_manual(values=c('lightblue2','steelblue','darkred','black'))+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        #legend.position = c(-100,100),
        legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=16), # c(0.72,0.90)
        legend.key.height  = unit(0.5,"cm"),legend.title = element_blank(),
        axis.text.x =element_blank(),axis.text.y =element_blank()) #element_text(face='bold',size=14)

pp4 <-ggplot(data = df11,aes(x=new.dose1,y=y,color=method)) +
  geom_line(size=1.3)+
  xlab('')+
  ylab('')+
  # xlim(2.5,10)+
  ylim(1,4.2)+
  scale_color_manual(values=c('lightblue2','steelblue','darkred','black'))+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        #legend.position = c(-100,100),
        legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=16), # c(0.72,0.90)
        legend.key.height  = unit(0.5,"cm"),legend.title = element_blank(),
        axis.text.x =element_text(face='bold',size=14),axis.text.y =element_blank()) #element_text(face='bold',size=14)

legend <- g_legend(pp4)

# 4. Produce the multiplot
#+++++++++++++++++++++++
grid.arrange(legend,arrangeGrob(pp1+ theme(legend.position="none"), pp3+ theme(legend.position="none"),pp2+ theme(legend.position="none"),pp4+ theme(legend.position="none"),
                                nrow=2), nrow = 2, heights = c(0.5,7))

#
x <- seq(0,10,l=1000)
par(mfrow=c(1,2))
plot(x,log(x+1),type = 'l',ylim = c(-1,3),ylab = 'log OR',xlab='Dose')
lines(x,log(log(x+1)),col=2)
lines(x,log(log(x)+1),col=3)

plot(x,x+1,ylim = c(-1,10),type = 'l',ylab = 'OR',xlab='Dose')
lines(x,log(x+1),col=2)
lines(x,log(x)+1,col=3)
#








rval <- t(S9ORspline2)
#  our x : doses
d <- seq(1,10,l=100)
knots<-unlist(round(quantile(d,c(0.25,0.5,0.75))))
new.dose<-rcspline.eval(d,knots,inclx = T)
new.dose1 <- c(new.dose[,1])
new.dose2 <- c(new.dose[,2])
n.d <- length(d)

# our beta1
beta1 <- rval[,c('BayesB1')]
#beta1 <- mybeta1[-c(1,5),]

# our beta2
beta2 <-rval[,c('BayesB2')]
#beta2 <- mybeta2[-c(1,5),]

y <- matrix(NA,100,100)
for (i in 1:length(beta1)) {
  y[i,] <- exp(beta1[i]*new.dose1+beta2[i]*new.dose2)
}
y_mean <- exp(mean(beta1)*new.dose1+mean(beta2)*new.dose2)
y_median <- exp(median(beta1)*new.dose1+median(beta2)*new.dose2)

par(mfrow = c(1, 1),las=1)
plot(new.dose1,y_mean,type = 'l',ylim = c(1,4),lwd=5)
cols <- rainbow(100)
for(i in 1:100){
  lines(new.dose1,y[i,], col = cols[i], lwd=1)
  }
lines(new.dose1,y_mean,lwd=5)
#lines(new.dose1,y_median,lwd=5,col=4)
lines(new.dose1,colMeans(y),lwd=5,col=2)
lines(new.dose1,log(new.dose1)+1,lwd=7,col='grey40')

exp(0.212)
































# # S1
# rval <- rbind(S8ORspline2.6$res1)
# #  our x : doses
# d <- seq(0,10,l=100)
# knots<-unlist(round(quantile(d,c(0.25,0.5,0.75))))
# new.dose<-rcspline.eval(d,knots,inclx = T)
# new.dose1 <- c(new.dose[,1])
# new.dose2 <- c(new.dose[,2])
# n.d <- length(d)
#
# # our beta1
# beta1 <- cbind(rval[c('true.beta1')],bayesB.beta1=rval[c('true.beta1')]+rval[c('BayesB1bias')],
#                  bayesN.beta1= rval[c('true.beta1')]+rval[c('BayesN1bias')],
#                  freq.beta1=rval[c('true.beta1')]+rval[c('Freq1bias')])
# colnames(beta1) <- c('true.beta1','bayesB.beta1','bayesN.beta1','freq.beta1')
# #beta1 <- mybeta1[-c(1,5),]
#
# # our beta2
# beta2 <- cbind(rval[c('true.beta2')],bayesB.beta1=rval[c('true.beta2')]+rval[c('BayesB2bias')],
#                  bayesN.beta1= rval[c('true.beta2')]+rval[c('BayesN2bias')],
#                  freq.beta1=rval[c('true.beta2')]+rval[c('Freq2bias')])
# colnames(beta2) <- c('true.beta2','bayesB.beta2','bayesN.beta2','freq.beta2')
# #beta2 <- mybeta2[-c(1,5),]
#
# # our y = OR
#
# y <- array(NA,dim=c(n.d,4,nrow(beta1)))
# for (k in 1:nrow(beta1)) {
#   for (j in 1:4) {
#     y[1:n.d,j,k] <- exp(beta1[k,j]*new.dose1+beta2[k,j]*new.dose2)
#   }
# }
#
# # plot x vs y = exp(beta1*x1+beta2*x2)
# df11 <- data.frame(new.dose1=new.dose1,y1 =c(y[,1,1],y[,2,1],
#                                              y[,3,1],y[,4,1]),method=rep(c('true curve','binomial Bayesian','normal Bayesian', 'one-stage (freq)'),each=100))
# df1 <- data.frame(new.dose1=new.dose1,y1 =cbind(y[,1,1]),y2=cbind(y[,2,1]),
#                   y3=cbind(y[,3,1]),y4=cbind(y[,4,1]))
# theme_set(
#   theme_minimal() +
#     theme(legend.position = "right",axis.text.x = element_text(size=12),axis.text.y = element_text(size=12))
# )
#
# # S1
# ggplot(data = df11,aes(x=new.dose1,y=y1,color=method)) +
#   geom_line(size=1.3)+
#   xlab('')+
#   ylab('')+
#   ylim(0.5,3.5)+
#   scale_color_manual(values=c('lightblue2','steelblue','darkred','black'))+
#   theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
#         legend.position = c(0.72,0.90), legend.key.size = unit(3, "cm"), legend.text = element_text(size=12),
#         legend.key.height  = unit(0.5,"cm"),legend.title = element_blank(), legend.background = element_blank(),
#         axis.text.x = element_text(face='bold',size=14),axis.text.y = element_text(face='bold',size=14))
#
# # S1 - zoom-in
# ggplot(data = df1,aes(x=new.dose1)) + # [9:11,]
#   geom_line(aes(y =y1),color='black' , size=1.3) +
#   geom_line(aes(y = y2),color='lightblue2' , size=1.3)+
#   geom_line(aes(y = y3),color='steelblue' , size=1.3)+
#   geom_line(aes(y = y4),color='darkred' , size=1.3)+
#   xlab('')+
#   ylab('')+
#   theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
#         legend.position = "right",axis.text.x = element_text(face='bold',size=14),axis.text.y =element_text(face='bold',size=14))
#
# # S2
# df2 <- data.frame(new.dose1=new.dose1,y1 =cbind(y[,1,2]),y2=cbind(y[,2,2]),
#                   y3=cbind(y[,3,2]),y4=cbind(y[,4,2]))
#
#
# ggplot(data = df2,aes(x=new.dose1)) +
#   geom_line(aes(y =y1),color='black' , size=1.3) +
#   geom_line(aes(y = y2),color='lightblue2' , size=1.3)+
#   geom_line(aes(y = y3),color='steelblue' , size=1.3)+
#   geom_line(aes(y = y4),color='darkred' , size=1.3)+
#   xlab('')+
#   ylab('')+
#   ylim(0.5,3.5)+
#   theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
#         legend.position = "right",axis.text.x = element_text(face='bold',size=14),axis.text.y = element_text(face='bold',size=14))
#
# # S2 - zoom-in
# ggplot(data = df2[9:11,],aes(x=new.dose1)) +
#   geom_line(aes(y =y1),color='black' , size=1.3) +
#   geom_line(aes(y = y2),color='lightblue2' , size=1.3)+
#   geom_line(aes(y = y3),color='steelblue' , size=1.3)+
#   geom_line(aes(y = y4),color='darkred' , size=1.3)+
#   xlab('')+
#   ylab('')+
#   theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
#         legend.position = "right",axis.text.x = element_text(face='bold',size=14),axis.text.y =element_text(face='bold',size=14))
#
# # S3
# df3 <- data.frame(new.dose1=new.dose1,y1 =cbind(y[,1,3]),y2=cbind(y[,2,3]),
#                   y3=cbind(y[,3,3]),y4=cbind(y[,4,3]))
#
#
# ggplot(data = df3,aes(x=new.dose1)) +
#   geom_line(aes(y =y1),color='black' , size=1.3) +
#   geom_line(aes(y = y2),color='lightblue2' , size=1.3)+
#   geom_line(aes(y = y3),color='steelblue' , size=1.3)+
#   geom_line(aes(y = y4),color='darkred' , size=1.3)+
#   xlab('')+
#   ylab('')+
#   ylim(0.5,3.5)+
#   theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
#         legend.position = "right",axis.text.x = element_text(face='bold',size=14),axis.text.y =element_text(face='bold',size=14))
#
#  # S3  zoom-in
#
# ggplot(data = df3[9:11,],aes(x=new.dose1)) +
#   geom_line(aes(y =y1),color='black' , size=1.3) +
#   geom_line(aes(y = y2),color='lightblue2' , size=1.3)+
#   geom_line(aes(y = y3),color='steelblue' , size=1.3)+
#   geom_line(aes(y = y4),color='darkred' , size=1.3)+
#   xlab('')+
#   ylab('')+
#   theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
#         legend.position = "right",axis.text.x = element_text(face='bold',size=14),axis.text.y =element_text(face='bold',size=14))
#









