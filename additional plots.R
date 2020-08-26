# load
load('S1ORspline2')
load('S2ORspline2')
load('S3ORspline2')


# S1
rval <- rbind(S4ORspline2$res1,S5ORspline2$res1,S3ORspline2$res1)
#  our x : doses
d <- 0:10
knots<-unlist(round(quantile(d,c(0.25,0.5,0.75))))
new.dose<-rcspline.eval(d,knots,inclx = T)
new.dose1 <- c(new.dose[,1])
new.dose2 <- c(new.dose[,2])
n.d <- length(d)

# our beta1
beta1 <- cbind(rval[c('true.beta1')],bayesB.beta1=rval[c('true.beta1')]+rval[c('BayesB1bias')],
                 bayesN.beta1= rval[c('true.beta1')]+rval[c('BayesN1bias')],
                 freq.beta1=rval[c('true.beta1')]+rval[c('Freq1bias')])
colnames(beta1) <- c('true.beta1','bayesB.beta1','bayesN.beta1','freq.beta1')
#beta1 <- mybeta1[-c(1,5),]

# our beta2
beta2 <- cbind(rval[c('true.beta2')],bayesB.beta1=rval[c('true.beta2')]+rval[c('BayesB2bias')],
                 bayesN.beta1= rval[c('true.beta2')]+rval[c('BayesN2bias')],
                 freq.beta1=rval[c('true.beta2')]+rval[c('Freq2bias')])
colnames(beta2) <- c('true.beta2','bayesB.beta2','bayesN.beta2','freq.beta2')
#beta2 <- mybeta2[-c(1,5),]

# our y = OR

y <- array(NA,dim=c(n.d,4,nrow(beta1)))
for (k in 1:nrow(beta1)) {
  for (j in 1:4) {
    y[1:n.d,j,k] <- exp(beta1[k,j]*new.dose1+beta2[k,j]*new.dose2)
  }
}

# plot x vs y = exp(beta1*x1+beta2*x2)
df11 <- data.frame(new.dose1=new.dose1,y1 =c(y[,1,1],y[,2,1],
                                             y[,3,1],y[,4,1]),method=rep(c('true curve','binomial Bayesian','normal Bayesian', 'one-stage (freq)'),each=11))
df1 <- data.frame(new.dose1=new.dose1,y1 =cbind(y[,1,1]),y2=cbind(y[,2,1]),
                  y3=cbind(y[,3,1]),y4=cbind(y[,4,1]))
theme_set(
  theme_minimal() +
    theme(legend.position = "right",axis.text.x = element_text(size=12),axis.text.y = element_text(size=12))
)

# S1
ggplot(data = df11,aes(x=new.dose1,y=y1,color=method)) +
  geom_line(size=1.3)+
  xlab('')+
  ylab('')+
  ylim(0.5,3.5)+
  scale_color_manual(values=c('lightblue2','steelblue','darkred','black'))+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        legend.position = c(0.72,0.90), legend.key.size = unit(3, "cm"), legend.text = element_text(size=12),
        legend.key.height  = unit(0.5,"cm"),legend.title = element_blank(), legend.background = element_blank(),
        axis.text.x = element_text(face='bold',size=14),axis.text.y = element_text(face='bold',size=14))

# S1 - zoom-in
ggplot(data = df1[9:11,],aes(x=new.dose1)) +
  geom_line(aes(y =y1),color='black' , size=1.3) +
  geom_line(aes(y = y2),color='lightblue2' , size=1.3)+
  geom_line(aes(y = y3),color='steelblue' , size=1.3)+
  geom_line(aes(y = y4),color='darkred' , size=1.3)+
  xlab('')+
  ylab('')+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        legend.position = "right",axis.text.x = element_text(face='bold',size=14),axis.text.y =element_text(face='bold',size=14))

# S2
df2 <- data.frame(new.dose1=new.dose1,y1 =cbind(y[,1,2]),y2=cbind(y[,2,2]),
                  y3=cbind(y[,3,2]),y4=cbind(y[,4,2]))


ggplot(data = df2,aes(x=new.dose1)) +
  geom_line(aes(y =y1),color='black' , size=1.3) +
  geom_line(aes(y = y2),color='lightblue2' , size=1.3)+
  geom_line(aes(y = y3),color='steelblue' , size=1.3)+
  geom_line(aes(y = y4),color='darkred' , size=1.3)+
  xlab('')+
  ylab('')+
  ylim(0.5,3.5)+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        legend.position = "right",axis.text.x = element_text(face='bold',size=14),axis.text.y = element_text(face='bold',size=14))

# S2 - zoom-in
ggplot(data = df2[9:11,],aes(x=new.dose1)) +
  geom_line(aes(y =y1),color='black' , size=1.3) +
  geom_line(aes(y = y2),color='lightblue2' , size=1.3)+
  geom_line(aes(y = y3),color='steelblue' , size=1.3)+
  geom_line(aes(y = y4),color='darkred' , size=1.3)+
  xlab('')+
  ylab('')+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        legend.position = "right",axis.text.x = element_text(face='bold',size=14),axis.text.y =element_text(face='bold',size=14))

# S3
df3 <- data.frame(new.dose1=new.dose1,y1 =cbind(y[,1,3]),y2=cbind(y[,2,3]),
                  y3=cbind(y[,3,3]),y4=cbind(y[,4,3]))


ggplot(data = df3,aes(x=new.dose1)) +
  geom_line(aes(y =y1),color='black' , size=1.3) +
  geom_line(aes(y = y2),color='lightblue2' , size=1.3)+
  geom_line(aes(y = y3),color='steelblue' , size=1.3)+
  geom_line(aes(y = y4),color='darkred' , size=1.3)+
  xlab('')+
  ylab('')+
  ylim(0.5,3.5)+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        legend.position = "right",axis.text.x = element_text(face='bold',size=14),axis.text.y =element_text(face='bold',size=14))

 # S3  zoom-in

ggplot(data = df3[9:11,],aes(x=new.dose1)) +
  geom_line(aes(y =y1),color='black' , size=1.3) +
  geom_line(aes(y = y2),color='lightblue2' , size=1.3)+
  geom_line(aes(y = y3),color='steelblue' , size=1.3)+
  geom_line(aes(y = y4),color='darkred' , size=1.3)+
  xlab('')+
  ylab('')+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        legend.position = "right",axis.text.x = element_text(face='bold',size=14),axis.text.y =element_text(face='bold',size=14))










