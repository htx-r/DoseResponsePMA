ORdose_singleplot <- function(rval,dmax=10, knots=NULL,ymin=0.5,ymax=3.5){

  #** x = dose
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


  #**  y = OR
  y <- matrix(NA,n.d,length(beta1))

  for (j in 1:length(beta1)) {
    y[,j] <- exp(beta1[j]*new.dose1+beta2[j]*new.dose2)
  }


  #** plot x = dose vs y= OR = exp(beta1*x1+beta2*x2)
   # data to plot
  df <- data.frame(dose=new.dose1,y=c(y[,1],y[,2],
                                      y[,3],y[,4]),method=rep(c('true curve','binomial Bayesian','normal Bayesian', 'one-stage (freq)'),each=100))
   # plot theme
  theme_set(
    theme_minimal() +
      theme(legend.position = "right",axis.text.x = element_text(size=12),axis.text.y = element_text(size=12))
  )

   # the final plot
  ggplot(data = df,aes(x=dose,y=y,color=method)) +
    geom_line(size=1.3)+
    xlab('')+
    ylab('')+
    ylim(ymin,ymax)+
    scale_color_manual(values=c('lightblue2','steelblue','darkred','black'))+
    theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
          legend.position = c(0.72,0.90), legend.key.size = unit(3, "cm"), legend.text = element_text(size=12),
          legend.key.height  = unit(0.5,"cm"),legend.title = element_blank(), legend.background = element_blank(),
          axis.text.x = element_text(face='bold',size=14),axis.text.y = element_text(face='bold',size=14))
}



ORdose_multiplot <- function(){

}


absoluteDose_singleplot <- function(){

}

absoluteDose_multiplot <- function(){

}


















