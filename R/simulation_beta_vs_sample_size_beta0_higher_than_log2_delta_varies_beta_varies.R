#------------------------------------------------------------------------------------
## Clear the environment in R ##
#------------------------------------------------------------------------------------
rm(list = ls())


#------------------------------------------------------------------------------------
## Install the required packages into ##
## R if it's not already installed.   ##
#------------------------------------------------------------------------------------
if(!require(latex2exp)) install.packages("latex2exp")
if(!require(ggplot2)) install.packages("ggplot2")


#------------------------------------------------------------------------------------
## Load the packages into the R environment ##
#------------------------------------------------------------------------------------
library(latex2exp)
library(ggplot2)


#------------------------------------------------------------------------------------
hbetap<-function(x,beta,p){
  return(-1*(beta*x^p - 0.5*((1+x)*log(1+x)+(1-x)*log(1-x))))
}


beta0 = 0.7#log(2)#
beta1 = 0.71
delta = 0.05
deltap = 0.1

plist = c(2,3,4)

beta1list = seq(beta1,beta1+deltap,deltap/10)

mstarbetapoptim = matrix(rep(0,length(beta1list)*length(plist)),nrow=length(plist))
bs = matrix(rep(0,length(beta1list)*length(plist)),nrow=length(plist))


for (j in 1:length(beta1list)){
  for (k in 1:length(plist)){
    denom2 = optim(par=0.5, fn = hbetap, beta = beta0, p = plist[k],
                   method="L-BFGS-B",lower = 0.8, upper = 1-exp(-20) )$value*(-1)
    
    mstarbetapoptim[k,j] = optim(par=0.51, fn = hbetap, beta = beta1list[j],
                                 p = plist[k], method="L-BFGS-B",lower = 0.8,
                                 upper = 1-exp(-20) )$par
    
    
    denom1 = hbetap(mstarbetapoptim[k,j],beta=beta0,p=plist[k])*(-1)
    print(denom1)
    
    bs[k,j] = log(delta)/(denom1-denom2)
  }
}


df = t(data.frame(bs))
df = data.frame(cbind(df,beta1list))
colnames(df) = c("p2","p3","p4","beta1")

dfdelta = data.frame("delta"=seq(0.05,0.95,0.01), "sample_size_p2"=(df[1,1]/log(delta))*log(seq(0.05,0.95,0.01)),
                     "sample_size_p3"=(df[1,2]/log(delta))*log(seq(0.05,0.95,0.01)),
                     "sample_size_p4"=(df[1,3]/log(delta))*log(seq(0.05,0.95,0.01)) )


## Plot 1
pdf("pbeta71deltavaries.pdf")
ggplot() +
  geom_line(data = dfdelta, aes(x=delta, y=sample_size_p2, linetype = "p=2")) +
  geom_line(data = dfdelta, aes(x=delta, y=sample_size_p3, linetype = "p=3")) +
  geom_line(data = dfdelta, aes(x=delta, y=sample_size_p4, linetype = "p=4")) +
  scale_y_continuous(trans = 'log10') +
  labs(x=TeX("$\\delta$"), y="Optimal sample size",
       title = TeX("Optimal sample size for $\\beta_0 = 0.7$ and $\\beta = 0.71$") ) +
  scale_linetype_manual(name="",values=c("p=2"=1,"p=3"=2,"p=4"=3))
dev.off()


## Plot 2
pdf("pbeta71deltavariesnotitle.pdf")
ggplot() +
  geom_line(data = dfdelta, aes(x=delta, y=sample_size_p2, linetype = "p=2")) +
  geom_line(data = dfdelta, aes(x=delta,y=sample_size_p3,linetype = "p=3")) +
  geom_line(data = dfdelta, aes(x=delta,y=sample_size_p4,linetype = "p=4")) +
  scale_y_continuous(trans='log10') +
  labs(x=TeX("$\\delta$"),y="Optimal sample size") +
  scale_linetype_manual(name="",values=c("p=2"=1,"p=3"=2,"p=4"=3))
dev.off()


## Plot 3
pdf("pdelta5betavaries.pdf")
ggplot() +
  geom_point(data = df, aes(x=beta1,y=p2,shape="p=2")) +
  geom_point(data = df, aes(x=beta1,y=p3,shape="p=3")) +
  geom_point(data = df, aes(x=beta1,y=p4,shape="p=4")) +
  scale_y_continuous(trans='log10') +
  labs(x=TeX("$\\beta$"), y="Optimal sample size", 
       title = TeX("Optimal sample size for $\\beta_0 = 0.7$ and $\\delta = 0.05$")) +
  scale_shape_manual(name="",values=c("p=2"=1,"p=3"=2,"p=4"=3))
dev.off()


## Plot 4
pdf("pdelta5betavariesnotitle.pdf")
ggplot() +
  geom_point(data = df, aes(x=beta1,y=p2,shape="p=2")) +
  geom_point(data = df, aes(x=beta1,y=p3,shape="p=3")) +
  geom_point(data = df, aes(x=beta1,y=p4,shape="p=4")) +
  scale_y_continuous(trans='log10') +
  labs(x=TeX("$\\beta$"), y="Optimal sample size") +
  scale_shape_manual(name="", values=c("p=2"=1,"p=3"=2,"p=4"=3))
dev.off()


