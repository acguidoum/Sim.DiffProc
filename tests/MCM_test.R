options(prompt="R> ",scipen=16,digits=4,warning=FALSE, message=FALSE)
library(Sim.DiffProc)

RNGkind(kind="L'Ecuyer-CMRG")
######

theta = 0.75; x0 = 1
fx <- expression( 0.5*theta^2*x )
gx <- expression( theta*x )
mod1 <- snssde1d(drift=fx,diffusion=gx,x0=x0,M=10,type="ito")
mod2 <- snssde1d(drift=fx,diffusion=gx,x0=x0,M=10,type="str")
## True values of means and variance for mod1 and mod2
E.mod1 <- function(t) x0*exp(0.5*theta^2*t)
E.mod2 <- function(t) x0
V.mod1 <- function(t) x0^2*exp(theta^2*t)*(exp(theta^2*t)-1)
V.mod2 <- function(t) x0^2 *(exp(theta^2*t)-1)
## function of the statistic(s) of interest.
sde.fun1d <- function(data, i){
     d <- data[i, ]
     return(c(mean(d),var(d)))
}
# Parallel MOnte Carlo for mod1
mcm.mod1 = MCM.sde(model=mod1,statistic=sde.fun1d,R=5, exact=list(m=E.mod1(1),S=V.mod1(1)))
mcm.mod1 = MCM.sde(model=mod1,statistic=sde.fun1d,R=5, exact=list(m=E.mod1(1),S=V.mod1(1)),parallel="snow",cl= parallel::makeCluster(getOption("cl.cores", 2)),ncpus=2)
mcm.mod1 = MCM.sde(model=mod1,statistic=sde.fun1d,R=5, exact=list(m=E.mod1(1),S=V.mod1(1)),parallel="snow",ncpus=2)
print(mcm.mod1)
# Parallel MOnte Carlo for mod2
mcm.mod2 = MCM.sde(model=mod2,statistic=sde.fun1d,R=5, exact=list(m=E.mod2(1),S=V.mod2(1)),parallel="snow",ncpus=2)
print(mcm.mod2)

## ----fig.cap=' MC output of mean and variance of `mod1`', fig.env='figure*'----
# plot(s) of Monte Carlo outputs of mod1
plot(mcm.mod1,index = 1)  # mean
plot(mcm.mod1,index = 2)  # variance

## ----fig.cap=' MC output of mean and variance of `mod2`', fig.env='figure*'----
# plot(s) of Monte Carlo outputs of mod2
plot(mcm.mod2,index = 1)  # mean
plot(mcm.mod2,index = 2)  # variance

## ------------------------------------------------------------------------
mu=1;sigma=0.5;theta=2
x0=0;y0=0;init=c(x0,y0)
f <- expression(1/mu*(theta-x), x)  
g <- expression(sqrt(sigma),0)
OUI <- snssde2d(drift=f,diffusion=g,M=10,Dt=0.015,x0=c(x=0,y=0))
## true values of first and second moment at time 10
Ex <- function(t) theta+(x0-theta)*exp(-t/mu)
Vx <- function(t) 0.5*sigma*mu *(1-exp(-2*(t/mu)))
Ey <- function(t) y0+theta*t+(x0-theta)*mu*(1-exp(-t/mu))
Vy <- function(t) sigma*mu^3*((t/mu)-2*(1-exp(-t/mu))+0.5*(1-exp(-2*(t/mu))))
covxy <- function(t) 0.5*sigma*mu^2 *(1-2*exp(-t/mu)+exp(-2*(t/mu)))
tvalue = list(m1=Ex(10),m2=Ey(10),S1=Vx(10),S2=Vy(10),C12=covxy(10))
## function of the statistic(s) of interest.
sde.fun2d <- function(data, i){
  d <- data[i,]
  return(c(mean(d$x),mean(d$y),var(d$x),var(d$y),cov(d$x,d$y)))
}
## Parallel Monte-Carlo of 'OUI' at time 10
mcm.mod2d = MCM.sde(OUI,statistic=sde.fun2d,time=10,R=5,exact=tvalue)
mcm.mod2d = MCM.sde(OUI,statistic=sde.fun2d,time=10,R=5,exact=tvalue,parallel="snow",cl= parallel::makeCluster(getOption("cl.cores", 2)),ncpus=2)
mcm.mod2d = MCM.sde(OUI,statistic=sde.fun2d,time=10,R=5,exact=tvalue,parallel="snow",ncpus=2)
print(mcm.mod2d)
plot(mcm.mod2d)


mu=1;sigma=0.5;theta=2
x0=0;y0=0;init=c(x0,y0)
f <- expression(1/mu*(theta-x), x)  
g <- expression(sqrt(sigma),0)
Sigma <- matrix(c(1, 0.75, 0.75, 1), nrow = 2, ncol = 2)
OUI <- snssde2d(drift=f,diffusion=g,corr=Sigma,M=10,Dt=0.015,x0=c(x=0,y=0))

## function of the statistic(s) of interest.
sde.fun2d <- function(data, i){
  d <- data[i,]
  return(c(mean(d$x),mean(d$y),var(d$x),var(d$y),cov(d$x,d$y)))
}
## Parallel Monte-Carlo of 'OUI' at time 10
mcm.mod2d = MCM.sde(OUI,statistic=sde.fun2d,time=10,R=5)
mcm.mod2d = MCM.sde(OUI,statistic=sde.fun2d,time=10,R=5,parallel="snow",cl= parallel::makeCluster(getOption("cl.cores", 2)),ncpus=2)
mcm.mod2d = MCM.sde(OUI,statistic=sde.fun2d,time=10,R=5,parallel="snow",ncpus=2)
print(mcm.mod2d)
plot(mcm.mod2d)

###

mu=1;sigma=0.5;theta=2
x0=0;y0=0;init=c(x0,y0)
f <- expression(1/mu*(theta-x), x)  
g <- expression(sqrt(sigma),0)
OUI <- snssde2d(drift=f,diffusion=g,M=10,Dt=0.015,x0=c(x=0,y=0),type="str")
## true values of first and second moment at time 10
Ex <- function(t) theta+(x0-theta)*exp(-t/mu)
Vx <- function(t) 0.5*sigma*mu *(1-exp(-2*(t/mu)))
Ey <- function(t) y0+theta*t+(x0-theta)*mu*(1-exp(-t/mu))
Vy <- function(t) sigma*mu^3*((t/mu)-2*(1-exp(-t/mu))+0.5*(1-exp(-2*(t/mu))))
covxy <- function(t) 0.5*sigma*mu^2 *(1-2*exp(-t/mu)+exp(-2*(t/mu)))
tvalue = list(m1=Ex(10),m2=Ey(10),S1=Vx(10),S2=Vy(10),C12=covxy(10))
## function of the statistic(s) of interest.
sde.fun2d <- function(data, i){
  d <- data[i,]
  return(c(mean(d$x),mean(d$y),var(d$x),var(d$y),cov(d$x,d$y)))
}
## Parallel Monte-Carlo of 'OUI' at time 10
mcm.mod2d = MCM.sde(OUI,statistic=sde.fun2d,time=10,R=5,exact=tvalue,parallel="snow",ncpus=2)
print(mcm.mod2d)
plot(mcm.mod2d)


mu=1;sigma=0.5;theta=2
x0=0;y0=0;init=c(x0,y0)
f <- expression(1/mu*(theta-x), x)  
g <- expression(sqrt(sigma),0)
Sigma <- matrix(c(1, 0.75, 0.75, 1), nrow = 2, ncol = 2)
OUI <- snssde2d(drift=f,diffusion=g,corr=Sigma,M=10,Dt=0.015,x0=c(x=0,y=0),type="str")

## function of the statistic(s) of interest.
sde.fun2d <- function(data, i){
  d <- data[i,]
  return(c(mean(d$x),mean(d$y),var(d$x),var(d$y),cov(d$x,d$y)))
}
## Parallel Monte-Carlo of 'OUI' at time 10
mcm.mod2d = MCM.sde(OUI,statistic=sde.fun2d,time=10,R=5)
mcm.mod2d = MCM.sde(OUI,statistic=sde.fun2d,time=10,R=5,parallel="snow",cl= parallel::makeCluster(getOption("cl.cores", 2)),ncpus=2)
mcm.mod2d = MCM.sde(OUI,statistic=sde.fun2d,time=10,R=5,parallel="snow",ncpus=2)
print(mcm.mod2d)
plot(mcm.mod2d)

## ------------------------------------------------------------------------
mu=0.5;sigma=0.25
fx <- expression(mu*y,0,0) 
gx <- expression(sigma*z,1,1)
modtra <- snssde3d(drift=fx,diffusion=gx,M=10)
## function of the statistic(s) of interest.
sde.fun3d <- function(data, i){
  d <- data[i,]
  return(c(mean(d$x),median(d$x),Mode(d$x),var(d$x),cov(d$x,d$y),cov(d$x,d$z)))
}
## Monte-Carlo at time = 10
mcm.mod3d = MCM.sde(modtra,statistic=sde.fun3d,R=5)
mcm.mod3d = MCM.sde(modtra,statistic=sde.fun3d,R=5,parallel="snow",cl= parallel::makeCluster(getOption("cl.cores", 2)),ncpus=2)
mcm.mod3d = MCM.sde(modtra,statistic=sde.fun3d,R=5,parallel="snow",ncpus=2)
print(mcm.mod3d)
plot(mcm.mod3d)

mu=0.5;sigma=0.25
fx <- expression(mu*y,0,0) 
gx <- expression(sigma*z,1,1)
Sigma <- matrix(c(1,-0.5,-0.25,-0.5,1,0.95,-0.25,0.95,1),nrow=3,ncol=3) 
modtra <- snssde3d(drift=fx,diffusion=gx,M=10,corr=Sigma)
## function of the statistic(s) of interest.
sde.fun3d <- function(data, i){
  d <- data[i,]
  return(c(mean(d$x),median(d$x),Mode(d$x),var(d$x),cov(d$x,d$y),cov(d$x,d$z)))
}
## Monte-Carlo at time = 10
mcm.mod3d = MCM.sde(modtra,statistic=sde.fun3d,R=5)
mcm.mod3d = MCM.sde(modtra,statistic=sde.fun3d,R=5,parallel="snow",cl= parallel::makeCluster(getOption("cl.cores", 2)),ncpus=2)
mcm.mod3d = MCM.sde(modtra,statistic=sde.fun3d,R=5,parallel="snow",ncpus=2)
print(mcm.mod3d)
plot(mcm.mod3d)

##

mu=0.5;sigma=0.25
fx <- expression(mu*y,0,0) 
gx <- expression(sigma*z,1,1)
modtra <- snssde3d(drift=fx,diffusion=gx,M=10,type="str")
## function of the statistic(s) of interest.
sde.fun3d <- function(data, i){
  d <- data[i,]
  return(c(mean(d$x),median(d$x),Mode(d$x),var(d$x),cov(d$x,d$y),cov(d$x,d$z)))
}
## Monte-Carlo at time = 10
mcm.mod3d = MCM.sde(modtra,statistic=sde.fun3d,R=5,parallel="snow",ncpus=2)
print(mcm.mod3d)
plot(mcm.mod3d)
plot(mcm.mod3d,index=2)

mu=0.5;sigma=0.25
fx <- expression(mu*y,0,0) 
gx <- expression(sigma*z,1,1)
Sigma <- matrix(c(1,-0.5,-0.25,-0.5,1,0.95,-0.25,0.95,1),nrow=3,ncol=3) 
modtra <- snssde3d(drift=fx,diffusion=gx,M=10,corr=Sigma,type="str")
## function of the statistic(s) of interest.
sde.fun3d <- function(data, i){
  d <- data[i,]
  return(c(mean(d$x),median(d$x),Mode(d$x),var(d$x),cov(d$x,d$y),cov(d$x,d$z)))
}
## Monte-Carlo at time = 10
mcm.mod3d = MCM.sde(modtra,statistic=sde.fun3d,R=5)
mcm.mod3d = MCM.sde(modtra,statistic=sde.fun3d,R=5,parallel="snow",cl= parallel::makeCluster(getOption("cl.cores", 2)),ncpus=2)
mcm.mod3d = MCM.sde(modtra,statistic=sde.fun3d,R=5,parallel="snow",ncpus=2)
print(mcm.mod3d)
plot(mcm.mod3d)
