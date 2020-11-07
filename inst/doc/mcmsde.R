## ----setup, echo = F, message = F, results = 'hide',screenshot.force=FALSE----
library(Sim.DiffProc)
library(knitr)
library(deSolve)
knitr::opts_chunk$set(comment="", prompt=TRUE, fig.show='hold',warning=FALSE, message=FALSE)
options(prompt="R> ",scipen=16,digits=5,warning=FALSE, message=FALSE,mc.cores=2)

## ----eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
#  MCM.sde(model, statistic, R = 1000, time, exact = NULL, names = NULL,level = 0.95,
#          parallel = c("no", "multicore", "snow"),ncpus = getOption("ncpus", 1L), cl = NULL, ...)

## ----eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
#  plot(x,index = 1,type=c("all","hist","qqplot","boxplot","CI"), ...)

## -------------------------------------------------------------------
set.seed(1234, kind = "L'Ecuyer-CMRG")
theta = 0.75; x0 = 1
fx <- expression( 0.5*theta^2*x )
gx <- expression( theta*x )
mod1 <- snssde1d(drift=fx,diffusion=gx,x0=x0,M=500,type="ito")
mod2 <- snssde1d(drift=fx,diffusion=gx,x0=x0,M=500,type="str")
## True values of means and variance for mod1 and mod2
E.mod1 <- function(t) x0 * exp(0.5 * theta^2 * t)
V.mod1 <- function(t) x0^2 * exp(theta^2 * t) * (exp(theta^2 * t) - 1)
E.mod2 <- function(t) x0 * exp(theta^2 * t)
V.mod2 <- function(t) x0^2 * exp(2 * theta^2 * t) * (exp(theta^2 * t) - 1)
## function of the statistic(s) of interest.
sde.fun1d <- function(data, i){
     d <- data[i, ]
     return(c(mean(d),var(d)))
}
# Parallel MOnte Carlo for mod1
mcm.mod1 = MCM.sde(model=mod1,statistic=sde.fun1d,R=20, exact=list(m=E.mod1(1),S=V.mod1(1)),parallel="snow",ncpus=2)
mcm.mod1
# Parallel MOnte Carlo for mod2
mcm.mod2 = MCM.sde(model=mod2,statistic=sde.fun1d,R=20, exact=list(m=E.mod2(1),S=V.mod2(1)),parallel="snow",ncpus=2)
mcm.mod2

## -------------------------------------------------------------------
set.seed(1234, kind = "L'Ecuyer-CMRG")
mu=1;sigma=0.5;theta=2
x0=0;y0=0;init=c(x0,y0)
f <- expression(1/mu*(theta-x), x)  
g <- expression(sqrt(sigma),0)
OUI <- snssde2d(drift=f,diffusion=g,M=500,Dt=0.015,x0=c(x=0,y=0))
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
mcm.mod2d = MCM.sde(OUI,statistic=sde.fun2d,time=10,R=20,exact=tvalue,parallel="snow",ncpus=2)
mcm.mod2d

## -------------------------------------------------------------------
set.seed(1234, kind = "L'Ecuyer-CMRG")
mu=0.5;sigma=0.25
fx <- expression(mu*y,0,0) 
gx <- expression(sigma*z,1,1)
Sigma <-matrix(c(1,0.3,-0.5,0.3,1,0.2,-0.5,0.2,1),nrow=3,ncol=3)
modtra <- snssde3d(drift=fx,diffusion=gx,M=500,type="str",corr=Sigma)
## function of the statistic(s) of interest.
sde.fun3d <- function(data, i){
  d <- data[i,]
  return(c(mean(d$x),median(d$x),Mode(d$x)))
}
## Monte-Carlo at time = 10
mcm.mod3d = MCM.sde(modtra,statistic=sde.fun3d,R=10,parallel="snow",ncpus=2)
mcm.mod3d

## ----eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
#  MEM.sde(drift, diffusion, corr = NULL, type = c("ito", "str"), solve = FALSE, parms = NULL, init = NULL, time = NULL, ...)

## -------------------------------------------------------------------
fx <- expression( 0.5*theta^2*x )
gx <- expression( theta*x )
start = c(m=1,S=0)
t = seq(0,1,by=0.001)
mem.mod1 = MEM.sde(drift=fx,diffusion=gx,type="ito",solve = TRUE,parms = c(theta=0.75), init = start, time = t)
mem.mod1
mem.mod2 = MEM.sde(drift=fx,diffusion=gx,type="str",solve = TRUE,parms = c(theta=0.75), init = start, time = t)
mem.mod2

## ----eval=FALSE, include=TRUE---------------------------------------
#  plot(mem.mod1$sol.ode, mem.mod2$sol.ode,ylab = c("m(t)"),select="m", xlab = "Time",main="",col = 2:3,lty=1)
#  legend("topleft",c(expression(m[mod1](t),m[mod2](t))),inset = .05, col=2:3,lty=1)
#  plot(mem.mod1$sol.ode, mem.mod2$sol.ode,ylab = c("S(t)"),select="S", xlab = "Time",main="",col = 2:3,lty=1)
#  legend("topleft",c(expression(S[mod1](t),S[mod2](t))),inset = .05, col=2:3,lty=1)

## -------------------------------------------------------------------
fx <- expression(1/mu*(theta-x), x)  
gx <- expression(sqrt(sigma),0)
start = c(m1=0,m2=0,S1=0,S2=0,C12=0)
t = seq(0,10,by=0.001)
mem.mod2d = MEM.sde(drift=fx,diffusion=gx,type="ito",solve = TRUE,parms = c(mu=1,sigma=0.5,theta=2), init = start, time = t)
mem.mod2d

## ----eval=FALSE, include=TRUE---------------------------------------
#  matplot.0D(mem.mod2d$sol.ode,main="")

## -------------------------------------------------------------------
fx <- expression(mu*y,0,0) 
gx <- expression(sigma*z,1,1)
RHO <- expression(0.75,0.5,-0.25)
start = c(m1=5,m2=0,m3=0,S1=0,S2=0,S3=0,C12=0,C13=0,C23=0)
t = seq(0,1,by=0.001)
mem.mod3d = MEM.sde(drift=fx,diffusion=gx,corr=RHO,type="ito",solve = TRUE,parms = c(mu=0.5,sigma=0.25), init = start, time = t)
mem.mod3d

## ----eval=FALSE, include=TRUE---------------------------------------
#  matplot.0D(mem.mod3d$sol.ode,main="",select=c("m1","m2","m3"))
#  matplot.0D(mem.mod3d$sol.ode,main="",select=c("S1","S2","S3"))
#  matplot.0D(mem.mod3d$sol.ode,main="",select=c("C12","C13","C23"))

