## ----setup, echo = F, message = F, results = 'hide',screenshot.force=FALSE----
library(Sim.DiffProc)
library(knitr)
knitr::opts_chunk$set(comment="",prompt=TRUE, fig.show='hold', warning=FALSE, message=FALSE)
options(prompt="R> ",scipen=16,digits=5,warning=FALSE, message=FALSE,
        width = 70)

## -------------------------------------------------------------------
set.seed(12345, kind = "L'Ecuyer-CMRG")
f <- expression( (1+2*x) ) ; g <- expression( 0.5*x^0.3 )
sim    <- snssde1d(drift=f,diffusion=g,x0=2,N=10^4,Dt=10^-4)
mydata <- sim$X

## ---- message=FALSE, warning=FALSE----------------------------------
fx <- expression( theta[1]+theta[2]*x ) ## drift coefficient of model
gx <- expression( theta[3]*x^theta[4] ) ## diffusion coefficient of model 
fitmod <- fitsde(data = mydata, drift = fx, diffusion = gx, start = list(theta1=1, theta2=1,
                theta3=1,theta4=1),pmle="euler")
fitmod

## -------------------------------------------------------------------
coef(fitmod)

## -------------------------------------------------------------------
summary(fitmod)

## -------------------------------------------------------------------
vcov(fitmod)
logLik(fitmod)
AIC(fitmod)
BIC(fitmod)

## -------------------------------------------------------------------
confint(fitmod, level=0.95)

## -------------------------------------------------------------------
set.seed(1234, kind = "L'Ecuyer-CMRG")
f <- expression( 3*(2-x) ) ; g <- expression( 0.5 )
sim <- snssde1d(drift=f,diffusion=g,x0=5,Dt=0.01)
HWV <- sim$X

## ---- message=FALSE, warning=FALSE----------------------------------
fx <- expression( theta[1]*(theta[2]- x) ) ## drift coefficient of model 
gx <- expression( theta[3] )           ## diffusion coefficient of model 
fitmod <- fitsde(data=HWV,drift=fx,diffusion=gx,start = list(theta1=1,theta2=1,
                  theta3=1),pmle="ozaki")
summary(fitmod)


## -------------------------------------------------------------------
confint(fitmod,parm=c("theta1","theta2"),level=0.95)

## -------------------------------------------------------------------
set.seed(1234, kind = "L'Ecuyer-CMRG")
f <- expression(-2*x*t) ; g <- expression(0.2*x)
sim <- snssde1d(drift=f,diffusion=g,N=1000,Dt=0.001,x0=10)
mydata <- sim$X

## ---- message=FALSE, warning=FALSE----------------------------------
fx <- expression( theta[1]*x*t ) ## drift coefficient of model 
gx <- expression( theta[2]*x )   ## diffusion coefficient of model 
fitmod <- fitsde(data=mydata,drift=fx,diffusion=gx,start = list(theta1=1,
                 theta2=1),pmle="shoji",lower=c(-3,0),upper=c(-1,1))
summary(fitmod)

## -------------------------------------------------------------------
set.seed(1234, kind = "L'Ecuyer-CMRG")
f <- expression(3*t*(sqrt(t)-x)) ; g <- expression(0.3*t)
sim <- snssde1d(drift=f,diffusion=g,M=1,N=1000,x0=2,Dt=0.001)
mydata <- sim$X

## ---- message=FALSE, warning=FALSE----------------------------------
fx <- expression( theta[1]*t* ( theta[2]*sqrt(t) - x ) ) ## drift coefficient of model 
gx <- expression( theta[3]*t ) ## diffusion coefficient of model 
fitmod <- fitsde(data=mydata,drift=fx,diffusion=gx,start = list(theta1=1,
                  theta2=1,theta3=1),pmle="kessler")
summary(fitmod)

## -------------------------------------------------------------------
set.seed(1234, kind = "L'Ecuyer-CMRG")
f <- expression( 2*x )
g <- expression( 0.3*x^0.5 )
sim <- snssde1d(drift=f,diffusion=g,M=1,N=10000,x0=2,Dt=0.0001)
mydata <- sim$X

## ---- message=FALSE, warning=FALSE----------------------------------
## True model
fx <- expression( theta[1]*x )
gx <- expression( theta[2]*x^theta[3] )
truemod <- fitsde(data=mydata,drift=fx,diffusion=gx,start = list(theta1=1,
                   theta2=1,theta3=1),pmle="euler")
## competing model 1
fx1 <- expression( theta[1]+theta[2]*x )
gx1 <- expression( theta[3]*x^theta[4] )
mod1 <- fitsde(data=mydata,drift=fx1,diffusion=gx1,start = list(theta1=1,
          theta2=1,theta3=1,theta4=1),pmle="euler")
## competing model 2
fx2 <- expression( theta[1]+theta[2]*x )
gx2 <- expression( theta[3]*sqrt(x) )
mod2 <- fitsde(data=mydata,drift=fx2,diffusion=gx2,start = list(theta1=1,
           theta2=1,theta3=1),pmle="euler")
## competing model 3
fx3 <- expression( theta[1] )
gx3 <- expression( theta[2]*x^theta[3] )
mod3 <- fitsde(data=mydata,drift=fx3,diffusion=gx3,start = list(theta1=1,
           theta2=1,theta3=1),pmle="euler")
## Computes AIC
AIC <- c(AIC(truemod),AIC(mod1),AIC(mod2),AIC(mod3))
Test <- data.frame(AIC,row.names = c("True mod","Comp mod1","Comp mod2","Comp mod3"))
Bestmod <- rownames(Test)[which.min(Test[,1])]
Bestmod

## -------------------------------------------------------------------
Theta1 <- c(coef(truemod)[[1]],coef(mod1)[[1]],coef(mod2)[[1]],coef(mod3)[[1]])
Theta2 <- c(coef(truemod)[[2]],coef(mod1)[[2]],coef(mod2)[[2]],coef(mod3)[[2]])
Theta3 <- c(coef(truemod)[[3]],coef(mod1)[[3]],coef(mod2)[[3]],coef(mod3)[[3]])
Theta4 <- c("",round(coef(mod1)[[4]],5),"","")
Parms  <- data.frame(Theta1,Theta2,Theta3,Theta4,row.names = c("True mod",
                      "Comp mod1","Comp mod2","Comp mod3"))
Parms

## ----01,fig.env='figure*', fig.cap=' The U.S. Interest Rates monthly form $06/1964$ to $12/1989$ ',fig.width=6,fig.height=4----
data(Irates)
rates <- Irates[, "r1"]
X <- window(rates, start = 1964.471, end = 1989.333)
plot(X)

## ---- message=FALSE, warning=FALSE----------------------------------
fx <- expression( theta[1]+theta[2]*x ) ## drift coefficient of CKLS model
gx <- expression( theta[3]*x^theta[4] ) ## diffusion coefficient of CKLS model
pmle <- eval(formals(fitsde.default)$pmle)
fitres <- lapply(1:4, function(i) fitsde(X,drift=fx,diffusion=gx,pmle=pmle[i],
                  start = list(theta1=1,theta2=1,theta3=1,theta4=1)))
Coef <- data.frame(do.call("cbind",lapply(1:4,function(i) coef(fitres[[i]]))))
Info <- data.frame(do.call("rbind",lapply(1:4,function(i) logLik(fitres[[i]]))),
                   do.call("rbind",lapply(1:4,function(i) AIC(fitres[[i]]))),
                   do.call("rbind",lapply(1:4,function(i) BIC(fitres[[i]]))),
                   row.names=pmle)
names(Coef) <- c(pmle)
names(Info) <- c("logLik","AIC","BIC")
Coef
Info

## ----02,fig.env='figure*', fig.cap='The path mean of the solution of the CKLS model with the estimated parameters and real data ',fig.width=6,fig.height=4----
set.seed(1234, kind = "L'Ecuyer-CMRG")
f <- expression( (2.076-0.263*x) )
g <- expression( 0.130*x^1.451 )
mod <- snssde1d(drift=f,diffusion=g,x0=X[1],M=500, N=length(X),t0=1964.471, T=1989.333)
mod
plot(mod,type="n",ylim=c(0,30))
lines(X,col=4,lwd=2)
lines(time(mod),apply(mod$X,1,mean),col=2,lwd=2)
lines(time(mod),apply(mod$X,1,bconfint,level=0.95)[1,],col=5,lwd=2)
lines(time(mod),apply(mod$X,1,bconfint,level=0.95)[2,],col=5,lwd=2)
legend("topleft",c("real data","mean path",paste("bound of", 95,"% confidence")),inset = .01,col=c(4,2,5),lwd=2,cex=0.8)


