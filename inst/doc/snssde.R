## ----setup, echo = F, message = F, results = 'hide',screenshot.force=FALSE----
library(Sim.DiffProc)
library(knitr)
knitr::opts_chunk$set(comment="",prompt=TRUE, fig.show='hold', warning=FALSE, message=FALSE)
options(prompt="R> ",scipen=16,digits=5,warning=FALSE, message=FALSE,
        width = 70)

## -------------------------------------------------------------------
theta = 0.5
f <- expression( (0.5*theta^2*x) )
g <- expression( theta*x )
mod1 <- snssde1d(drift=f,diffusion=g,x0=10,M=1000,type="ito") # Using It?
mod2 <- snssde1d(drift=f,diffusion=g,x0=10,M=1000,type="str") # Using Stratonovich 
mod1
mod2

## -------------------------------------------------------------------
s = 1
mean(mod1, at = s)
moment(mod1, at = s , center = TRUE , order = 2) ## variance
Median(mod1, at = s)
Mode(mod1, at =s)
quantile(mod1 , at = s)
kurtosis(mod1 , at = s)
skewness(mod1 , at = s)
cv(mod1 , at = s )
min(mod1 , at = s)
max(mod1 , at = s)
moment(mod1, at = s , center= TRUE , order = 4)
moment(mod1, at = s , center= FALSE , order = 4)

## -------------------------------------------------------------------
summary(mod1, at = 1)
summary(mod2, at = 1)

## -------------------------------------------------------------------
x1 <- rsde1d(object = mod1, at = 1)  # X(t=1) | X(0)=x0 (It? SDE)
x2 <- rsde1d(object = mod2, at = 1)  # X(t=1) | X(0)=x0 (Stratonovich SDE)
head(x1,n=10)
head(x2,n=10)
summary(data.frame(x1,x2))

## ----01,fig.env='figure*', fig.cap='  '-----------------------------
mu1 = log(10); sigma1= sqrt(theta^2)  # log mean and log variance for mod1 
mu2 = log(10)-0.5*theta^2 ; sigma2 = sqrt(theta^2) # log mean and log variance for mod2
AppdensI <- dsde1d(mod1, at = 1)
AppdensS <- dsde1d(mod2, at = 1)
plot(AppdensI , dens = function(x) dlnorm(x,meanlog=mu1,sdlog = sigma1))
plot(AppdensS , dens = function(x) dlnorm(x,meanlog=mu2,sdlog = sigma2))

## ----02,fig.env='figure*', fig.cap=' ',eval=FALSE, include=TRUE-----
#  ## It?
#  plot(mod1,ylab=expression(X^mod1))
#  lines(time(mod1),apply(mod1$X,1,mean),col=2,lwd=2)
#  lines(time(mod1),apply(mod1$X,1,bconfint,level=0.95)[1,],col=4,lwd=2)
#  lines(time(mod1),apply(mod1$X,1,bconfint,level=0.95)[2,],col=4,lwd=2)
#  legend("topleft",c("mean path",paste("bound of", 95,"% confidence")),inset = .01,col=c(2,4),lwd=2,cex=0.8)
#  ## Stratonovich
#  plot(mod2,ylab=expression(X^mod2))
#  lines(time(mod2),apply(mod2$X,1,mean),col=2,lwd=2)
#  lines(time(mod2),apply(mod2$X,1,bconfint,level=0.95)[1,],col=4,lwd=2)
#  lines(time(mod2),apply(mod2$X,1,bconfint,level=0.95)[2,],col=4,lwd=2)
#  legend("topleft",c("mean path",paste("bound of",95,"% confidence")),col=c(2,4),inset =.01,lwd=2,cex=0.8)

## ----100, echo=FALSE, fig.cap=' ', fig.env='figure*',fig.width=7,fig.height=7----
knitr::include_graphics("Figures/fig07.png")

## ----101, echo=FALSE, fig.cap=' mod1: It? and mod2: Stratonovich ', fig.env='figure*',fig.width=7,fig.height=7----
knitr::include_graphics("Figures/fig08.png")

## -------------------------------------------------------------------
x0=5;y0=0
mu=3;sigma=0.5
fx <- expression(-(x/mu),x)  
gx <- expression(sqrt(sigma),0)
mod2d <- snssde2d(drift=fx,diffusion=gx,Dt=0.01,M=1000,x0=c(x0,y0),method="smilstein")
mod2d

## -------------------------------------------------------------------
s = 5
mean(mod2d, at = s)
moment(mod2d, at = s , center = TRUE , order = 2) ## variance
Median(mod2d, at = s)
Mode(mod2d, at = s)
quantile(mod2d , at = s)
kurtosis(mod2d , at = s)
skewness(mod2d , at = s)
cv(mod2d , at = s )
min(mod2d , at = s)
max(mod2d , at = s)
moment(mod2d, at = s , center= TRUE , order = 4)
moment(mod2d, at = s , center= FALSE , order = 4)

## -------------------------------------------------------------------
summary(mod2d, at = s)

## ----03,fig.env='figure*', fig.cap=' ',eval=FALSE, include=TRUE-----
#  plot(mod2d)

## ----102, echo=FALSE, fig.cap=' Ornstein-Uhlenbeck process and its integral ', fig.env='figure*',fig.width=7,fig.height=7----
knitr::include_graphics("Figures/fig09.png")

## -------------------------------------------------------------------
out <- rsde2d(object = mod2d, at = 10) 
head(out,n=10)
summary(out)
cov(out)

## ----04,fig.env='figure*', fig.cap='  '-----------------------------
mx <- apply(mod2d$X,1,mean)
my <- apply(mod2d$Y,1,mean)
Sx <- apply(mod2d$X,1,var)
Sy <- apply(mod2d$Y,1,var)
Cxy <- sapply(1:1001,function(i) cov(mod2d$X[i,],mod2d$Y[i,]))
out_b <- data.frame(mx,my,Sx,Sy,Cxy)
matplot(time(mod2d), out_b, type = "l", xlab = "time", ylab = "",col=2:6,lwd=2,lty=2:6,las=1)
legend("topleft",c(expression(hat(E)(X[t]),hat(E)(Y[t]),hat(Var)(X[t]),hat(Var)(Y[t]),hat(Cov)(X[t],Y[t]))),inset = .05,col=2:6,lty=2:6,lwd=2,cex=0.9)

## ----05,fig.env='figure*', fig.cap='  '-----------------------------
denM <- dsde2d(mod2d,pdf="M",at =10)
denM
plot(denM, main="Marginal Density")

## ----06,fig.env='figure*', fig.cap='  '-----------------------------
denJ <- dsde2d(mod2d, pdf="J", n=100,at =10)
denJ
plot(denJ,display="contour",main="Bivariate Transition Density at time t=10")
plot(denJ,display="image",main="Bivariate Transition Density at time t=10")

## ----07,fig.env='figure*', fig.cap='  '-----------------------------
plot(denJ,main="Bivariate Transition Density at time t=10")

## ----eval=FALSE, include=TRUE---------------------------------------
#  for (i in seq(1,10,by=0.005)){
#  plot(dsde2d(mod2d, at = i,n=100),display="contour",main=paste0('Transition Density \n t = ',i))
#  }

## -------------------------------------------------------------------
mu = 4; sigma=0.1
fx <- expression( y ,  (mu*( 1-x^2 )* y - x)) 
gx <- expression( 0 ,2*sigma)
mod2d <- snssde2d(drift=fx,diffusion=gx,N=10000,Dt=0.01,type="str",method="rk1")
mod2d

## ----9,fig.env='figure*', fig.cap='  '------------------------------
plot2d(mod2d) ## in plane (O,X,Y)
plot(mod2d)   ## back in time

## -------------------------------------------------------------------
fx <- expression(4*(-1-x)*y , 4*(1-y)*x , 4*(1-z)*y) 
gx <- rep(expression(0.2),3)
mod3d <- snssde3d(x0=c(x=2,y=-2,z=-2),drift=fx,diffusion=gx,M=1000)
mod3d

## -------------------------------------------------------------------
s = 1
mean(mod3d, at = s)
moment(mod3d, at = s , center = TRUE , order = 2) ## variance
Median(mod3d, at = s)
Mode(mod3d, at = s)
quantile(mod3d , at = s)
kurtosis(mod3d , at = s)
skewness(mod3d , at = s)
cv(mod3d , at = s )
min(mod3d , at = s)
max(mod3d , at = s)
moment(mod3d, at = s , center= TRUE , order = 4)
moment(mod3d, at = s , center= FALSE , order = 4)

## -------------------------------------------------------------------
summary(mod3d, at = s)

## ----10,fig.env='figure*', fig.cap='  ',eval=FALSE, include=TRUE----
#  plot(mod3d,union = TRUE)         ## back in time
#  plot3D(mod3d,display="persp")    ## in space (O,X,Y,Z)

## ----103, echo=FALSE, fig.cap=' ', fig.env='figure*',fig.width=7,fig.height=7----
knitr::include_graphics("Figures/fig10.png")

## ----104, echo=FALSE, fig.cap=' 3D SDEs ', fig.env='figure*',fig.width=7,fig.height=7----
knitr::include_graphics("Figures/fig11.png")

## -------------------------------------------------------------------
out <- rsde3d(object = mod3d, at = s) 
head(out,n=10)
summary(out)
cov(out)

## ----11,fig.env='figure*', fig.cap='  '-----------------------------
den <- dsde3d(mod3d,pdf="M",at =1)
den
plot(den, main="Marginal Density") 

## ----111,fig.env='figure*', fig.cap='  ',eval=FALSE, include=TRUE----
#  denJ <- dsde3d(mod3d,pdf="J")
#  plot(denJ,display="rgl")

## -------------------------------------------------------------------
K = 4; s = 1; sigma = 0.2
fx <- expression( (-K*x/sqrt(x^2+y^2+z^2)) , (-K*y/sqrt(x^2+y^2+z^2)) , (-K*z/sqrt(x^2+y^2+z^2)) ) 
gx <- rep(expression(sigma),3)
mod3d <- snssde3d(drift=fx,diffusion=gx,N=10000,x0=c(x=1,y=1,z=1))
mod3d

## ----12,fig.env='figure*', fig.cap='  '-----------------------------
plot3D(mod3d,display="persp",col="blue")

## -------------------------------------------------------------------
fx <- expression(y,0,0) 
gx <- expression(z,1,1)
modtra <- snssde3d(drift=fx,diffusion=gx,M=1000,type="str")
modtra
summary(modtra)

## ----13,fig.env='figure*', fig.cap='  ',eval=FALSE, include=TRUE----
#  plot(modtra$X,plot.type="single",ylab=expression(X[t]))
#  lines(time(modtra),apply(modtra$X,1,mean),col=2,lwd=2)
#  legend("topleft",c("mean path"),col=2,lwd=2,cex=0.8)

## ----105, echo=FALSE, fig.cap=' Simulation of $X_t$ ', fig.env='figure*',fig.width=7,fig.height=7----
knitr::include_graphics("Figures/fig12.png")

## ----14, fig.cap='  ', fig.env='figure*'----------------------------
den <- dsde3d(modtra,pdf="Marginal",at=1)
den$resx
MASS::truehist(den$ech$x,xlab = expression(X[t==1]));box()
lines(den$resx,col="red",lwd=2)
legend("topleft",c("Distribution histogram","Kernel Density"),inset =.01,pch=c(15,NA),lty=c(NA,1),col=c("cyan","red"), lwd=2,cex=0.8)

## ----15, fig.cap='  ', fig.env='figure*'----------------------------
m1  <- apply(modtra$X,1,mean) ## m1(t)
S1  <- apply(modtra$X,1,var)  ## s1(t)
out_a <- data.frame(m1,S1)
matplot(time(modtra), out_a, type = "l", xlab = "time", ylab = "", col=2:3,lwd=2,lty=2:3,las=1)
legend("topleft",c(expression(m[1](t),S[1](t))),inset = .09,col=2:3,lty=2:3,lwd=2,cex=0.9)

## ----16, fig.cap='  ', fig.env='figure*'----------------------------
color.palette=colorRampPalette(c('white','green','blue','red'))
filled.contour(time(modtra), time(modtra), cov(t(modtra$X)), color.palette=color.palette,plot.title = title(main = expression(paste("Covariance empirique:",cov(X[s],X[t]))),xlab = "time", ylab = "time"),key.title = title(main = ""))

