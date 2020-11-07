## ----setup, echo = F, message = F, results = 'hide',screenshot.force=FALSE----
library(Sim.DiffProc)
library(knitr)
knitr::opts_chunk$set(comment="",prompt=TRUE, fig.show='hold', warning=FALSE, message=FALSE)
options(prompt="R> ",scipen=16,digits=5,warning=FALSE, message=FALSE,
        width = 70)

## -------------------------------------------------------------------
set.seed(1234, kind = "L'Ecuyer-CMRG")
theta = 0.5
f <- expression( (0.5*theta^2*x) )
g <- expression( theta*x )
mod1 <- snssde1d(drift=f,diffusion=g,x0=10,M=1000,type="ito") # Using Ito
mod2 <- snssde1d(drift=f,diffusion=g,x0=10,M=1000,type="str") # Using Stratonovich 
mod1
mod2

## -------------------------------------------------------------------
summary(mod1, at = 1)
summary(mod2, at = 1)

## -------------------------------------------------------------------
x1 <- rsde1d(object = mod1, at = 1)  # X(t=1) | X(0)=x0 (Ito SDE)
x2 <- rsde1d(object = mod2, at = 1)  # X(t=1) | X(0)=x0 (Stratonovich SDE)
head(data.frame(x1,x2),n=5)

## ----01,fig.env='figure*', fig.cap=' ',eval=FALSE, include=TRUE-----
#  mu1 = log(10); sigma1= sqrt(theta^2)  # log mean and log variance for mod1
#  mu2 = log(10)-0.5*theta^2 ; sigma2 = sqrt(theta^2) # log mean and log variance for mod2
#  AppdensI <- dsde1d(mod1, at = 1)
#  AppdensS <- dsde1d(mod2, at = 1)
#  plot(AppdensI , dens = function(x) dlnorm(x,meanlog=mu1,sdlog = sigma1))
#  plot(AppdensS , dens = function(x) dlnorm(x,meanlog=mu2,sdlog = sigma2))

## ----001, echo=FALSE, fig.cap='Approximate transitional density for $X_{t}|X_{0}$ at time $t-s=1$ with log-normal curves. mod1: Ito and mod2: Stratonovich  ', fig.env='figure*',fig.width=10,fig.height=10----
knitr::include_graphics(c("Figures/fig007.png","Figures/fig008.png"))

## ----02,fig.env='figure*', fig.cap=' ',eval=FALSE, include=TRUE-----
#  ## Ito
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

## ----100, echo=FALSE, fig.cap='mod1: Ito and mod2: Stratonovich  ', fig.env='figure*',fig.width=7,fig.height=7----
knitr::include_graphics(c("Figures/fig07.png","Figures/fig08.png"))

## -------------------------------------------------------------------
set.seed(1234, kind = "L'Ecuyer-CMRG")
x0=5;y0=0
mu=3;sigma=0.5
fx <- expression(-(x/mu),x)  
gx <- expression(sqrt(sigma),0)
mod2d <- snssde2d(drift=fx,diffusion=gx,Dt=0.01,M=1000,x0=c(x0,y0),method="smilstein")
mod2d

## ----eval=FALSE, include=TRUE---------------------------------------
#  summary(mod2d, at = 10)

## ----03,fig.env='figure*', fig.cap=' ',eval=FALSE, include=TRUE-----
#  ## in time
#  plot(mod2d)
#  ## in plane (O,X,Y)
#  plot2d(mod2d,type="n")
#  points2d(mod2d,col=rgb(0,100,0,50,maxColorValue=255), pch=16)

## ----102, echo=FALSE, fig.cap=' Ornstein-Uhlenbeck process and its integral ', fig.env='figure*',fig.width=10,fig.height=10----
knitr::include_graphics(c("Figures/fig09.png","Figures/fig009.png"))

## -------------------------------------------------------------------
out <- rsde2d(object = mod2d, at = 10) 
head(out,n=3)

## ----04,fig.env='figure*', fig.cap=' ',eval=FALSE, include=TRUE-----
#  ## the marginal density
#  denM <- dsde2d(mod2d,pdf="M",at =10)
#  plot(denM, main="Marginal Density")
#  ## the Joint density
#  denJ <- dsde2d(mod2d, pdf="J", n=100,at =10)
#  plot(denJ,display="contour",main="Bivariate Transition Density at time t=10")

## ----1002, echo=FALSE, fig.cap='Marginal and Joint density at time t=10 ', fig.env='figure*',fig.width=10,fig.height=10----
knitr::include_graphics(c("Figures/fig1001.png","Figures/fig1002.png"))

## ----07,fig.env='figure*', fig.cap=' ',eval=FALSE, include=TRUE-----
#  plot(denJ,display="persp",main="Bivariate Transition Density at time t=10")

## ----10002, echo=FALSE, fig.cap='Marginal and Joint density at time t=10 ', fig.env='figure*',fig.width=10,fig.height=10----
knitr::include_graphics(c("Figures/fig1003.png"))

## ----eval=FALSE, include=TRUE---------------------------------------
#  for (i in seq(1,10,by=0.005)){
#  plot(dsde2d(mod2d, at = i,n=100),display="contour",main=paste0('Transition Density \n t = ',i))
#  }

## -------------------------------------------------------------------
set.seed(1234, kind = "L'Ecuyer-CMRG")
mu = 4; sigma=0.1
fx <- expression( y ,  (mu*( 1-x^2 )* y - x)) 
gx <- expression( 0 ,2*sigma)
mod2d <- snssde2d(drift=fx,diffusion=gx,N=10000,Dt=0.01,type="str",method="rk1")

## ----9,fig.env='figure*', fig.cap=' ',eval=FALSE, include=TRUE------
#  plot(mod2d,ylim=c(-8,8))   ## back in time
#  plot2d(mod2d)              ## in plane (O,X,Y)

## ----100002, echo=FALSE, fig.cap='The stochastic Van-der-Pol equation', fig.env='figure*',fig.width=10,fig.height=10----
knitr::include_graphics(c("Figures/fig1004.png","Figures/fig1005.png"))

## -------------------------------------------------------------------
set.seed(1234, kind = "L'Ecuyer-CMRG")
mu = 1.2; sigma=0.1;nu=2;theta=0.5
fx <- expression( mu*x ,nu*(theta-y)) 
gx <- expression( x*sqrt(y) ,sigma*sqrt(y))
Sigma <- matrix(c(1,0.3,0.3,2),nrow=2,ncol=2) # correlation matrix
HM <- snssde2d(drift=fx,diffusion=gx,Dt=0.001,x0=c(100,1),corr=Sigma,M=1000)
HM

## -------------------------------------------------------------------
out <- rsde2d(object = HM, at = 1) 
head(out,n=3)

## ----eval=FALSE, include=TRUE---------------------------------------
#  denJ <- dsde2d(HM,pdf="J",at =1,lims=c(-100,900,0.4,0.75))
#  plot(denJ,display="contour",main="Bivariate Transition Density at time t=10")
#  plot(denJ,display="persp",main="Bivariate Transition Density at time t=10")

## -------------------------------------------------------------------
set.seed(1234, kind = "L'Ecuyer-CMRG")
fx <- expression(4*(-1-x)*y , 4*(1-y)*x , 4*(1-z)*y) 
gx <- rep(expression(0.2),3)
mod3d <- snssde3d(x0=c(x=2,y=-2,z=-2),drift=fx,diffusion=gx,M=1000)
mod3d

## ----eval=FALSE, include=TRUE---------------------------------------
#  s = 1
#  mean(mod3d, at = s)
#  moment(mod3d, at = s , center = TRUE , order = 2) ## variance
#  Median(mod3d, at = s)
#  Mode(mod3d, at = s)
#  quantile(mod3d , at = s)
#  kurtosis(mod3d , at = s)
#  skewness(mod3d , at = s)
#  cv(mod3d , at = s )
#  min(mod3d , at = s)
#  max(mod3d , at = s)
#  moment(mod3d, at = s , center= TRUE , order = 4)
#  moment(mod3d, at = s , center= FALSE , order = 4)

## ----eval=FALSE, include=TRUE---------------------------------------
#  summary(mod3d, at = t)

## ----10,fig.env='figure*', fig.cap=' ',eval=FALSE, include=TRUE-----
#  plot(mod3d,union = TRUE)         ## back in time
#  plot3D(mod3d,display="persp")    ## in space (O,X,Y,Z)

## ----103, echo=FALSE, fig.cap=' Flow of $1000$ trajectories of $(X_t ,Y_t ,Z_t)$ ', fig.env='figure*',fig.width=7,fig.height=7----
knitr::include_graphics(c("Figures/fig10.png","Figures/fig11.png"))

## -------------------------------------------------------------------
out <- rsde3d(object = mod3d, at = 1) 
head(out,n=3)

## ----11,fig.env='figure*', fig.cap=' Marginal density of $X_t$, $Y_t$ and $Z_t$ at time $t=1$ ',fig.width=3.5,fig.height=3.5----
den <- dsde3d(mod3d,pdf="M",at =1)
plot(den, main="Marginal Density") 

## ----111,fig.env='figure*', fig.cap='  ',eval=FALSE, include=TRUE,fig.width=5,fig.height=5----
#  denJ <- dsde3d(mod3d,pdf="J")
#  plot(denJ,display="rgl")

## -------------------------------------------------------------------
set.seed(1234, kind = "L'Ecuyer-CMRG")
K = 4; s = 1; sigma = 0.2
fx <- expression( (-K*x/sqrt(x^2+y^2+z^2)) , (-K*y/sqrt(x^2+y^2+z^2)) , (-K*z/sqrt(x^2+y^2+z^2)) ) 
gx <- rep(expression(sigma),3)
mod3d <- snssde3d(drift=fx,diffusion=gx,N=10000,x0=c(x=1,y=1,z=1))

## ----12,fig.env='figure*',fig.width=3.5,fig.height=3.5, fig.cap=' Attractive model for 3D diffusion processes '----
plot3D(mod3d,display="persp",col="blue")

## -------------------------------------------------------------------
set.seed(1234, kind = "L'Ecuyer-CMRG")
fx <- expression(y,0,0) 
gx <- expression(z,1,1)
Sigma <-matrix(c(1,0.2,0.5,0.2,1,-0.7,0.5,-0.7,1),nrow=3,ncol=3)
modtra <- snssde3d(drift=fx,diffusion=gx,M=1000,corr=Sigma)
modtra

## ----14, fig.cap='  ', fig.env='figure*',eval=FALSE, include=TRUE----
#  X <- rsde3d(modtra,at=1)$x
#  MASS::truehist(X,xlab = expression(X[t==1]));box()
#  lines(density(X),col="red",lwd=2)
#  legend("topleft",c("Distribution histogram","Kernel Density"),inset =.01,pch=c(15,NA),lty=c(NA,1),col=c("cyan","red"), lwd=2,cex=0.8)
#  ## Cov-Matrix
#  color.palette=colorRampPalette(c('white','green','blue','red'))
#  filled.contour(time(modtra), time(modtra), cov(t(modtra$X)), color.palette=color.palette,plot.title = title(main = expression(paste("Covariance empirique:",cov(X[s],X[t]))),xlab = "time", ylab = "time"),key.title = title(main = ""))

## ----1000002, echo=FALSE, fig.cap='The histogram and kernel density of $X_t$ at time $t=1$. Emprical variance-covariance matrix', fig.env='figure*',fig.width=10,fig.height=10----
knitr::include_graphics(c("Figures/fig1007.png","Figures/fig1006.png"))

