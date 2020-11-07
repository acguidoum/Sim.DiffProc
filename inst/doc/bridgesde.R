## ----setup, echo = F, message = F, results = 'hide',screenshot.force=FALSE----
library(Sim.DiffProc)
library(knitr)
knitr::opts_chunk$set(comment="",prompt=TRUE, fig.show='hold', warning=FALSE, message=FALSE)
options(prompt="R> ",scipen=16,digits=5,warning=FALSE, message=FALSE,
        width = 70)

## -------------------------------------------------------------------
set.seed(1234, kind = "L'Ecuyer-CMRG")
f <- expression((1-x)/(1-t))
g <- expression(x)
mod <- bridgesde1d(drift=f,diffusion=g,x0=3,y=1,M=1000)
mod

## ----01,fig.env='figure*', fig.cap=' ',eval=FALSE, include=TRUE-----
#  plot(mod,ylab=expression(X[t]))
#  lines(time(mod),apply(mod$X,1,mean),col=2,lwd=2)
#  legend("topleft","mean path",inset = .01,col=2,lwd=2,cex=0.8,bty="n")

## -------------------------------------------------------------------
x <- rsde1d(object = mod, at = 0.55) 
head(x, n = 3)

## ----04,fig.env='figure*', fig.cap='  ',eval=FALSE, include=TRUE----
#  dens <- dsde1d(mod, at = 0.55)
#  plot(dens,hist=TRUE) ## histgramme
#  plot(dens,add=TRUE)  ## kernel density

## ----33, echo=FALSE, fig.cap='Bridge sde 1D. Histgramme and kernel density estimation for $X_{t}|X_{0}=3,X_{T}=1$', fig.env='figure*',fig.width=7,fig.height=7----
knitr::include_graphics(c("Figures/fig03.png","Figures/fig1008.png"))

## -------------------------------------------------------------------
set.seed(1234, kind = "L'Ecuyer-CMRG")
fx <- expression(-(1+y)*x , -(1+x)*y)
gx <- expression(0.2*(1-y),0.1*(1-x))
Sigma <-matrix(c(1,0.3,0.3,1),nrow=2,ncol=2)
mod2 <- bridgesde2d(drift=fx,diffusion=gx,x0=c(1,-0.5),y=c(1,0.5),Dt=0.01,M=1000,type="str",corr=Sigma)
mod2

## ----06,fig.env='figure*', fig.cap=' ',eval=FALSE, include=TRUE-----
#  plot(mod2,col=c('#FF00004B','#0000FF82'))

## ----333, echo=FALSE, fig.cap='  Bridge sde 2D ', fig.env='figure*',fig.width=5,fig.height=5----
knitr::include_graphics("Figures/fig04.png")

## -------------------------------------------------------------------
x2 <- rsde2d(object = mod2, at = 5) 
head(x2, n = 3)

## ----09,fig.env='figure*', fig.cap='  ',eval=FALSE, include=TRUE----
#  ## Marginal
#  denM <- dsde2d(mod2,pdf="M",at = 5)
#  plot(denM, main="Marginal Density")
#  ## Joint
#  denJ <- dsde2d(mod2, pdf="J", n=100,at = 5)
#  plot(denJ,display="contour",main="Bivariate Transition Density at time t=5")

## ----103, echo=FALSE, fig.cap='The marginal and joint density of $X_{t}|X_{0}=1,X_{T}=1$ and $Y_{t}|Y_{0}=-0.5,Y_{T}=0.5$ at time $t=5$', fig.env='figure*',fig.width=7,fig.height=7----
knitr::include_graphics(c("Figures/fig1009.png","Figures/fig1010.png"))

## ----11,fig.env='figure*', fig.cap='  ',eval=FALSE, include=TRUE----
#  plot(denJ,main="Bivariate Transition Density at time t=5")

## ----33303, echo=FALSE, fig.cap='$3$D plot of the transition density of $X_{t}|X_{0}=1,X_{T}=1$ and $Y_{t}|Y_{0}=-0.5,Y_{T}=0.5$ at time $t=5$  ', fig.env='figure*',fig.width=5,fig.height=5----
knitr::include_graphics("Figures/fig1011.png")

## ----eval=FALSE, include=TRUE---------------------------------------
#  for (i in seq(1,9,by=0.005)){
#  plot(dsde2d(mod2, at = i,n=100),display="contour",main=paste0('Transition Density \n t = ',i))
#  }

## -------------------------------------------------------------------
set.seed(1234, kind = "L'Ecuyer-CMRG")
fx <- expression(-4*(1+x)*y, 4*(1-y)*x, 4*(1-z)*y)
gx <- rep(expression(0.2),3)
mod3 <- bridgesde3d(x0=c(0,-1,0.5),y=c(0,-2,0.5),drift=fx,diffusion=gx,M=1000)
mod3

## ----12,fig.env='figure*', fig.cap=' ',eval=FALSE, include=TRUE-----
#  plot(mod3) ## in time
#  plot3D(mod3,display = "persp",main="3D Bridge SDE's") ## in space

## ----3333, echo=FALSE, fig.cap=' Bridge sde 3D ', fig.env='figure*',fig.width=7,fig.height=7----
knitr::include_graphics(c("Figures/fig05.png","Figures/fig06.png"))

## -------------------------------------------------------------------
x3 <- rsde3d(object = mod3, at = 0.75) 
head(x3, n = 3)

## ----15,fig.env='figure*', fig.cap='  ',eval=FALSE, include=TRUE----
#  ## Marginal
#  denM <- dsde3d(mod3,pdf="M",at =0.75)
#  plot(denM, main="Marginal Density")
#  ## Joint
#  denJ <- dsde3d(mod3,pdf="J",at=0.75)
#  plot(denJ,display="rgl")

