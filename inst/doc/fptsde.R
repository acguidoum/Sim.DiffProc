## ----setup, echo = F, message = F, results = 'hide',screenshot.force=FALSE----
library(Sim.DiffProc)
library(knitr)
knitr::opts_chunk$set(comment="", prompt=TRUE, fig.show='hold',warning=FALSE, message=FALSE)
options(prompt="R> ",scipen=16,digits=5,warning=FALSE, message=FALSE,
        width = 70)

## -------------------------------------------------------------------
f <- expression( (1-0.5*x) )
g <- expression( 1 )
mod1d <- snssde1d(drift=f,diffusion=g,x0=1.7,M=1000,method="taylor")

## -------------------------------------------------------------------
St  <- expression(2*(1-sinh(0.5*t)) )
fpt1d <- fptsde1d(mod1d, boundary = St)
fpt1d
head(fpt1d$fpt, n = 5)

## ----eval=FALSE, include=TRUE---------------------------------------
#  mean(fpt1d)
#  moment(fpt1d , center = TRUE , order = 2) ## variance
#  Median(fpt1d)
#  Mode(fpt1d)
#  quantile(fpt1d)
#  kurtosis(fpt1d)
#  skewness(fpt1d)
#  cv(fpt1d)
#  min(fpt1d)
#  max(fpt1d)
#  moment(fpt1d , center= TRUE , order = 4)
#  moment(fpt1d , center= FALSE , order = 4)

## ----2,fig.env='figure*', fig.cap='  ',eval=FALSE, include=TRUE-----
#  plot(dfptsde1d(fpt1d),hist=TRUE,nbins="FD")  ## histogramm
#  plot(dfptsde1d(fpt1d))              ## kernel density

## ----eval=FALSE, include=TRUE---------------------------------------
#  require(fptdApprox)
#  x <- character(4)
#  x[1] <- "m * x"
#  x[2] <- "(sigma^2) * x^2"
#  x[3] <- "dnorm((log(x) - (log(y) + (m - sigma^2/2) * (t- s)))/(sigma * sqrt(t - s)),0,1)/(sigma * sqrt(t - s) * x)"
#  x[4] <- "plnorm(x,log(y) + (m - sigma^2/2) * (t - s),sigma * sqrt(t - s))"
#  Lognormal <- diffproc(x)
#  res1 <- Approx.fpt.density(Lognormal, 0, 10, 1, "7 + 3.2 * t + 1.4 * t * sin(1.75 * t)",list(m = 0.48,sigma = 0.07))

## -------------------------------------------------------------------
## Set the model X(t)
f <- expression( 0.48*x )
g <- expression( 0.07*x )
mod1 <- snssde1d(drift=f,diffusion=g,x0=1,T=10,M=1000)
## Set the boundary S(t)
St  <- expression( 7 + 3.2 * t + 1.4 * t * sin(1.75 * t) )
## Generate the fpt
fpt1 <- fptsde1d(mod1, boundary = St)
head(fpt1$fpt, n = 5)
summary(fpt1)

## ----3,fig.env='figure*', fig.cap='  ',eval=FALSE, include=TRUE-----
#  plot(res1$y ~ res1$x, type = 'l',main = 'Approximation First-Passage-Time Density', ylab = 'Density', xlab = expression(tau[S(t)]),cex.main = 0.95,lwd=2)
#  plot(dfptsde1d(fpt1,bw="bcv"),add=TRUE)
#  legend('topright', lty = c(1, NA), col = c(1,'#BBCCEE'),pch=c(NA,15),legend = c('Approx.fpt.density()', 'fptsde1d()'), lwd = 2, bty = 'n')

## ----33, echo=FALSE, fig.cap=' `fptsde1d()` vs `Approx.fpt.density()` ', fig.env='figure*',fig.width=7,fig.height=7----
knitr::include_graphics("Figures/fig01.png")

## ----eval=FALSE, include=TRUE---------------------------------------
#  require(DiffusionRgqd)
#  G1 <- function(t)
#       {
#   theta[1] * (10+0.2 * sin(2 * pi * t) + 0.3 * prod(sqrt(t),
#   1+cos(3 * pi * t)))
#   }
#  G2 <- function(t){-theta[1]}
#  Q2 <- function(t){0.1}
#  res2 = GQD.TIpassage(8, 12, 1, 4, 1 / 100, theta = c(0.5))

## -------------------------------------------------------------------
## Set the model X(t)
theta1=0.5
f <- expression( theta1*x*(10+0.2*sin(2*pi*t)+0.3*sqrt(t)*(1+cos(3*pi*t))-x) )
g <- expression( sqrt(0.1)*x )
mod2 <- snssde1d(drift=f,diffusion=g,x0=8,t0=1,T=4,M=1000)
## Set the boundary S(t)
St  <- expression( 12 )
## Generate the fpt
fpt2 <- fptsde1d(mod2, boundary = St)
head(fpt2$fpt, n = 5)
summary(fpt2)

## ----4,fig.env='figure*', fig.cap='  ',eval=FALSE, include=TRUE-----
#  plot(dfptsde1d(fpt2),hist=TRUE,nbins = "Scott",main = 'Approximation First-Passage-Time Density', ylab = 'Density', xlab = expression(tau[S(t)]), cex.main = 0.95)
#  lines(res2$density ~ res2$time, type = 'l',lwd=2)
#  legend('topright', lty = c(1, NA), col = c(1,'#FF00004B'),pch=c(NA,15),legend = c('GQD.TIpassage()', 'fptsde1d()'), lwd = 2, bty = 'n')

## ----44, echo=FALSE, fig.cap='`fptsde1d()` vs `GQD.TIpassage()` ', fig.env='figure*',fig.width=7,fig.height=7----
knitr::include_graphics("Figures/fig02.png")

## -------------------------------------------------------------------
fx <- expression(5*(-1-y)*x , 5*(-1-x)*y)
gx <- expression(0.5*y,0.5*x)
mod2d <- snssde2d(drift=fx,diffusion=gx,x0=c(x=1,y=-1),M=1000,type="str")

## -------------------------------------------------------------------
St <- expression(sin(2*pi*t))
fpt2d <- fptsde2d(mod2d, boundary = St)
head(fpt2d$fpt, n = 5)

## ----eval=FALSE, include=TRUE---------------------------------------
#  mean(fpt2d)
#  moment(fpt2d , center = TRUE , order = 2) ## variance
#  Median(fpt2d)
#  Mode(fpt2d)
#  quantile(fpt2d)
#  kurtosis(fpt2d)
#  skewness(fpt2d)
#  cv(fpt2d)
#  min(fpt2d)
#  max(fpt2d)
#  moment(fpt2d , center= TRUE , order = 4)
#  moment(fpt2d , center= FALSE , order = 4)

## -------------------------------------------------------------------
summary(fpt2d)

## ----6,fig.env='figure*', fig.cap='  ',eval=FALSE, include=TRUE-----
#  denM <- dfptsde2d(fpt2d, pdf = 'M')
#  plot(denM)

## ----7,fig.env='figure*', fig.cap='  ',eval=FALSE, include=TRUE-----
#  denJ <- dfptsde2d(fpt2d, pdf = 'J',n=100)
#  plot(denJ,display="contour",main="Bivariate Density of F.P.T",xlab=expression(tau[x]),ylab=expression(tau[y]))
#  plot(denJ,display="image",main="Bivariate Density of F.P.T",xlab=expression(tau[x]),ylab=expression(tau[y]))

## ----8, echo=TRUE, fig.cap='  ', fig.env='figure*', message=FALSE, warning=FALSE,eval=FALSE, include=TRUE----
#  plot(denJ,display="persp",main="Bivariate Density of F.P.T",xlab=expression(tau[x]),ylab=expression(tau[y]))

## -------------------------------------------------------------------
fx <- expression(4*(-1-x)*y , 4*(1-y)*x , 4*(1-z)*y) 
gx <- rep(expression(0.2),3)
Sigma <-matrix(c(1,0.3,-0.5,0.3,1,0.2,-0.5,0.2,1),nrow=3,ncol=3)
mod3d <- snssde3d(drift=fx,diffusion=gx,x0=c(x=2,y=-2,z=0),M=1000,corr=Sigma)

## -------------------------------------------------------------------
St <- expression(-1.5+3*t)
fpt3d <- fptsde3d(mod3d, boundary = St)
head(fpt3d$fpt, n = 5)

## ----eval=FALSE, include=TRUE---------------------------------------
#  mean(fpt3d)
#  moment(fpt3d , center = TRUE , order = 2) ## variance
#  Median(fpt3d)
#  Mode(fpt3d)
#  quantile(fpt3d)
#  kurtosis(fpt3d)
#  skewness(fpt3d)
#  cv(fpt3d)
#  min(fpt3d)
#  max(fpt3d)
#  moment(fpt3d , center= TRUE , order = 4)
#  moment(fpt3d , center= FALSE , order = 4)

## -------------------------------------------------------------------
summary(fpt3d)

## ----10,fig.env='figure*', fig.cap='  ',eval=FALSE, include=TRUE----
#  denM <- dfptsde3d(fpt3d, pdf = "M")
#  plot(denM)

## ----111,fig.env='figure*', fig.cap='  ',eval=FALSE, include=TRUE----
#  denJ <- dfptsde3d(fpt3d,pdf="J")
#  plot(denJ,display="rgl")

