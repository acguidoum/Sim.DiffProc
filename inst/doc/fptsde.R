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
head(fpt1d$fpt, n = 10)

## -------------------------------------------------------------------
mean(fpt1d)
moment(fpt1d , center = TRUE , order = 2) ## variance
Median(fpt1d)
Mode(fpt1d)
quantile(fpt1d)
kurtosis(fpt1d)
skewness(fpt1d)
cv(fpt1d)
min(fpt1d)
max(fpt1d)
moment(fpt1d , center= TRUE , order = 4)
moment(fpt1d , center= FALSE , order = 4)

## -------------------------------------------------------------------
summary(fpt1d)

## ----1 ,fig.env='figure*', fig.cap='  '-----------------------------
plot(time(mod1d),mod1d$X[,1],type="l",lty=3,ylab="X(t)",xlab="time",axes=F)
curve(2*(1-sinh(0.5*x)),add=TRUE,col=2)
points(fpt1d$fpt[1],2*(1-sinh(0.5*fpt1d$fpt[1])),pch=19,col=4,cex=0.5)
lines(c(fpt1d$fpt[1],fpt1d$fpt[1]),c(0,2*(1-sinh(0.5*fpt1d$fpt[1]))),lty=2,col=4)
axis(1, fpt1d$fpt[1], bquote(tau[S(t)]==.(fpt1d$fpt[1])),col=4,col.ticks=4)
legend('topleft',col=c(1,2,4),lty=c(1,1,NA),pch=c(NA,NA,19),legend=c(expression(X[t]),expression(S(t)),expression(tau[S(t)])),cex=0.8,bty = 'n')
box()

## ----2,fig.env='figure*', fig.cap='  '------------------------------
plot(dfptsde1d(fpt1d),hist=TRUE,nbins="FD")  ## histogramm
plot(dfptsde1d(fpt1d))              ## kernel density

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
fpt1
head(fpt1$fpt, n = 10)
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
fpt2
head(fpt2$fpt, n = 10)
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
fpt2d
head(fpt2d$fpt, n = 10)

## -------------------------------------------------------------------
mean(fpt2d)
moment(fpt2d , center = TRUE , order = 2) ## variance
Median(fpt2d)
Mode(fpt2d)
quantile(fpt2d)
kurtosis(fpt2d)
skewness(fpt2d)
cv(fpt2d)
min(fpt2d)
max(fpt2d)
moment(fpt2d , center= TRUE , order = 4)
moment(fpt2d , center= FALSE , order = 4)

## -------------------------------------------------------------------
summary(fpt2d)

## ----5 ,fig.env='figure*', fig.cap='  '-----------------------------
plot(ts.union(mod2d$X[,1],mod2d$Y[,1]),col=1:2,lty=3,plot.type="single",type="l",ylab= "",xlab="time",axes=F)
curve(sin(2*pi*x),add=TRUE,col=3)
points(fpt2d$fpt$x[1],sin(2*pi*fpt2d$fpt$x[1]),pch=19,col=4,cex=0.5)
lines(c(fpt2d$fpt$x[1],fpt2d$fpt$x[1]),c(sin(2*pi*fpt2d$fpt$x[1]),-10),lty=2,col=4)
axis(1, fpt2d$fpt$x[1], bquote(tau[X[S(t)]]==.(fpt2d$fpt$x[1])),col=4,col.ticks=4)
points(fpt2d$fpt$y[1],sin(2*pi*fpt2d$fpt$y[1]),pch=19,col=5,cex=0.5)
lines(c(fpt2d$fpt$y[1],fpt2d$fpt$y[1]),c(sin(2*pi*fpt2d$fpt$y[1]),-10),lty=2,col=5)
axis(1, fpt2d$fpt$y[1], bquote(tau[Y[S(t)]]==.(fpt2d$fpt$y[1])),col=5,col.ticks=5)
legend('topright',col=1:5,lty=c(1,1,1,NA,NA),pch=c(NA,NA,NA,19,19),legend=c(expression(X[t]),expression(Y[t]),expression(S(t)),expression(tau[X[S(t)]]),expression(tau[Y[S(t)]])),cex=0.8,inset = .01)
box()

## ----6,fig.env='figure*', fig.cap='  '------------------------------
denM <- dfptsde2d(fpt2d, pdf = 'M')
plot(denM)

## ----7,fig.env='figure*', fig.cap='  '------------------------------
denJ <- dfptsde2d(fpt2d, pdf = 'J',n=100)
plot(denJ,display="contour",main="Bivariate Density of F.P.T",xlab=expression(tau[x]),ylab=expression(tau[y]))
plot(denJ,display="image",main="Bivariate Density of F.P.T",xlab=expression(tau[x]),ylab=expression(tau[y]))

## ----8, echo=TRUE, fig.cap='  ', fig.env='figure*', message=FALSE, warning=FALSE----
plot(denJ,display="persp",main="Bivariate Density of F.P.T",xlab=expression(tau[x]),ylab=expression(tau[y]))

## -------------------------------------------------------------------
fx <- expression(4*(-1-x)*y , 4*(1-y)*x , 4*(1-z)*y) 
gx <- rep(expression(0.2),3)
mod3d <- snssde3d(drift=fx,diffusion=gx,x0=c(x=2,y=-2,z=0),M=1000)

## -------------------------------------------------------------------
St <- expression(-1.5+3*t)
fpt3d <- fptsde3d(mod3d, boundary = St)
fpt3d
head(fpt3d$fpt, n = 10)

## -------------------------------------------------------------------
mean(fpt3d)
moment(fpt3d , center = TRUE , order = 2) ## variance
Median(fpt3d)
Mode(fpt3d)
quantile(fpt3d)
kurtosis(fpt3d)
skewness(fpt3d)
cv(fpt3d)
min(fpt3d)
max(fpt3d)
moment(fpt3d , center= TRUE , order = 4)
moment(fpt3d , center= FALSE , order = 4)

## -------------------------------------------------------------------
summary(fpt3d)

## ----9 ,fig.env='figure*', fig.cap='  '-----------------------------
plot(ts.union(mod3d$X[,1],mod3d$Y[,1],mod3d$Z[,1]),col=1:3,lty=3,plot.type="single",type="l",ylab="",xlab="time",axes=F)
curve(-1.5+3*x,add=TRUE,col=4)
points(fpt3d$fpt$x[1],-1.5+3*fpt3d$fpt$x[1],pch=19,col=5,cex=0.5)
lines(c(fpt3d$fpt$x[1],fpt3d$fpt$x[1]),c(-1.5+3*fpt3d$fpt$x[1],-10),lty=2,col=5)
axis(1, fpt3d$fpt$x[1], bquote(tau[X[S(t)]]==.(fpt3d$fpt$x[1])),col=5,col.ticks=5)
points(fpt3d$fpt$y[1],-1.5+3*fpt3d$fpt$y[1],pch=19,col=6,cex=0.5)
lines(c(fpt3d$fpt$y[1],fpt3d$fpt$y[1]),c(-1.5+3*fpt3d$fpt$y[1],-10),lty=2,col=6)
axis(1, fpt3d$fpt$y[1], bquote(tau[Y[S(t)]]==.(fpt3d$fpt$y[1])),col=6,col.ticks=6)
points(fpt3d$fpt$z[1],-1.5+3*fpt3d$fpt$z[1],pch=19,col=7,cex=0.5)
lines(c(fpt3d$fpt$z[1],fpt3d$fpt$z[1]),c(-1.5+3*fpt3d$fpt$z[1],-10),lty=2,col=7)
axis(1, fpt3d$fpt$z[1], bquote(tau[Z[S(t)]]==.(fpt3d$fpt$z[1])),col=7,col.ticks=7)
legend('topright',col=1:7,lty=c(1,1,1,1,NA,NA,NA),pch=c(NA,NA,NA,NA,19,19,19),legend=c(expression(X[t]),expression(Y[t]),expression(Z[t]),expression(S(t)),expression(tau[X[S(t)]]),expression(tau[Y[S(t)]]),expression(tau[Z[S(t)]])),cex=0.8,inset = .01)
box()

## ----10,fig.env='figure*', fig.cap='  '-----------------------------
denM <- dfptsde3d(fpt3d, pdf = "M")
denM
plot(denM)

## ----111,fig.env='figure*', fig.cap='  ',eval=FALSE, include=TRUE----
#  denJ <- dfptsde3d(fpt3d,pdf="J")
#  plot(denJ,display="rgl")

