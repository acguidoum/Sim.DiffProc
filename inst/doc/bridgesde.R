## ----setup, echo = F, message = F, results = 'hide',screenshot.force=FALSE----
library(Sim.DiffProc)
library(knitr)
knitr::opts_chunk$set(comment="",prompt=TRUE, fig.show='hold', warning=FALSE, message=FALSE)
options(prompt="R> ",scipen=16,digits=5,warning=FALSE, message=FALSE,
        width = 70)

## -------------------------------------------------------------------
f <- expression((1-x)/(1-t))
g <- expression(x)
mod <- bridgesde1d(drift=f,diffusion=g,x0=3,y=1,M=1000,method="milstein")
mod
summary(mod) ## default: summary at time = (T-t0)/2

## ----01,fig.env='figure*', fig.cap=' ',eval=FALSE, include=TRUE-----
#  plot(mod,ylab=expression(X[t]))
#  lines(time(mod),apply(mod$X,1,mean),col=2,lwd=2)
#  legend("topleft","mean path",inset = .01,col=2,lwd=2,cex=0.8,bty="n")

## ----33, echo=FALSE, fig.cap=' Bridge sde 1D ', fig.env='figure*',fig.width=7,fig.height=7----
knitr::include_graphics("Figures/fig03.png")

## ----02, fig.cap='  ', fig.env='figure*'----------------------------
m  <- apply(mod$X,1,mean) 
S  <- apply(mod$X,1,var)
out <- data.frame(m,S)
matplot(time(mod), out, type = "l", xlab = "time", ylab = "", col=2:3,lwd=2,lty=2:3,las=1)
legend("topright",c(expression(m(t),S(t))),col=2:3,lty=2:3,lwd=2,bty="n")

## -------------------------------------------------------------------
s = 0.55
mean(mod, at = s)
moment(mod, at = s , center = TRUE , order = 2) ## variance
Median(mod, at = s)
Mode(mod, at = s)
quantile(mod , at = s)
kurtosis(mod , at = s)
skewness(mod , at = s)
cv(mod , at = s )
min(mod , at = s)
max(mod , at = s)
moment(mod, at = s , center= TRUE , order = 4)
moment(mod, at = s , center= FALSE , order = 4)

## -------------------------------------------------------------------
summary(mod, at = 0.55)

## -------------------------------------------------------------------
x <- rsde1d(object = mod, at = s) 
head(x, n = 10)
summary(x)

## ----03 ,fig.env='figure*', fig.cap='  '----------------------------
plot(time(mod),mod$X[,1],type="l",ylab="X(t)",xlab="time",axes=F,lty=3)
points(s,x[1],pch=19,col=2,cex=0.5)
lines(c(s,s),c(0,x[1]),lty=2,col=2)
lines(c(0,s),c(x[1],x[1]),lty=2,col=2)
axis(1, s, bquote(at==.(s)),col=2,col.ticks=2)
axis(2, x[1], bquote(X[t==.(s)]),col=2,col.ticks=2)
legend('topright',col=2,pch=19,legend=bquote(X[t==.(s)]==.(x[1])),bty = 'n')
box()

## ----04,fig.env='figure*', fig.cap='  '-----------------------------
dens <- dsde1d(mod, at = s)
dens
plot(dens,hist=TRUE) ## histgramme
plot(dens,add=TRUE)  ## kernel density

## ----05,fig.env='figure*', fig.cap=' Transitional densitie at time $t-s = 0.25,0.75$ '----
plot(dsde1d(mod,at=0.75))
plot(dsde1d(mod,at=0.25),add=TRUE)
legend('topright',col=c('#0000FF4B','#FF00004B'),pch=15,legend=c("t-s=0.25","t-s=0.75"),bty = 'n')

## -------------------------------------------------------------------
fx <- expression(-(1+y)*x , -(1+x)*y)
gx <- expression(0.2*(1-y),0.1*(1-x))
mod2 <- bridgesde2d(drift=fx,diffusion=gx,x0=c(1,-0.5),y=c(1,0.5),Dt=0.01,M=1000,type="str",method="rk1")
mod2
summary(mod2) ## default: summary at time = (T-t0)/2

## ----06,fig.env='figure*', fig.cap=' ',eval=FALSE, include=TRUE-----
#  plot(mod2,col=c('#FF00004B','#0000FF82'))

## ----333, echo=FALSE, fig.cap='  Bridge sde 2D ', fig.env='figure*',fig.width=7,fig.height=7----
knitr::include_graphics("Figures/fig04.png")

## ----07, fig.cap='  ', fig.env='figure*'----------------------------
m1  <- apply(mod2$X,1,mean) 
m2  <- apply(mod2$Y,1,mean) 
S1  <- apply(mod2$X,1,var)
S2  <- apply(mod2$Y,1,var)
C12 <- sapply(1:dim(mod2$X)[1],function(i) cov(mod2$X[i,],mod2$Y[i,]))
out2 <- data.frame(m1,m2,S1,S2,C12)
matplot(time(mod2), out2, type = "l", xlab = "time", ylab = "", col=2:6,lwd=2,lty=2:6,las=1)
legend("top",c(expression(m[1](t),m[2](t),S[1](t),S[2](t),C[12](t))),col=2:6,lty=2:6,lwd=2,bty="n")

## -------------------------------------------------------------------
s = 6.75
mean(mod2, at = s)
moment(mod2, at = s , center = TRUE , order = 2) ## variance
Median(mod2, at = s)
Mode(mod2, at = s)
quantile(mod2 , at = s)
kurtosis(mod2 , at = s)
skewness(mod2 , at = s)
cv(mod2 , at = s )
min(mod2 , at = s)
max(mod2 , at = s)
moment(mod2 , at = s , center= TRUE , order = 4)
moment(mod2 , at = s , center= FALSE , order = 4)

## -------------------------------------------------------------------
summary(mod2, at = 6.75)

## -------------------------------------------------------------------
x2 <- rsde2d(object = mod2, at = s) 
head(x2, n = 10)
summary(x2)

## ----08,fig.env='figure*', fig.cap=' '------------------------------
plot(ts.union(mod2$X[,1],mod2$Y[,1]),col=1:2,lty=3,plot.type="single",type="l",ylab= "",xlab="time",axes=F)
points(s,x2$x[1],pch=19,col=3,cex=0.8)
points(s,x2$y[1],pch=19,col=4,cex=0.8)
lines(c(s,s),c(-10,x2$x[1]),lty=2,col=6)
lines(c(0,s),c(x2$x[1],x2$x[1]),lty=2,col=3)
lines(c(0,s),c(x2$y[1],x2$y[1]),lty=2,col=4)
axis(1, s, bquote(at==.(s)),col=6,col.ticks=6)
axis(2, x2$x[1], bquote(X[t==.(s)]),col=3,col.ticks=3)
axis(2, x2$y[1], bquote(Y[t==.(s)]),col=4,col.ticks=4)
legend('topright',legend=bquote(c(X[t==.(s)]==.(x2$x[1]),Y[t==.(s)]==.(x2$y[1]))),bty = 'n')
box()

## ----09,fig.env='figure*', fig.cap='  '-----------------------------
denM <- dsde2d(mod2,pdf="M",at =s)
denM
plot(denM, main="Marginal Density")

## ----10,fig.env='figure*', fig.cap='  '-----------------------------
denJ <- dsde2d(mod2, pdf="J", n=100,at =s)
denJ
plot(denJ,display="contour",main="Bivariate Transition Density at time t=6.755")
plot(denJ,display="image",main="Bivariate Transition Density at time t=6.755")

## ----11,fig.env='figure*', fig.cap='  '-----------------------------
plot(denJ,main="Bivariate Transition Density at time t=6.75")

## ----eval=FALSE, include=TRUE---------------------------------------
#  for (i in seq(1,9,by=0.005)){
#  plot(dsde2d(mod2, at = i,n=100),display="contour",main=paste0('Transition Density \n t = ',i))
#  }

## -------------------------------------------------------------------
fx <- expression(-4*(1+x)*y, 4*(1-y)*x, 4*(1-z)*y)
gx <- rep(expression(0.2),3)
mod3 <- bridgesde3d(x0=c(0,-1,0.5),y=c(0,-2,0.5),drift=fx,diffusion=gx,M=1000)
mod3
summary(mod3) ## default: summary at time = (T-t0)/2

## ----12,fig.env='figure*', fig.cap=' ',eval=FALSE, include=TRUE-----
#  plot(mod3) ## in time
#  plot3D(mod3,display = "persp",main="3D Bridge SDE's") ## in space

## ----3333, echo=FALSE, fig.cap=' ', fig.env='figure*',fig.width=7,fig.height=7----
knitr::include_graphics("Figures/fig05.png")

## ----33333, echo=FALSE, fig.cap='  Bridge sde 3D ', fig.env='figure*',fig.width=7,fig.height=7----
knitr::include_graphics("Figures/fig06.png")

## ----13, fig.cap='  ', fig.env='figure*'----------------------------
m1  <- apply(mod3$X,1,mean) 
m2  <- apply(mod3$Y,1,mean) 
m3  <- apply(mod3$Z,1,mean) 
S1  <- apply(mod3$X,1,var)
S2  <- apply(mod3$Y,1,var)
S3  <- apply(mod3$Z,1,var)
out3 <- data.frame(m1,m2,m3,S1,S2,S3)
matplot(time(mod3), out3, type = "l", xlab = "time", ylab = "", col=2:7,lwd=2,lty=2:7,las=1)
legend("bottom",c(expression(m[1](t),m[2](t),m[3](t),S[1](t),S[2](t),S[3](t))),col=2:7,lty=2:7,lwd=2,bty="n")

## -------------------------------------------------------------------
s = 0.75
mean(mod3, at = s)
moment(mod3, at = s , center = TRUE , order = 2) ## variance
Median(mod3, at = s)
Mode(mod3, at = s)
quantile(mod3 , at = s)
kurtosis(mod3 , at = s)
skewness(mod3 , at = s)
cv(mod3 , at = s )
min(mod3 , at = s)
max(mod3 , at = s)
moment(mod3 , at = s , center= TRUE , order = 4)
moment(mod3 , at = s , center= FALSE , order = 4)

## -------------------------------------------------------------------
summary(mod3, at = 0.75)

## -------------------------------------------------------------------
x3 <- rsde3d(object = mod3, at = s) 
head(x3, n = 10)
summary(x3)

## ----14,fig.env='figure*', fig.cap=' '------------------------------
plot(ts.union(mod3$X[,1],mod3$Y[,1],mod3$Z[,1]),col=1:3,lty=3,plot.type="single",type="l",ylab= "",xlab="time",axes=F)
points(s,x3$x[1],pch=19,col=4,cex=0.8)
points(s,x3$y[1],pch=19,col=5,cex=0.8)
points(s,x3$z[1],pch=19,col=6,cex=0.8)
lines(c(s,s),c(-10,x3$x[1]),lty=2,col=7)
lines(c(0,s),c(x3$x[1],x3$x[1]),lty=2,col=4)
lines(c(0,s),c(x3$y[1],x3$y[1]),lty=2,col=5)
lines(c(0,s),c(x3$z[1],x3$z[1]),lty=2,col=6)
axis(1, s, bquote(at==.(s)),col=7,col.ticks=7)
axis(2, x3$x[1], bquote(X[t==.(s)]),col=4,col.ticks=4)
axis(2, x3$y[1], bquote(Y[t==.(s)]),col=5,col.ticks=5)
axis(2, x3$z[1], bquote(Z[t==.(s)]),col=6,col.ticks=6)
legend("bottomleft",legend=bquote(c(X[t==.(s)]==.(x3$x[1]),Y[t==.(s)]==.(x3$y[1]),Z[t==.(s)]==.(x3$z[1]))),bty = 'n',cex=0.75)
box()

## ----15,fig.env='figure*', fig.cap='  '-----------------------------
denM <- dsde3d(mod3,pdf="M",at =s)
denM
plot(denM, main="Marginal Density")

## ----111,fig.env='figure*', fig.cap='  ',eval=FALSE, include=TRUE----
#  denJ <- dsde3d(mod3,pdf="J",at=0.75)
#  plot(denJ,display="rgl")

