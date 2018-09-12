options(prompt="R> ",scipen=16,digits=4,warning=FALSE, message=FALSE)
library(Sim.DiffProc) 

## -------------------------------------------------------------------

f <- expression( (1-0.5*x) )
g <- expression( 1 )
mod1d <- snssde1d(drift=f,diffusion=g,x0=1.7,M=20)

St  <- expression(2*(1-sinh(0.5*t)) )
fpt1d <- fptsde1d(mod1d, boundary = St)
print(fpt1d)
head(fpt1d$fpt, n = 4)


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

summary(fpt1d)


plot(time(mod1d),mod1d$X[,1],type="l",lty=3,ylab="X(t)",xlab="time",axes=F)
curve(2*(1-sinh(0.5*x)),add=TRUE,col=2)
points(fpt1d$fpt[1],2*(1-sinh(0.5*fpt1d$fpt[1])),pch=19,col=4,cex=0.5)
lines(c(fpt1d$fpt[1],fpt1d$fpt[1]),c(0,2*(1-sinh(0.5*fpt1d$fpt[1]))),lty=2,col=4)
axis(1, fpt1d$fpt[1], bquote(tau[S(t)]==.(fpt1d$fpt[1])),col=4,col.ticks=4)
legend('topleft',col=c(1,2,4),lty=c(1,1,NA),pch=c(NA,NA,19),legend=c(expression(X[t]),expression(S(t)),expression(tau[S(t)])),cex=0.8,bty = 'n')
box()

print(dfptsde1d(fpt1d))
plot(dfptsde1d(fpt1d),hist=TRUE,nbins="FD")  ## histogramm
plot(dfptsde1d(fpt1d))              ## kernel density

## -------------------------------------------------------------------

f <- expression( (1-0.5*x) )
g <- expression( 1 )
mod1d <- snssde1d(drift=f,diffusion=g,x0=2.1,M=20,type="str")

St  <- expression(2*(1-sinh(0.5*t)) )
fpt1d <- fptsde1d(mod1d, boundary = St)
print(fpt1d)
head(fpt1d$fpt, n = 4)


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

summary(fpt1d)


plot(time(mod1d),mod1d$X[,1],type="l",lty=3,ylab="X(t)",xlab="time",axes=F)
curve(2*(1-sinh(0.5*x)),add=TRUE,col=2)
points(fpt1d$fpt[1],2*(1-sinh(0.5*fpt1d$fpt[1])),pch=19,col=4,cex=0.5)
lines(c(fpt1d$fpt[1],fpt1d$fpt[1]),c(0,2*(1-sinh(0.5*fpt1d$fpt[1]))),lty=2,col=4)
axis(1, fpt1d$fpt[1], bquote(tau[S(t)]==.(fpt1d$fpt[1])),col=4,col.ticks=4)
legend('topleft',col=c(1,2,4),lty=c(1,1,NA),pch=c(NA,NA,19),legend=c(expression(X[t]),expression(S(t)),expression(tau[S(t)])),cex=0.8,bty = 'n')
box()

print(dfptsde1d(fpt1d))
plot(dfptsde1d(fpt1d),hist=TRUE,nbins="FD")  ## histogramm
plot(dfptsde1d(fpt1d))


## -------------------------------------------------------------------


fx <- expression(5*(-1-y)*x , 5*(-1-x)*y)
gx <- expression(0.5*y,0.5*x)
mod2d <- snssde2d(drift=fx,diffusion=gx,x0=c(x=1,y=-1),M=20)

St <- expression(sin(2*pi*t))
fpt2d <- fptsde2d(mod2d, boundary = St)
print(fpt2d)
head(fpt2d$fpt, n = 4)

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

summary(fpt2d)

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

denM <- dfptsde2d(fpt2d, pdf = 'M')
print(denM)
plot(denM)
plot(denM, hist=TRUE)

denJ <- dfptsde2d(fpt2d, pdf = 'J',n=100)
print(denJ)
plot(denJ,display="contour",main="Bivariate Density of F.P.T",xlab=expression(tau[x]),ylab=expression(tau[y]))
plot(denJ,display="image",main="Bivariate Density of F.P.T",xlab=expression(tau[x]),ylab=expression(tau[y]))
plot(denJ,display="persp",main="Bivariate Density of F.P.T",xlab=expression(tau[x]),ylab=expression(tau[y]))

## -------------------------------------------------------------------

fx <- expression(5*(-1-y)*x , 5*(-1-x)*y)
gx <- expression(0.5*y,0.5*x)
mod2d <- snssde2d(drift=fx,diffusion=gx,x0=c(x=-1,y=1),M=20,type="str")

St <- expression(sin(2*pi*t))
fpt2d <- fptsde2d(mod2d, boundary = St)
print(fpt2d)
head(fpt2d$fpt, n = 4)

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

summary(fpt2d)

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

denM <- dfptsde2d(fpt2d, pdf = 'M')
print(denM)
plot(denM)
plot(denM, hist=TRUE)

denJ <- dfptsde2d(fpt2d, pdf = 'J',n=100)
print(denJ)
plot(denJ,display="contour",main="Bivariate Density of F.P.T",xlab=expression(tau[x]),ylab=expression(tau[y]))
plot(denJ,display="image",main="Bivariate Density of F.P.T",xlab=expression(tau[x]),ylab=expression(tau[y]))
plot(denJ,display="persp",main="Bivariate Density of F.P.T",xlab=expression(tau[x]),ylab=expression(tau[y]))

## -------------------------------------------------------------------

fx <- expression(4*(-1-x)*y , 4*(1-y)*x , 4*(1-z)*y) 
gx <- rep(expression(0.2),3)
mod3d <- snssde3d(drift=fx,diffusion=gx,x0=c(x=2,y=-2,z=0),M=20)

St <- expression(-1.5+3*t)
fpt3d <- fptsde3d(mod3d, boundary = St)
print(fpt3d)
head(fpt3d$fpt, n = 4)

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

summary(fpt3d)

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

denM <- dfptsde3d(fpt3d, pdf = "M")
print(denM)
plot(denM)
plot(denM,hist=TRUE)

denJ <- dfptsde3d(fpt3d,pdf="J")
print(denJ)
plot(denJ,display="rgl")

## -------------------------------------------------------------------

fx <- expression(4*(1-x)*y , 4*(1-y)*x , 4*(1-z)*y) 
gx <- rep(expression(0.2),3)
mod3d <- snssde3d(drift=fx,diffusion=gx,x0=c(x=-2,y=2,z=-2),M=20,type="str")

St <- expression(-1.5+3*t)
fpt3d <- fptsde3d(mod3d, boundary = St)
print(fpt3d)
head(fpt3d$fpt, n = 4)

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

summary(fpt3d)

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

denM <- dfptsde3d(fpt3d, pdf = "M")
print(denM)
plot(denM)
plot(denM,hist=TRUE)

denJ <- dfptsde3d(fpt3d,pdf="J")
print(denJ)
plot(denJ,display="rgl")