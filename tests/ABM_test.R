options(prompt="R> ",scipen=16,digits=4,warning=FALSE, message=FALSE)
library(Sim.DiffProc)

op <- par(mfrow = c(2, 2))

## Brownian motion
set.seed(1234)
X <- BM()
X <- BM(N =1000,M=10,x0=0,t0=0,T=1,Dt=0.001)
plot(X,plot.type="single")
lines(as.vector(time(X)),rowMeans(X),col="red")

## Brownian bridge
set.seed(1234)
X <- BB()
X <- BB(N =1000,M=10,x0=0,y=0,t0=0,T=1,Dt=0.001)
plot(X,plot.type="single")
lines(as.vector(time(X)),rowMeans(X),col="red")

## Geometric Brownian motion
set.seed(1234)
X <- GBM()
X <- GBM(N =1000,M=10,x0=1,t0=0,T=1,Dt=0.001,theta=1,sigma=1)
plot(X,plot.type="single")
lines(as.vector(time(X)),rowMeans(X),col="red")

## Arithmetic Brownian motion
set.seed(1234)
X <- ABM()
X <- ABM(N =1000,M=10,x0=0,t0=0,T=1,Dt=0.001,theta=1,sigma=1)
plot(X,plot.type="single")
lines(as.vector(time(X)),rowMeans(X),col="red")

par(op)