options(prompt="R> ",scipen=16,digits=4,warning=FALSE, message=FALSE)
library(Sim.DiffProc)

set.seed(1234)
X <- HWV()

set.seed(1234)
X <- HWV(Dt=0.08)
plot(X,plot.type="single")


set.seed(1234)
X <- HWV(N=1000,M=10,mu = 4, theta = 2.5,sigma = 1,x0=10)
plot(X,plot.type="single")
lines(as.vector(time(X)),rowMeans(X),col="red")

#####

set.seed(1234)
X <- OU()

set.seed(1234)
X <- OU(N=1000,M=10,mu = 4,sigma = 1,x0=10)
plot(X,plot.type="single")
lines(as.vector(time(X)),rowMeans(X),col="red")
