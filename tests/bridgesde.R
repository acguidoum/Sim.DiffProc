library(Sim.DiffProc)


### 2-dim Ito bridge sde
set.seed(1234)

fx <- expression(4*(-1-x) , x)
gx <- expression(0.2 , 0)
res <- bridgesde2d(drift=fx,diffusion=gx,Dt=0.005,M=2000)
res
summary(res) ## Monte-Carlo statistics at time T/2=2.5
summary(res,at=1) ## Monte-Carlo statistics at time 1
summary(res,at=4) ## Monte-Carlo statistics at time 4
##
plot(res,type="n")
lines(time(res),apply(res$X,1,mean),col=3,lwd=2)
lines(time(res),apply(res$Y,1,mean),col=4,lwd=2)
legend("topright",c(expression(E(X[t])),expression(E(Y[t]))),lty=1,inset = .7,col=c(3,4))
##
plot2d(res)

