library(Sim.DiffProc)

## sde 2-dim
set.seed(1234)

fx <- expression( y , (4*( 1-x^2 )* y - x))
gx <- expression( 0 , 0.2)

mod2d2 <- snssde2d(drift=fx,diffusion=gx,type="str",T=100,N=10000)
mod2d2
plot(mod2d2,pos=2)
dev.new()
plot(mod2d2,union = FALSE)
dev.new()
plot2d(mod2d2,type="n") ## in plane (O,X,Y)
points2d(mod2d2,col=rgb(0,100,0,50,maxColorValue=255), pch=16)


