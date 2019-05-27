library(Sim.DiffProc)


## Example 1: 

# SDE's 2d
fx <- expression(5*(-1-y)*x , 5*(-1-x)*y)
gx <- expression(0.5 , 0.5)
mod2d <- snssde2d(drift=fx,diffusion=gx,x0=c(2,-2),M=2000)

# boundary

St <- expression(-1+5*t)

# random fpt

out <- fptsde2d(mod2d,boundary=St)
out
summary(out)

# Marginal density 

denM <- dfptsde2d(out,pdf="M")
denM
plot(denM)

# Joint density

denJ <- dfptsde2d(out,pdf="J",n=200,lims=c(0.28,0.4,0.04,0.13))
denJ
plot(denJ)
plot(denJ,display="image")
plot(denJ,display="image",drawpoints=TRUE,cex=0.5,pch=19,col.pt='green')
plot(denJ,display="contour")
plot(denJ,display="contour",color.palette=colorRampPalette(c('white','green','blue','red')))