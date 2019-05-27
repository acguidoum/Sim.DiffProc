library(Sim.DiffProc)

## 2-dim SDE
set.seed(1234)

# SDE's 2d
fx <- expression(3*(2-y),2*x)
gx <- expression(1,y)
mod2d <- snssde2d(drift=fx,diffusion=gx,x0=c(1,2),M=5000)

# random 
r2d <- rsde2d(mod2d,at=0.5)
summary(r2d)

# Marginal density 

denM <- dsde2d(mod2d,pdf="M", at=0.5)
denM
plot(denM)

# Joint density
denJ <- dsde2d(mod2d,pdf="J", at= 0.5)
denJ
plot(denJ)
plot(denJ,display="contour")
