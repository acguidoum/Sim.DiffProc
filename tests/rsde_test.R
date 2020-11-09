options(prompt="R> ",scipen=16,digits=4,warning=FALSE, message=FALSE)
library(Sim.DiffProc)


f <- expression(0)
g <- expression(1)
res1 <- snssde1d(drift=f,diffusion=g,M=1)
x <- rsde1d(res1, at = 1)
x <- rsde1d(res1, at = 0.02154)

#####1d

f <- expression(0)
g <- expression(1)
res1 <- snssde1d(drift=f,diffusion=g,M=10)


x <- rsde1d(res1, at = 1)
x <- rsde1d(res1, at = 0.02154)
print(dsde1d(res1))
print(dsde1d(res1, at = 0.02154))
plot(dsde1d(res1))
plot(dsde1d(res1),dens=function(x) dnorm(x,0,1))
plot(dsde1d(res1, at = 0.02154),dens=function(x) dnorm(x,0,sqrt(0.02154)),hist=TRUE,xlim=c(-0.5,1))



## 2-dim SDE
set.seed(1234)

fx <- expression(3*(2-y),2*x)
gx <- expression(1,y)
mod2d <- snssde2d(drift=fx,diffusion=gx,x0=c(1,2),M=1)

# random 
r2d <- rsde2d(mod2d,at=0.5)
r2d <- rsde2d(mod2d,at=0.51475)

# SDE's 2d
fx <- expression(3*(2-y),2*x)
gx <- expression(1,y)
mod2d <- snssde2d(drift=fx,diffusion=gx,x0=c(1,2),M=10)

# random 
r2d <- rsde2d(mod2d,at=0.5)
r2d <- rsde2d(mod2d,at=0.51475)
summary(r2d)

# Marginal density 

denM <- dsde2d(mod2d,pdf="M", at=0.5)
print(denM)
plot(denM)

denM <- dsde2d(mod2d,pdf="M", at=0.51475)
print(denM)
plot(denM,hist=TRUE)
plot(denM)

# Joint density
denJ <- dsde2d(mod2d,pdf="J", at= 0.5)
print(denJ)
denJ <- dsde2d(mod2d,pdf="J", at= 0.51475)
print(denJ)
plot(denJ,display="persp")
plot(denJ,display="contour")
plot(denJ,display="image")

######## 3d

fx <- expression(4*(-1-x)*y, 4*(1-y)*x, 4*(1-z)*y)
gx <- rep(expression(0.2),3)
mod3d <- snssde3d(drift=fx,diffusion=gx,x0=c(1,2,1),M=1)

# random 
r3d <- rsde3d(mod3d ,at=0.5)
r3d <- rsde3d(mod3d ,at=0.51475)

fx <- expression(4*(-1-x)*y, 4*(1-y)*x, 4*(1-z)*y)
gx <- rep(expression(0.2),3)
mod3d <- snssde3d(drift=fx,diffusion=gx,x0=c(1,2,1),M=10)

# random 
r3d <- rsde3d(mod3d ,at=0.5)
r3d <- rsde3d(mod3d ,at=0.51475)
summary(r3d)

# Marginal density 

denM <- dsde3d(mod3d,pdf="M", at=0.5)
print(denM)
plot(denM)

denM <- dsde3d(mod3d,pdf="M", at=0.51475)
print(denM)
plot(denM)
plot(denM,hist=TRUE)

# Joint density
#denJ <- dsde3d(mod3d,pdf="J", at= 0.5)
#print(denJ)
#denJ <- dsde3d(mod3d,pdf="J", at= 0.51475)
#print(denJ)

###

#####1d

f <- expression(-2*x)
g <- expression(0.5)
res1 <- bridgesde1d(drift=f,diffusion=g,M=10,Dt=0.001,x0=0.5,y=0.5)


x <- rsde1d(res1, at = 1)
x <- rsde1d(res1, at = 0.02154)
print(dsde1d(res1))
print(dsde1d(res1, at = 0.02154))
plot(dsde1d(res1))
plot(dsde1d(res1, at = 0.02154),hist=TRUE,xlim=c(-0.5,1))



## 2-dim SDE
set.seed(1234)

# SDE's 2d
fx <- expression(4*(-1-x) , x)
gx <- expression(0.2 , 0)
mod2d <- bridgesde2d(drift=fx,diffusion=gx,M=10)

# random 
r2d <- rsde2d(mod2d,at=0.5)
r2d <- rsde2d(mod2d,at=0.51475)
summary(r2d)

# Marginal density 

denM <- dsde2d(mod2d,pdf="M", at=0.5)
print(denM)
plot(denM)

denM <- dsde2d(mod2d,pdf="M", at=0.51475)
print(denM)
plot(denM,hist=TRUE)
plot(denM)

# Joint density
denJ <- dsde2d(mod2d,pdf="J", at= 0.5)
print(denJ)
denJ <- dsde2d(mod2d,pdf="J", at= 0.51475)
print(denJ)
plot(denJ,display="persp")
plot(denJ,display="contour")
plot(denJ,display="image")

######## 3d

fx <- expression(4*(-1-x)*y, 4*(1-y)*x, 4*(1-z)*y)
gx <- rep(expression(0.2),3)
mod3d <- bridgesde3d(drift=fx,diffusion=gx,x0=c(0,-1,0.5),y=c(0,-2,0.5),M=10)

# random 
r3d <- rsde3d(mod3d ,at=0.5)
r3d <- rsde3d(mod3d ,at=0.51475)
summary(r3d)

# Marginal density 

denM <- dsde3d(mod3d,pdf="M", at=0.5)
print(denM)
plot(denM)

denM <- dsde3d(mod3d,pdf="M", at=0.51475)
print(denM)
plot(denM)
plot(denM,hist=TRUE)

# Joint density
# denJ <- dsde3d(mod3d,pdf="J", at= 0.5)
# print(denJ)
# denJ <- dsde3d(mod3d,pdf="J", at= 0.51475)
# print(denJ)
