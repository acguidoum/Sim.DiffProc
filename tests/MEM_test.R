options(prompt="R> ",scipen=16,digits=4,warning=FALSE, message=FALSE)
library(Sim.DiffProc)


###################

f <- expression(mu*x)
g <- expression(sigma*x)
para <- c(mu=2,sigma=0.5)
t    <- seq(0,1,by=0.001)
init <- c(m=1,S=0)

# Ito
res1 <- MEM.sde(drift=f,diffusion=g,solve=TRUE,init=init,parms=para,time=t)
print(res1)
summary(res1,at = 0.5)

# Str

res2 <- MEM.sde(drift=f,diffusion=g,type = "str",solve=TRUE,init=init,parms=para,time=t)
print(res2)
summary(res2,at = 0.5)

###################

f <- expression(1/mu*(theta-x), x)  
g <- expression(sqrt(sigma),0)
para=c(mu=0.75,sigma=0.1,theta=2)
init=c(m1=0,m2=0,S1=0,S2=0,C12=0)
t <- seq(0,10,by=0.001)

# Ito
res2d1 <- MEM.sde(drift=f,diffusion=g,solve=TRUE,init=init,parms=para,time=t)
print(res2d1)
summary(res2d1)

# Str

res2d2 <- MEM.sde(drift=f,diffusion=g,type = "str",solve=TRUE,init=init,parms=para,time=t)
print(res2d2)
summary(res2d2)

###################

f <- expression(sigma*(y-x),rho*x-y-x*z,x*y-bet*z)
g <- expression(0.1,0.1,0.1)
para=c(sigma=10,rho=28,bet=8/3)
ini=c(m1=1,m2=1,m3=1,S1=0,S2=0,S3=0,C12=0,C13=0,C23=0)

# Ito

res3d1 = MEM.sde(drift=f,diffusion=g,solve=TRUE,parms=para,init=ini,time=seq(0,1,by=0.01))
print(res3d1)
summary(res3d1)

# Str

res3d2 = MEM.sde(drift=f,diffusion=g,type = "str",solve=TRUE,parms=para,init=ini,time=seq(0,1,by=0.01))
print(res3d2)
summary(res3d2)