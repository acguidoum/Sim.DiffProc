## Thu Apr 30 12:03:41 2020
options(prompt="R> ",scipen=20,digits=5,scientific = 5,width = 70,warning=FALSE, message=FALSE)

############################################################################
#                               Demo 1                                     # 
#                              1-dim SDE                                   #
############################################################################        
set.seed(1234)
f <- expression(2*(3-x) )
g <- expression(1)
mod1 <- snssde1d(drift=f,diffusion=g,M=4000,x0=10,Dt=0.01)
mod1
summary(mod1)
plot(mod1)
lines(time(mod1),apply(mod1$X,1,mean),col=2,lwd=2)
lines(time(mod1),apply(mod1$X,1,bconfint,level=0.95)[1,],col=4,lwd=2)
lines(time(mod1),apply(mod1$X,1,bconfint,level=0.95)[2,],col=4,lwd=2)
legend("topright",c("mean path",paste("bound of", 95,"% percent confidence")),
       inset = .01,col=c(2,4),lwd=2,cex=0.8)

############################################################################
#                               Demo 2                                     # 
#                              2-dim SDE                                   #
############################################################################        
set.seed(1234)
fx <- expression( y ,(4*( 1-x^2 )* y - x))
gx <- expression( 0 , 0.2)

res1 <- snssde2d(drift=fx,diffusion=gx,type="str",T=100,N=10000)
res1
plot(res1,pos=2)
plot(res1,union=FALSE)
dev.new()
plot2d(res1,type="n") ## in plane (O,X,Y)
points2d(res1,col=rgb(0,100,0,50,maxColorValue=255), pch=16)

############################################################################
#                               Demo 3                                     # 
#                              3-dim SDE                                   #
############################################################################        
set.seed(1234)

fx <- expression(4*(-1-x)*y, 4*(1-y)*x, 4*(1-z)*y)
gx <- rep(expression(0.2),3)
res <- snssde3d(x0=c(2,-2,-2),drift=fx,diffusion=gx,M=1000)
res
summary(res)
plot(res,pos=2)
dev.new()
plot3D(res,display="persp")

############################################################################
#                               Demo 4                                     # 
#                          2-dim Bridge SDE                                #
############################################################################
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

############################################################################
#                               Demo 5                                     # 
#                            Bridge 3-dim SDE                              #
############################################################################  
set.seed(1234)

fx <- expression(4*(-1-x)*y, 4*(1-y)*x, 4*(1-z)*y)
gx <- rep(expression(0.2),3)

res <- bridgesde3d(x0=c(0,-1,0.5),y=c(0,-2,0.5),drift=fx,diffusion=gx,M=1000)
res
summary(res) ## Monte-Carlo statistics at time T/2=0.5
summary(res,at=0.25) ## Monte-Carlo statistics at time 0.25
summary(res,at=0.75) ## Monte-Carlo statistics at time 0.75
##
plot(res)
plot(res,union=FALSE)
##
plot(res,type="n")
lines(time(res),apply(res$X,1,mean),col=3,lwd=2)
lines(time(res),apply(res$Y,1,mean),col=4,lwd=2)
lines(time(res),apply(res$Z,1,mean),col=5,lwd=2)
legend("topleft",c(expression(E(X[t])),expression(E(Y[t])),expression(E(Z[t]))),lty=1,inset = .01,col=c(3,4,5))
##
plot3D(res,display = "persp",main="3-dim bridge sde")

############################################################################
#                               Demo 6                                     # 
#                              1-dim FPT                                   #
############################################################################        

## X(t) Brownian motion
## S(t) = 0 (time-dependent boundary)
set.seed(1234)

f <- expression( -4*x )
g <- expression( 0.5 )
mod <- snssde1d(drift=f,diffusion=g,x0=2,M=5000)

# boundary
St <- expression(0) 

# random
out <- fptsde1d(mod, boundary=St)
out
summary(out)
# density approximate
den <- dfptsde1d(out)
den
plot(den)
plot(den,hist=TRUE)



############################################################################
#                               Demo 7                                     # 
#                            Fiting 1-dim SDE                              #
############################################################################  
set.seed(1234)

## Application to real data
## CKLS modele vs CIR modele 
## CKLS (mod1):  dX(t) = (theta1+theta2* X(t))* dt + theta3 * X(t)^theta4 * dW(t)
## CIR  (mod2):  dX(t) = (theta1+theta2* X(t))* dt + theta3 * sqrt(X(t))  * dW(t)
set.seed(1234)

data(Irates)
rates <- Irates[,"r1"]
rates <- window(rates, start=1964.471, end=1989.333)

fx1 <- expression(theta[1]+theta[2]*x)
gx1 <- expression(theta[3]*x^theta[4])
gx2 <- expression(theta[3]*sqrt(x))

fitmod1 <- fitsde(rates,drift=fx1,diffusion=gx1,pmle="euler",start = list(theta1=1,theta2=1,
                  theta3=1,theta4=1),optim.method = "L-BFGS-B")
fitmod2 <- fitsde(rates,drift=fx1,diffusion=gx2,pmle="euler",start = list(theta1=1,theta2=1,
                  theta3=1),optim.method = "L-BFGS-B")	
summary(fitmod1)
summary(fitmod2)
coef(fitmod1)
coef(fitmod2)
confint(fitmod1,parm=c('theta2','theta3'))
confint(fitmod2,parm=c('theta2','theta3'))
AIC(fitmod1)
AIC(fitmod2)	

############################################################################
#                               Demo 8                                     # 
#                 Transition Density of Brownian motion                    #
############################################################################  

f <- expression(0); g <- expression(1)
B <- snssde1d(drift=f,diffusion=g,M=1000)
for (i in seq(B$t0,B$T,by=B$Dt)){
plot(dsde1d(B, at = i),main=paste0('Transition Density \n t = ',i))
}

############################################################################
#                               Demo 9                                     # 
#     Bivariate Transition Density of 2 Brownian motion                    #
############################################################################  

fx <- expression(0,0)
gx <- expression(1,1)
B <- snssde2d(drift=fx,diffusion=gx,M=1000)
B
for (i in seq(B$Dt,B$T,by=B$Dt)){
plot(dsde2d(B, at = i),display="contour",main=paste0('Transition Density \n t = ',i))
}

############################################################################
#                               Demo 10                                    # 
#     Bivariate Transition Density of 3 Brownian motion                    #
############################################################################  

fx <- rep(expression(0),3)
gx <- rep(expression(1),3)
B <- snssde3d(drift=fx,diffusion=gx,M=1000)
plot(dsde3d(B, at = 1),display="rgl")
plot(dsde3d(B, at = 1),display="persp")

############################################################################
#                               Demo 11                                    # 
#                         Parallel Monte Carlo                             #
############################################################################  

theta = 0.75; x0 = 1
fx <- expression( (0.5*theta^2*x) )
gx <- expression( theta*x )
mod1 <- snssde1d(drift=fx,diffusion=gx,x0=x0,M=5000,type="ito")
mod2 <- snssde1d(drift=fx,diffusion=gx,x0=x0,M=5000,type="str")	
E.mod1 <- function(t) x0*exp(0.5*theta^2*t)
E.mod2 <- function(t) x0
V.mod1 <- function(t) x0^2*exp(theta^2*t)*(exp(theta^2*t)-1)
V.mod2 <- function(t) x0^2 *(exp(theta^2*t)-1)
sde.fun1d <- function(data, i){
   d <- data[i, ]
   return(c(mean(d),var(d)))
 }
mcm.mod1 = MCM.sde(model=mod1,statistic=sde.fun1d,R=100,
           exact=list(m=E.mod1(1),S=V.mod1(1)),parallel="snow",ncpus=4)
mcm.mod1
mcm.mod2 = MCM.sde(model=mod2,statistic=sde.fun1d,R=100,
            exact=list(m=E.mod2(1),S=V.mod2(1)),parallel="snow",ncpus=4)
mcm.mod2

############################################################################
#                               Demo 12                                    # 
#                            Moment equations                              #
############################################################################  
fx <- expression( (0.5*theta^2*x) )
gx <- expression( theta*x )
start = c(m=1,S=0)
t = seq(0,1,by=0.001)	
mem.mod1 = MEM.sde(drift=fx,diffusion=gx,type="ito",solve = TRUE, 
                     parms = c(theta=0.75), init = start, time = t) 
mem.mod1
mem.mod2 = MEM.sde(drift=fx,diffusion=gx,type="str",solve = TRUE, 
                     parms = c(theta=0.75), init = start, time = t)
mem.mod2

plot(mem.mod1$sol.ode, mem.mod2$sol.ode,ylab = c("m(t)"),select="m",
       xlab = "Time",main="",col = 2:3,lty=1)
legend("topleft",c(expression(m[mod1](t),m[mod2](t))),inset = .05,
         col=2:3,lty=1)
plot(mem.mod1$sol.ode, mem.mod2$sol.ode,ylab = c("S(t)"),select="S",
       xlab = "Time",main="",col = 2:3,lty=1)
legend("topleft",c(expression(S[mod1](t),S[mod2](t))),inset = .05,
         col=2:3,lty=1)
		 
############################################################################
#                               Demo 13                                    # 
#                  Converting Sim.DiffProc Objects to LaTeX                #
############################################################################  

f <- expression(-mu1 * x) 
g <- expression(mu2 * sqrt(x)) 
TEX.sde(object = c(drift = f, diffusion = g))

# Example 2

f <- expression(mu1*cos(mu2+z),mu1*sin(mu2+z),0) 
g <- expression(sigma,sigma,alpha) 
TEX.sde(object = c(drift = f, diffusion = g))

## LaTeX mathematic for object of class 'MEM.sde'
## Copy and paste the following output in your LaTeX file

# Example 3
f <- expression(mu1*cos(mu2+z),mu1*sin(mu2+z),0) 
g <- expression(sigma,sigma,alpha) 
mem.mod3d <- MEM.sde(drift = f, diffusion = g)
TEX.sde(object = mem.mod3d)
