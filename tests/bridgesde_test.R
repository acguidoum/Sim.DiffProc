options(prompt="R> ",scipen=16,digits=4,warning=FALSE, message=FALSE)
library(Sim.DiffProc) 



####### 1d

f <- expression( 2*(1-x) )
g <- expression( 1 )


mod1 <- bridgesde1d(drift=f,diffusion=g,x0=2,y=1,M=10,N=50,Dt=0.05)
print(mod1)
summary(mod1)				
plot(mod1)
bconfint(mod1)

mean(mod1)
moment(mod1, center = TRUE , order = 2) ## variance
Median(mod1)
Mode(mod1)
quantile(mod1)
kurtosis(mod1)
skewness(mod1)
cv(mod1)
min(mod1)
max(mod1)
moment(mod1, center= TRUE )
moment(mod1, center= FALSE , order = 4)

##

mod1 <- bridgesde1d(drift=f,diffusion=g,x0=2,y=1,M=10,N=50,T=2,type="str")
print(mod1)
summary(mod1)		
plot(mod1,type="n",plot.type="single")
lines(mod1,col=2,lwd=2)
points(mod1,pch=21,col=5,cex=0.1)

s = 0.510547
mean(mod1, at = s)
moment(mod1, at = s , center = TRUE , order = 2) ## variance
Median(mod1, at = s)
Mode(mod1, at = s)
quantile(mod1 , at = s)
kurtosis(mod1 , at = s)
skewness(mod1 , at = s)
cv(mod1 , at = s )
min(mod1 , at = s)
max(mod1 , at = s)
moment(mod1, at = s , center= TRUE , order = 4)
moment(mod1, at = s , center= FALSE , order = 4)
bconfint(mod1, at =s)

######## 2d

fx <- expression(4*(-1-x) , x)
gx <- expression(0.2 , 0)

mod2 <- bridgesde2d(drift=fx,diffusion=gx,Dt=0.005,M=10,N=50)
print(mod2)				
summary(mod2)
mean(mod2)
moment(mod2, center = TRUE , order = 2) ## variance
Median(mod2)
Mode(mod2)
quantile(mod2)
kurtosis(mod2)
skewness(mod2)
cv(mod2)
min(mod2)
max(mod2)
moment(mod2, center= TRUE , order = 4)
moment(mod2, center= FALSE , order = 4)

				
plot(mod2,type="n")
lines(mod2,col=2,lwd=2)
points(mod2,pch=21,col=5,cex=0.1)
bconfint(mod2)
plot2d(mod2,type="n")
lines2d(mod2,col=4)
points2d(mod2,pch=19,cex=0.1)

Sigma <- matrix(c(1, 0.75, 0.75, 1), nrow = 2, ncol = 2)
mod2 <- bridgesde2d(drift=fx,diffusion=gx,corr=Sigma,Dt=0.005,M=10,N=50)
print(mod2)				
summary(mod2)

##

mod2 <- bridgesde2d(drift=fx,diffusion=gx,T=5,M=10,N=50,type="str")
print(mod2)				
summary(mod2)				
plot(mod2)
bconfint(mod2)
plot2d(mod2,type="n")
lines2d(mod2,col=4)
points2d(mod2,pch=19,cex=0.1)

mean(mod2, at = s)
moment(mod2, at = s , center = TRUE , order = 2) ## variance
Median(mod2, at = s)
Mode(mod2, at = s)
quantile(mod2 , at = s)
kurtosis(mod2 , at = s)
skewness(mod2 , at = s)
cv(mod2 , at = s )
min(mod2 , at = s)
max(mod2 , at = s)
moment(mod2, at = s , center= TRUE , order = 4)
moment(mod2, at = s , center= FALSE , order = 4)
bconfint(mod2, at =s)

Sigma <- matrix(c(1, 0.75, 0.75, 1), nrow = 2, ncol = 2)
mod2 <- bridgesde2d(drift=fx,diffusion=gx,corr=Sigma,Dt=0.005,M=10,N=50,type="str")
print(mod2)				
summary(mod2)

######## 3d

fx <- expression(4*(-1-x)*y, 4*(1-y)*x, 4*(1-z)*y)
gx <- rep(expression(0.2),3)

mod3 <- bridgesde3d(drift=fx,diffusion=gx,x0=c(0,-1,0.5),y=c(0,-2,0.5),M=10,Dt=0.01,N=50)
print(mod3)				
summary(mod3)
mean(mod3)
moment(mod3, center = TRUE , order = 2) ## variance
Median(mod3)
Mode(mod3)
quantile(mod3)
kurtosis(mod3)
skewness(mod3)
cv(mod3)
min(mod3)
max(mod3)
moment(mod3, center= TRUE , order = 4)
moment(mod3, center= FALSE , order = 4)
				
plot(mod3)
bconfint(mod3)
plot3D(mod3,display="persp",main="3-dim bridge sde")

Sigma <- matrix(c(1,-0.5,-0.25,-0.5,1,0.95,-0.25,0.95,1),nrow=3,ncol=3) 
mod3 <- bridgesde3d(drift=fx,diffusion=gx,corr=Sigma,x0=c(0,-1,0.5),y=c(0,-2,0.5),M=10,Dt=0.01,N=50)
print(mod3)				
summary(mod3)

##


mod3 <- bridgesde3d(drift=fx,diffusion=gx,x0=c(0,-1,0.5),y=c(0,-2,0.5),M=10,N=50,type="str")
print(mod3)
summary(mod3)				
bconfint(mod3)
plot(mod3,type="n")
lines(mod3,col=2,lwd=2)
points(mod3,pch=21,col=5,cex=0.1)

mean(mod3, at = s)
moment(mod3, at = s , center = TRUE , order = 2) ## variance
Median(mod3, at = s)
Mode(mod3, at = s)
quantile(mod3 , at = s)
kurtosis(mod3 , at = s)
skewness(mod3 , at = s)
cv(mod3 , at = s )
min(mod3 , at = s)
max(mod3 , at = s)
moment(mod3, at = s , center= TRUE , order = 4)
moment(mod3, at = s , center= FALSE , order = 4)
bconfint(mod3, at =s)

Sigma <- matrix(c(1,-0.5,-0.25,-0.5,1,0.95,-0.25,0.95,1),nrow=3,ncol=3) 
mod3 <- bridgesde3d(drift=fx,diffusion=gx,corr=Sigma,x0=c(0,-1,0.5),y=c(0,-2,0.5),M=10,Dt=0.01,N=50,type="str")
print(mod3)				
summary(mod3)

#############################

f <- expression( 2*(1-x) )
g <- expression( 1 )

mod1 <- bridgesde1d(drift=f,diffusion=g,x0=2,y=2,M=1,N=50,Dt=0.05,method="predcorr")			
plot(mod1,type="n")
lines(mod1,col=2)
points(mod1,cex=0.1,pch=19)

mean(mod1)
moment(mod1, center = TRUE , order = 2) ## variance
Median(mod1)
quantile(mod1)
kurtosis(mod1)
skewness(mod1)
cv(mod1)
min(mod1)
max(mod1)

####

fx <- expression(4*(-1-x) , x)
gx <- expression(0.2 , 0)

mod2 <- bridgesde2d(drift=fx,diffusion=gx,M=1,Dt=0.05,N=50)
plot(mod2,union = FALSE)
plot(mod2,type="n")
lines(mod2,col=2)
points(mod2,cex=0.1,pch=19)
plot2d(mod2)

mean(mod2)
moment(mod2, center = TRUE , order = 2) ## variance
Median(mod2)
quantile(mod2)
kurtosis(mod2)
skewness(mod2)
cv(mod2)
min(mod2)
max(mod2)

####

fx <- expression(4*(-1-x), 4*(1-y), 4*(1-z))
gx <- rep(expression(0.2),3)

mod3 <- bridgesde3d(drift=fx,diffusion=gx,M=10,Dt=0.05,N=50,method="predcorr")
plot(mod3,union = FALSE)
plot(mod3,type="n")
lines(mod3,col=2)
points(mod3,cex=0.1,pch=19)
plot3D(mod3)

mean(mod3)
moment(mod3, center = TRUE , order = 2) ## variance
Median(mod3)
quantile(mod3)
kurtosis(mod3)
skewness(mod3)
cv(mod3)
min(mod3)
max(mod3)
