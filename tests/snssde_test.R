options(prompt="R> ",scipen=16,digits=8,warning=FALSE, message=FALSE)
library(Sim.DiffProc) 


meth <- eval(formals(snssde1d.default)$method)


######## 1d


f <- expression( 2*(1-x) )
g <- expression( 1 )


mod1 <- lapply(2:length(meth), function(i) 
                snssde1d(N=50,drift=f,diffusion=g,x0=2,M=2,method=meth[i]))
print(mod1[[1]])
summary(mod1[[1]])				
plot(mod1[[1]],plot.type="single")
bconfint(mod1[[1]])

mean(mod1[[1]])
moment(mod1[[1]], center = TRUE , order = 2) ## variance
Median(mod1[[1]])
Mode(mod1[[1]])
quantile(mod1[[1]])
kurtosis(mod1[[1]])
skewness(mod1[[1]])
cv(mod1[[1]])
min(mod1[[1]])
max(mod1[[1]])
moment(mod1[[1]], center= TRUE , order = 4)
moment(mod1[[1]], center= FALSE , order = 4)

##

mod1 <- lapply(2:length(meth), function(i) 
                snssde1d(N=50,drift=f,diffusion=g,x0=2,M=2,method=meth[i],type="str"))
print(mod1[[1]])
summary(mod1[[1]])		
plot(mod1[[1]],type="n",plot.type="single")
lines(mod1[[1]],col=2,lwd=2)
points(mod1[[1]],pch=21,col=5,cex=0.5)

s = 0.0554747
mean(mod1[[1]], at = s)
moment(mod1[[1]], at = s , center = TRUE , order = 2) ## variance
Median(mod1[[1]], at = s)
Mode(mod1[[1]], at = s)
quantile(mod1[[1]] , at = s)
kurtosis(mod1[[1]] , at = s)
skewness(mod1[[1]] , at = s)
cv(mod1[[1]] , at = s )
min(mod1[[1]] , at = s)
max(mod1[[1]] , at = s)
moment(mod1[[1]], at = s , center= TRUE , order = 4)
moment(mod1[[1]], at = s , center= FALSE , order = 4)
bconfint(mod1[[1]], at =s)

######## 2d

fx <- expression(4*(-1-x) , x)
gx <- expression(0.2 , 0)

mod2 <- lapply(2:length(meth), function(i) 
                snssde2d(N=50,drift=fx,diffusion=gx,Dt=0.005,M=2,method=meth[i]))
print(mod2[[1]])				
summary(mod2[[1]])
mean(mod2[[1]])
moment(mod2[[1]], center = TRUE , order = 2) ## variance
Median(mod2[[1]])
Mode(mod2[[1]])
quantile(mod2[[1]])
kurtosis(mod2[[1]])
skewness(mod2[[1]])
cv(mod2[[1]])
min(mod2[[1]])
max(mod2[[1]])
moment(mod2[[1]], center= TRUE , order = 4)
moment(mod2[[1]], center= FALSE , order = 4)

				
plot(mod2[[1]],type="n")
lines(mod2[[1]],col=2,lwd=2)
points(mod2[[1]],pch=21,col=5,cex=0.1)
bconfint(mod2[[1]])
plot2d(mod2[[1]],type="n")
lines2d(mod2[[1]],col=4)
points2d(mod2[[1]],pch=19,cex=0.1)

Sigma <- matrix(c(1, 0.75, 0.75, 1), nrow = 2, ncol = 2)
mod2 <- lapply(1:2, function(i) 
                snssde2d(N=50,drift=fx,diffusion=gx,corr=Sigma,Dt=0.005,M=2,method=meth[i]))
print(mod2[[1]])				
summary(mod2[[1]])
##

mod2 <- lapply(2:length(meth), function(i) 
                snssde2d(N=50,drift=fx,diffusion=gx,Dt=0.005,M=2,method=meth[i],type="str"))
print(mod2[[1]])				
summary(mod2[[1]])				
plot(mod2[[1]])
bconfint(mod2[[1]])
plot2d(mod2[[1]],type="n")
lines2d(mod2[[1]],col=4)
points2d(mod2[[1]],pch=19,cex=0.1)

mean(mod2[[1]], at = s)
moment(mod2[[1]], at = s , center = TRUE , order = 2) ## variance
Median(mod2[[1]], at = s)
Mode(mod2[[1]], at = s)
quantile(mod2[[1]] , at = s)
kurtosis(mod2[[1]] , at = s)
skewness(mod2[[1]] , at = s)
cv(mod2[[1]] , at = s )
min(mod2[[1]] , at = s)
max(mod2[[1]] , at = s)
moment(mod2[[1]], at = s , center= TRUE , order = 4)
moment(mod2[[1]], at = s , center= FALSE , order = 4)
bconfint(mod2[[1]], at =s)

Sigma <- matrix(c(1, 0.75, 0.75, 1), nrow = 2, ncol = 2)
mod2 <- lapply(1:2, function(i) 
                snssde2d(N=50,drift=fx,diffusion=gx,corr=Sigma,Dt=0.005,M=2,method=meth[i],type="str"))
print(mod2[[1]])				
summary(mod2[[1]])
	       
######## 3d

fx <- expression(4*(-1-x)*y, 4*(1-y)*x, 4*(1-z)*y)
gx <- rep(expression(0.2),3)

mod3 <- lapply(2:length(meth), function(i) 
                snssde3d(N=50,drift=fx,diffusion=gx,M=2,method=meth[i]))
print(mod3[[1]])				
summary(mod3[[1]])
mean(mod3[[1]])
moment(mod3[[1]], center = TRUE , order = 2) ## variance
Median(mod3[[1]])
Mode(mod3[[1]])
quantile(mod3[[1]])
kurtosis(mod3[[1]])
skewness(mod3[[1]])
cv(mod3[[1]])
min(mod3[[1]])
max(mod3[[1]])
moment(mod3[[1]], center= TRUE , order = 4)
moment(mod3[[1]], center= FALSE , order = 4)
				
plot(mod3[[1]])
bconfint(mod3[[1]])
plot3D(mod3[[1]],display="persp",main="3-dim bridge sde")

Sigma <- matrix(c(1,-0.5,-0.25,-0.5,1,0.95,-0.25,0.95,1),nrow=3,ncol=3) 
mod3 <- lapply(1:2, function(i) 
                snssde3d(N=50,drift=fx,diffusion=gx,corr=Sigma,M=2,method=meth[i]))
print(mod3[[1]])				
summary(mod3[[1]])
	       
##


mod3 <- lapply(2:length(meth), function(i) 
                snssde3d(N=50,drift=fx,diffusion=gx,M=2,method=meth[i],type="str"))
print(mod3[[1]])
summary(mod3[[1]])				
bconfint(mod3[[1]])
plot(mod3[[1]],type="n")
lines(mod3[[1]],col=2,lwd=2)
points(mod3[[1]],pch=21,col=5,cex=0.1)

mean(mod3[[1]], at = s)
moment(mod3[[1]], at = s , center = TRUE , order = 2) ## variance
Median(mod3[[1]], at = s)
Mode(mod3[[1]], at = s)
quantile(mod3[[1]] , at = s)
kurtosis(mod3[[1]] , at = s)
skewness(mod3[[1]] , at = s)
cv(mod3[[1]] , at = s )
min(mod3[[1]] , at = s)
max(mod3[[1]] , at = s)
moment(mod3[[1]], at = s , center= TRUE , order = 4)
moment(mod3[[1]], at = s , center= FALSE , order = 4)
bconfint(mod3[[1]], at =s)

Sigma <- matrix(c(1,-0.5,-0.25,-0.5,1,0.95,-0.25,0.95,1),nrow=3,ncol=3) 
mod3 <- lapply(1:2, function(i) 
                snssde3d(N=50,drift=fx,diffusion=gx,corr=Sigma,M=2,method=meth[i],type="str"))
print(mod3[[1]])				
summary(mod3[[1]])
#############################
s= 0.00458
f <- expression( 2*(1-x) )
g <- expression( 1 )

mod1 <- snssde1d(drift=f,diffusion=g,x0=2,N=50,M=1,Dt=0.001)			
plot(mod1,type="n")
lines(mod1,col=2)
points(mod1,cex=0.1,pch=19)

mean(mod1,at = s)
moment(mod1,at = s, center = TRUE , order = 2) ## variance
Median(mod1,at = s)
quantile(mod1,at = s)
kurtosis(mod1,at = s)
skewness(mod1,at = s)
cv(mod1,at = s)
min(mod1,at = s)
max(mod1,at = s)


####

fx <- expression(4*(-1-x) , x)
gx <- expression(0.2 , 0)

mod2 <- snssde2d(drift=fx,diffusion=gx,M=5,N=50,Dt=0.001)
plot(mod2,union = FALSE)
plot(mod2,type="n")
lines(mod2,col=2)
points(mod2,cex=0.1,pch=19)
plot2d(mod2)

mean(mod2,at = s)
moment(mod2,at = s, center = TRUE , order = 2) ## variance
Median(mod2,at = s)
quantile(mod2,at = s)
kurtosis(mod2,at = s)
skewness(mod2,at = s)
cv(mod2,at = s)
min(mod2,at = s)
max(mod2,at = s)

####

fx <- expression(4*(-1-x), 4*(1-y), 4*(1-z))
gx <- rep(expression(0.2),3)

mod3 <- snssde3d(drift=fx,diffusion=gx,N=50,M=5,Dt=0.001)
plot(mod3,type="n")
lines(mod3,col=2)
points(mod3,cex=0.1,pch=19)

mean(mod3,at = s)
moment(mod3,at = s, center = TRUE , order = 2) ## variance
Median(mod3,at = s)
quantile(mod3,at = s)
kurtosis(mod3,at = s)
skewness(mod3,at = s)
cv(mod3,at = s)
min(mod3,at = s)
max(mod3,at = s)


