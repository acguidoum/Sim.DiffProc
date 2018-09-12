options(prompt="R> ",scipen=16,digits=4,warning=FALSE, message=FALSE)
library(Sim.DiffProc)




X <- GBM(N =1000,theta=4,sigma=1)
## Estimation: true theta=c(4,1)
fx <- expression(theta[1]*x)
gx <- expression(theta[2]*x)
#gx2 <- expression(theta[3]*sqrt(x))


pmle <- eval(formals(fitsde.default)$pmle)

fres <- lapply(1:4, function(i) fitsde(data=X,drift=fx,diffusion=gx,
	             pmle=pmle[i],start = list(theta1=1,theta2=1),
				 optim.method = "L-BFGS-B"))
Coef <- data.frame(do.call("cbind",lapply(1:4,function(i) coef(fres[[i]]))))
names(Coef) <- c(pmle)
Summary <- data.frame(do.call("rbind",lapply(1:4,function(i) logLik(fres[[i]]))),
                      do.call("rbind",lapply(1:4,function(i) AIC(fres[[i]]))),
                      do.call("rbind",lapply(1:4,function(i) BIC(fres[[i]]))),
                      row.names=pmle)
names(Summary) <- c("logLik","AIC","BIC")
Coef	
Summary

print(fres[[2]])
summary(fres[[2]])
vcov(fres[[2]])
confint(fres[[2]])