## Thu Apr 30 10:50:58 2020
## Original file Copyright © 2020 A.C. Guidoum, K. Boukhetala
## This file is part of the R package Sim.DiffProc
## Department of Probabilities & Statistics
## Faculty of Mathematics
## University of Science and Technology Houari Boumediene
## BP 32 El-Alia, U.S.T.H.B, Algiers
## Algeria

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## A copy of the GNU General Public License is available at
## http://www.r-project.org/Licenses/
## Unlimited use and distribution (see LICENCE).
###################################################################################################



#####
##### bridgesde1d

bridgesde1d <- function(N, ...)  UseMethod("bridgesde1d")

bridgesde1d.default <- function(N =1000,M=1,x0=0,y=0,t0=0,T=1,Dt,drift,diffusion,
                              alpha=0.5,mu=0.5,type=c("ito","str"), method=c(
                              "euler","milstein","predcorr","smilstein","taylor",
                              "heun","rk1","rk2","rk3"),...)
                     {
    if (!is.numeric(x0) || !is.numeric(y)) 
	             stop("'x0' and 'y' must be numeric")
    if (any(!is.numeric(t0) || !is.numeric(T))) 
	             stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) 
	             stop(" 'N' must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0))  
	             stop(" 'M' must be a positive integer ")
    if (any(!is.expression(drift) || !is.expression(diffusion) )) 
	             stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't' and 'x'")
    if (missing(type)) type <- "ito"
    method <- match.arg(method)
    if (method =="predcorr"){
    if (any(alpha > 1 || alpha < 0)) 
	      stop("please use '0 <= alpha <= 1' ")
    if (any(mu > 1 || mu < 0))       
	      stop("please use '0 <= mu <= 1' ")
                            }
    if (t0 < 0 || T < 0) 
	          stop(" please use positive times! (0 <= t0 < T) ")
    if (missing(Dt)) {
        t <- seq(t0, T, by=Dt)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
		T <- t[N + 1]
    }
    Dt <- (T - t0)/N
	X1 <- snssde1d(N,M,x0,t0,T,Dt,drift,diffusion,alpha,mu,type, method,...)$X
    if (M > 1){X2 <- apply(data.frame(snssde1d(N,M,x0=y,t0,T,Dt,drift,diffusion,alpha,mu,type, method,...)$X),2,rev)
	}else{
	X2 <- rev(snssde1d(N,M,x0=y,t0,T,Dt,drift,diffusion,alpha,mu,type, method,...)$X)}        
    G <- rep(NA,M)
    if (M > 1){
        for(j in 1:M){
              if (X1[1,M] >= X2[1,M]){
                    if (!all(X1[,j] > X2[,j]))
                        G[j] <- min(which((X1[,j]-X2[,j]) <= 0)) - 1
             }else{ if (!all(X1[,j] < X2[,j])) 
                        G[j] <- min(which((X1[,j]-X2[,j]) >= 0)) - 1
						}
                   }
    }else{
             if (X1[1] >= X2[1]){
                    if (!all(X1 > X2))
                        G <- min(which((X1-X2) <= 0)) - 1
            }else{  if (!all(X1 < X2)) 
                        G <- min(which((X1-X2) >= 0)) - 1
						}       
   }
   G[which(G==0)]=NA
   name <- "X"
   name <- if(M > 1) 
             paste("X",1:length(which(!is.na(G))),sep="")
   if (M == 1){ 
       if (is.na(G) ){stop( "A crossing has been no realized,trying again (Repeat)..." )
	   }else{X <- ts(c(X1[1:G],X2[-(1:G)]),start=t0,deltat=deltat(X1),names=name)}
   }else if (M > 1){ if(length(which(is.na(G))) == M){stop( "A crossing has been no realized,trying again (Repeat)..." )
                 }else if (length(which(!is.na(G))) == 1){X <- ts(c(X1[,-which(is.na(G))][1:G[which(!is.na(G))]],X2[,-which(is.na(G))][-(1:G[which(!is.na(G))])]),start=t0,deltat=deltat(X1),names=name)
                 }else if (length(which( is.na(G))) == 0){X <- ts(sapply(1:length(which(!is.na(G))), function(j) c(X1[,j][1:G[!is.na(G)][j]],X2[,j][-(1:G[!is.na(G)][j])])),start=t0,deltat=deltat(X1),names=name)
                 }else{ X1 <- X1[,-c(which(is.na(G)))]
				        X2 <- X2[,-c(which(is.na(G)))]
                        X <- ts(sapply(1:length(which(!is.na(G))), function(j) c(X1[,j][1:G[!is.na(G)][j]],X2[,j][-(1:G[!is.na(G)][j])])),
                                 start=t0,deltat=deltat(X1),names=name)
                       }
   }
    structure(list(X=X,drift=drift[[1]], diffusion=diffusion[[1]],type=type,method=method, 
                   x0=x0,y=y, N=N,M=M,Dt=Dt,t0=t0,T=T,C=G,dim="1d",call=match.call()),class="bridgesde1d")
}

###

print.bridgesde1d <- function(x, digits=NULL, ...)
           {
    class(x) <- "bridgesde1d"
	Ito = "It\xf4"
    Encoding(Ito) <- "latin1"
    if (x$method=="euler")         
	         sch <- "Euler scheme with order 0.5"
    else if (x$method=="milstein") 
	         sch <- "First-order Milstein scheme"
    else if (x$method=="predcorr") 
	         sch <- "Predictor-corrector method with order 1"
    else if (x$method=="smilstein")
	         sch <- "Second-order Milstein scheme"
    else if (x$method=="taylor")   
	         sch <- "Taylor scheme with order 1.5"
    else if (x$method=="heun")     
	         sch <- "Heun scheme with order 2"
    else if (x$method=="rk1")      
	         sch <- "Runge-Kutta method with order 1"
    else if (x$method=="rk2")      
	         sch <- "Runge-Kutta method with order 2"
    else if (x$method=="rk3")      
	         sch <- "Runge-Kutta method with order 3"
	##Dr <- gsub(pattern = 'x', replacement = 'X(t)', x = as.expression(x$drift), ignore.case = F,fixed = T)
    ##DD <- gsub(pattern = 'x', replacement = 'X(t)', x = as.expression(x$diffusion), ignore.case = F,fixed = T)
    Dr <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)))), list(e = x$drift))))   
    DD <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)))), list(e = x$diffusion))))
    if(x$type=="ito"){
    cat(Ito," Bridge Sde 1D:","\n",        
        "\t| dX(t) = ", Dr," * dt + ", DD," * dW(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N = ",format(x$N+1,digits=digits),".","\n",
        "\t| Crossing realized","\t| C = ",format(length(which(!is.na(x$C))),digits=digits)," among ",format(x$M,digits=digits),".","\n",
        "\t| Initial value","\t\t| x0 = ",format(x$x0,digits=digits),".","\n",
        "\t| Ending value","\t\t| y = ",format(x$y,digits=digits),".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",       
        sep="")}else{
    cat("Stratonovich Bridge Sde 1D:","\n",
        "\t| dX(t) = ", Dr," * dt + ", DD," o dW(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N = ",format(x$N+1,digits=digits),".","\n",
        "\t| Crossing realized","\t| C = ",format(length(which(!is.na(x$C))),digits=digits)," among ",format(x$M,digits=digits),".","\n",
        "\t| Initial value","\t\t| x0 = ",format(x$x0,digits=digits),".","\n",
        "\t| Ending value","\t\t| y = ",format(x$y,digits=digits),".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}
    invisible(x)
}

time.bridgesde1d <- function(x,...)
                    {
    class(x) <- "bridgesde1d"
    as.vector(time(x$X))
}

summary.bridgesde1d  <- function(object, at,digits=NULL, ...)
           {   
    class(object) <- "bridgesde1d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
	if (is.null(digits)){digits = base::options()$digits}
    if (any(object$T < at || object$t0 > at) )  
	       stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$C)))
    if (C == 1){ X = matrix(object$X,nrow=length(object$X),ncol=1)}else{X = object$X}
    xx   <- as.vector(X[which(time(object)==as.character(at)),])
    if (length(xx) == 0){
	if (C == 1){ F  <- lapply(1:C,function(i) approxfun(time(object),object$X))}else{
                       F  <- lapply(1:C,function(i) approxfun(time(object),object$X[,i]))}
                       xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    cat("\nMonte-Carlo Statistics for X(t) at time t = ",at,"\n",
	    "\t| Crossing realized ",C," among ",object$M,"\n",
        sep="")
    res <- as.data.frame(matrix(c(mean(xx,na.rm = TRUE),var(xx,na.rm = TRUE),Median(xx),Mode(xx),
                               quantile(xx,0.25,na.rm = TRUE),quantile(xx,0.75,na.rm = TRUE),min(xx,na.rm=TRUE),
                               max(xx,na.rm=TRUE),skewness(xx),kurtosis(xx),cv(xx),moment(xx,order=3,center=TRUE),
                               moment(xx,order=4,center=TRUE),moment(xx,order=5,center=TRUE),moment(xx,order=6,center=TRUE)),
                               ncol=1))
    dimnames(res) <- list(c("Mean","Variance","Median","Mode","First quartile","Third quartile","Minimum","Maximum","Skewness","Kurtosis",
	                         "Coef-variation","3th-order moment","4th-order moment","5th-order moment","6th-order moment"),c(""))
    print(round(res,digits=digits), quote = FALSE, right = TRUE,...)
    invisible(object)
}


mean.bridgesde1d <- function(x,at,...)
                    {
    class(x) <- "bridgesde1d"
    if (missing(at)) {at = (as.numeric(x$T)-as.numeric(x$t0))/2}
    if (any(x$T < at || x$t0 > at) )  
	         stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(x$C)))
    if (C == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (C==1){ F  <- lapply(1:C,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:C,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    mean(xx,na.rm = TRUE,...)
}

cv.bridgesde1d <- function(x,at,...)
                    {
    class(x) <- "bridgesde1d"
    if (missing(at)) {at = (as.numeric(x$T)-as.numeric(x$t0))/2}
    if (any(x$T < at || x$t0 > at) )  
	      stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(x$C)))
    if (C == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (C==1){ F  <- lapply(1:C,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:C,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    cv(xx,...)
}

min.bridgesde1d <- function(x,at,...)
                    {
    class(x) <- "bridgesde1d"
    if (missing(at)) {at = (as.numeric(x$T)-as.numeric(x$t0))/2}
    if (any(x$T < at || x$t0 > at) )  
	     stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(x$C)))
    if (C == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (C==1){ F  <- lapply(1:C,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:C,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    min(xx,na.rm = TRUE)
}

max.bridgesde1d <- function(x,at,...)
                    {
    class(x) <- "bridgesde1d"
    if (missing(at)) {at = (as.numeric(x$T)-as.numeric(x$t0))/2}
    if (any(x$T < at || x$t0 > at) )  
	     stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(x$C)))
    if (C == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (C==1){ F  <- lapply(1:C,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:C,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    max(xx,na.rm = TRUE)
}

skewness.bridgesde1d <- function(x,at,...)
                    {
    class(x) <- "bridgesde1d"
    if (missing(at)) {at = (as.numeric(x$T)-as.numeric(x$t0))/2}
    if (any(x$T < at || x$t0 > at) )  
	      stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(x$C)))
    if (C == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (C==1){ F  <- lapply(1:C,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:C,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
	skewness(xx)
}

kurtosis.bridgesde1d <- function(x,at,...)
                    {
    class(x) <- "bridgesde1d"
    if (missing(at)) {at = (as.numeric(x$T)-as.numeric(x$t0))/2}
    if (any(x$T < at || x$t0 > at) )  
	      stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(x$C)))
    if (C == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (C==1){ F  <- lapply(1:C,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:C,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
	kurtosis(xx)
}

Median.bridgesde1d <- function(x,at,...)
                    {
    class(x) <- "bridgesde1d"
    if (missing(at)) {at = (as.numeric(x$T)-as.numeric(x$t0))/2}
    if (any(x$T < at || x$t0 > at) )  
	       stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(x$C)))
    if (C == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (C==1){ F  <- lapply(1:C,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:C,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
	Median(xx)
}

Mode.bridgesde1d <- function(x,at,...)
                    {
    class(x) <- "bridgesde1d"
    if (missing(at)) {at = (as.numeric(x$T)-as.numeric(x$t0))/2}
    if (any(x$T < at || x$t0 > at) )  
	      stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(x$C)))
    if (C == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (C==1){ F  <- lapply(1:C,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:C,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
	Mode(xx)
}

quantile.bridgesde1d <- function(x,at,...)
                    {
    class(x) <- "bridgesde1d"
    if (missing(at)) {at = (as.numeric(x$T)-as.numeric(x$t0))/2}
    if (any(x$T < at || x$t0 > at) )  
	      stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(x$C)))
    if (C == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (C==1){ F  <- lapply(1:C,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:C,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    quantile(xx,...)
}

moment.bridgesde1d <- function(x,at,...)
                    {
    class(x) <- "bridgesde1d"
    if (missing(at)) {at = (as.numeric(x$T)-as.numeric(x$t0))/2}
    if (any(x$T < at || x$t0 > at) )  
	       stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(x$C)))
    if (C == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (C==1){ F  <- lapply(1:C,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:C,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    moment(xx,...)
}

bconfint.bridgesde1d <- function(x,at,...)
                    {
    class(x) <- "bridgesde1d"
    if (missing(at)) {at = (as.numeric(x$T)-as.numeric(x$t0))/2}
    if (any(x$T < at || x$t0 > at) )  
	        stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(x$C)))
    if (C == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (C==1){ F  <- lapply(1:C,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:C,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    bconfint(xx,...)
}

##
## Plot

plot.bridgesde1d <- function(x,...)
                 {
    class(x) <- "bridgesde1d"
    X <- x$X
	if (x$M > 10) {plot(X,plot.type="single",...)}else{plot(X,...)}
}

lines.bridgesde1d <- function(x,...)
                 {
    class(x) <- "bridgesde1d"
    X <- x$X
	if (length(which(!is.na(x$C))) >=2){
    for (i in 1:dim(X)[2]){
    lines(time(x),X[,i],...)}}else{
	lines(time(x),X,...)}
}


points.bridgesde1d <- function(x,...)
                 {
    class(x) <- "bridgesde1d"
    X <- x$X
	if (length(which(!is.na(x$C))) >=2){
    for (i in 1:dim(X)[2]){
    points(time(x),X[,i],...)}}else{
	points(time(x),X,...)}
}


#####
##### bridgesde2d

bridgesde2d <- function(N, ...)  UseMethod("bridgesde2d")

bridgesde2d.default <- function(N =1000,M=1,x0=c(0,0),y=c(0,0),t0=0,T=1,Dt,drift,diffusion,
                              corr = NULL,alpha=0.5,mu=0.5,type=c("ito","str"), method=c(
                              "euler","milstein","predcorr","smilstein","taylor",
                              "heun","rk1","rk2","rk3"),...)
                     {
    if (!is.numeric(x0) || length(x0) != 2) 
	                  stop("'x0' must be numeric, and length(x0) = 2 ")
    if (!is.numeric(y)  || length(y) != 2) 
	                  stop("'y' must be numeric, and length(y) = 2 ")
    if (any(!is.numeric(t0) || !is.numeric(T))) 
	                  stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) 
	                  stop(" 'N' must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0))  
	                  stop(" 'M' must be a positive integer ")
	if (length(drift) !=2 ) 
	                  stop("drift must be expression 2d (vector of 2 expression)")
	if (length(diffusion) !=2 ) 
	                  stop("diffusion must be expression 2d (vector of 2 expression)")
    if (any(!is.expression(drift) || !is.expression(diffusion) )) 
	                  stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't', 'x' and 'y'")
	if (!is.null(corr)) {
        if (any(!is.matrix(corr) || det(corr) <= 0 || nrow(corr)!=2 || ncol(corr)!=2 || !isSymmetric(corr)))
                     stop("the correlation structure of W1(t) and W2(t) must be a real symmetric\n  positive-definite square matrix of dimension 2") 		   
     }
    if (missing(type)) type <- "ito"
    method <- match.arg(method)
	if (!is.null(corr) && method != "euler" && method != "milstein") {method="euler"}
    if (method =="predcorr"){
    if (any(alpha > 1 || alpha < 0)) 
	      stop("please use '0 <= alpha <= 1' ")
    if (any(mu > 1 || mu < 0))       
	      stop("please use '0 <= mu <= 1' ")
                            }
    if (t0 < 0 || T < 0) 
	         stop(" please use positive times! (0 <= t0 < T) ")
    if (missing(Dt)) {
        t <- seq(t0, T, by=Dt)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
		T <- t[N + 1]
    }
    Dt <- (T - t0)/N
    X1 <- snssde2d(N,M,x0=c(x0[1],x0[2]),t0,T,Dt,drift,diffusion,corr,alpha,mu,type, method,...)$X
    Y1 <- snssde2d(N,M,x0=c(x0[1],x0[2]),t0,T,Dt,drift,diffusion,corr,alpha,mu,type, method,...)$Y
    if (M > 1){X2 <- apply(data.frame(snssde2d(N,M,x0=c(y[1],y[2]),t0,T,Dt,drift,diffusion,corr,alpha,mu,type, method,...)$X),2,rev)
	           Y2 <- apply(data.frame(snssde2d(N,M,x0=c(y[1],y[2]),t0,T,Dt,drift,diffusion,corr,alpha,mu,type, method,...)$Y),2,rev)
    }else{X2 <- rev(snssde2d(N,M,x0=c(y[1],y[2]),t0,T,Dt,drift,diffusion,corr,alpha,mu,type, method,...)$X)
	      Y2 <- rev(snssde2d(N,M,x0=c(y[1],y[2]),t0,T,Dt,drift,diffusion,corr,alpha,mu,type, method,...)$Y)
	    }        
    Gx = Gy <- rep(NA,M)
    if (M > 1){
        for(j in 1:M){
              if (X1[1,M] >= X2[1,M]){
                    if (!all(X1[,j] > X2[,j]))
                        Gx[j] <- min(which((X1[,j]-X2[,j]) <= 0)) - 1
             }else{ if (!all(X1[,j] < X2[,j])) 
                        Gx[j] <- min(which((X1[,j]-X2[,j]) >= 0)) - 1
						}
                   }
    }else{
             if (X1[1] >= X2[1]){
                    if (!all(X1 > X2))
                        Gx <- min(which((X1-X2) <= 0)) - 1
            }else{  if (!all(X1 < X2)) 
                        Gx <- min(which((X1-X2) >= 0)) - 1
						}       
   }
   if (M > 1){
        for(j in 1:M){
              if (Y1[1,M] >= Y2[1,M]){
                    if (!all(Y1[,j] > Y2[,j]))
                        Gy[j] <- min(which((Y1[,j]-Y2[,j]) <= 0)) - 1
             }else{ if (!all(Y1[,j] < Y2[,j])) 
                        Gy[j] <- min(which((Y1[,j]-Y2[,j]) >= 0)) - 1
						}
                   }
    }else{
             if (Y1[1] >= Y2[1]){
                    if (!all(Y1 > Y2))
                        Gy <- min(which((Y1-Y2) <= 0)) - 1
            }else{  if (!all(Y1 < Y2)) 
                        Gy <- min(which((Y1-Y2) >= 0)) - 1
						}       
   }
   Gx[which(Gx==0)]=NA
   Gy[which(Gy==0)]=NA
   G <- na.omit(data.frame(Gx,Gy))
   Gx <- G$Gx;Gy <- G$Gy
   namex <- "X"
   namey <- "Y"
   namex <- if(M > 1) paste("X",1:length(which(!is.na(Gx))),sep="")
   namey <- if(M > 1) paste("Y",1:length(which(!is.na(Gy))),sep="")
   if (M == 1){ if (is.na(Gx) ){stop( "A crossing has been no realized,trying again (Repeat)..." )}else{X <- ts(c(X1[1:Gx],X2[-(1:Gx)]),start=t0,deltat=deltat(X1),names=namex)}
   }else if (M > 1){ if(length(which(is.na(Gx))) == M){stop( "A crossing has been no realized,trying again (Repeat)..." )
                 }else if (length(which(!is.na(Gx))) == 1){X <- ts(c(X1[,-which(is.na(Gx))][1:Gx[which(!is.na(Gx))]],X2[,-which(is.na(Gx))][-(1:Gx[which(!is.na(Gx))])]),start=t0,deltat=deltat(X1),names=namex)
                 }else if (length(which( is.na(Gx))) == 0){X <- ts(sapply(1:length(which(!is.na(Gx))), function(j) c(X1[,j][1:Gx[!is.na(Gx)][j]],X2[,j][-(1:Gx[!is.na(Gx)][j])])),start=t0,deltat=deltat(X1),names=namex)
                 }else{ X1 <- X1[,-c(which(is.na(Gx)))]
				        X2 <- X2[,-c(which(is.na(Gx)))]
                        X <- ts(sapply(1:length(which(!is.na(Gx))), function(j) c(X1[,j][1:Gx[!is.na(Gx)][j]],X2[,j][-(1:Gx[!is.na(Gx)][j])])),
                                 start=t0,deltat=deltat(X1),names=namex)
                       }
   }
   if (M == 1){ if (is.na(Gy) ){stop( "A crossing has been no realized,trying again (Repeat)..." )}else{Y <- ts(c(Y1[1:Gy],Y2[-(1:Gy)]),start=t0,deltat=deltat(Y1),names=namey)}
   }else if (M > 1){ if(length(which(is.na(Gy))) == M){stop( "A crossing has been no realized,trying again (Repeat)..." )
                 }else if (length(which(!is.na(Gy))) == 1){Y <- ts(c(Y1[,-which(is.na(Gy))][1:Gy[which(!is.na(Gy))]],Y2[,-which(is.na(Gy))][-(1:Gy[which(!is.na(Gy))])]),start=t0,deltat=deltat(Y1),names=namey)
                 }else if (length(which( is.na(Gy))) == 0){Y <- ts(sapply(1:length(which(!is.na(Gy))), function(j) c(Y1[,j][1:Gy[!is.na(Gy)][j]],Y2[,j][-(1:Gy[!is.na(Gy)][j])])),start=t0,deltat=deltat(Y1),names=namey)
                 }else{ Y1 <- Y1[,-c(which(is.na(Gy)))]
				        Y2 <- Y2[,-c(which(is.na(Gy)))]
                        Y <- ts(sapply(1:length(which(!is.na(Gy))), function(j) c(Y1[,j][1:Gy[!is.na(Gy)][j]],Y2[,j][-(1:Gy[!is.na(Gy)][j])])),
                                 start=t0,deltat=deltat(Y1),names=namey)
                       }
   }
    structure(list(X=X,Y=Y, driftx=drift[[1]], diffx=diffusion[[1]],drifty=drift[[2]], diffy=diffusion[[2]],corrmat=corr,type=type,method=method, 
                   x0=x0,y=y, N=N,M=M,Dt=Dt,t0=t0,T=T,Cx=Gx,Cy=Gy,dim="2d",call=match.call()),class="bridgesde2d")
}


###

print.bridgesde2d <- function(x, digits=NULL, ...)
           {
    class(x) <- "bridgesde2d"
	Ito = "It\xf4"
    Encoding(Ito) <- "latin1"
    if (x$method=="euler")         
	         sch <- "Euler scheme with order 0.5"
    else if (x$method=="milstein") 
	         sch <- "First-order Milstein scheme"
    else if (x$method=="predcorr") 
	         sch <- "Predictor-corrector method with order 1"
    else if (x$method=="smilstein")
	         sch <- "Second-order Milstein scheme"
    else if (x$method=="taylor")   
	         sch <- "Taylor scheme with order 1.5"
    else if (x$method=="heun")     
	         sch <- "Heun scheme with order 2"
    else if (x$method=="rk1")      
	         sch <- "Runge-Kutta method with order 1"
    else if (x$method=="rk2")      
	         sch <- "Runge-Kutta method with order 2"
    else if (x$method=="rk3")      
	         sch <- "Runge-Kutta method with order 3"
	Drx <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)))), list(e = x$driftx))))   
    DDx <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)))), list(e = x$diffx))))
	Dry <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)))), list(e = x$drifty))))   
    DDy <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)))), list(e = x$diffy)))) 
    if (is.null(x$corrmat)){ 
    if(x$type=="ito"){
    cat(Ito," Bridge Sde 2D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," * dW1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," * dW2(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N = ",format(x$N+1,digits=digits),".","\n",
        "\t| Crossing realized","\t| C = ",format(length(which(!is.na(x$Cx))),digits=digits)," among ",format(x$M,digits=digits),".","\n",
        "\t| Initial values","\t| x0 = ","(",format(x$x0[1],digits=digits),",",format(x$x0[2],digits=digits),")",".","\n",
        "\t| Ending values","\t\t| y = ","(",format(x$y[1],digits=digits),",",format(x$y[2],digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}else{
    cat("Stratonovich Bridge Sde 2D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," o dW1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," o dW2(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N = ",format(x$N+1,digits=digits),".","\n",
        "\t| Crossing realized","\t| C = ",format(length(which(!is.na(x$Cx))),digits=digits)," among ",format(x$M,digits=digits),".","\n",
        "\t| Initial values","\t| x0 = ","(",format(x$x0[1],digits=digits),",",format(x$x0[2],digits=digits),")",".","\n",
        "\t| Ending values","\t\t| y = ","(",format(x$y[1],digits=digits),",",format(x$y[2],digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}} else {
    if(x$type=="ito"){
    cat(Ito," Bridge Sde 2D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," * dB1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," * dB2(t)","\n",
		"\t| Correlation structure:",sep="")		 
	    prmatrix(x$corrmat,rowlab = rep("            ", 2), collab = rep("", 2),digits=digits)
    cat("Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N = ",format(x$N+1,digits=digits),".","\n",
        "\t| Crossing realized","\t| C = ",format(length(which(!is.na(x$Cx))),digits=digits)," among ",format(x$M,digits=digits),".","\n",
        "\t| Initial values","\t| x0 = ","(",format(x$x0[1],digits=digits),",",format(x$x0[2],digits=digits),")",".","\n",
        "\t| Ending values","\t\t| y = ","(",format(x$y[1],digits=digits),",",format(x$y[2],digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}else{
    cat("Stratonovich Bridge Sde 2D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," o dB1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," o dB2(t)","\n",
		"\t| Correlation structure:",sep="")		 
	    prmatrix(x$corrmat,rowlab = rep("            ", 2), collab = rep("", 2),digits=digits)
    cat("Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N = ",format(x$N+1,digits=digits),".","\n",
        "\t| Crossing realized","\t| C = ",format(length(which(!is.na(x$Cx))),digits=digits)," among ",format(x$M,digits=digits),".","\n",
        "\t| Initial values","\t| x0 = ","(",format(x$x0[1],digits=digits),",",format(x$x0[2],digits=digits),")",".","\n",
        "\t| Ending values","\t\t| y = ","(",format(x$y[1],digits=digits),",",format(x$y[2],digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}		
		}
    invisible(x)
}

##
## Plot

plot.bridgesde2d <- function(x,...) .plot.bridgesde2d(x,...)

lines.bridgesde2d <- function(x,...)
                 {
    class(x) <- "bridgesde2d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}				  
    for (i in 1:length(which(!is.na(x$Cx)))){lines(time(x),X[,i],...)}
    for (i in 1:length(which(!is.na(x$Cy)))){lines(time(x),Y[,i],...)}
}

points.bridgesde2d <- function(x,...)
                 {
    class(x) <- "bridgesde2d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}				  
    for (i in 1:length(which(!is.na(x$Cx)))){points(time(x),X[,i],...)}
    for (i in 1:length(which(!is.na(x$Cy)))){points(time(x),Y[,i],...)}
}

plot2d.bridgesde2d <- function(x,...) .plot2d.bridgesde2d(x,...)

lines2d.bridgesde2d <- function(x,...)
                 {
    class(x) <- "bridgesde2d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    X <- X[,1]
    Y <- Y[,1]
    lines2d(X,Y,...)
}

points2d.bridgesde2d <- function(x,...)
                 {
    class(x) <- "bridgesde2d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    X <- X[,1]
    Y <- Y[,1]
    points2d(X,Y,...)
}

##
## summary

summary.bridgesde2d <- function(object,at,digits=NULL,...)
                    {			
    class(object) <- "bridgesde2d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
	if (is.null(digits)){digits = base::options()$digits}
    if (any(object$T < at || object$t0 > at) )  
	     stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
	if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    cat("\nMonte-Carlo Statistics for (X(t),Y(t)) at time t = ",at,"\n",
	    "\t| Crossing realized ",C," among ",object$M,"\n",
        sep="")
    res <- as.data.frame(matrix(c(mean(x,na.rm = TRUE),var(x,na.rm = TRUE),Median(x),Mode(x),
                               quantile(x,0.25,na.rm = TRUE),quantile(x,0.75,na.rm = TRUE),min(x,na.rm=TRUE),
                               max(x,na.rm=TRUE),skewness(x),kurtosis(x),cv(x),moment(x,order=3,center=TRUE),
                               moment(x,order=4,center=TRUE),moment(x,order=5,center=TRUE),moment(x,order=6,center=TRUE),
                               mean(y,na.rm = TRUE),var(y,na.rm = TRUE),Median(y),Mode(y),
                               quantile(y,0.25,na.rm = TRUE),quantile(y,0.75,na.rm = TRUE),min(y,na.rm=TRUE),
                               max(y,na.rm=TRUE),skewness(y),kurtosis(y),cv(y),moment(y,order=3,center=TRUE),
                               moment(y,order=4,center=TRUE),moment(y,order=5,center=TRUE),moment(y,order=6,center=TRUE)),
                               ncol=2))
	dimnames(res) <- list(c("Mean","Variance","Median","Mode","First quartile","Third quartile","Minimum","Maximum","Skewness","Kurtosis",
	                         "Coef-variation","3th-order moment","4th-order moment","5th-order moment","6th-order moment"),c("X","Y"))
    print(round(res,digits=digits), quote = FALSE, right = TRUE,...)
    invisible(object)
}



mean.bridgesde2d <- function(x,at,...)
                    {
	object <- x
    class(object) <- "bridgesde2d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
    if (any(object$T < at || object$t0 > at) )  
	     stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
	if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
	return(c(mean(x,na.rm = TRUE,...),mean(y,na.rm = TRUE,...)))
}


cv.bridgesde2d  <- function(x,at,...)
                    {
	object <- x
    class(object) <- "bridgesde2d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
    if (any(object$T < at || object$t0 > at) )  
	      stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
	if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
	return(c(cv(x,...),cv(y,...)))
}


max.bridgesde2d  <- function(x,at,...)
                    {
	object <- x
    class(object) <- "bridgesde2d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
    if (any(object$T < at || object$t0 > at) )  
	      stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
	if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
	return(c(max(x,na.rm = TRUE),max(y,na.rm = TRUE)))
}

min.bridgesde2d  <- function(x,at,...)
                    {
	object <- x
    class(object) <- "bridgesde2d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
    if (any(object$T < at || object$t0 > at) )  
	         stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
	if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
	return(c(min(x,na.rm = TRUE),min(y,na.rm = TRUE)))
}

skewness.bridgesde2d <- function(x,at,...)
                    {
	object <- x
    class(object) <- "bridgesde2d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
    if (any(object$T < at || object$t0 > at) )  
	      stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
	if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
	return(c(skewness(x),skewness(y)))
}

kurtosis.bridgesde2d <- function(x,at,...)
                    {
	object <- x
    class(object) <- "bridgesde2d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
    if (any(object$T < at || object$t0 > at) )  
	     stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
	if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
	return(c(kurtosis(x),kurtosis(y)))
}

Median.bridgesde2d <- function(x, at,...)
                    {
	object <- x
    class(object) <- "bridgesde2d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
    if (any(object$T < at || object$t0 > at) )  
	     stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
	if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
	return(c(Median(x),Median(y)))
}

Mode.bridgesde2d <- function(x, at,...)
                    {
	object <- x
    class(object) <- "bridgesde2d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
    if (any(object$T < at || object$t0 > at) )  
	     stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
	if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
	return(c(Mode(x),Mode(y)))
}

quantile.bridgesde2d <- function(x,at,...)
                    {
	object <- x
    class(object) <- "bridgesde2d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
    if (any(object$T < at || object$t0 > at) )  
	     stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
	if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
	return(list(x=quantile(x,...),y=quantile(y,...)))
}

moment.bridgesde2d <- function(x,at,...)
                    {
	object <- x
    class(object) <- "bridgesde2d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
    if (any(object$T < at || object$t0 > at) )  
	     stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
	if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
	return(c(moment(x,...),moment(y,...)))
}

bconfint.bridgesde2d <- function(x,at,...)
                    {
	object <- x
    class(object) <- "bridgesde2d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
    if (any(object$T < at || object$t0 > at) )  
	     stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
	if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
	return(list(x=bconfint(x,...),y=bconfint(y,...)))
}


time.bridgesde2d <- function(x,...)
                    {
    class(x) <- "bridgesde2d"
    as.vector(time(x$X))
}

#####
##### bridgesde3d

bridgesde3d <- function(N, ...)  UseMethod("bridgesde3d")

bridgesde3d.default <- function(N =1000,M=1,x0=c(0,0,0),y=c(0,0,0),t0=0,T=1,Dt,drift,diffusion,
                               corr=NULL, alpha=0.5,mu=0.5,type=c("ito","str"), method=c(
                              "euler","milstein","predcorr","smilstein","taylor",
                              "heun","rk1","rk2","rk3"),...)
                     {
    if (!is.numeric(x0) || length(x0) != 3) 
	                            stop("'x0' must be numeric, and length(x0) = 3 ")
    if (!is.numeric(y)  || length(y) != 3) 
	                            stop("'y' must be numeric, and length(y) = 3 ")
    if (any(!is.numeric(t0) || !is.numeric(T))) 
	                            stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) 
	                            stop(" 'N' must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0))  
	                            stop(" 'M' must be a positive integer ")
	if (length(drift) !=3 ) 
	               stop("drift must be expression 3d (vector of 3 expression)")
	if (length(diffusion) !=3 ) 
	               stop("diffusion must be expression 3d (vector of 3 expression)")
    if (any(!is.expression(drift) || !is.expression(diffusion) )) 
	               stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't', 'x', 'y' and 'z'")
    if (missing(type)) type <- "ito"
	if (!is.null(corr)) {
        if (any(!is.matrix(corr) || det(corr) <= 0 || nrow(corr)!=3 || ncol(corr)!=3 || !isSymmetric(corr)))
                     stop("the correlation structure of W1(t), W2(t) and W3(t) must be a real symmetric\n  positive-definite square matrix of dimension 3") 	   
     }
    method <- match.arg(method)
	if (!is.null(corr) && method != "euler" && method != "milstein") {method="euler"}
    if (method =="predcorr"){
    if (any(alpha > 1 || alpha < 0)) 
	      stop("please use '0 <= alpha <= 1' ")
    if (any(mu > 1 || mu < 0))       
	      stop("please use '0 <= mu <= 1' ")
                            }
    if ( t0 < 0 || T < 0 ) 
	      stop(" please use positive times! (0 <= t0 < T) ")
    if (missing(Dt)) {
        t <- seq(t0, T, by=Dt)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
		T <- t[N + 1]
    } 
    Dt <- (T - t0)/N	
	X1 <- snssde3d(N,M,x0=c(x0[1],x0[2],x0[3]),t0,T,Dt,drift,diffusion,corr,alpha,mu,type, method,...)$X
    Y1 <- snssde3d(N,M,x0=c(x0[1],x0[2],x0[3]),t0,T,Dt,drift,diffusion,corr,alpha,mu,type, method,...)$Y
	Z1 <- snssde3d(N,M,x0=c(x0[1],x0[2],x0[3]),t0,T,Dt,drift,diffusion,corr,alpha,mu,type, method,...)$Z
    if (M > 1){X2 <- apply(data.frame(snssde3d(N,M,x0=c(y[1],y[2],y[3]),t0,T,Dt,drift,diffusion,corr,alpha,mu,type, method,...)$X),2,rev)
	           Y2 <- apply(data.frame(snssde3d(N,M,x0=c(y[1],y[2],y[3]),t0,T,Dt,drift,diffusion,corr,alpha,mu,type, method,...)$Y),2,rev)
			   Z2 <- apply(data.frame(snssde3d(N,M,x0=c(y[1],y[2],y[3]),t0,T,Dt,drift,diffusion,corr,alpha,mu,type, method,...)$Z),2,rev)
    }else{X2 <- rev(snssde3d(N,M,x0=c(y[1],y[2],y[3]),t0,T,Dt,drift,diffusion,corr,alpha,mu,type, method,...)$X)
	      Y2 <- rev(snssde3d(N,M,x0=c(y[1],y[2],y[3]),t0,T,Dt,drift,diffusion,corr,alpha,mu,type, method,...)$Y)
	      Z2 <- rev(snssde3d(N,M,x0=c(y[1],y[2],y[3]),t0,T,Dt,drift,diffusion,corr,alpha,mu,type, method,...)$Z)		  
	    }        
    Gx = Gy = Gz <- rep(NA,M)
    if (M > 1){
        for(j in 1:M){
              if (X1[1,M] >= X2[1,M]){
                    if (!all(X1[,j] > X2[,j]))
                        Gx[j] <- min(which((X1[,j]-X2[,j]) <= 0)) - 1
             }else{ if (!all(X1[,j] < X2[,j])) 
                        Gx[j] <- min(which((X1[,j]-X2[,j]) >= 0)) - 1
						}
                   }
    }else{
             if (X1[1] >= X2[1]){
                    if (!all(X1 > X2))
                        Gx <- min(which((X1-X2) <= 0)) - 1
            }else{  if (!all(X1 < X2)) 
                        Gx <- min(which((X1-X2) >= 0)) - 1
						}       
   }
   if (M > 1){
        for(j in 1:M){
              if (Y1[1,M] >= Y2[1,M]){
                    if (!all(Y1[,j] > Y2[,j]))
                        Gy[j] <- min(which((Y1[,j]-Y2[,j]) <= 0)) - 1
             }else{ if (!all(Y1[,j] < Y2[,j])) 
                        Gy[j] <- min(which((Y1[,j]-Y2[,j]) >= 0)) - 1
						}
                   }
    }else{
             if (Y1[1] >= Y2[1]){
                    if (!all(Y1 > Y2))
                        Gy <- min(which((Y1-Y2) <= 0)) - 1
            }else{  if (!all(Y1 < Y2)) 
                        Gy <- min(which((Y1-Y2) >= 0)) - 1
						}       
   }
   if (M > 1){
        for(j in 1:M){
              if (Z1[1,M] >= Z2[1,M]){
                    if (!all(Z1[,j] > Z2[,j]))
                        Gz[j] <- min(which((Z1[,j]-Z2[,j]) <= 0)) - 1
             }else{ if (!all(Z1[,j] < Z2[,j])) 
                        Gz[j] <- min(which((Z1[,j]-Z2[,j]) >= 0)) - 1
						}
                   }
    }else{
             if (Z1[1] >= Z2[1]){
                    if (!all(Z1 > Z2))
                        Gz <- min(which((Z1-Z2) <= 0)) - 1
            }else{  if (!all(Z1 < Z2)) 
                        Gz <- min(which((Z1-Z2) >= 0)) - 1
						}       
   }   
   Gx[which(Gx==0)]=NA
   Gy[which(Gy==0)]=NA
   Gz[which(Gz==0)]=NA
   G <- na.omit(data.frame(Gx,Gy,Gz))
   Gx <- G$Gx;Gy <- G$Gy; Gz <- G$Gz
   namex <- "X"
   namey <- "Y"
   namey <- "Z"
   namex <- if(M > 1) paste("X",1:length(which(!is.na(Gx))),sep="")
   namey <- if(M > 1) paste("Y",1:length(which(!is.na(Gy))),sep="")
   namez <- if(M > 1) paste("Z",1:length(which(!is.na(Gz))),sep="")
   if (M == 1){ if (is.na(Gx) ){stop( "A crossing has been no realized,trying again (Repeat)..." )}else{X <- ts(c(X1[1:Gx],X2[-(1:Gx)]),start=t0,deltat=deltat(X1),names=namex)}
   }else if (M > 1){ if(length(which(is.na(Gx))) == M){stop( "A crossing has been no realized,trying again (Repeat)..." )
                 }else if (length(which(!is.na(Gx))) == 1){X <- ts(c(X1[,-which(is.na(Gx))][1:Gx[which(!is.na(Gx))]],X2[,-which(is.na(Gx))][-(1:Gx[which(!is.na(Gx))])]),start=t0,deltat=deltat(X1),names=namex)
                 }else if (length(which( is.na(Gx))) == 0){X <- ts(sapply(1:length(which(!is.na(Gx))), function(j) c(X1[,j][1:Gx[!is.na(Gx)][j]],X2[,j][-(1:Gx[!is.na(Gx)][j])])),start=t0,deltat=deltat(X1),names=namex)
                 }else{ X1 <- X1[,-c(which(is.na(Gx)))]
				        X2 <- X2[,-c(which(is.na(Gx)))]
                        X <- ts(sapply(1:length(which(!is.na(Gx))), function(j) c(X1[,j][1:Gx[!is.na(Gx)][j]],X2[,j][-(1:Gx[!is.na(Gx)][j])])),
                                 start=t0,deltat=deltat(X1),names=namex)
                       }
   }
   if (M == 1){ if (is.na(Gy) ){stop( "A crossing has been no realized,trying again (Repeat)..." )}else{Y <- ts(c(Y1[1:Gy],Y2[-(1:Gy)]),start=t0,deltat=deltat(Y1),names=namey)}
   }else if (M > 1){ if(length(which(is.na(Gy))) == M){stop( "A crossing has been no realized,trying again (Repeat)..." )
                 }else if (length(which(!is.na(Gy))) == 1){Y <- ts(c(Y1[,-which(is.na(Gy))][1:Gy[which(!is.na(Gy))]],Y2[,-which(is.na(Gy))][-(1:Gy[which(!is.na(Gy))])]),start=t0,deltat=deltat(Y1),names=namey)
                 }else if (length(which( is.na(Gy))) == 0){Y <- ts(sapply(1:length(which(!is.na(Gy))), function(j) c(Y1[,j][1:Gy[!is.na(Gy)][j]],Y2[,j][-(1:Gy[!is.na(Gy)][j])])),start=t0,deltat=deltat(Y1),names=namey)
                 }else{ Y1 <- Y1[,-c(which(is.na(Gy)))]
				        Y2 <- Y2[,-c(which(is.na(Gy)))]
                        Y <- ts(sapply(1:length(which(!is.na(Gy))), function(j) c(Y1[,j][1:Gy[!is.na(Gy)][j]],Y2[,j][-(1:Gy[!is.na(Gy)][j])])),
                                 start=t0,deltat=deltat(Y1),names=namey)
                       }
   }
   if (M == 1){ if (is.na(Gz) ){stop( "A crossing has been no realized,trying again (Repeat)..." )}else{Z <- ts(c(Z1[1:Gz],Z2[-(1:Gz)]),start=t0,deltat=deltat(Z1),names=namez)}
   }else if (M > 1){ if(length(which(is.na(Gz))) == M){stop( "A crossing has been no realized,trying again (Repeat)..." )
                 }else if (length(which(!is.na(Gz))) == 1){Z <- ts(c(Z1[,-which(is.na(Gz))][1:Gz[which(!is.na(Gz))]],Z2[,-which(is.na(Gz))][-(1:Gz[which(!is.na(Gz))])]),start=t0,deltat=deltat(Z1),names=namez)
                 }else if (length(which( is.na(Gz))) == 0){Z <- ts(sapply(1:length(which(!is.na(Gz))), function(j) c(Z1[,j][1:Gz[!is.na(Gz)][j]],Z2[,j][-(1:Gz[!is.na(Gz)][j])])),start=t0,deltat=deltat(Z1),names=namez)
                 }else{ Z1 <- Z1[,-c(which(is.na(Gz)))]
				        Z2 <- Z2[,-c(which(is.na(Gz)))]
                        Z <- ts(sapply(1:length(which(!is.na(Gz))), function(j) c(Z1[,j][1:Gz[!is.na(Gz)][j]],Z2[,j][-(1:Gz[!is.na(Gz)][j])])),
                                 start=t0,deltat=deltat(Z1),names=namez)
                       }
   }
   structure(list(X=X,Y=Y,Z=Z, driftx=drift[[1]], diffx=diffusion[[1]],drifty=drift[[2]], diffy=diffusion[[2]],driftz=drift[[3]],diffz=diffusion[[3]],
                   type=type,corrmat=corr,method=method, x0=x0,y=y, N=N,M=M,Dt=Dt,t0=t0,T=T,Cx=Gx,Cy=Gy,Cz=Gz,dim="3d",call=match.call()),class="bridgesde3d")
}

###

print.bridgesde3d <- function(x, digits=NULL, ...)
           {
    class(x) <- "bridgesde3d"
	Ito = "It\xf4"
    Encoding(Ito) <- "latin1"
    if (x$method=="euler")         
	         sch <- "Euler scheme with order 0.5"
    else if (x$method=="milstein") 
	         sch <- "First-order Milstein scheme"
    else if (x$method=="predcorr") 
	         sch <- "Predictor-corrector method with order 1"
    else if (x$method=="smilstein")
	         sch <- "Second-order Milstein scheme"
    else if (x$method=="taylor")   
	         sch <- "Taylor scheme with order 1.5"
    else if (x$method=="heun")     
	         sch <- "Heun scheme with order 2"
    else if (x$method=="rk1")      
	         sch <- "Runge-Kutta method with order 1"
    else if (x$method=="rk2")      
	         sch <- "Runge-Kutta method with order 2"
    else if (x$method=="rk3")      
	         sch <- "Runge-Kutta method with order 3"
    Drx <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$driftx))))   
    DDx <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$diffx))))
	Dry <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$drifty))))   
    DDy <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$diffy))))
    Drz <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$driftz))))   
    DDz <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$diffz))))
    if (is.null(x$corrmat)){ 
    if(x$type=="ito"){
    cat(Ito," Bridge Sde 3D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," * dW1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," * dW2(t)","\n",
        "\t| dZ(t) = ", Drz," * dt + ", DDz," * dW3(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N = ",format(x$N+1,digits=digits),".","\n",
        "\t| Crossing realized","\t| C = ",format(length(which(!is.na(x$Cx))),digits=digits)," among ",format(x$M,digits=digits),".","\n",
        "\t| Initial values","\t| x0 = ","(",format(x$x0[1],digits=digits),",",format(x$x0[2],digits=digits),",",format(x$x0[3],digits=digits),")",".","\n",
        "\t| Ending values","\t\t| y  = ","(",format(x$y[1],digits=digits),",",format(x$y[2],digits=digits),",",format(x$y[3],digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}else{
    cat("Stratonovich Bridge Sde 3D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," o dW1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," o dW2(t)","\n",
        "\t| dZ(t) = ", Drz," * dt + ", DDz," o dW3(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N = ",format(x$N+1,digits=digits),".","\n",
        "\t| Crossing realized","\t| C = ",format(length(which(!is.na(x$Cx))),digits=digits)," among ",format(x$M,digits=digits),".","\n",
        "\t| Initial values","\t| x0 = ","(",format(x$x0[1],digits=digits),",",format(x$x0[2],digits=digits),",",format(x$x0[3],digits=digits),")",".","\n",
        "\t| Ending values","\t\t| y  = ","(",format(x$y[1],digits=digits),",",format(x$y[2],digits=digits),",",format(x$y[3],digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}} else {
    if(x$type=="ito"){
    cat(Ito," Bridge Sde 3D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," * dB1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," * dB2(t)","\n",
        "\t| dZ(t) = ", Drz," * dt + ", DDz," * dB3(t)","\n",
		"\t| Correlation structure:",sep="")		 
		prmatrix(x$corrmat,rowlab = rep("            ", 3), collab = rep("", 3),digits=digits)
    cat("Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N = ",format(x$N+1,digits=digits),".","\n",
        "\t| Crossing realized","\t| C = ",format(length(which(!is.na(x$Cx))),digits=digits)," among ",format(x$M,digits=digits),".","\n",
        "\t| Initial values","\t| x0 = ","(",format(x$x0[1],digits=digits),",",format(x$x0[2],digits=digits),",",format(x$x0[3],digits=digits),")",".","\n",
        "\t| Ending values","\t\t| y  = ","(",format(x$y[1],digits=digits),",",format(x$y[2],digits=digits),",",format(x$y[3],digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}else{
    cat("Stratonovich Bridge Sde 3D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," o dB1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," o dB2(t)","\n",
        "\t| dZ(t) = ", Drz," * dt + ", DDz," o dB3(t)","\n",
		"\t| Correlation structure:",sep="")		 
		prmatrix(x$corrmat,rowlab = rep("            ", 3), collab = rep("", 3),digits=digits)
    cat("Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N = ",format(x$N+1,digits=digits),".","\n",
        "\t| Crossing realized","\t| C = ",format(length(which(!is.na(x$Cx))),digits=digits)," among ",format(x$M,digits=digits),".","\n",
        "\t| Initial values","\t| x0 = ","(",format(x$x0[1],digits=digits),",",format(x$x0[2],digits=digits),",",format(x$x0[3],digits=digits),")",".","\n",
        "\t| Ending values","\t\t| y  = ","(",format(x$y[1],digits=digits),",",format(x$y[2],digits=digits),",",format(x$y[3],digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}		
		}
    invisible(x)
}

##
## summary

summary.bridgesde3d <- function(object,at,digits=NULL,...)
                    {
    class(object) <- "bridgesde3d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
	if (is.null(digits)){digits = base::options()$digits}
    if (any(object$T < at || object$t0 > at) )  
	     stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
    if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(C)==1){ Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
    cat("\nMonte-Carlo Statistics for (X(t),Y(t),Z(t)) at time t = ",at,"\n",
	    "\t| Crossing realized ",C," among ",object$M,"\n",
        sep="")
    res <- as.data.frame(matrix(c(mean(x,na.rm = TRUE),var(x,na.rm = TRUE),Median(x),Mode(x),
                               quantile(x,0.25,na.rm = TRUE),quantile(x,0.75,na.rm = TRUE),min(x,na.rm=TRUE),max(x,na.rm=TRUE),
                               skewness(x),kurtosis(x),cv(x),moment(x,order=3,center=TRUE),
                               moment(x,order=4,center=TRUE),moment(x,order=5,center=TRUE),moment(x,order=6,center=TRUE),
                               mean(y,na.rm = TRUE),var(y,na.rm = TRUE),Median(y),Mode(y),
                               quantile(y,0.25,na.rm = TRUE),quantile(y,0.75,na.rm = TRUE),min(y,na.rm=TRUE),max(y,na.rm=TRUE),
                               skewness(y),kurtosis(y),cv(y),moment(y,order=3,center=TRUE),
                               moment(y,order=4,center=TRUE),moment(y,order=5,center=TRUE),moment(y,order=6,center=TRUE),
                               mean(z,na.rm = TRUE),var(z,na.rm = TRUE),Median(z),Mode(z),
                               quantile(z,0.25,na.rm = TRUE),quantile(z,0.75,na.rm = TRUE),min(z,na.rm=TRUE),max(z,na.rm=TRUE),
                               skewness(z),kurtosis(z),cv(z),moment(z,order=3,center=TRUE),
                               moment(z,order=4,center=TRUE),moment(z,order=5,center=TRUE),moment(z,order=6,center=TRUE)),
                               ncol=3))
	dimnames(res) <- list(c("Mean","Variance","Median","Mode","First quartile","Third quartile","Minimum","Maximum","Skewness","Kurtosis",
	                         "Coef-variation","3th-order moment","4th-order moment","5th-order moment","6th-order moment"),c("X","Y","Z"))
    print(round(res,digits=digits), quote = FALSE, right = TRUE,...)
    invisible(object)
}



mean.bridgesde3d <- function(x,at,...)
                    {
	object <- x
    class(object) <- "bridgesde3d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
    if (any(object$T < at || object$t0 > at) )  
	     stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
    if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(C)==1){ Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
	return(c(mean(x,na.rm = TRUE,...),mean(y,na.rm = TRUE,...),mean(z,na.rm = TRUE,...)))
}


cv.bridgesde3d  <- function(x,at,...)
                    {
	object <- x
    class(object) <- "bridgesde3d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
    if (any(object$T < at || object$t0 > at) )  
	     stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
    if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(C)==1){ Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
	return(c(cv(x,...),cv(y,...),cv(z,...)))
}


max.bridgesde3d  <- function(x,at,...)
                    {
	object <- x
    class(object) <- "bridgesde3d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
    if (any(object$T < at || object$t0 > at) )  
	     stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
    if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(C)==1){ Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
	return(c(max(x,na.rm = TRUE),max(y,na.rm = TRUE),max(z,na.rm=TRUE)))
}

min.bridgesde3d  <- function(x,at,...)
                    {
	object <- x
    class(object) <- "bridgesde3d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
    if (any(object$T < at || object$t0 > at) )  
	     stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
    if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(C)==1){ Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
	return(c(min(x,na.rm = TRUE),min(y,na.rm = TRUE),min(z,na.rm=TRUE)))
}

skewness.bridgesde3d <- function(x,at,...)
                    {
	object <- x
    class(object) <- "bridgesde3d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
    if (any(object$T < at || object$t0 > at) ) 
        	stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
    if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(C)==1){ Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
	return(c(skewness(x),skewness(y),skewness(z)))
}

kurtosis.bridgesde3d <- function(x,at,...)
                    {
	object <- x
    class(object) <- "bridgesde3d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
    if (any(object$T < at || object$t0 > at) )  
	     stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
    if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(C)==1){ Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
	return(c(kurtosis(x),kurtosis(y),kurtosis(z)))
}

Median.bridgesde3d <- function(x,at,...)
                    {
	object <- x
    class(object) <- "bridgesde3d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
    if (any(object$T < at || object$t0 > at) )  
	     stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
    if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(C)==1){ Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
	return(c(Median(x),Median(y),Median(z)))
}

Mode.bridgesde3d <- function(x,at,...)
                    {
	object <- x
    class(object) <- "bridgesde3d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
    if (any(object$T < at || object$t0 > at) )  
	     stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
    if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(C)==1){ Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
	return(c(Mode(x),Mode(y),Mode(z)))
}

quantile.bridgesde3d <- function(x,at,...)
                    {
	object <- x
    class(object) <- "bridgesde3d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
    if (any(object$T < at || object$t0 > at) )  
	     stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
    if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(C)==1){ Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
	return(list(x=quantile(x,...),y=quantile(y,...),z=quantile(z,...)))
}

moment.bridgesde3d <- function(x,at,...)
                    {
	object <- x
    class(object) <- "bridgesde3d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
    if (any(object$T < at || object$t0 > at) )  
	     stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
    if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(C)==1){ Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
	return(c(moment(x,...),moment(y,...),moment(z,...)))
}

bconfint.bridgesde3d <- function(x,at,...)
                    {
	object <- x
    class(object) <- "bridgesde3d"
    if (missing(at)) {at = (as.numeric(object$T)-as.numeric(object$t0))/2}
    if (any(object$T < at || object$t0 > at) )  
	     stop( " please use 't0 <= at <= T'")
	C = length(which(!is.na(object$Cx)))
    if (as.numeric(C) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(C)==1){ Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(C)==1){ Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(C)==1){ Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(C),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
	return(list(x=bconfint(x,...),y=bconfint(y,...),z=bconfint(z,...)))
}



time.bridgesde3d <- function(x,...)
                    {
    class(x) <- "bridgesde3d"
    as.vector(time(x$X))
}

##
## Plot

plot.bridgesde3d <- function(x,...) .plot.bridgesde3d(x,...)

lines.bridgesde3d <- function(x,...)
                 {
    class(x) <- "bridgesde3d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    if (length(which(!is.na(x$Cz))) == 1){Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{Z = x$Z}		
    for (i in 1:length(which(!is.na(x$Cx)))){lines(time(x),X[,i],...)}
    for (i in 1:length(which(!is.na(x$Cy)))){lines(time(x),Y[,i],...)} 
    for (i in 1:length(which(!is.na(x$Cz)))){lines(time(x),Z[,i],...)} 	
}

points.bridgesde3d <- function(x,...)
                 {
    class(x) <- "bridgesde3d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    if (length(which(!is.na(x$Cz))) == 1){Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{Z = x$Z}		
    for (i in 1:length(which(!is.na(x$Cx)))){points(time(x),X[,i],...)}
    for (i in 1:length(which(!is.na(x$Cy)))){points(time(x),Y[,i],...)} 
    for (i in 1:length(which(!is.na(x$Cz)))){points(time(x),Z[,i],...)} 
}


plot3D.bridgesde3d <- function(x,display = c("persp","rgl"),...)
                 {
	class(x) <- "bridgesde3d"		 
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    if (length(which(!is.na(x$Cz))) == 1){Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{Z = x$Z}	
    plot3D(X[,1],Y[,1],Z[,1],display,...)
}

