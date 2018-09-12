## Mon Sep 03 23:51:29 2018
## Original file Copyright © 2018 A.C. Guidoum, K. Boukhetala
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
##### snssde1d

snssde1d <- function(N, ...)  UseMethod("snssde1d")

snssde1d.default <- function(N =1000,M=1,x0=0,t0=0,T=1,Dt=NULL,drift,diffusion,alpha=0.5,mu=0.5,
                     type=c("ito","str"), method=c("euler","milstein","predcorr",
                     "smilstein","taylor","heun","rk1","rk2","rk3"),...)
        {
    if (!is.numeric(x0)) 
	     stop("'x0' must be numeric")
    if (any(!is.numeric(t0) || !is.numeric(T))) 
	     stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) 
	     stop(" 'N' must be an integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0))  
	     stop(" 'M' must be an integer ")
    if (!is.expression(drift) || !is.expression(diffusion)) 
	     stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't' and 'x'")
    if (missing(drift))     
	     drift     <- expression(0)
    if (missing(diffusion)) 
	     diffusion <- expression(1)
    #if (length(all.vars(drift)) > 2 & all.vars(drift) != "t"  & all.vars(drift) != "x" ) stop("coefficient of 'drift' must be expressions in 't' and 'x'")
    #if (length(all.vars(diffusion)) > 2 & all.vars(diffusion) != "t"  & all.vars(diffusion) != "x" ) stop("coefficient of 'diffusion' must be expressions in 't' and 'x'")
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
    if (is.null(Dt)) {
        Dt <- (T - t0)/N
        t <- seq(t0, T, by=Dt)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
		T <- t[N + 1]
    }
    if (method=="euler")         
	      res <- .Euler1D(N,M,x0,t0,T,Dt,drift,diffusion,type)
    else if (method=="predcorr") 
	      res <- .PredCorr1D(N,M,x0,t0,T,Dt,alpha,mu,drift,diffusion,type)
    else if (method=="milstein") 
	      res <- .Milstein1D(N,M,x0,t0,T,Dt,drift,diffusion,type)
    else if (method=="smilstein")
	      res <- .SMilstein1D(N,M,x0,t0,T,Dt,drift,diffusion,type)
    else if (method=="taylor")   
	      res <- .STS1D(N,M,x0,t0,T,Dt,drift,diffusion,type)
    else if (method=="heun")     
	      res <- .Heun1D(N,M,x0,t0,T,Dt,drift,diffusion,type)
    else if (method=="rk1")      
	      res <- .RK1D(N,M,x0,t0,T,Dt,drift,diffusion,type,order=1)
    else if (method=="rk2")      
	      res <- .RK1D(N,M,x0,t0,T,Dt,drift,diffusion,type,order=2)
    else if (method=="rk3")      
	      res <- .RK1D(N,M,x0,t0,T,Dt,drift,diffusion,type,order=3)
		  
    structure(list(X=res$X,drift=drift[[1]], diffusion=diffusion[[1]],type=type,method=method, 
                   x0=as.numeric(format(x0)), N=as.numeric(format(N)), M=as.numeric(format(M)),Dt=as.numeric(format(Dt)),t0=as.numeric(format(t0)),T=as.numeric(format(T)),dim="1d",call=match.call()),class="snssde1d")
}

###

print.snssde1d <- function(x, digits=NULL, ...)
           {
    class(x) <- "snssde1d"
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
    Dr <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)))), list(e = x$drift))))   
    DD <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)))), list(e = x$diffusion))))
    if(x$type=="ito"){
    cat(Ito," Sde 1D:","\n",        
        "\t| dX(t) = ", Dr," * dt + ", DD," * dW(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N+1,digits=digits),".","\n",
        "\t| Number of simulation","\t| M  = ",format(x$M,digits=digits),".","\n",
        "\t| Initial value","\t\t| x0 = ",format(x$x0,digits=digits),".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",       
        sep="")}else{
    cat("Stratonovich Sde 1D:","\n",
        "\t| dX(t) = ", Dr," * dt + ", DD," o dW(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N+1,digits=digits),".","\n",
        "\t| Number of simulation","\t| M  = ",format(x$M,digits=digits),".","\n",
        "\t| Initial value","\t\t| x0 = ",format(x$x0,digits=digits),".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}
    invisible(x)
}

##
## Plot

plot.snssde1d <- function(x,...)
                 {
    class(x) <- "snssde1d"
    X <- x$X
	if (x$M > 10) {plot(X,plot.type="single",...)}else{plot(X,...)}
}

lines.snssde1d <- function(x,...)
                 {
    class(x) <- "snssde1d"
    X <- x$X
	if (as.numeric(x$M) >=2){
    for (i in 1:dim(X)[2]){
    lines(time(x),X[,i],...)}}else{
	lines(time(x),X,...)}
}


points.snssde1d <- function(x,...)
                 {
    class(x) <- "snssde1d"
    X <- x$X
	if (as.numeric(x$M) >=2){
    for (i in 1:dim(X)[2]){
    points(time(x),X[,i],...)}}else{
	points(time(x),X,...)}
}

# add.mean <- function(x,lty=NULL,lwd=NULL,col=NULL,cex=NULL,...)
                 # {
    # class(x) <- "snssde1d"
    # X <- x$X
    # if (is.null(lty)) {lty = 1}
    # if (is.null(lwd)) {lwd = 1}
    # if (is.null(col)) {col = 2}
    # if (is.null(cex)) {cex = 0.8}
	# if (as.numeric(x$M) >=2){
    # lines(time(x),rowMeans(X,na.rm = TRUE),lwd=lwd,lty=lty,col=col,...)}else{
	# lines(time(x),X,lwd=lwd,lty=lty,col=col,...)}
    # legend("topright",c("mean path"),inset = .01,lty=lty,col=col,lwd=lwd,cex=cex,...)
# }


# add.bconfint.snssde1d <- function(x,level=0.95,lty=NULL,lwd=NULL,col=NULL,cex=NULL,...)
                 # {
    # class(x) <- "snssde1d"
    # if (is.null(lty)) {lty = 1}
    # if (is.null(lwd)) {lwd = 1}
    # if (is.null(col)) {col = 4}
    # if (is.null(cex)) {cex = 0.8}
    # lines(time(x),apply(x$X,1,bconfint,level)[1,],lwd=lwd,lty=lty,col=col,...)
    # lines(time(x),apply(x$X,1,bconfint,level)[2,],lwd=lwd,lty=lty,col=col,...)
    # legend("topleft",c(paste("bound of",level*100,"% confidence")),inset = .01,lty=lty,col=col,lwd=lwd,cex=cex,...)
# }

##
## summary

summary.snssde1d  <- function(object, at,digits=NULL, ...)
           {   
    class(object) <- "snssde1d"
    if (missing(at)) {at = as.numeric(object$T)}
	if (is.null(digits)){digits = base::options()$digits}
    if (any(object$T < at || object$t0 > at) )  
	      stop( " please use 't0 <= at <= T'")
    if (object$M == 1){ X = matrix(object$X,nrow=length(object$X),ncol=1)}else{X = object$X}
    xx   <- as.vector(X[which(time(object)==as.character(at)),])
    if (length(xx) == 0){
	 if (object$M==1){ F  <- lapply(1:object$M,function(i) approxfun(time(object),object$X))}else{
                       F  <- lapply(1:object$M,function(i) approxfun(time(object),object$X[,i]))}
                       xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    cat("\n\tMonte-Carlo Statistics for X(t) at time t = ",at,"\n",
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

mean.snssde1d <- function(x,at,...)
                    {
    class(x) <- "snssde1d"
    if (missing(at)) {at = as.numeric(x$T)}
    if (any(x$T < at || x$t0 > at) )  
	       stop( " please use 't0 <= at <= T'")
    if (x$M == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (x$M==1){ F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    mean(xx,na.rm = TRUE,...)
}

cv.snssde1d <- function(x,at,...)
                    {
    class(x) <- "snssde1d"
    if (missing(at)) {at = as.numeric(x$T)}
    if (any(x$T < at || x$t0 > at) )  
	      stop( " please use 't0 <= at <= T'")
    if (x$M == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (x$M==1){ F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    cv(xx,...)
}

max.snssde1d <- function(x,at,...)
                    {
    class(x) <- "snssde1d"
    if (missing(at)) {at = as.numeric(x$T)}
    if (any(x$T < at || x$t0 > at) )  
	        stop( " please use 't0 <= at <= T'")
    if (x$M == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (x$M==1){ F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    max(xx,na.rm = TRUE)
}

min.snssde1d <- function(x,at,...)
                    {
    class(x) <- "snssde1d"
    if (missing(at)) {at = as.numeric(x$T)}
    if (any(x$T < at || x$t0 > at) )  
	       stop( " please use 't0 <= at <= T'")
    if (x$M == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (x$M==1){ F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    min(xx,na.rm = TRUE)
}

skewness.snssde1d <- function(x,at,...)
                    {
    class(x) <- "snssde1d"
    if (missing(at)) {at = as.numeric(x$T)}
    if (any(x$T < at || x$t0 > at) )  
	        stop( " please use 't0 <= at <= T'")
    if (x$M == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (x$M==1){ F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    skewness(xx)
}

kurtosis.snssde1d <- function(x,at,...)
                    {
    class(x) <- "snssde1d"
    if (missing(at)) {at = as.numeric(x$T)}
    if (any(x$T < at || x$t0 > at) )  
	         stop( " please use 't0 <= at <= T'")
    if (x$M == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (x$M==1){ F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    kurtosis(xx)
}

Median.snssde1d <- function(x,at,...)
                    {
    class(x) <- "snssde1d"
    if (missing(at)) {at = as.numeric(x$T)}
    if (any(x$T < at || x$t0 > at) )  
	     stop( " please use 't0 <= at <= T'")
    if (x$M == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (x$M==1){ F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    Median(xx)
}

Mode.snssde1d <- function(x,at,...)
                    {
    class(x) <- "snssde1d"
    if (missing(at)) {at = as.numeric(x$T)}
    if (any(x$T < at || x$t0 > at) )  
	       stop( " please use 't0 <= at <= T'")
    if (x$M == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (x$M==1){ F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    Mode(xx)
}

quantile.snssde1d <- function(x,at,...)
                    {
    class(x) <- "snssde1d"
    if (missing(at)) {at = as.numeric(x$T)}
    if (any(x$T < at || x$t0 > at) )  
	         stop( " please use 't0 <= at <= T'")
    if (x$M == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (x$M==1){ F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    quantile(xx,...)
}

moment.snssde1d <- function(x,at,...)
                    {
    class(x) <- "snssde1d"
    if (missing(at)) {at = as.numeric(x$T)}
    if (any(x$T < at || x$t0 > at) )  
	         stop( " please use 't0 <= at <= T'")
    if (x$M == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (x$M==1){ F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    moment(xx,...)
}

bconfint.snssde1d <- function(x,at,...)
                    {
    class(x) <- "snssde1d"
    if (missing(at)) {at = as.numeric(x$T)}
    if (any(x$T < at || x$t0 > at) )  
	             stop( " please use 't0 <= at <= T'")
    if (x$M == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (x$M==1){ F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    bconfint(xx,...)
}


time.snssde1d <- function(x,...)
                    {
    class(x) <- "snssde1d"
    as.vector(time(x$X))
}


################################################################################
################################################################################
##### snssde2d

snssde2d <- function(N, ...)  UseMethod("snssde2d")

snssde2d.default <- function(N =1000,M=1,x0=c(0,0),t0=0,T=1,Dt=NULL,drift,diffusion,alpha=0.5,mu=0.5,
                     type=c("ito","str"), method=c("euler","milstein","predcorr",
                     "smilstein","taylor","heun","rk1","rk2","rk3"),...)
        {
    if (any(!is.numeric(x0) || length(x0) !=2)) 
	           stop("'x0' must be numeric, and length(x0) = 2")
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
	if (missing(drift)) 
	          drift <- expression(0,0)
    if (missing(diffusion))    
	          diffusion  <- expression(1,1)
    # if (length(all.vars(drift)) == 3 && all.vars(drift) != "t"  && all.vars(drift) != "x" && all.vars(drift) !="y") stop("coefficient of 'driftx' must be expressions in 't', 'x' and 'y'")
    # if (length(all.vars(drifty)) == 3 && all.vars(drifty) != "t"  && all.vars(drifty) != "x" && all.vars(drifty) !="y") stop("coefficient of 'drifty' must be expressions in 't', 'x' and 'y'")
    # if (length(all.vars(diffx)) == 3 && all.vars(diffx) != "t"  && all.vars(diffx) != "x" && all.vars(diffx) != "x" ) stop("coefficient of 'diffx' must be expressions in 't', 'x' and 'y'")
    # if (length(all.vars(diffy)) == 3 && all.vars(diffy) != "t"  && all.vars(diffy) != "x" && all.vars(diffy) != "x" ) stop("coefficient of 'diffy' must be expressions in 't', 'x' and 'y'")
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
    if (is.null(Dt)) {
        Dt <- (T - t0)/N
        t <- seq(t0, T, by=Dt)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
		T <- t[N + 1]
    }   
    if (method=="euler")         
	      res <- .Euler2D(N,M,x0=x0[1],y0=x0[2],t0,T,Dt,driftx=drift[1],diffx=diffusion[1],drifty=drift[2],diffy=diffusion[2],type)
    else if (method=="predcorr") 
	      res <- .PredCorr2D(N,M,x0=x0[1],y0=x0[2],t0,T,Dt,alpha,mu,driftx=drift[1],diffx=diffusion[1],drifty=drift[2],diffy=diffusion[2],type)
    else if (method=="milstein") 
	      res <- .Milstein2D(N,M,x0=x0[1],y0=x0[2],t0,T,Dt,driftx=drift[1],diffx=diffusion[1],drifty=drift[2],diffy=diffusion[2],type)
    else if (method=="smilstein")
	      res <- .SMilstein2D(N,M,x0=x0[1],y0=x0[2],t0,T,Dt,driftx=drift[1],diffx=diffusion[1],drifty=drift[2],diffy=diffusion[2],type)
    else if (method=="taylor")   
	      res <- .STS2D(N,M,x0=x0[1],y0=x0[2],t0,T,Dt,driftx=drift[1],diffx=diffusion[1],drifty=drift[2],diffy=diffusion[2],type)
    else if (method=="heun")     
	      res <- .Heun2D(N,M,x0=x0[1],y0=x0[2],t0,T,Dt,driftx=drift[1],diffx=diffusion[1],drifty=drift[2],diffy=diffusion[2],type)
    else if (method=="rk1")      
	      res <- .RK2D(N,M,x0=x0[1],y0=x0[2],t0,T,Dt,driftx=drift[1],diffx=diffusion[1],drifty=drift[2],diffy=diffusion[2],type,order=1)
    else if (method=="rk2")      
	      res <- .RK2D(N,M,x0=x0[1],y0=x0[2],t0,T,Dt,driftx=drift[1],diffx=diffusion[1],drifty=drift[2],diffy=diffusion[2],type,order=2)
    else if (method=="rk3")      
	      res <- .RK2D(N,M,x0=x0[1],y0=x0[2],t0,T,Dt,driftx=drift[1],diffx=diffusion[1],drifty=drift[2],diffy=diffusion[2],type,order=3)
		  
    structure(list(X=res$X,Y=res$Y, driftx=drift[[1]], diffx=diffusion[[1]],drifty=drift[[2]], diffy=diffusion[[2]],type=type,method=method, 
                   x0=as.numeric(format(x0[1])),y0=as.numeric(format(x0[2])), N=as.numeric(format(N)),M=as.numeric(format(M)),Dt=as.numeric(format(Dt)),t0=as.numeric(format(t0)),T=as.numeric(format(T)),dim="2d",call=match.call()),class="snssde2d")
}


###

print.snssde2d <- function(x, digits=NULL, ...)
           {
    class(x) <- "snssde2d"
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
    if(x$type=="ito"){
    cat(Ito," Sde 2D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," * dW1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," * dW2(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N+1,digits=digits),".","\n",
        "\t| Number of simulation","\t| M  = ",format(x$M,digits=digits),".","\n",
        "\t| Initial values","\t| (x0,y0) = ","(",format(x$x0,digits=digits),",",format(x$y0,digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}else{
    cat("Stratonovich Sde 2D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," o dW1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," o dW2(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N+1,digits=digits),".","\n",
        "\t| Number of simulation","\t| M  = ",format(x$M,digits=digits),".","\n",
        "\t| Initial values","\t| (x0,y0) = ","(",format(x$x0,digits=digits),",",format(x$y0,digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}
    invisible(x)
}

##
## Plot

plot.snssde2d <- function(x,...) .plot.snssde2d(x,...)

lines.snssde2d <- function(x,...)
                 {
    class(x) <- "snssde2d"
    if (as.numeric(x$M) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{
                  X = x$X
                  Y = x$Y}    
    for (i in 1:as.numeric(x$M)){
    lines(time(x),X[,i],...)
    lines(time(x),Y[,i],...)}
}

points.snssde2d <- function(x,...)
                 {
    class(x) <- "snssde2d"
    if (as.numeric(x$M) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{
                  X = x$X
                  Y = x$Y} 
    for (i in 1:as.numeric(x$M)){
    points(time(x),X[,i],...)
    points(time(x),Y[,i],...)}
}

plot2d.snssde2d <- function(x,...) .plot2d.snssde2d(x,...)

lines2d.snssde2d <- function(x,...)
        {
    class(x) <- "snssde2d"
    if (as.numeric(x$M) == 1){
    lines(as.vector(x$X),as.vector(x$Y),...)}else{
    lines(as.vector(x$X[,1]),as.vector(x$Y[,1]),...)}
}

points2d.snssde2d <- function(x,...)
        {
    class(x) <- "snssde2d"
    if (as.numeric(x$M) == 1){
    points(as.vector(x$X),as.vector(x$Y),...)}else{
    points(as.vector(x$X[,1]),as.vector(x$Y[,1]),...)}
}

##
## summary

summary.snssde2d  <- function(object,at,digits=NULL, ...)
           {   
    class(object) <- "snssde2d"
    if (missing(at)) {at = as.numeric(object$T)}
	if (is.null(digits)){digits = base::options()$digits}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	             stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    cat("\n\tMonte-Carlo Statistics for (X(t),Y(t)) at time t = ",at,"\n",
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

mean.snssde2d <- function(x,at,...)
                    {
	object <- x				
    class(object) <- "snssde2d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	            stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
	return(c(mean(x,na.rm = TRUE,...),mean(y,na.rm = TRUE,...)))
}

cv.snssde2d <- function(x,at,...)
                    {
	object <- x				
    class(object) <- "snssde2d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	          stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
	return(c(cv(x,...),cv(y,...)))
}


max.snssde2d <- function(x,at,...)
                    {
	object <- x				
    class(object) <- "snssde2d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	             stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
	return(c(max(x,na.rm = TRUE),max(y,na.rm = TRUE)))
}

min.snssde2d <- function(x,at,...)
                    {
	object <- x				
    class(object) <- "snssde2d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	           stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
	return(c(min(x,na.rm = TRUE),min(y,na.rm = TRUE)))
}

skewness.snssde2d <- function(x,at,...)
                    {
	object <- x				
    class(object) <- "snssde2d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	             stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
	return(c(skewness(x),skewness(y)))
}

kurtosis.snssde2d <- function(x,at,...)
                    {
	object <- x				
    class(object) <- "snssde2d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	           stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
	return(c(kurtosis(x),kurtosis(y)))
}

Median.snssde2d <- function(x,at,...)
                    {
	object <- x				
    class(object) <- "snssde2d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	          stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
	return(c(Median(x),Median(y)))
}

Mode.snssde2d <- function(x,at,...)
                    {
	object <- x				
    class(object) <- "snssde2d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	        stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
	return(c(Mode(x),Mode(y)))
}


quantile.snssde2d <- function(x,at,...)
                    {
	object <- x				
    class(object) <- "snssde2d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) ) 
           	stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
	return(list(x=quantile(x,...),y=quantile(y,...)))
}

moment.snssde2d <- function(x,at,...)
                    {
	object <- x				
    class(object) <- "snssde2d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	        stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
	return(c(moment(x,...),moment(y,...)))
}

bconfint.snssde2d <- function(x,at,...)
                    {
	object <- x				
    class(object) <- "snssde2d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) ) 
           	stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
	return(list(x=bconfint(x,...),y=bconfint(y,...)))
}

time.snssde2d <- function(x,...)
                    {
    class(x) <- "snssde2d"
    as.vector(time(x$X))
}


################################################################################
################################################################################
##### snssde3d


snssde3d <- function(N, ...)  UseMethod("snssde3d")

snssde3d.default <- function(N =1000,M=1,x0=c(0,0,0),t0=0,T=1,Dt=NULL,drift,diffusion,
                             alpha=0.5,mu=0.5,type=c("ito","str"), method=c("euler","milstein",
                             "predcorr","smilstein","taylor","heun","rk1","rk2","rk3"),...)
        {
    if (any(!is.numeric(x0) || length(x0) !=3)) 
	         stop("'x0' must be numeric, and length(x0) = 3")
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
	if (missing(drift)) 
	         drift <- expression(0,0,0)
    if (missing(diffusion))    
	         diffusion  <- expression(1,1,1)
    if (missing(type)) 
	         type <- "ito"
    method <- match.arg(method)
    if (method =="predcorr"){
    if (any(alpha > 1 || alpha < 0)) 
	       stop("please use '0 <= alpha <= 1' ")
    if (any(mu > 1 || mu < 0))       
	       stop("please use '0 <= mu <= 1' ")
                            }
    if (t0 < 0 || T < 0) 
	        stop(" please use positive times! (0 <= t0 < T) ")
    if (is.null(Dt)) {
        Dt <- (T - t0)/N
        t <- seq(t0, T, by=Dt)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
		T <- t[N + 1]
    }   
    if (method=="euler")         
	     res <- .Euler3D(N,M,x0=x0[1],y0=x0[2],z0=x0[3],t0,T,Dt,driftx=drift[1],diffx=diffusion[1],drifty=drift[2],diffy=diffusion[2],driftz=drift[3],diffz=diffusion[3],type)
    else if (method=="predcorr") 
	     res <- .PredCorr3D(N,M,x0=x0[1],y0=x0[2],z0=x0[3],t0,T,Dt,alpha,mu,driftx=drift[1],diffx=diffusion[1],drifty=drift[2],diffy=diffusion[2],driftz=drift[3],diffz=diffusion[3],type)
    else if (method=="milstein") 
	     res <- .Milstein3D(N,M,x0=x0[1],y0=x0[2],z0=x0[3],t0,T,Dt,driftx=drift[1],diffx=diffusion[1],drifty=drift[2],diffy=diffusion[2],driftz=drift[3],diffz=diffusion[3],type)
    else if (method=="smilstein")
	     res <- .SMilstein3D(N,M,x0=x0[1],y0=x0[2],z0=x0[3],t0,T,Dt,driftx=drift[1],diffx=diffusion[1],drifty=drift[2],diffy=diffusion[2],driftz=drift[3],diffz=diffusion[3],type)
    else if (method=="taylor")   
	     res <- .STS3D(N,M,x0=x0[1],y0=x0[2],z0=x0[3],t0,T,Dt,driftx=drift[1],diffx=diffusion[1],drifty=drift[2],diffy=diffusion[2],driftz=drift[3],diffz=diffusion[3],type)
    else if (method=="heun")     
	     res <- .Heun3D(N,M,x0=x0[1],y0=x0[2],z0=x0[3],t0,T,Dt,driftx=drift[1],diffx=diffusion[1],drifty=drift[2],diffy=diffusion[2],driftz=drift[3],diffz=diffusion[3],type)
    else if (method=="rk1")      
	     res <- .RK3D(N,M,x0=x0[1],y0=x0[2],z0=x0[3],t0,T,Dt,driftx=drift[1],diffx=diffusion[1],drifty=drift[2],diffy=diffusion[2],driftz=drift[3],diffz=diffusion[3],type,order=1)
    else if (method=="rk2")      
	     res <- .RK3D(N,M,x0=x0[1],y0=x0[2],z0=x0[3],t0,T,Dt,driftx=drift[1],diffx=diffusion[1],drifty=drift[2],diffy=diffusion[2],driftz=drift[3],diffz=diffusion[3],type,order=2)
    else if (method=="rk3")      
	     res <- .RK3D(N,M,x0=x0[1],y0=x0[2],z0=x0[3],t0,T,Dt,driftx=drift[1],diffx=diffusion[1],drifty=drift[2],diffy=diffusion[2],driftz=drift[3],diffz=diffusion[3],type,order=3)
		 
    structure(list(X=res$X,Y=res$Y,Z=res$Z,driftx=drift[[1]], diffx=diffusion[[1]],drifty=drift[[2]], diffy=diffusion[[2]],driftz=drift[[3]], 
                   diffz=diffusion[[3]],type=type,method=method,x0=as.numeric(format(x0[1])),y0=as.numeric(format(x0[2])),z0=as.numeric(format(x0[3])),N=as.numeric(format(N)),M=as.numeric(format(M)),Dt=as.numeric(format(Dt)),t0=as.numeric(format(t0)),T=as.numeric(format(T)),dim="3d",call=match.call()),class="snssde3d")
}


###

print.snssde3d <- function(x, digits=NULL, ...)
           {
    class(x) <- "snssde3d"
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
    if(x$type=="ito"){
    cat(Ito," Sde 3D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," * dW1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," * dW2(t)","\n",
        "\t| dZ(t) = ", Drz," * dt + ", DDz," * dW3(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N+1,digits=digits),".","\n",
        "\t| Number of simulation","\t| M  = ",format(x$M,digits=digits),".","\n",
        "\t| Initial values","\t| (x0,y0,z0) = ","(",format(x$x0,digits=digits),",",format(x$y0,digits=digits),",",format(x$z0,digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}else{
    cat("Stratonovich Sde 3D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," o dW1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," o dW2(t)","\n",
        "\t| dZ(t) = ", Drz," * dt + ", DDz," o dW3(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N+1,digits=digits),".","\n",
        "\t| Number of simulation","\t| M  = ",format(x$M,digits=digits),".","\n",
        "\t| Initial values","\t| (x0,y0,z0) = ","(",format(x$x0,digits=digits),",",format(x$y0,digits=digits),",",format(x$z0,digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}
    invisible(x)
}

##
## Plot

plot.snssde3d <- function(x,...) .plot.snssde3d(x,...)

lines.snssde3d <- function(x,...)
                 {
    class(x) <- "snssde3d"
    if (as.numeric(x$M) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)
                  Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{
                  X = x$X
                  Y = x$Y
                  Z = x$Z}    
    for (i in 1:as.numeric(x$M)){
    lines(time(x),X[,i],...)
    lines(time(x),Y[,i],...)
    lines(time(x),Z[,i],...)}
}

points.snssde3d <- function(x,...)
                 {
    class(x) <- "snssde3d"
    if (as.numeric(x$M) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)
                  Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{
                  X = x$X
                  Y = x$Y
                  Z = x$Z} 
    for (i in 1:as.numeric(x$M)){
    points(time(x),X[,i],...)
    points(time(x),Y[,i],...)
    points(time(x),Z[,i],...)}
}


plot3D.snssde3d <- function(x,display = c("persp","rgl"),...)
                 {
	class(x) <- "snssde3d"		 
    if (as.numeric(x$M) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)
                  Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{
                  X = x$X
                  Y = x$Y
                  Z = x$Z} 
    plot3D(X[,1],Y[,1],Z[,1],display,...)
}

##
## summary


summary.snssde3d  <- function(object,at,digits=NULL, ...)
           {   
    class(object) <- "snssde3d"
    if (missing(at)) {at = as.numeric(object$T)}
	if (is.null(digits)){digits = base::options()$digits}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	        stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(object$M)==1){ Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
    cat("\n  Monte-Carlo Statistics for (X(t),Y(t),Z(t)) at time t = ",at,"\n",
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

mean.snssde3d <- function(x,at,...)
                    {
	object <- x				
    class(object) <- "snssde3d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	         stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(object$M)==1){ Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
	return(c(mean(x,na.rm = TRUE,...),mean(y,na.rm = TRUE,...),mean(z,na.rm = TRUE,...)))
}

cv.snssde3d <- function(x,at,...)
                    {
	object <- x				
    class(object) <- "snssde3d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	         stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(object$M)==1){ Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
	return(c(cv(x,...),cv(y,...),cv(z,...)))
}


skewness.snssde3d <- function(x,at,...)
                    {
	object <- x				
    class(object) <- "snssde3d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	          stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(object$M)==1){ Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
	return(c(skewness(x),skewness(y),skewness(z)))
}

kurtosis.snssde3d <- function(x,at,...)
                    {
	object <- x				
    class(object) <- "snssde3d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	           stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(object$M)==1){ Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
	return(c(kurtosis(x),kurtosis(y),kurtosis(z)))
}

Median.snssde3d <- function(x,at,...)
                    {
	object <- x				
    class(object) <- "snssde3d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	          stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(object$M)==1){ Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
	return(c(Median(x),Median(y),Median(z)))
}

Mode.snssde3d <- function(x,at,...)
                    {
	object <- x				
    class(object) <- "snssde3d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) ) 
            	stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(object$M)==1){ Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
	return(c(Mode(x),Mode(y),Mode(z)))
}

quantile.snssde3d <- function(x,at,...)
                    {
	object <- x				
    class(object) <- "snssde3d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	             stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(object$M)==1){ Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
	return(list(x=quantile(x,...),y=quantile(y,...),z=quantile(z,...)))
}

moment.snssde3d <- function(x,at,...)
                    {
	object <- x				
    class(object) <- "snssde3d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	          stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(object$M)==1){ Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
	return(c(moment(x,...),moment(y,...),moment(z,...)))
}

min.snssde3d <- function(x,at,...)
                    {
	object <- x				
    class(object) <- "snssde3d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	            stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(object$M)==1){ Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
	return(c(min(x,na.rm = TRUE),min(y,na.rm = TRUE),min(z,na.rm = TRUE)))
}

max.snssde3d <- function(x,at,...)
                    {
	object <- x				
    class(object) <- "snssde3d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	             stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(object$M)==1){ Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
	return(c(max(x,na.rm = TRUE),max(y,na.rm = TRUE),max(z,na.rm = TRUE)))
}

bconfint.snssde3d <- function(x,at,...)
                    {
	object <- x				
    class(object) <- "snssde3d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	              stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(object$M)==1){ Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
	return(list(x=bconfint(x,...),y=bconfint(y,...),z=bconfint(z,...)))
}

time.snssde3d <- function(x,...)
                    {
    class(x) <- "snssde3d"
    as.vector(time(x$X))
}
