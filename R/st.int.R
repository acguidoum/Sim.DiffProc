## Fri Apr 03 08:54:53 2020
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
##### st.int

st.int <- function(expr, ...)  UseMethod("st.int")

st.int.default <- function(expr, lower = 0, upper = 1, M = 1, subdivisions = 1000L,
                          type=c("ito","str"),...)
          {
    if (any(!is.numeric(subdivisions) || (subdivisions - floor(subdivisions) > 0) || subdivisions <= 1L)) 
	       stop(" subdivisions must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0)) 
	       stop(" 'M' must be a positive integer ")	
    if (!is.expression(expr))  
	       stop(" 'expr' must be expression in t and w 'expr(t,w)'")
    if (missing(type)) 
	       type <- "ito"
    if (any(!is.finite(lower) && !is.finite(upper)))    
	       stop("a limit is missing")
    if (any(is.na(lower) && is.na(upper)))              
	       stop("a limit is missing")
    if (any(lower < 0 || upper < 0 || upper <= lower) ) 
	       stop(" limit of integration. please use positive bound 'upper > lower >= 0' ") 
    t <- seq(lower ,upper, by=(upper-lower)/subdivisions)
    fun <- function(t,w) eval(expr)
    Ito <- function()  {
            w = c(0,cumsum(rnorm(subdivisions+1,mean=0,sd=sqrt((upper-lower)/subdivisions))))
            dw   <- diff(w)
	       St <- cumsum(sapply(1:(subdivisions+1), function(i) fun(t[i],w[i])*dw[i]))
	       St    }
    Str <- function()  {
            w = c(0,cumsum(rnorm(subdivisions+1,mean=0,sd=sqrt((upper-lower)/subdivisions))))
            dw   <- diff(w)
	       St <- cumsum(sapply(1:(subdivisions+1), function(i) 0.5*(fun(t[i],w[i])+fun(t[i+1],w[i+1]))*dw[i]))
	       St    }
    if (type=="ito"){res <- data.frame(sapply(1:M,function(i) Ito()))}
    else { res <- data.frame(sapply(1:M,function(i) Str()))}
    names(res) <- paste("X",1:M,sep="")
    X <- ts(res, start = lower, deltat = (upper-lower)/subdivisions )
    structure(list(X=X, fun=expr[[1]], type=type, subdivisions=subdivisions, M = M, 
                   Dt=(upper-lower)/subdivisions,t0=lower,T=upper,call=match.call()),class="st.int")
}

###

print.st.int <- function(x, digits=NULL, ...)
           {
    class(x) <- "st.int"
	Ito = "It\xf4"
    Encoding(Ito) <- "latin1"
    if(x$type=="ito"){
    cat(Ito," integral:","\n",        
        "\t| X(t)   = integral (f(s,w) * dw(s))","\n", 
        "\t| f(t,w) = ",deparse(x$fun),"\n",
        "Summary:","\n",
        "\t| Number of subinterval","\t| N = ",format(x$subdivisions+1,digits=digits),".","\n",
        "\t| Number of simulation","\t| M = ",format(x$M,digits=digits),".","\n",
        "\t| Limits of integration","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",    
        sep="")}else{
    cat("Stratonovich integral:","\n",
        "\t| X(t)   = integral (f(s,w) o dw(s))","\n", 
        "\t| f(t,w) = ",deparse(x$fun),"\n",
        "Summary:","\n",
        "\t| Number of subinterval","\t| N = ",format(x$subdivisions+1,digits=digits),".","\n",
        "\t| Number of simulation","\t| M = ",format(x$M,digits=digits),".","\n",
        "\t| Limits of integration","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        sep="")}
    invisible(x)
}

##
## Plot

plot.st.int <- function(x,...)
                 {
    class(x) <- "st.int"
    X <- x$X
	if (x$M > 10) {plot(X,plot.type="single",...)}else{plot(X,...)}
}

lines.st.int <- function(x,...)
                 {
    class(x) <- "st.int"
    X <- x$X
    for (i in 1:dim(X)[2]){
    lines(time(x),X[,i],...)}
}

points.st.int <- function(x,...)
                 {
    class(x) <- "st.int"
    X <- x$X
    for (i in 1:dim(X)[2]){
    points(time(x),X[,i],...)}
}


##
## summary

summary.st.int  <- function(object, at,digits=NULL, ...)
           {   
    class(object) <- "st.int"
    if (missing(at)) {at = as.numeric(object$T)}
	if (is.null(digits)){digits = base::options()$digits}
    if (any(object$T < at || object$t0 > at) )  
	         stop( " please use 'lower <= at <= upper'")
    if (object$M == 1){ X = matrix(object$X,nrow=length(object$X),ncol=1)}else{X = object$X}
    xx   <- as.vector(X[which(time(object)==as.character(at)),])
    if (length(xx) == 0){
	 if (object$M==1){ F  <- lapply(1:object$M,function(i) approxfun(time(object),object$X))}else{
                       F  <- lapply(1:object$M,function(i) approxfun(time(object),object$X[,i]))}
                       xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
	if(object$type=="ito"){int=paste("integral(f(s,w) * dw(s))")}else{int=paste("integral (f(s,w) o dw(s))")}
    cat("\nMonte-Carlo Statistics for ", int," at time t = ",at,"\n",
        "| f(t,w) = ",deparse(object$fun),"\n",
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

mean.st.int <- function(x,at,...)
                    {
    class(x) <- "st.int"
    if (missing(at)) {at = as.numeric(x$T)}
    if (any(x$T < at || x$t0 > at) )  
	         stop( " please use 'lower <= at <= upper'")
    if (x$M == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (x$M==1){ F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    mean(xx,na.rm = TRUE,...)
}

cv.st.int <- function(x,at,...)
                    {
    class(x) <- "st.int"
    if (missing(at)) {at = as.numeric(x$T)}
    if (any(x$T < at || x$t0 > at) )  
	       stop( " please use 'lower <= at <= upper'")
    if (x$M == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (x$M==1){ F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    cv(xx,...)
}

max.st.int <- function(x,at,...)
                    {
    class(x) <- "st.int"
    if (missing(at)) {at = as.numeric(x$T)}
    if (any(x$T < at || x$t0 > at) )  
	        stop( " please use 'lower <= at <= upper'")
    if (x$M == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (x$M==1){ F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    max(xx,na.rm = TRUE)
}

min.st.int <- function(x,at,...)
                    {
    class(x) <- "st.int"
    if (missing(at)) {at = as.numeric(x$T)}
    if (any(x$T < at || x$t0 > at) )  
	         stop( " please use 'lower <= at <= upper'")
    if (x$M == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (x$M==1){ F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    min(xx,na.rm = TRUE)
}

skewness.st.int <- function(x,at,...)
                    {
    class(x) <- "st.int"
    if (missing(at)) {at = as.numeric(x$T)}
    if (any(x$T < at || x$t0 > at) )  
	        stop( " please use 'lower <= at <= upper'")
    if (x$M == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (x$M==1){ F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    skewness(xx)
}

kurtosis.st.int <- function(x,at,...)
                    {
    class(x) <- "st.int"
    if (missing(at)) {at = as.numeric(x$T)}
    if (any(x$T < at || x$t0 > at) )  
	         stop( " please use 'lower <= at <= upper'")
    if (x$M == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (x$M==1){ F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    kurtosis(xx)
}

Median.st.int <- function(x,at,...)
                    {
    class(x) <- "st.int"
    if (missing(at)) {at = as.numeric(x$T)}
    if (any(x$T < at || x$t0 > at) )  
	        stop( " please use 'lower <= at <= upper'")
    if (x$M == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (x$M==1){ F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    Median(xx)
}

Mode.st.int <- function(x,at,...)
                    {
    class(x) <- "st.int"
    if (missing(at)) {at = as.numeric(x$T)}
    if (any(x$T < at || x$t0 > at) )  stop( " please use 'lower <= at <= upper'")
    if (x$M == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (x$M==1){ F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    Mode(xx)
}

quantile.st.int <- function(x,at,...)
                    {
    class(x) <- "st.int"
    if (missing(at)) {at = as.numeric(x$T)}
    if (any(x$T < at || x$t0 > at) )  
	         stop( " please use 'lower <= at <= upper'")
    if (x$M == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (x$M==1){ F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    quantile(xx,...)
}

moment.st.int <- function(x,at,...)
                    {
    class(x) <- "st.int"
    if (missing(at)) {at = as.numeric(x$T)}
    if (any(x$T < at || x$t0 > at) )  
	         stop( " please use 'lower <= at <= upper'")
    if (x$M == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (x$M==1){ F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    moment(xx,...)
}

bconfint.st.int <- function(x,at,...)
                    {
    class(x) <- "st.int"
    if (missing(at)) {at = as.numeric(x$T)}
    if (any(x$T < at || x$t0 > at) )  
	         stop( " please use 'lower <= at <= upper'")
    if (x$M == 1){ X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}
    xx   <- as.vector(X[which(time(x)==as.character(at)),])
    if (length(xx) == 0){
	 if (x$M==1){ F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X))}else{
                  F  <- lapply(1:x$M,function(i) approxfun(time(x),x$X[,i]))}
                  xx <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    bconfint(xx,...)
}

time.st.int <- function(x,...)
                    {
    class(x) <- "st.int"
    as.vector(time(x$X))
}
