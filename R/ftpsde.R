## Mon May 27 03:39:04 2019
## Original file Copyright © 2019 A.C. Guidoum, K. Boukhetala
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
##### fptsde1d

fptsde1d <- function(object, ...)  UseMethod("fptsde1d")

fptsde1d.default <- function(object,boundary,...)
                     {
    class(object) <- "snssde1d"
    if (any(!is.expression(boundary))) 
	          stop(" must be expressions of a constant or time-dependent boundary ")
    R   <- data.frame(object$X)
    Bn  <- function(t)  eval(boundary)
    if(object$x0==Bn(object$t0)) 
	          warning(paste("x0 = S(t0) ==> crossing realized at time",object$t0))
    F   <- lapply(1:as.numeric(object$M),function(i) stats::approxfun(time(object),as.vector(R[,i])) )
      if (as.numeric(object$x0) > Bn(as.numeric(object$t0))){
           v1  <- sapply(1:length(F),function(i) ifelse(length(which(R[,i] <= Bn(time(object))))==0, NA,min(which(R[,i] <= Bn(time(object))))))
           fpt <- sapply(1:length(F),function(i) ifelse(is.na(v1[i]),NA,stats::uniroot(f = function(x) (F[[i]](x)- Bn(x)),lower=time(object)[v1[i]-1],upper=time(object)[v1[i]],tol= .Machine$double.eps)$root))
      }else if (as.numeric(object$x0) < Bn(as.numeric(object$t0))){
           v2  <- sapply(1:length(F),function(i) ifelse(length(which(R[,i] <= Bn(time(object))))==as.numeric(object$N)+1, NA,min(which(R[,i] >= Bn(time(object))))))
           fpt <- sapply(1:length(F),function(i) ifelse(is.na(v2[i]),NA,stats::uniroot(f = function(x) (F[[i]](x)- Bn(x)),lower=time(object)[v2[i]-1],upper=time(object)[v2[i]],tol= .Machine$double.eps)$root))
      }else{
           fpt <- rep(0,as.numeric(object$M))
                     }
	fpt <- fpt[which(!is.na(fpt))]
    #structure(list(SDE=object,boundary=boundary[[1]],fpt=fpt),class="rfptsde1d")
	#return(fpt)
	#if (length(fpt != object$M)) {cat("missing output are removed\n")}
    structure(list(fpt=fpt,boundary=boundary[[1]],obj=object,call=match.call()),class="fptsde1d")
}

print.fptsde1d <- function(x, digits=NULL, ...)
           {
    class(x) <- "fptsde1d"
	Ito = "It\xf4"
    Encoding(Ito) <- "latin1"	
    S <- function(t) eval(x$boundary)
    # Dr <- gsub(pattern = 'x', replacement = 'X(t)', x = as.expression(x$obj$drift), ignore.case = F,fixed = T)
    # DD <- gsub(pattern = 'x', replacement = 'X(t)', x = as.expression(x$obj$diffusion), ignore.case = F,fixed = T)
    Dr <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)))), list(e = x$obj$drift))))   
    DD <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)))), list(e = x$obj$diffusion))))
    if (as.numeric(x$obj$x0) > S(as.numeric(x$obj$t0))){
    fpt_x <- paste("T(S(t),X(t)) = inf{t >= ",as.numeric(x$obj$t0) ,": X(t) <= ", deparse(x$boundary),"}")
    }else{
    fpt_x <- paste("T(S(t),X(t)) = inf{t >= ",as.numeric(x$obj$t0) ,": X(t) >= ", deparse(x$boundary),"}")
    }
    if(x$obj$type=="ito"){
    cat(Ito," Sde 1D:","\n",        
        "\t| dX(t) = ", Dr," * dt + ", DD," * dW(t)","\n",
		"\t| t in [",format(x$obj$t0,digits=digits),",",format(x$obj$T,digits=digits),"].","\n",
        "Boundary:","\n",
        "\t| S(t) = ",deparse(x$boundary),"\n",
        "F.P.T:","\n",
        "\t| ",fpt_x,"\n",
        "\t| Crossing realized ",format(length(x$fpt),digits=digits)," among ",format(x$obj$M,digits=digits),".","\n",      
        sep="")}else{
    cat("Stratonovich Sde 1D:","\n",
        "\t| dX(t) = ", Dr," * dt + ", DD," o dW(t)","\n",
		"\t| t in [",format(x$obj$t0,digits=digits),",",format(x$obj$T,digits=digits),"].","\n",
        "Boundary:","\n",
        "\t| S(t) = ",deparse(x$boundary),"\n",
        "F.P.T:","\n",
        "\t| ",fpt_x,"\n",
        "\t| Crossing realized ",format(length(x$fpt),digits=digits)," among ",format(x$obj$M,digits=digits),".","\n",
        sep="")}
    invisible(x)
}

summary.fptsde1d  <- function(object,digits=NULL, ...)
           {   
    class(object) <- "fptsde1d"
	if (is.null(digits)){digits = base::options()$digits}
    S <- function(t) eval(object$boundary)
    if (as.numeric(object$obj$x0) > S(as.numeric(object$obj$t0))){
    fpt_x <- paste("T(S(t),X(t)) = inf{t >= ",as.numeric(object$obj$t0) ,": X(t) <= ", deparse(object$boundary),"}")
    }else{
    fpt_x <- paste("T(S(t),X(t)) = inf{t >= ",as.numeric(object$obj$t0) ,": X(t) >= ", deparse(object$boundary),"}")
    }
    xx <- object$fpt
    cat("\nMonte-Carlo Statistics of F.P.T:","\n",
        "|",fpt_x,"\n",
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

mean.fptsde1d <- function(x,...)
                    {
    class(x) <- "fptsde1d"
    mean(x$fpt,na.rm = TRUE,...)
}

cv.fptsde1d <- function(x,...)
                    {
    class(x) <- "fptsde1d"
    cv(x$fpt,...)
}

min.fptsde1d <- function(x,...)
                    {
    class(x) <- "fptsde1d"
    min(x$fpt,na.rm = TRUE)
}

max.fptsde1d <- function(x,...)
                    {
    class(x) <- "fptsde1d"
    max(x$fpt,na.rm = TRUE)
}

skewness.fptsde1d <- function(x,...)
                    {
    class(x) <- "fptsde1d"
	skewness(x$fpt)
}

kurtosis.fptsde1d <- function(x,...)
                    {
    class(x) <- "fptsde1d"
	kurtosis(x$fpt)
}

Median.fptsde1d <- function(x,...)
                    {
    class(x) <- "fptsde1d"
	Median(x$fpt)
}

Mode.fptsde1d <- function(x,...)
                    {
    class(x) <- "fptsde1d"
	Mode(x$fpt)
}

quantile.fptsde1d <- function(x,...)
                    {
    class(x) <- "fptsde1d"
    quantile(x$fpt,...)
}

moment.fptsde1d <- function(x,...)
                    {
    class(x) <- "fptsde1d"
    moment(x$fpt,...)
}

## density fpt

dfptsde1d <- function(object, ...)  UseMethod("dfptsde1d")

dfptsde1d.default <- function(object,...)
                     {
    class(object) <- "fptsde1d"
    # if (any(!is.expression(boundary))) stop(" must be expressions of a constant or time-dependent boundary ")
    # R   <- data.frame(object$X)
    # Bn  <- function(t)  eval(boundary)
    # F   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),as.vector(R[,i])) )
      # if (as.numeric(object$x0) > Bn(as.numeric(object$t0))){
           # v1  <- sapply(1:length(F),function(i) ifelse(length(which(R[,i] <= Bn(time(object))))==0, NA,min(which(R[,i] <= Bn(time(object))))))
           # fpt <- sapply(1:length(F),function(i) ifelse(is.na(v1[i]),NA,uniroot(f = function(x) (F[[i]](x)- Bn(x)),lower=time(object)[v1[i]-1],upper=time(object)[v1[i]],tol= .Machine$double.eps)$root))
      # }else if (as.numeric(object$x0) < Bn(as.numeric(object$t0))){
           # v2  <- sapply(1:length(F),function(i) ifelse(length(which(R[,i] <= Bn(time(object))))==as.numeric(object$N)+1, NA,min(which(R[,i] >= Bn(time(object))))))
           # fpt <- sapply(1:length(F),function(i) ifelse(is.na(v2[i]),NA,uniroot(f = function(x) (F[[i]](x)- Bn(x)),lower=time(object)[v2[i]-1],upper=time(object)[v2[i]],tol= .Machine$double.eps)$root))
      # }else{
           # fpt <- rep(0,as.numeric(object$M))
                     # }
      # structure(list(SDE=object,ech=fpt,res=density.default(fpt,na.rm = TRUE,...),boundary=boundary[[1]]),class="dfptsde1d")
	  structure(list(SDE=object$obj,ech=object$fpt,res=density.default(object$fpt,na.rm = TRUE,...),boundary=object$boundary),class="dfptsde1d")
}

print.dfptsde1d <- function(x, digits=NULL, ...)
           {
    class(x) <- "dfptsde1d"
	if (is.null(digits)){digits = base::options()$digits}
    S <- function(t) eval(x$boundary)
    if (as.numeric(x$SDE$x0) > S(as.numeric(x$SDE$t0))){
    cat("\nKernel density of F.P.T:","\n", 
	"|T(S,X) = inf{t >= ",x$SDE$t0 ," : X(t) <= ", deparse(x$boundary),"}","\n",
        sep="")
    }else{
    cat("\nKernel density of F.P.T:","\n", 
	"|T(S,X) = inf{t >= ",x$SDE$t0 ," : X(t) >= ", deparse(x$boundary),"}","\n",
        sep="")}
    cat(
	"\nData: ", x$res$data.name, " (", x$res$n, " obs.);",
	"\tBandwidth 'bw' = ", formatC(x$res$bw, digits = digits), "\n\n", sep = "")
    out <- as.data.frame(x$res[c("x","y")])
    names(out) <- c("x","f(x)")
    print(summary(out, digits = digits), ...)
    invisible(x)
}

plot.dfptsde1d <- function(x,hist=FALSE,...) .plot.dfptsde1d(x,hist,...)


################################################################################
################################################################################
#####
##### fptsde2d

fptsde2d <- function(object, ...)  UseMethod("fptsde2d")


fptsde2d.default <- function(object,boundary,...)
                     {
    class(object) <- "snssde2d"
    if (any(!is.expression(boundary))) 
	         stop(" must be expression of a constant or time-dependent boundary ")
    R1   <- data.frame(object$X)
    R2   <- data.frame(object$Y)
    Bn  <- function(t)  eval(boundary) 
    if(object$x0==Bn(object$t0) | object$y0==Bn(object$t0) ) 
	         warning(paste("x0 = S(t0) or y0 =S(t0) ==> crossing realized at time",object$t0))
    F1   <- lapply(1:as.numeric(object$M),function(i) stats::approxfun(time(object),as.vector(R1[,i])) )
    F2   <- lapply(1:as.numeric(object$M),function(i) stats::approxfun(time(object),as.vector(R2[,i])) )
    if (as.numeric(object$x0) > Bn(as.numeric(object$t0))){
           v11  <- sapply(1:length(F1),function(i) ifelse(length(which(R1[,i] <= Bn(time(object))))==0, NA,min(which(R1[,i] <= Bn(time(object))))))
           fptx <- sapply(1:length(F1),function(i) ifelse(is.na(v11[i]),NA,stats::uniroot(f = function(x) (F1[[i]](x)- Bn(x)),lower=time(object)[v11[i]-1],upper=time(object)[v11[i]],tol= .Machine$double.eps)$root))
      }else if (as.numeric(object$x0) < Bn(as.numeric(object$t0))){
           v12  <- sapply(1:length(F1),function(i) ifelse(length(which(R1[,i] <= Bn(time(object))))==as.numeric(object$N)+1, NA,min(which(R1[,i] >= Bn(time(object))))))
           fptx <- sapply(1:length(F1),function(i) ifelse(is.na(v12[i]),NA,stats::uniroot(f = function(x) (F1[[i]](x)- Bn(x)),lower=time(object)[v12[i]-1],upper=time(object)[v12[i]],tol= .Machine$double.eps)$root))
      }else{
           fptx <- rep(0,as.numeric(object$M))
                     }
    if (as.numeric(object$y0) > Bn(as.numeric(object$t0))){
           v21  <- sapply(1:length(F2),function(i) ifelse(length(which(R2[,i] <= Bn(time(object))))==0, NA,min(which(R2[,i] <= Bn(time(object))))))
           fpty <- sapply(1:length(F2),function(i) ifelse(is.na(v21[i]),NA,stats::uniroot(f = function(x) (F2[[i]](x)- Bn(x)),lower=time(object)[v21[i]-1],upper=time(object)[v21[i]],tol= .Machine$double.eps)$root))
      }else if (as.numeric(object$y0) < Bn(as.numeric(object$t0))){
           v22  <- sapply(1:length(F2),function(i) ifelse(length(which(R2[,i] <= Bn(time(object))))==as.numeric(object$N)+1, NA,min(which(R2[,i] >= Bn(time(object))))))
           fpty <- sapply(1:length(F2),function(i) ifelse(is.na(v22[i]),NA,stats::uniroot(f = function(x) (F2[[i]](x)- Bn(x)),lower=time(object)[v22[i]-1],upper=time(object)[v22[i]],tol= .Machine$double.eps)$root))
      }else{
           fpty <- rep(0,as.numeric(object$M))
                     }
    #	
	#if (length(which(is.na(fptx))) != length(which(is.na(fpty)))) cat("missing output are removed\n")
	fpt <- na.omit(data.frame(fptx,fpty))
	out <- data.frame(fpt$fptx,fpt$fpty)
	names(out) <- c("x","y")
	#return(out)
    structure(list(obj=object,boundary=boundary[[1]],fpt=out,call=match.call()),class="fptsde2d")
}

print.fptsde2d <- function(x, digits=NULL, ...)
           {
    class(x) <- "fptsde2d"
	Ito = "It\xf4"
    Encoding(Ito) <- "latin1"
    S <- function(t) eval(x$boundary)
	Drx <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)))), list(e = x$obj$driftx))))   
    DDx <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)))), list(e = x$obj$diffx))))
	Dry <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)))), list(e = x$obj$drifty))))   
    DDy <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)))), list(e = x$obj$diffy)))) 
	# Drx <- gsub(pattern = 'x', replacement = 'X(t)', x = gsub(pattern = 'y', replacement = 'Y(t)', x = as.expression(x$obj$driftx), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
	# DDx <- gsub(pattern = 'x', replacement = 'X(t)', x = gsub(pattern = 'y', replacement = 'Y(t)', x = as.expression(x$obj$diffx), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
	# Dry <- gsub(pattern = 'x', replacement = 'X(t)', x = gsub(pattern = 'y', replacement = 'Y(t)', x = as.expression(x$obj$drifty), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
	# DDy <- gsub(pattern = 'x', replacement = 'X(t)', x = gsub(pattern = 'y', replacement = 'Y(t)', x = as.expression(x$obj$diffy), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
    if (as.numeric(x$obj$x0) > S(as.numeric(x$obj$t0))){
    fpt_x <- paste("T(S(t),X(t)) = inf{t >= ",as.numeric(x$obj$t0) ,": X(t) <= ", deparse(x$boundary),"}")
    }else{
    fpt_x <- paste("T(S(t),X(t)) = inf{t >= ",as.numeric(x$obj$t0) ,": X(t) >= ", deparse(x$boundary),"}")
    }
    if (as.numeric(x$obj$y0) > S(as.numeric(x$obj$t0))){
    fpt_y <- paste("T(S(t),Y(t)) = inf{t >= ",as.numeric(x$obj$t0) ,": Y(t) <= ", deparse(x$boundary),"}")
    }else{
    fpt_y <- paste("T(S(t),Y(t)) = inf{t >= ",as.numeric(x$obj$t0) ,": Y(t) >= ", deparse(x$boundary),"}")
    }
    if(x$obj$type=="ito"){
    cat(Ito," Sde 2D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," * dW1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," * dW2(t)","\n",
		"\t| t in [",format(x$obj$t0,digits=digits),",",format(x$obj$T,digits=digits),"].","\n",
        "Boundary:","\n",
        "\t| S(t) = ",deparse(x$boundary),"\n",
        "F.P.T:","\n",
        "\t| ",fpt_x,"\n",
        "\t| \tAnd ","\n",
        "\t| ",fpt_y,"\n",
        "\t| Crossing realized ",format(dim(x$fpt)[1],digits=digits)," among ",format(x$obj$M,digits=digits),".","\n",      
        sep="")}else{
    cat("Stratonovich Sde 2D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," o dW1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," o dW2(t)","\n",
		"\t| t in [",format(x$obj$t0,digits=digits),",",format(x$obj$T,digits=digits),"].","\n",
        "Boundary:","\n",
        "\t| S(t) = ",deparse(x$boundary),"\n",
        "F.P.T:","\n",
        "\t| ",fpt_x,"\n",
        "\t| \tAnd ","\n",
        "\t| ",fpt_y,"\n",
        "\t| Crossing realized ",format(dim(x$fpt)[1],digits=digits)," among ",format(x$obj$M,digits=digits),".","\n",      
        sep="")}	
    invisible(x)
}


###

summary.fptsde2d <- function(object,digits=NULL, ...)
           {  
    class(object) <- "fptsde2d"
	if (is.null(digits)){digits = base::options()$digits}
    S <- function(t) eval(object$boundary)
    if (as.numeric(object$obj$x0) > S(as.numeric(object$obj$t0))){
    fpt_x <- paste("T(S(t),X(t)) = inf{t >= ",as.numeric(object$obj$t0) ,": X(t) <= ", deparse(object$boundary),"}")
    }else{
    fpt_x <- paste("T(S(t),X(t)) = inf{t >= ",as.numeric(object$obj$t0) ,": X(t) >= ", deparse(object$boundary),"}")
    }
    if (as.numeric(object$obj$y0) > S(as.numeric(object$obj$t0))){
    fpt_y <- paste("T(S(t),Y(t)) = inf{t >= ",as.numeric(object$obj$t0) ,": Y(t) <= ", deparse(object$boundary),"}")
    }else{
    fpt_y <- paste("T(S(t),Y(t)) = inf{t >= ",as.numeric(object$obj$t0) ,": Y(t) >= ", deparse(object$boundary),"}")
    }
    cat("\nMonte-Carlo Statistics for the F.P.T of (X(t),Y(t))","\n",
        "\t| ",fpt_x,"\n",
        "\t| ","\t And","\n",
        "\t| ",fpt_y,"\n",
       sep="")
    x <- object$fpt[,1]
    y <- object$fpt[,2]
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
	                         "Coef-variation","3th-order moment","4th-order moment","5th-order moment","6th-order moment"),c("T(S,X)","T(S,Y)"))
    print(round(res,digits=digits), quote = FALSE, right = TRUE,...)
    invisible(object)
}

mean.fptsde2d <- function(x,...)
                    {
    class(x) <- "fptsde2d"
	return(c(mean(x$fpt[,1],na.rm = TRUE,...),mean(x$fpt[,2],na.rm = TRUE,...)))
}

cv.fptsde2d <- function(x,...)
                    {
    class(x) <- "fptsde2d"
	return(c(cv(x$fpt[,1],...),cv(x$fpt[,2],...)))
}

min.fptsde2d <- function(x,...)
                    {
    class(x) <- "fptsde2d"
	return(c(min(x$fpt[,1],na.rm = TRUE),min(x$fpt[,2],na.rm = TRUE)))
}

max.fptsde2d <- function(x,...)
                    {
    class(x) <- "fptsde2d"
	return(c(max(x$fpt[,1],na.rm = TRUE),max(x$fpt[,2],na.rm = TRUE)))
}

skewness.fptsde2d <- function(x,...)
                    {
    class(x) <- "fptsde2d"
	return(c(skewness(x$fpt[,1]),skewness(x$fpt[,2])))
}

kurtosis.fptsde2d <- function(x,...)
                    {
    class(x) <- "fptsde2d"
	return(c(kurtosis(x$fpt[,1]),kurtosis(x$fpt[,2])))
}

Median.fptsde2d <- function(x,...)
                    {
    class(x) <- "fptsde2d"
	return(c(Median(x$fpt[,1]),Median(x$fpt[,2])))
}

Mode.fptsde2d <- function(x,...)
                    {
    class(x) <- "fptsde2d"
	return(c(Mode(x$fpt[,1]),Mode(x$fpt[,2])))
}

quantile.fptsde2d <- function(x,...)
                    {
    class(x) <- "fptsde2d"
	return(list(x=quantile(x$fpt[,1],...),y=quantile(x$fpt[,2],...)))
}

moment.fptsde2d <- function(x,...)
                    {
    class(x) <- "fptsde2d"
	return(c(moment(x$fpt[,1],...),moment(x$fpt[,2],...)))
}

## density fpt2d

dfptsde2d <- function(object, ...)  UseMethod("dfptsde2d")

dfptsde2d.default <- function(object,pdf=c("Joint","Marginal"),...)
                     {
    class(object) <- "fptsde2d"
    # if (any(!is.expression(boundary))) stop(" must be expression of a constant or time-dependent boundary ")
    # R1   <- data.frame(object$X)
    # R2   <- data.frame(object$Y)
    # Bn  <- function(t)  eval(boundary) 
    # F1   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),as.vector(R1[,i])) )
    # F2   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),as.vector(R2[,i])) )
    # if (as.numeric(object$x0) > Bn(as.numeric(object$t0))){
           # v11  <- sapply(1:length(F1),function(i) ifelse(length(which(R1[,i] <= Bn(time(object))))==0, NA,min(which(R1[,i] <= Bn(time(object))))))
           # fptx <- sapply(1:length(F1),function(i) ifelse(is.na(v11[i]),NA,uniroot(f = function(x) (F1[[i]](x)- Bn(x)),lower=time(object)[v11[i]-1],upper=time(object)[v11[i]],tol= .Machine$double.eps)$root))
      # }else if (as.numeric(object$x0) < Bn(as.numeric(object$t0))){
           # v12  <- sapply(1:length(F1),function(i) ifelse(length(which(R1[,i] <= Bn(time(object))))==as.numeric(object$N)+1, NA,min(which(R1[,i] >= Bn(time(object))))))
           # fptx <- sapply(1:length(F1),function(i) ifelse(is.na(v12[i]),NA,uniroot(f = function(x) (F1[[i]](x)- Bn(x)),lower=time(object)[v12[i]-1],upper=time(object)[v12[i]],tol= .Machine$double.eps)$root))
      # }else{
           # fptx <- rep(0,as.numeric(object$M))
                     # }
    # if (as.numeric(object$y0) > Bn(as.numeric(object$t0))){
           # v21  <- sapply(1:length(F2),function(i) ifelse(length(which(R2[,i] <= Bn(time(object))))==0, NA,min(which(R2[,i] <= Bn(time(object))))))
           # fpty <- sapply(1:length(F2),function(i) ifelse(is.na(v21[i]),NA,uniroot(f = function(x) (F2[[i]](x)- Bn(x)),lower=time(object)[v21[i]-1],upper=time(object)[v21[i]],tol= .Machine$double.eps)$root))
      # }else if (as.numeric(object$y0) < Bn(as.numeric(object$t0))){
           # v22  <- sapply(1:length(F2),function(i) ifelse(length(which(R2[,i] <= Bn(time(object))))==as.numeric(object$N)+1, NA,min(which(R2[,i] >= Bn(time(object))))))
           # fpty <- sapply(1:length(F2),function(i) ifelse(is.na(v22[i]),NA,uniroot(f = function(x) (F2[[i]](x)- Bn(x)),lower=time(object)[v22[i]-1],upper=time(object)[v22[i]],tol= .Machine$double.eps)$root))
      # }else{
           # fpty <- rep(0,as.numeric(object$M))
                     # }
    # if (length(which(is.na(fptx))) != length(which(is.na(fpty)))) cat("missing output are removed\n")
	# fpt <- na.omit(data.frame(fptx,fpty))
	# out <- data.frame(fpt$fptx,fpt$fpty)
	out <- object$fpt
	#names(out) <- c("x","y")
    pdf <- match.arg(pdf)
    if (pdf=="Marginal"){
    structure(list(fpt=out,resx=density.default(out[,"x"],na.rm = TRUE,...),resy=density.default(out[,"y"],na.rm = TRUE,...),boundary=object$boundary,pdf=pdf,SDE=object),class="dfptsde2d")
    }else{
    structure(list(fpt=out,res=MASS::kde2d(out$x, out$y, ...),boundary=object$boundary,pdf=pdf,SDE=object),class="dfptsde2d")
    }
}

print.dfptsde2d <- function(x, digits=NULL, ...)
           {
    class(x) <- "dfptsde2d"
	if (is.null(digits)){digits = base::options()$digits}
    S <- function(t) eval(x$boundary)
    if (x$pdf=="Joint") {
    if (as.numeric(x$SDE$obj$x0) > S(as.numeric(x$SDE$obj$t0))){
     fpt_x <- paste("T(S,X) = inf{t >= 0 : X(t) <= ", deparse(x$boundary),"} ")
     }else{
     fpt_x <- paste("T(S,X) = inf{t >= 0 : X(t) >= ", deparse(x$boundary)," }")
     }
     if (as.numeric(x$SDE$obj$y0) > S(as.numeric(x$SDE$obj$t0))){
     fpt_y <- paste("T(S,Y) = inf{t >= 0 : Y(t) <= ", deparse(x$boundary),"}")
     }else{
     fpt_y <- paste("T(S,Y) = inf{t >= 0 : Y(t) >= ", deparse(x$boundary),"}")
     }
     cat("\nJoint density for the F.P.T of (X(t),Y(t))","\n",
        "\t| ",fpt_x,"\n",
        "\t| \tAnd ","\n",
        "\t| ",fpt_y,"\n",
        sep="")
    cat(
	"\nData: (x,y)", " (2 x ", dim(x$fpt)[1], " obs.);", "\n\n", sep = "")
	##"\tBandwidth 'bw' = ", formatC(MASS::bandwidth.nrd(x$res$x), digits = digits),"~~", formatC(MASS::bandwidth.nrd(x$res$y), digits = digits), "\n\n", sep = "")
    out3 <- list(x=x$res$x,y=x$res$y,z=as.vector(x$res$z))
    names(out3) <- c("x","y","f(x,y)")
    print(summary.data.frame(out3, digits = digits), ...)}else{
    if (as.numeric(x$SDE$obj$x0) > S(as.numeric(x$SDE$obj$t0))){
    cat("\nMarginal density for the F.P.T of X(t)","\n", 
	    "\t| T(S,X) = inf{t >= ",x$SDE$obj$t0 ," : X(t) <= ", deparse(x$boundary),"}","\n",
        sep="")
    }else{
    cat("\nMarginal density for the F.P.T of X(t)","\n", 
	    "\t| T(S,X) = inf{t >= ",x$SDE$obj$t0 ," : X(t) >= ", deparse(x$boundary),"}","\n",
        sep="")}
    cat(
	"\nData: ", x$resx$data.name, " (", x$resx$n, " obs.);",
	"\tBandwidth 'bw' = ", formatC(x$resx$bw, digits = digits), "\n\n", sep = "")
    out1 <- as.data.frame(x$resx[c("x","y")])
    names(out1) <- c("x","f(x)")
    print(summary(out1, digits = digits), ...)

    if (as.numeric(x$SDE$obj$y0) > S(as.numeric(x$SDE$obj$t0))){
    cat("\nMarginal density for the F.P.T of Y(t)","\n", 
	    "\t| T(S,Y) = inf{t >= ",x$SDE$obj$t0 ," : Y(t) <= ", deparse(x$boundary),"}","\n",
        sep="")
    }else{
    cat("\nMarginal density for the F.P.T of Y(t)","\n", 
	    "\t| T(S,Y) = inf{t >= ",x$SDE$obj$t0 ," : Y(t) >= ", deparse(x$boundary),"}","\n",
        sep="")}
    cat(
	"\nData: ", x$resy$data.name, " (", x$resy$n, " obs.);",
	"\tBandwidth 'bw' = ", formatC(x$resy$bw, digits = digits), "\n\n", sep = "")
    out2 <- as.data.frame(x$resy[c("x","y")])
    names(out2) <- c("y","f(y)")
    print(summary(out2, digits = digits), ...)
   }
    invisible(x)
}

plot.dfptsde2d <- function(x,display=c("persp","rgl","image","contour"),hist=FALSE,...) .plot.dfptsde2d(x,display,hist,...)

################################################################################
################################################################################
#####
##### fptsde3d

fptsde3d <- function(object, ...)  UseMethod("fptsde3d")

fptsde3d.default <- function(object,boundary,...)
                     {
    class(object) <- "snssde3d"
    if (any(!is.expression(boundary))) 
	         stop(" must be expression of a constant or time-dependent boundary ")
    R1   <- data.frame(object$X)
    R2   <- data.frame(object$Y)
    R3   <- data.frame(object$Z)	
    Bn  <- function(t)  eval(boundary)
    if(object$x0==Bn(object$t0) | object$y0==Bn(object$t0) | object$z0==Bn(object$t0) ) 
	         warning(paste("x0 = S(t0) or y0 =S(t0) or z0 =S(t0) ==> crossing realized at time",object$t0))
    F1   <- lapply(1:as.numeric(object$M),function(i) stats::approxfun(time(object),as.vector(R1[,i])) )
    F2   <- lapply(1:as.numeric(object$M),function(i) stats::approxfun(time(object),as.vector(R2[,i])) )
    F3   <- lapply(1:as.numeric(object$M),function(i) stats::approxfun(time(object),as.vector(R3[,i])) )
    if (as.numeric(object$x0) > Bn(as.numeric(object$t0))){
           v11  <- sapply(1:length(F1),function(i) ifelse(length(which(R1[,i] <= Bn(time(object))))==0, NA,min(which(R1[,i] <= Bn(time(object))))))
           fptx <- sapply(1:length(F1),function(i) ifelse(is.na(v11[i]),NA,stats::uniroot(f = function(x) (F1[[i]](x)- Bn(x)),lower=time(object)[v11[i]-1],upper=time(object)[v11[i]],tol= .Machine$double.eps)$root))
      }else if (as.numeric(object$x0) < Bn(as.numeric(object$t0))){
           v12  <- sapply(1:length(F1),function(i) ifelse(length(which(R1[,i] <= Bn(time(object))))==as.numeric(object$N)+1, NA,min(which(R1[,i] >= Bn(time(object))))))
           fptx <- sapply(1:length(F1),function(i) ifelse(is.na(v12[i]),NA,stats::uniroot(f = function(x) (F1[[i]](x)- Bn(x)),lower=time(object)[v12[i]-1],upper=time(object)[v12[i]],tol= .Machine$double.eps)$root))
      }else{
           fptx <- rep(0,as.numeric(object$M))
                     }
    if (as.numeric(object$y0) > Bn(as.numeric(object$t0))){
           v21  <- sapply(1:length(F2),function(i) ifelse(length(which(R2[,i] <= Bn(time(object))))==0, NA,min(which(R2[,i] <= Bn(time(object))))))
           fpty <- sapply(1:length(F2),function(i) ifelse(is.na(v21[i]),NA,stats::uniroot(f = function(x) (F2[[i]](x)- Bn(x)),lower=time(object)[v21[i]-1],upper=time(object)[v21[i]],tol= .Machine$double.eps)$root))
      }else if (as.numeric(object$y0) < Bn(as.numeric(object$t0))){
           v22  <- sapply(1:length(F2),function(i) ifelse(length(which(R2[,i] <= Bn(time(object))))==as.numeric(object$N)+1, NA,min(which(R2[,i] >= Bn(time(object))))))
           fpty <- sapply(1:length(F2),function(i) ifelse(is.na(v22[i]),NA,stats::uniroot(f = function(x) (F2[[i]](x)- Bn(x)),lower=time(object)[v22[i]-1],upper=time(object)[v22[i]],tol= .Machine$double.eps)$root))
      }else{
           fpty <- rep(0,as.numeric(object$M))
                     }
    if (as.numeric(object$z0) > Bn(as.numeric(object$t0))){
           v31  <- sapply(1:length(F3),function(i) ifelse(length(which(R3[,i] <= Bn(time(object))))==0, NA,min(which(R3[,i] <= Bn(time(object))))))
           fptz <- sapply(1:length(F3),function(i) ifelse(is.na(v31[i]),NA,stats::uniroot(f = function(x) (F3[[i]](x)- Bn(x)),lower=time(object)[v31[i]-1],upper=time(object)[v31[i]],tol= .Machine$double.eps)$root))
      }else if (object$z0 < Bn(as.numeric(object$t0))){
           v32  <- sapply(1:length(F3),function(i) ifelse(length(which(R3[,i] <= Bn(time(object))))==as.numeric(object$N)+1, NA,min(which(R3[,i] >= Bn(time(object))))))
           fptz <- sapply(1:length(F3),function(i) ifelse(is.na(v32[i]),NA,stats::uniroot(f = function(x) (F3[[i]](x)- Bn(x)),lower=time(object)[v32[i]-1],upper=time(object)[v32[i]],tol= .Machine$double.eps)$root))
      }else{
           fptz <- rep(0,as.numeric(object$M))
                     }
    #structure(list(SDE=object,boundary=boundary[[1]],fptx=fptx,fpty=fpty,fptz=fptz),class="fptsde3d")
	#if (length(which(is.na(fptx))) != length(which(is.na(fpty))) | length(which(is.na(fptx))) != length(which(is.na(fptz))) | length(which(is.na(fpty))) != length(which(is.na(fptz)))) cat("missing output are removed\n")
	fpt <- na.omit(data.frame(fptx,fpty,fptz))
	out <- data.frame(fpt$fptx,fpt$fpty,fpt$fptz)
	names(out) <- c("x","y","z")
    structure(list(obj=object,boundary=boundary[[1]],fpt=out,call=match.call()),class="fptsde3d")
}

print.fptsde3d <- function(x, digits=NULL, ...)
           {
    class(x) <- "fptsde3d"
	Ito = "It\xf4"
    Encoding(Ito) <- "latin1"
    S <- function(t) eval(x$boundary)
    Drx <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$obj$driftx))))   
    DDx <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$obj$diffx))))
	Dry <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$obj$drifty))))   
    DDy <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$obj$diffy))))
    Drz <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$obj$driftz))))   
    DDz <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$obj$diffz))))
    # Drx <- gsub(pattern = 'x', replacement = 'X(t)', x = gsub(pattern = 'y', replacement = 'Y(t)', x = gsub(pattern = 'z', replacement = 'Z(t)', x = as.expression(x$obj$driftx), ignore.case = F,fixed = T), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
	# DDx <- gsub(pattern = 'x', replacement = 'X(t)', x = gsub(pattern = 'y', replacement = 'Y(t)', x = gsub(pattern = 'z', replacement = 'Z(t)', x = as.expression(x$obj$diffx), ignore.case = F,fixed = T), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
    # Dry <- gsub(pattern = 'x', replacement = 'X(t)', x = gsub(pattern = 'y', replacement = 'Y(t)', x = gsub(pattern = 'z', replacement = 'Z(t)', x = as.expression(x$obj$drifty), ignore.case = F,fixed = T), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
	# DDy <- gsub(pattern = 'x', replacement = 'X(t)', x = gsub(pattern = 'y', replacement = 'Y(t)', x = gsub(pattern = 'z', replacement = 'Z(t)', x = as.expression(x$obj$diffy), ignore.case = F,fixed = T), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
	# Drz <- gsub(pattern = 'x', replacement = 'X(t)', x = gsub(pattern = 'y', replacement = 'Y(t)', x = gsub(pattern = 'z', replacement = 'Z(t)', x = as.expression(x$obj$driftz), ignore.case = F,fixed = T), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
	# DDz <- gsub(pattern = 'x', replacement = 'X(t)', x = gsub(pattern = 'y', replacement = 'Y(t)', x = gsub(pattern = 'z', replacement = 'Z(t)', x = as.expression(x$obj$diffz), ignore.case = F,fixed = T), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
    if (as.numeric(x$obj$x0) > S(as.numeric(x$obj$t0))){
    fpt_x <- paste("T(S(t),X(t)) = inf{t >= ",as.numeric(x$obj$t0) ,": X(t) <= ", deparse(x$boundary),"}")
    }else{
    fpt_x <- paste("T(S(t),X(t)) = inf{t >= ",as.numeric(x$obj$t0) ,": X(t) >= ", deparse(x$boundary),"}")
    }
    if (as.numeric(x$obj$y0) > S(as.numeric(x$obj$t0))){
    fpt_y <- paste("T(S(t),Y(t)) = inf{t >= ",as.numeric(x$obj$t0) ,": Y(t) <= ", deparse(x$boundary),"}")
    }else{
    fpt_y <- paste("T(S(t),Y(t)) = inf{t >= ",as.numeric(x$obj$t0) ,": Y(t) >= ", deparse(x$boundary),"}")
    }
    if (as.numeric(x$obj$z0) > S(as.numeric(x$obj$t0))){
    fpt_z <- paste("T(S(t),Z(t)) = inf{t >= ",as.numeric(x$obj$t0) ,": Z(t) <= ", deparse(x$boundary),"}")
    }else{
    fpt_z <- paste("T(S(t),Z(t)) = inf{t >= ",as.numeric(x$obj$t0) ,": Z(t) >= ", deparse(x$boundary),"}")
    }
    if(x$obj$type=="ito"){
    cat(Ito," Sde 3D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," * dW1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," * dW2(t)","\n",
        "\t| dZ(t) = ", Drz," * dt + ", DDz," * dW3(t)","\n",
		"\t| t in [",format(x$obj$t0,digits=digits),",",format(x$obj$T,digits=digits),"].","\n",
        "Boundary:","\n",
        "\t| S(t) = ",deparse(x$boundary),"\n",
        "F.P.T:","\n",
        "\t| ",fpt_x,"\n",
        "\t| \tAnd ","\n",
        "\t| ",fpt_y,"\n",
        "\t| \tAnd ","\n",
        "\t| ",fpt_z,"\n",
        "\t| Crossing realized ",format(dim(x$fpt)[1],digits=digits)," among ",format(x$obj$M,digits=digits),".","\n",      
        sep="")}else{
    cat("Stratonovich Sde 3D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," o dW1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," o dW2(t)","\n",
        "\t| dZ(t) = ", Drz," * dt + ", DDz," o dW3(t)","\n",
		"\t| t in [",format(x$obj$t0,digits=digits),",",format(x$obj$T,digits=digits),"].","\n",
        "Boundary:","\n",
        "\t| S(t) = ",deparse(x$boundary),"\n",
        "F.P.T:","\n",
        "\t| ",fpt_x,"\n",
        "\t| \tAnd ","\n",
        "\t| ",fpt_y,"\n",
        "\t| \tAnd ","\n",
        "\t| ",fpt_z,"\n",
        "\t| Crossing realized ",format(dim(x$fpt)[1],digits=digits)," among ",format(x$obj$M,digits=digits),".","\n",      
        sep="")}	
    invisible(x)
}


###

summary.fptsde3d  <- function(object,digits=NULL, ...)
           {   
    class(object) <- "fptsde3d"
	if (is.null(digits)){digits = base::options()$digits}
    S <- function(t) eval(object$boundary)
    if (as.numeric(object$obj$x0) > S(as.numeric(object$obj$t0))){
    fpt_x <- paste("T(S(t),X(t)) = inf{t >= ",as.numeric(object$obj$t0) ,": X(t) <= ", deparse(object$boundary),"}")
    }else{
    fpt_x <- paste("T(S(t),X(t)) = inf{t >= ",as.numeric(object$obj$t0) ,": X(t) >= ", deparse(object$boundary),"}")
    }
    if (as.numeric(object$obj$y0) > S(as.numeric(object$obj$t0))){
    fpt_y <- paste("T(S(t),Y(t)) = inf{t >= ",as.numeric(object$obj$t0) ,": Y(t) <= ", deparse(object$boundary),"}")
    }else{
    fpt_y <- paste("T(S(t),Y(t)) = inf{t >= ",as.numeric(object$obj$t0) ,": Y(t) >= ", deparse(object$boundary),"}")
    }
    if (as.numeric(object$obj$z0) > S(as.numeric(object$obj$t0))){
    fpt_z <- paste("T(S(t),Z(t)) = inf{t >= ",as.numeric(object$obj$t0) ,": Z(t) <= ", deparse(object$boundary),"}")
    }else{
    fpt_z <- paste("T(S(t),Z(t)) = inf{t >= ",as.numeric(object$obj$t0) ,": Z(t) >= ", deparse(object$boundary),"}")
    }
    cat("\nMonte-Carlo Statistics for the F.P.T of (X(t),Y(t),Z(t))","\n",
        "\t| ",fpt_x,"\n",
        "\t| ","\t And","\n",
        "\t| ",fpt_y,"\n",
        "\t| ","\t And","\n",
        "\t| ",fpt_z,"\n",
       sep="")
    x <- object$fpt[,1]
    y <- object$fpt[,2]
    z <- object$fpt[,3]	
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
	                         "Coef-variation","3th-order moment","4th-order moment","5th-order moment","6th-order moment"),c("T(S,X)","T(S,Y)","T(S,Z)"))
    print(round(res,digits=digits), quote = FALSE, right = TRUE,...)
    invisible(object)
}

mean.fptsde3d <- function(x,...)
                    {
    class(x) <- "fptsde3d"
	return(c(mean(x$fpt[,1],na.rm = TRUE,...),mean(x$fpt[,2],na.rm = TRUE,...),mean(x$fpt[,3],na.rm = TRUE,...)))
}

cv.fptsde3d <- function(x,...)
                    {
    class(x) <- "fptsde3d"
	return(c(cv(x$fpt[,1],...),cv(x$fpt[,2],...),cv(x$fpt[,3],...)))
}

min.fptsde3d <- function(x,...)
                    {
    class(x) <- "fptsde3d"
	return(c(min(x$fpt[,1],na.rm = TRUE),min(x$fpt[,2],na.rm = TRUE),min(x$fpt[,3],na.rm = TRUE)))
}

max.fptsde3d <- function(x,...)
                    {
    class(x) <- "fptsde3d"
	return(c(max(x$fpt[,1],na.rm = TRUE),max(x$fpt[,2],na.rm = TRUE),max(x$fpt[,3],na.rm = TRUE)))
}

skewness.fptsde3d <- function(x,...)
                    {
    class(x) <- "fptsde3d"
	return(c(skewness(x$fpt[,1]),skewness(x$fpt[,2]),skewness(x$fpt[,3])))
}

kurtosis.fptsde3d <- function(x,...)
                    {
    class(x) <- "fptsde3d"
	return(c(kurtosis(x$fpt[,1]),kurtosis(x$fpt[,2]),kurtosis(x$fpt[,3])))
}

Median.fptsde3d <- function(x,...)
                    {
    class(x) <- "fptsde3d"
	return(c(Median(x$fpt[,1]),Median(x$fpt[,2]),Median(x$fpt[,3])))
}

Mode.fptsde3d <- function(x,...)
                    {
    class(x) <- "fptsde3d"
	return(c(Mode(x$fpt[,1]),Mode(x$fpt[,2]),Mode(x$fpt[,3])))
}

quantile.fptsde3d <- function(x,...)
                    {
    class(x) <- "fptsde3d"
	return(list(x=quantile(x$fpt[,1],...),y=quantile(x$fpt[,2],...),z=quantile(x$fpt[,3],...)))
}

moment.fptsde3d <- function(x,...)
                    {
    class(x) <- "fptsde3d"
	return(c(moment(x$fpt[,1],...),moment(x$fpt[,2],...),moment(x$fpt[,3],...)))
}

## density fpt3d

dfptsde3d <- function(object, ...)  UseMethod("dfptsde3d")

dfptsde3d.default <- function(object,pdf=c("Joint","Marginal"),...)
                     {
    class(object) <- "fptsde3d"	
	out <- object$fpt
    pdf <- match.arg(pdf)
    if (pdf=="Joint"){
        if (!requireNamespace("sm", quietly = TRUE)) {
                cat("The sm package is not available, you need to install the package.\\n")
                pdf <- "Marginal" }
      }
    if (pdf=="Marginal"){
    structure(list(fpt=out,resx=density.default(out[,"x"],na.rm = TRUE,...),resy=density.default(out[,"y"],na.rm = TRUE,...),resz=density.default(out[,"z"],na.rm = TRUE,...),boundary=object$boundary,pdf=pdf,SDE=object),class="dfptsde3d")
    }else{
	structure(list(fpt=out,res=sm::sm.density(out,display="none", ...),boundary=object$boundary,pdf=pdf,SDE=object),class="dfptsde3d")
    }
}

print.dfptsde3d <- function(x, digits=NULL, ...)
           {
    class(x) <- "dfptsde2d"
	if (is.null(digits)){digits = base::options()$digits}
    S <- function(t) eval(x$boundary)
    if (x$pdf=="Joint") {
    if (as.numeric(x$SDE$obj$x0) > S(as.numeric(x$SDE$obj$t0))){
     fpt_x <- paste("T(S,X) = inf{t >= 0 : X(t) <= ", deparse(x$boundary),"} ")
     }else{
     fpt_x <- paste("T(S,X) = inf{t >= 0 : X(t) >= ", deparse(x$boundary)," }")
     }
     if (as.numeric(x$SDE$obj$y0) > S(as.numeric(x$SDE$obj$t0))){
     fpt_y <- paste("T(S,Y) = inf{t >= 0 : Y(t) <= ", deparse(x$boundary),"}")
     }else{
     fpt_y <- paste("T(S,Y) = inf{t >= 0 : Y(t) >= ", deparse(x$boundary),"}")
     }
     if (as.numeric(x$SDE$obj$z0) > S(as.numeric(x$SDE$obj$t0))){
     fpt_z <- paste("T(S,Z) = inf{t >= 0 : Z(t) <= ", deparse(x$boundary),"}")
     }else{
     fpt_z <- paste("T(S,Z) = inf{t >= 0 : Z(t) >= ", deparse(x$boundary),"}")
     }
     cat("\nJoint density for the F.P.T of (X(t),Y(t),Z(t))","\n",
        "\t| ",fpt_x,"\n",
        "\t| \tAnd ","\n",
        "\t| ",fpt_y,"\n",
        "\t| \tAnd ","\n",
        "\t| ",fpt_z,"\n",
        sep="")
     cat(
	 "\nData: (x,y,z)", " (3 x ", dim(x$fpt)[1], " obs.);", 
	 "\tBandwidth 'bw' = c(", round(x$res$h[1], digits = 3),",",formatC(x$res$h[2], digits = 3),",",formatC(x$res$h[3], digits = 3),")\n\n", sep = "")
    out3 <- list(x=x$res$eval.points[,1],y=x$res$eval.points[,2],z=x$res$eval.points[,3],d=as.vector(x$res$estimate))
    names(out3) <- c("x","y","z","f(x,y,z)")
    print(summary.data.frame(out3,digits = digits), ...)}else{
    if (as.numeric(x$SDE$obj$x0) > S(as.numeric(x$SDE$obj$t0))){
    cat("\nMarginal density for the F.P.T of X(t)","\n", 
	    "\t| T(S,X) = inf{t >= ",x$SDE$obj$t0 ," : X(t) <= ", deparse(x$boundary),"}","\n",
        sep="")
    }else{
    cat("\nMarginal density for the F.P.T of X(t)","\n", 
	    "\t| T(S,X) = inf{t >= ",x$SDE$obj$t0 ," : X(t) >= ", deparse(x$boundary),"}","\n",
        sep="")}
    cat(
	"\nData: ", x$resx$data.name, " (", x$resx$n, " obs.);",
	"\tBandwidth 'bw' = ", formatC(x$resx$bw, digits = digits), "\n\n", sep = "")
    out1 <- as.data.frame(x$resx[c("x","y")])
    names(out1) <- c("x","f(x)")
    print(summary(out1, digits = digits), ...)

    if (as.numeric(x$SDE$obj$y0) > S(as.numeric(x$SDE$obj$t0))){
    cat("\nMarginal density for the F.P.T of Y(t)","\n", 
	    "\t| T(S,Y) = inf{t >= ",x$SDE$obj$t0 ," : Y(t) <= ", deparse(x$boundary),"}","\n",
        sep="")
    }else{
    cat("\nMarginal density for the F.P.T of Y(t)","\n", 
	    "\t| T(S,Y) = inf{t >= ",x$SDE$obj$t0 ," : Y(t) >= ", deparse(x$boundary),"}","\n",
        sep="")}
    cat(
	"\nData: ", x$resy$data.name, " (", x$resy$n, " obs.);",
	"\tBandwidth 'bw' = ", formatC(x$resy$bw, digits = digits), "\n\n", sep = "")
    out2 <- as.data.frame(x$resy[c("x","y")])
    names(out2) <- c("y","f(y)")
    print(summary(out2, digits = digits), ...)

    if (as.numeric(x$SDE$obj$z0) > S(as.numeric(x$SDE$obj$t0))){
    cat("\nMarginal density for the F.P.T of Z(t)","\n", 
	    "\t| T(S,Z) = inf{t >= ",x$SDE$obj$t0 ," : Z(t) <= ", deparse(x$boundary),"}","\n",
        sep="")
    }else{
    cat("\nMarginal density for the F.P.T of Z(t)","\n", 
	    "\t| T(S,Z) = inf{t >= ",x$SDE$obj$t0 ," : Z(t) >= ", deparse(x$boundary),"}","\n",
        sep="")}
    cat(
	"\nData: ", x$resz$data.name, " (", x$resz$n, " obs.);",
	"\tBandwidth 'bw' = ", formatC(x$resz$bw, digits = digits), "\n\n", sep = "")
    out3 <- as.data.frame(x$resy[c("x","y")])
    names(out3) <- c("z","f(z)")
    print(summary(out3, digits = digits), ...)
	}
    invisible(x)
}

plot.dfptsde3d <- function(x,display="rgl",hist=FALSE,...) .plot.dfptsde3d(x,display,hist,...)
