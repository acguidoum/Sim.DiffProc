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
##### rsde1d



rsde1d <- function(object, ...)  UseMethod("rsde1d")

rsde1d.default <- function(object,at,...)
                     {
    if (class(object) == "snssde1d") {M=object$M}
	else if (class(object) == "bridgesde1d") {M=length(na.omit(object$C))}
    if (missing(at)) 
	          at = as.numeric(object$T)
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	      stop( " please use 't0 <= at <= T'")	
	if (M == 1){ X = matrix(object$X,nrow=length(object$X),ncol=1)}else{X = object$X}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (M==1){ F  <- lapply(1:M,function(i) stats::approxfun(time(object),object$X))}else{
                       F  <- lapply(1:M,function(i) stats::approxfun(time(object),object$X[,i]))}
                       x <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    return(x)
}

###


#####
##### rsde2d

rsde2d <- function(object, ...)  UseMethod("rsde2d")

rsde2d.default <- function(object,at,...)
                    { 			
	if (class(object) == "snssde2d") {M=object$M}
	else if (class(object) == "bridgesde2d") {M=length(na.omit(object$Cx))}
    if (missing(at)) 
	        at = as.numeric(object$T)
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	         stop( " please use 't0 <= at <= T'")
    if (as.numeric(M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(M)==1){ Fx   <- lapply(1:as.numeric(M),function(i) stats::approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(M),function(i) stats::approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(M)==1){ Fy   <- lapply(1:as.numeric(M),function(i) stats::approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(M),function(i) stats::approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    return(data.frame(x,y))
}

##
##


#####
##### rsde3d

rsde3d <- function(object, ...)  UseMethod("rsde3d")

rsde3d.default <- function(object,at,...)
                    { 
	if (class(object) == "snssde3d") {M=object$M}
	else if (class(object) == "bridgesde3d") {M=length(na.omit(object$Cx))}
    if (missing(at)) 
	         at = as.numeric(object$T)
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  
	         stop( " please use 't0 <= at <= T'")
    if (as.numeric(M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==as.character(at)),])
    y   <- as.vector(Y[which(time(object)==as.character(at)),])
    z   <- as.vector(Z[which(time(object)==as.character(at)),])
    if (length(x) == 0){
	if (as.numeric(M)==1){ Fx   <- lapply(1:as.numeric(M),function(i) stats::approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(M),function(i) stats::approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(M)==1){ Fy   <- lapply(1:as.numeric(M),function(i) stats::approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(M),function(i) stats::approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(M)==1){ Fz   <- lapply(1:as.numeric(M),function(i) stats::approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(M),function(i) stats::approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }
    return(data.frame(x,y,z))
}

###

