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
##### Milstein1D
 
.Milstein1D <- function(N =1000,M=1,x0=0,t0=0,T=1,Dt=NULL,drift,diffusion,
                          type=c("ito","str"),...)
                       {
    DSx  <- Simplify(Deriv(diffusion,"x",cache.exp=FALSE))  
    if (type=="ito"){A    <- function(t,x)  eval(drift)}else{
	driftstr <- eval(Simplify(substitute(expression(e1 - 0.5 * e2 * de2), list(e1 = drift[[1]], e2 = diffusion[[1]], de2 = Deriv(diffusion,"x",cache.exp=FALSE)[[1]]))))
	A  <- function(t,x) eval(driftstr)
    # A  <- function(t,x)  eval(drift) - 0.5 * eval(diffusion) * eval(DSx)
	}
    S  <- function(t,x)  eval(diffusion)
    Sx <- function(t,x)  eval(DSx)
    if (is.null(Dt)) {
        Dt <- (T - t0)/N
        t <- seq(t0, T, by=Dt)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
		T <- t[N + 1]
    }
    W <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    X <- matrix(x0, N+1, M)
    for (i in 1L:N) {X[i + 1L,] <- X[i,] + A(t[i],X[i,]) * Dt + S(t[i],X[i,]) * W[i,]+0.5* S(t[i],X[i,])*Sx(t[i],X[i,]) *(W[i,]^2 - Dt) }
    name <- "X"
    name <- if(M > 1) paste("X",1:M,sep="")
    X <- ts(X, start = t0, deltat = Dt, names=name)
    return(list(X=X))
}   
      

#####
##### Milstein2D

.Milstein2D <- function(N =1000,M=1,x0=0,y0=0,t0=0,T=1,Dt=NULL,driftx,diffx,drifty,diffy,
                          type=c("ito","str"),...)
                       {					     
    DSx  <- Simplify(Deriv(diffx,"x",cache.exp=FALSE))
    DSy  <- Simplify(Deriv(diffy,"y",cache.exp=FALSE))  
    if (type=="ito"){Ax    <- function(t,x,y) eval(driftx)
                     Ay    <- function(t,x,y) eval(drifty)}else{
	driftstrx <- eval(Simplify(substitute(expression(e1 - 0.5 * e2 * de2), list(e1 = driftx[[1]], e2 = diffx[[1]], de2 = Deriv(diffx,"x",cache.exp=FALSE)[[1]]))))
	driftstry <- eval(Simplify(substitute(expression(e1 - 0.5 * e2 * de2), list(e1 = drifty[[1]], e2 = diffy[[1]], de2 = Deriv(diffy,"y",cache.exp=FALSE)[[1]]))))
    Ax <- function(t,x,y) eval(driftstrx) 
    Ay <- function(t,x,y) eval(driftstry)
    # Ax  <- function(t,x,y) eval(driftx) - 0.5 * eval(diffx) * eval(DSx)
    # Ay  <- function(t,x,y) eval(drifty) - 0.5 * eval(diffy) * eval(DSy)
	}
    Sx  <- function(t,x,y) eval(diffx)
    Sy  <- function(t,x,y) eval(diffy) 
    dSx <- function(t,x,y) eval(DSx)
    dSy <- function(t,x,y) eval(DSy) 
    if (is.null(Dt)) {
        Dt <- (T - t0)/N
        t <- seq(t0, T, by=Dt)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
		T <- t[N + 1]
    }
    W1 <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    W2 <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    X <- matrix(x0, N+1, M)
    Y <- matrix(y0, N+1, M)	
    for (i in 1L:N) {
	X[i + 1L,] <- X[i,] + Ax(t[i],X[i,],Y[i,]) * Dt + Sx(t[i],X[i,],Y[i,]) * W1[i,]+0.5* Sx(t[i],X[i,],Y[i,])*dSx(t[i],X[i,],Y[i,]) *(W1[i,]^2 - Dt) 
	Y[i + 1L,] <- Y[i,] + Ay(t[i],X[i,],Y[i,]) * Dt + Sy(t[i],X[i,],Y[i,]) * W2[i,]+0.5* Sy(t[i],X[i,],Y[i,])*dSy(t[i],X[i,],Y[i,]) *(W2[i,]^2 - Dt)
	}
    name <- c("X","Y")
    name <- if(M > 1) c(paste(name[1],1:M,sep=""),paste(name[2],1:M,sep=""))
    X <- ts(X, start = t0, deltat = Dt, names=name[1:M])
    Y <- ts(Y, start = t0, deltat = Dt, names=name[(M+1):(2*M)])
    return(list(X=X,Y=Y))
}


#####
##### Milstein3D

.Milstein3D <- function(N =1000,M=1,x0=0,y0=0,z0=0,t0=0,T=1,Dt=NULL,driftx,diffx,drifty,diffy,
                     driftz,diffz,type=c("ito","str"),...)
                       {
    DSx  <- Simplify(Deriv(diffx,"x",cache.exp=FALSE))
    DSy  <- Simplify(Deriv(diffy,"y",cache.exp=FALSE))
    DSz  <- Simplify(Deriv(diffz,"z",cache.exp=FALSE))  
    if (type=="ito"){
        Ax    <- function(t,x,y,z) eval(driftx)
        Ay    <- function(t,x,y,z) eval(drifty)
        Az    <- function(t,x,y,z) eval(driftz)}else{
		driftstrx <- eval(Simplify(substitute(expression(e1 - 0.5 * e2 * de2), list(e1 = driftx[[1]], e2 = diffx[[1]], de2 = Deriv(diffx,"x",cache.exp=FALSE)[[1]]))))
		driftstry <- eval(Simplify(substitute(expression(e1 - 0.5 * e2 * de2), list(e1 = drifty[[1]], e2 = diffy[[1]], de2 = Deriv(diffy,"y",cache.exp=FALSE)[[1]]))))
	    driftstrz <- eval(Simplify(substitute(expression(e1 - 0.5 * e2 * de2), list(e1 = driftz[[1]], e2 = diffz[[1]], de2 = Deriv(diffz,"z",cache.exp=FALSE)[[1]]))))
        Ax <- function(t,x,y,z) eval(driftstrx) 
        Ay <- function(t,x,y,z) eval(driftstry)
        Az <- function(t,x,y,z) eval(driftstrz)
    # Ax <- function(t,x,y,z) eval(driftx) - 0.5 * eval(diffx) * eval(DSx)
    # Ay <- function(t,x,y,z) eval(drifty) - 0.5 * eval(diffy) * eval(DSy)
    # Az <- function(t,x,y,z) eval(driftz) - 0.5 * eval(diffz) * eval(DSz)
	}
    Sx  <- function(t,x,y,z) eval(diffx)
    Sy  <- function(t,x,y,z) eval(diffy)
    Sz  <- function(t,x,y,z) eval(diffz) 
    dSx <- function(t,x,y,z) eval(DSx)
    dSy <- function(t,x,y,z) eval(DSy)
    dSz <- function(t,x,y,z) eval(DSz)
    if (is.null(Dt)) {
        Dt <- (T - t0)/N
        t <- seq(t0, T, by=Dt)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
		T <- t[N + 1]
    }
    W1 <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    W2 <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
	W3 <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    X <- matrix(x0, N+1, M)
    Y <- matrix(y0, N+1, M)	
	Z <- matrix(z0, N+1, M)	
    for (i in 1L:N) {
	X[i + 1L,] <- X[i,] + Ax(t[i],X[i,],Y[i,],Z[i,]) * Dt + Sx(t[i],X[i,],Y[i,],Z[i,]) * W1[i,]+0.5* Sx(t[i],X[i,],Y[i,],Z[i,])*dSx(t[i],X[i,],Y[i,],Z[i,]) *(W1[i,]^2 - Dt) 
	Y[i + 1L,] <- Y[i,] + Ay(t[i],X[i,],Y[i,],Z[i,]) * Dt + Sy(t[i],X[i,],Y[i,],Z[i,]) * W2[i,]+0.5* Sy(t[i],X[i,],Y[i,],Z[i,])*dSy(t[i],X[i,],Y[i,],Z[i,]) *(W2[i,]^2 - Dt)
	Z[i + 1L,] <- Z[i,] + Az(t[i],X[i,],Y[i,],Z[i,]) * Dt + Sz(t[i],X[i,],Y[i,],Z[i,]) * W3[i,]+0.5* Sz(t[i],X[i,],Y[i,],Z[i,])*dSz(t[i],X[i,],Y[i,],Z[i,]) *(W3[i,]^2 - Dt)
	}    
    name <- c("X","Y","Z")
    name <- if(M > 1) c(paste(name[1],1:M,sep=""),paste(name[2],1:M,sep=""),paste(name[3],1:M,sep=""))
    X <- ts(X, start = t0, deltat = Dt, names=name[1:M])
    Y <- ts(Y, start = t0, deltat = Dt, names=name[(M+1):(2*M)])
    Z <- ts(Z, start = t0, deltat = Dt, names=name[(2*M+1):(3*M)])
    return(list(X=X,Y=Y,Z=Z))
} 

