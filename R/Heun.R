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
##### Heun1D

.Heun1D <- function(N =1000,M=1,x0=0,t0=0,T=1,Dt=NULL,drift,diffusion,
                          type=c("ito","str"),...)
                       { 					   
    if (type=="ito") {A    <- function(t,x)  eval(drift)}else{
    driftstr <- eval(Simplify(substitute(expression(e1 - 0.5 * e2 * de2), list(e1 = drift[[1]], e2 = diffusion[[1]], de2 = Deriv(diffusion,"x",cache.exp=FALSE)[[1]]))))
	A  <- function(t,x) eval(driftstr)
    #A  <- function(t,x) eval(drift) - 0.5 * eval(diffusion) * eval(Deriv(diffusion,"x"))
	}
    S  <- function(t,x) eval(diffusion)
    if (is.null(Dt)) {
        Dt <- (T - t0)/N
        t <- seq(t0, T, by=Dt)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
		T <- t[N + 1]
    }
    Wu <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
	Wd <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    X <- matrix(x0, N+1, M)	
	XX <- matrix(x0, N+1, M)	
    for (i in 1L:N) {
	XX[i + 1L,]= X[i,]+A(t[i],X[i,])*Dt+S(t[i],X[i,])*Wu[i,]
	X[i + 1L,] <- X[i,] + 0.5*(A(t[i],X[i,])+A(t[i],XX[i,])) * Dt +0.5* (S(t[i],X[i,])+S(t[i],XX[i,])) * Wd[i,]  }
    name <- "X"
    name <- if(M > 1) paste("X",1:M,sep="")
    X <- ts(X, start = t0, deltat = Dt, names=name)
    return(list(X=X))
}   

  
#####
##### Heun2D

.Heun2D <- function(N =1000,M=1,x0=0,y0=0,t0=0,T=1,Dt=NULL,driftx,diffx,drifty,diffy,
                          type=c("ito","str"),order=c(1,2,3),...)
                       {
    if (type=="ito"){
    Ax <- function(t,x,y)  eval(driftx)
    Ay <- function(t,x,y)  eval(drifty) }else{
	driftstrx <- eval(Simplify(substitute(expression(e1 - 0.5 * e2 * de2), list(e1 = driftx[[1]], e2 = diffx[[1]], de2 = Deriv(diffx,"x",cache.exp=FALSE)[[1]]))))
	driftstry <- eval(Simplify(substitute(expression(e1 - 0.5 * e2 * de2), list(e1 = drifty[[1]], e2 = diffy[[1]], de2 = Deriv(diffy,"y",cache.exp=FALSE)[[1]]))))
    Ax <- function(t,x,y) eval(driftstrx) 
    Ay <- function(t,x,y) eval(driftstry)
    # Ax <- function(t,x,y) eval(driftx) - 0.5 * eval(diffx) * eval(Deriv(diffx,"x"))
    # Ay <- function(t,x,y) eval(drifty) - 0.5 * eval(diffy) * eval(Deriv(diffy,"y"))
                         }
    Sx <- function(t,x,y) eval(diffx)
    Sy <- function(t,x,y) eval(diffy) 
    if (is.null(Dt)) {
        Dt <- (T - t0)/N
        t <- seq(t0, T, by=Dt)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
		T <- t[N + 1]
    }
    Wux <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
	Wdx <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    Wuy <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
	Wdy <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    X = XX <- matrix(x0, N+1, M)	
	Y = YY <- matrix(y0, N+1, M)		
    for (i in 1L:N) {
	XX[i + 1L,] <- X[i,]+Ax(t[i],X[i,],Y[i,])*Dt+Sx(t[i],X[i,],Y[i,])*Wux[i,]
	YY[i + 1L,] <- Y[i,]+Ay(t[i],X[i,],Y[i,])*Dt+Sy(t[i],X[i,],Y[i,])*Wuy[i,]
	X[i + 1L,] <- X[i,] + 0.5*(Ax(t[i],X[i,],Y[i,])+Ax(t[i],XX[i,],Y[i,])) * Dt +0.5* (Sx(t[i],X[i,],Y[i,])+Sx(t[i],XX[i,],Y[i,])) * Wdx[i,]
    Y[i + 1L,] <- Y[i,] + 0.5*(Ay(t[i],X[i,],Y[i,])+Ay(t[i],X[i,],YY[i,])) * Dt +0.5* (Sy(t[i],X[i,],Y[i,])+Sy(t[i],X[i,],YY[i,])) * Wdy[i,]	
	}	
    name <- c("X","Y")
    name <- if(M > 1) c(paste(name[1],1:M,sep=""),paste(name[2],1:M,sep=""))
    X <- ts(X, start = t0, deltat = Dt, names=name[1:M])
    Y <- ts(Y, start = t0, deltat = Dt, names=name[(M+1):(2*M)])
    return(list(X=X,Y=Y))
}

#####
##### Heun3D

.Heun3D <- function(N =1000,M=1,x0=0,y0=0,z0=0,t0=0,T=1,Dt=NULL,driftx,diffx,drifty,diffy,
                     driftz,diffz,type=c("ito","str"),...)
                       {
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
    # Ax <- function(t,x,y,z) eval(driftx) - 0.5 * eval(diffx) * eval(Deriv(diffx,"x"))
    # Ay <- function(t,x,y,z) eval(drifty) - 0.5 * eval(diffy) * eval(Deriv(diffy,"y"))
    # Az <- function(t,x,y,z) eval(driftz) - 0.5 * eval(diffz) * eval(Deriv(diffz,"z"))
	}
    Sx <- function(t,x,y,z) eval(diffx)
    Sy <- function(t,x,y,z) eval(diffy) 
    Sz <- function(t,x,y,z) eval(diffz)
    if (is.null(Dt)) {
        Dt <- (T - t0)/N
        t <- seq(t0, T, by=Dt)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
		T <- t[N + 1]
    }
    Wux <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
	Wdx <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    Wuy <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
	Wdy <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    Wuz <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
	Wdz <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    X = XX <- matrix(x0, N+1, M)	
	Y = YY <- matrix(y0, N+1, M)
	Z = ZZ <- matrix(z0, N+1, M)	
    for (i in 1L:N) {
	XX[i + 1L,] <- X[i,]+Ax(t[i],X[i,],Y[i,],Z[i,])*Dt+Sx(t[i],X[i,],Y[i,],Z[i,])*Wux[i,]
	YY[i + 1L,] <- Y[i,]+Ay(t[i],X[i,],Y[i,],Z[i,])*Dt+Sy(t[i],X[i,],Y[i,],Z[i,])*Wuy[i,]
	ZZ[i + 1L,] <- Z[i,]+Az(t[i],X[i,],Y[i,],Z[i,])*Dt+Sz(t[i],X[i,],Y[i,],Z[i,])*Wuz[i,]
	X[i + 1L,] <- X[i,] + 0.5*(Ax(t[i],X[i,],Y[i,],Z[i,])+Ax(t[i],XX[i,],Y[i,],Z[i,])) * Dt +0.5* (Sx(t[i],X[i,],Y[i,],Z[i,])+Sx(t[i],XX[i,],Y[i,],Z[i,])) * Wdx[i,]
    Y[i + 1L,] <- Y[i,] + 0.5*(Ay(t[i],X[i,],Y[i,],Z[i,])+Ay(t[i],X[i,],YY[i,],Z[i,])) * Dt +0.5* (Sy(t[i],X[i,],Y[i,],Z[i,])+Sy(t[i],X[i,],YY[i,],Z[i,])) * Wdy[i,]
    Z[i + 1L,] <- Z[i,] + 0.5*(Az(t[i],X[i,],Y[i,],Z[i,])+Az(t[i],X[i,],Y[i,],ZZ[i,])) * Dt +0.5* (Sz(t[i],X[i,],Y[i,],Z[i,])+Sy(t[i],X[i,],Y[i,],ZZ[i,])) * Wdz[i,]	
	}	
    name <- c("X","Y","Z")
    name <- if(M > 1) c(paste(name[1],1:M,sep=""),paste(name[2],1:M,sep=""),paste(name[3],1:M,sep=""))
    X <- ts(X, start = t0, deltat = Dt, names=name[1:M])
    Y <- ts(Y, start = t0, deltat = Dt, names=name[(M+1):(2*M)])
    Z <- ts(Z, start = t0, deltat = Dt, names=name[(2*M+1):(3*M)])
    return(list(X=X,Y=Y,Z=Z))
} 

