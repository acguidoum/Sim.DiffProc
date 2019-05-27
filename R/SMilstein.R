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
##### SMilstein1D
 
.SMilstein1D <- function(N =1000,M=1,x0=0,t0=0,T=1,Dt=NULL,drift,diffusion,
                          type=c("ito","str"),...)
                       {
	DAx   <- Simplify(Deriv(drift,"x",cache.exp=FALSE))
    DAxx  <- Simplify(Deriv(drift,"x",nderiv=2,cache.exp=FALSE))	
    DSx   <- Simplify(Deriv(diffusion,"x",cache.exp=FALSE))  
    DSxx  <- Simplify(Deriv(diffusion,"x",nderiv=2,cache.exp=FALSE))
    DSxxx <- Simplify(Deriv(diffusion,"x",nderiv=3,cache.exp=FALSE))					     					   	
    if (type=="ito"){
    A    <- function(t,x)  eval(drift)
    Ax   <- function(t,x)  eval(DAx)
    Axx  <- function(t,x)  eval(DAxx)
    }else{
    A    <- function(t,x)  eval(drift) - 0.5 * eval(diffusion) * eval(DSx)
    Ax   <- function(t,x)  eval(DAx) - 0.5 * (eval(DSx) * eval(DSx)+ eval(diffusion) * eval(DSxx))
    Axx  <- function(t,x)  eval(DAxx) - 0.5 * ( eval(DSxx) * eval(DSx)+ eval(DSx) * eval(DSxx)+ eval(DSx) * eval(DSxx) + eval(diffusion) * eval(DSxxx) )
                  }
    S    <- function(t,x)  eval(diffusion)
    Sx   <- function(t,x)  eval(DSx)
    Sxx  <- function(t,x)  eval(DSxx)
    if (is.null(Dt)) {
        Dt <- (T - t0)/N
        t <- seq(t0, T, by=Dt)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
		T <- t[N + 1]
    }
    W <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    X <- matrix(x0, N+1, M)
	for (i in 1L:N){
	X[i + 1L,] = X[i,] + A(t[i],X[i,])*Dt + S(t[i],X[i,])*W[i,] +0.5 *S(t[i],X[i,]) *
              Sx(t[i],X[i,])*(W[i,]^2-Dt)+ Dt^(3 /2)*(0.5 *A(t[i],X[i,])*Sx(t[i],X[i,]) + 
              0.5 * Ax(t[i],X[i,])*S(t[i],X[i,])+ 0.25 *(S(t[i] ,X[i,])^2) *Sxx(t[i] ,X[i,]))*
              W[i,]+ (Dt^2) * (0.5*A(t[i],X[i,])*Ax(t[i],X[i,])+ 0.25 *Axx(t[i],X[i,])*(S(t[i],X[i,])^2))
			  }
    name <- "X"
    name <- if(M > 1) paste("X",1:M,sep="")
    X <- ts(X, start = t0, deltat = Dt, names=name)
    return(list(X=X))
}   
      
        

#####
##### SMilstein2D

.SMilstein2D <- function(N =1000,M=1,x0=0,y0=0,t0=0,T=1,Dt=NULL,driftx,diffx,drifty,diffy,
                          type=c("ito","str"),...)
                       {
	DAx   <- Simplify(Deriv(driftx,"x",cache.exp=FALSE))
    DAxx  <- Simplify(Deriv(driftx,"x",nderiv=2,cache.exp=FALSE))	
    DSx   <- Simplify(Deriv(diffx,"x",cache.exp=FALSE))  
    DSxx  <- Simplify(Deriv(diffx,"x",nderiv=2,cache.exp=FALSE))
    DSxxx <- Simplify(Deriv(diffx,"x",nderiv=3,cache.exp=FALSE))	
	DAy   <- Simplify(Deriv(drifty,"y",cache.exp=FALSE))
    DAyy  <- Simplify(Deriv(drifty,"y",nderiv=2,cache.exp=FALSE))	
    DSy   <- Simplify(Deriv(diffy,"y",cache.exp=FALSE))  
    DSyy  <- Simplify(Deriv(diffy,"y",nderiv=2,cache.exp=FALSE))
    DSyyy <- Simplify(Deriv(diffy,"y",nderiv=3,cache.exp=FALSE))					   
    if (type=="ito"){
    Ax    <- function(t,x,y)  eval(driftx)
    dAx   <- function(t,x,y)  eval(DAx)
    dAxx  <- function(t,x,y)  eval(DAxx)
    Ay    <- function(t,x,y)  eval(drifty)
    dAy   <- function(t,x,y)  eval(DAy)
    dAyy  <- function(t,x,y)  eval(DAyy)}else{
    Ax    <- function(t,x,y)  eval(driftx) - 0.5 * eval(diffx) * eval(DSx)
    dAx   <- function(t,x,y)  eval(DAx) - 0.5 * (eval(DSx) * eval(DSx)+ eval(diffx) * eval(DSxx))
    dAxx  <- function(t,x,y)  eval(DAxx) - 0.5 * ( eval(DSxx) * eval(DSx)+ eval(DSx) * eval(DSxx)+ eval(DSx) * eval(DSxx) + eval(diffx) * eval(DSxxx) )
    Ay    <- function(t,x,y)  eval(drifty) - 0.5 * eval(diffy) * eval(DSy)
    dAy   <- function(t,x,y)  eval(DAy) - 0.5 * (eval(DSy) * eval(DSy)+ eval(diffy) * eval(DSyy))
    dAyy  <- function(t,x,y)  eval(DAyy) - 0.5 * ( eval(DSyy) * eval(DSy)+ eval(DSy) * eval(DSyy)+ eval(DSy) * eval(DSyy) + eval(diffy) * eval(DSyyy) )
                  }
    Sx    <- function(t,x,y)  eval(diffx)
    dSx   <- function(t,x,y)  eval(DSx)
    dSxx  <- function(t,x,y)  eval(DSxx)
    Sy    <- function(t,x,y)  eval(diffy)
    dSy   <- function(t,x,y)  eval(DSy)
    dSyy  <- function(t,x,y)  eval(DSyy)
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
	for (i in 1L:N){
	X[i + 1L,] = X[i,] + Ax(t[i],X[i,],Y[i,])*Dt + Sx(t[i],X[i,],Y[i,])*W1[i,] +0.5 *Sx(t[i],X[i,],Y[i,]) *
              dSx(t[i],X[i,],Y[i,])*(W1[i,]^2-Dt)+ Dt^(3 /2)*(0.5 *Ax(t[i],X[i,],Y[i,])*dSx(t[i],X[i,],Y[i,]) + 
              0.5 * dAx(t[i],X[i,],Y[i,])*Sx(t[i],X[i,],Y[i,])+ 0.25 *(Sx(t[i] ,X[i,],Y[i,])^2) *dSxx(t[i] ,X[i,],Y[i,]))*
              W1[i,]+ (Dt^2) * (0.5*Ax(t[i],X[i,],Y[i,])*dAx(t[i],X[i,],Y[i,])+ 0.25 *dAxx(t[i],X[i,],Y[i,])*(Sx(t[i],X[i,],Y[i,])^2))
	Y[i + 1L,] = Y[i,] + Ay(t[i],X[i,],Y[i,])*Dt + Sy(t[i],X[i,],Y[i,])*W2[i,] +0.5 *Sy(t[i],X[i,],Y[i,]) *
              dSy(t[i],X[i,],Y[i,])*(W2[i,]^2-Dt)+ Dt^(3 /2)*(0.5 *Ay(t[i],X[i,],Y[i,])*dSy(t[i],X[i,],Y[i,]) + 
              0.5 * dAy(t[i],X[i,],Y[i,])*Sy(t[i],X[i,],Y[i,])+ 0.25 *(Sy(t[i] ,X[i,],Y[i,])^2) *dSyy(t[i] ,X[i,],Y[i,]))*
              W2[i,]+ (Dt^2) * (0.5*Ay(t[i],X[i,],Y[i,])*dAy(t[i],X[i,],Y[i,])+ 0.25 *dAyy(t[i],X[i,],Y[i,])*(Sy(t[i],X[i,],Y[i,])^2))
			  }
    name <- c("X","Y")
    name <- if(M > 1) c(paste(name[1],1:M,sep=""),paste(name[2],1:M,sep=""))
    X <- ts(X, start = t0, deltat = Dt, names=name[1:M])
    Y <- ts(Y, start = t0, deltat = Dt, names=name[(M+1):(2*M)])
    return(list(X=X,Y=Y))
} 

#####
##### SMilstein3D

.SMilstein3D <- function(N =1000,M=1,x0=0,y0=0,z0=0,t0=0,T=1,Dt=NULL,driftx,diffx,drifty,diffy,
                     driftz,diffz,type=c("ito","str"),...)
                       {
	DAx   <- Simplify(Deriv(driftx,"x",cache.exp=FALSE))
    DAxx  <- Simplify(Deriv(driftx,"x",nderiv=2,cache.exp=FALSE))	
    DSx   <- Simplify(Deriv(diffx,"x",cache.exp=FALSE))  
    DSxx  <- Simplify(Deriv(diffx,"x",nderiv=2,cache.exp=FALSE))
    DSxxx <- Simplify(Deriv(diffx,"x",nderiv=3,cache.exp=FALSE))	
	DAy   <- Simplify(Deriv(drifty,"y",cache.exp=FALSE))
    DAyy  <- Simplify(Deriv(drifty,"y",nderiv=2,cache.exp=FALSE))	
    DSy   <- Simplify(Deriv(diffy,"y",cache.exp=FALSE))  
    DSyy  <- Simplify(Deriv(diffy,"y",nderiv=2,cache.exp=FALSE))
    DSyyy <- Simplify(Deriv(diffy,"y",nderiv=3,cache.exp=FALSE))
	DAz   <- Simplify(Deriv(driftz,"z",cache.exp=FALSE))
    DAzz  <- Simplify(Deriv(driftz,"z",nderiv=2,cache.exp=FALSE))	
    DSz   <- Simplify(Deriv(diffz,"z",cache.exp=FALSE))  
    DSzz  <- Simplify(Deriv(diffz,"z",nderiv=2,cache.exp=FALSE))
    DSzzz <- Simplify(Deriv(diffz,"z",nderiv=3,cache.exp=FALSE))
    if (type=="ito"){
    Ax    <- function(t,x,y,z)  eval(driftx)
    dAx   <- function(t,x,y,z)  eval(DAx)
    dAxx  <- function(t,x,y,z)  eval(DAxx)
    Ay    <- function(t,x,y,z)  eval(drifty)
    dAy   <- function(t,x,y,z)  eval(DAy)
    dAyy  <- function(t,x,y,z)  eval(DAyy)
    Az    <- function(t,x,y,z)  eval(driftz)
    dAz   <- function(t,x,y,z)  eval(DAz)
    dAzz  <- function(t,x,y,z)  eval(DAzz)}else{
    Ax    <- function(t,x,y,z)  eval(driftx) - 0.5 * eval(diffx) * eval(DSx)
    dAx   <- function(t,x,y,z)  eval(DAx) - 0.5 * (eval(DSx) * eval(DSx)+ eval(diffx) * eval(DSxx))
    dAxx  <- function(t,x,y,z)  eval(DAxx) - 0.5 * ( eval(DSxx) * eval(DSx)+ eval(DSx) * eval(DSxx)+ eval(DSx) * eval(DSxx) + eval(diffx) * eval(DSxxx) )
    Ay    <- function(t,x,y,z)  eval(drifty) - 0.5 * eval(diffy) * eval(DSy)
    dAy   <- function(t,x,y,z)  eval(DAy) - 0.5 * (eval(DSy) * eval(DSy)+ eval(diffy) * eval(DSyy))
    dAyy  <- function(t,x,y,z)  eval(DAyy) - 0.5 * ( eval(DSyy) * eval(DSy)+ eval(DSy) * eval(DSyy)+ eval(DSy) * eval(DSyy) + eval(diffy) * eval(DSyyy) )
    Az    <- function(t,x,y,z)  eval(driftz) - 0.5 * eval(diffz) * eval(DSz)
    dAz   <- function(t,x,y,z)  eval(DAz) - 0.5 * (eval(DSz) * eval(DSz)+ eval(diffz) * eval(DSzz))
    dAzz  <- function(t,x,y,z)  eval(DAzz) - 0.5 * ( eval(DSzz) * eval(DSz)+ eval(DSz) * eval(DSzz)+ eval(DSz) * eval(DSzz) + eval(diffz) * eval(DSzzz) )
                  }
    Sx    <- function(t,x,y,z)  eval(diffx)
    dSx   <- function(t,x,y,z)  eval(DSx)
    dSxx  <- function(t,x,y,z)  eval(DSxx)
    Sy    <- function(t,x,y,z)  eval(diffy)
    dSy   <- function(t,x,y,z)  eval(DSy)
    dSyy  <- function(t,x,y,z)  eval(DSyy)
    Sz    <- function(t,x,y,z)  eval(diffz)
    dSz   <- function(t,x,y,z)  eval(DSz)
    dSzz  <- function(t,x,y,z)  eval(DSzz)
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
	for (i in 1L:N){
	X[i + 1L,] = X[i,] + Ax(t[i],X[i,],Y[i,],Z[i,])*Dt + Sx(t[i],X[i,],Y[i,],Z[i,])*W1[i,] +0.5 *Sx(t[i],X[i,],Y[i,],Z[i,]) *
              dSx(t[i],X[i,],Y[i,],Z[i,])*(W1[i,]^2-Dt)+ Dt^(3 /2)*(0.5 *Ax(t[i],X[i,],Y[i,],Z[i,])*dSx(t[i],X[i,],Y[i,],Z[i,]) + 
              0.5 * dAx(t[i],X[i,],Y[i,],Z[i,])*Sx(t[i],X[i,],Y[i,],Z[i,])+ 0.25 *(Sx(t[i] ,X[i,],Y[i,],Z[i,])^2) *dSxx(t[i] ,X[i,],Y[i,],Z[i,]))*
              W1[i,]+ (Dt^2) * (0.5*Ax(t[i],X[i,],Y[i,],Z[i,])*dAx(t[i],X[i,],Y[i,],Z[i,])+ 0.25 *dAxx(t[i],X[i,],Y[i,],Z[i,])*(Sx(t[i],X[i,],Y[i,],Z[i,])^2))
	Y[i + 1L,] = Y[i,] + Ay(t[i],X[i,],Y[i,],Z[i,])*Dt + Sy(t[i],X[i,],Y[i,],Z[i,])*W2[i,] +0.5 *Sy(t[i],X[i,],Y[i,],Z[i,]) *
              dSy(t[i],X[i,],Y[i,],Z[i,])*(W2[i,]^2-Dt)+ Dt^(3 /2)*(0.5 *Ay(t[i],X[i,],Y[i,],Z[i,])*dSy(t[i],X[i,],Y[i,],Z[i,]) + 
              0.5 * dAy(t[i],X[i,],Y[i,],Z[i,])*Sy(t[i],X[i,],Y[i,],Z[i,])+ 0.25 *(Sy(t[i] ,X[i,],Y[i,],Z[i,])^2) *dSyy(t[i] ,X[i,],Y[i,],Z[i,]))*
              W2[i,]+ (Dt^2) * (0.5*Ay(t[i],X[i,],Y[i,],Z[i,])*dAy(t[i],X[i,],Y[i,],Z[i,])+ 0.25 *dAyy(t[i],X[i,],Y[i,],Z[i,])*(Sy(t[i],X[i,],Y[i,],Z[i,])^2))
	Z[i + 1L,] = Z[i,] + Az(t[i],X[i,],Y[i,],Z[i,])*Dt + Sz(t[i],X[i,],Y[i,],Z[i,])*W3[i,] +0.5 *Sz(t[i],X[i,],Y[i,],Z[i,]) *
              dSz(t[i],X[i,],Y[i,],Z[i,])*(W3[i,]^2-Dt)+ Dt^(3 /2)*(0.5 *Az(t[i],X[i,],Y[i,],Z[i,])*dSz(t[i],X[i,],Y[i,],Z[i,]) + 
              0.5 * dAz(t[i],X[i,],Y[i,],Z[i,])*Sz(t[i],X[i,],Y[i,],Z[i,])+ 0.25 *(Sz(t[i] ,X[i,],Y[i,],Z[i,])^2) *dSzz(t[i] ,X[i,],Y[i,],Z[i,]))*
              W3[i,]+ (Dt^2) * (0.5*Az(t[i],X[i,],Y[i,],Z[i,])*dAz(t[i],X[i,],Y[i,],Z[i,])+ 0.25 *dAzz(t[i],X[i,],Y[i,],Z[i,])*(Sz(t[i],X[i,],Y[i,],Z[i,])^2))
			  }
    name <- c("X","Y","Z")
    name <- if(M > 1) c(paste(name[1],1:M,sep=""),paste(name[2],1:M,sep=""),paste(name[3],1:M,sep=""))
    X <- ts(X, start = t0, deltat = Dt, names=name[1:M])
    Y <- ts(Y, start = t0, deltat = Dt, names=name[(M+1):(2*M)])
    Z <- ts(Z, start = t0, deltat = Dt, names=name[(2*M+1):(3*M)])
    return(list(X=X,Y=Y,Z=Z))
}




