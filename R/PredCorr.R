## Thu Mar 30 21:48:29 2017
## Original file Copyright © 2017 A.C. Guidoum, K. Boukhetala
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
##### PredCorr1D

.PredCorr1D <- function(N =1000,M=1,x0=0,t0=0,T=1,Dt=NULL,alpha=0.5,mu=0.5,
                               drift,diffusion,type=c("ito","str"),...)
                       {
    DSx  <- Simplify(Deriv(diffusion,"x",cache.exp=FALSE)) 
    if (type=="ito"){A    <- function(t,x)  eval(drift)}else{
	driftstr <- eval(Simplify(substitute(expression(e1 - 0.5 * e2 * de2), list(e1 = drift[[1]], e2 = diffusion[[1]], de2 = Deriv(diffusion,"x",cache.exp=FALSE)[[1]]))))
	A  <- function(t,x) eval(driftstr)
    # A <- function(t,x) eval(drift) - 0.5 * eval(diffusion) * eval(DSx)
	}
    S    <- function(t,x)  eval(diffusion)
    Sx   <- function(t,x)  eval(DSx)
    SS   <- function(t,x)  A(t,x) - mu * S(t,x) * Sx(t,x)
    if (is.null(Dt)) {
        Dt <- (T - t0)/N
        t <- seq(t0, T, by=Dt)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
		T <- t[N + 1]
    }
    Wu <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
	Wd <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    X = XX <- matrix(x0, N+1, M)	
    for (i in 1L:N) {
	XX[i + 1L,] = XX[i,] + A(t[i],XX[i,])*Dt + S(t[i],XX[i,])*Wu[i,]
	X[i + 1L,] <- X[i,] + (alpha*SS(t[i+1],XX[i+1,])+(1-alpha)*SS(t[i],X[i,]))*Dt+(mu*S(t[i+1],XX[i+1,])+(1-mu)*S(t[i],X[i,]))*Wd[i,] 
	}	   
    name <- "X"
    name <- if(M > 1) paste("X",1:M,sep="")
    X <- ts(X, start = t0, deltat = Dt, names=name)
    return(list(X=X))
}   
      
  


#####
##### PredCorr2D

.PredCorr2D <- function(N =1000,M=1,x0=0,y0=0,t0=0,T=1,Dt=NULL,alpha=0.5,mu=0.5,driftx,
                        diffx,drifty,diffy,type=c("ito","str"),...)
                       {
    DSx  <- Simplify(Deriv(diffx,"x",cache.exp=FALSE))
    DSy  <- Simplify(Deriv(diffy,"y",cache.exp=FALSE))
    if (type=="ito"){
    Ax    <- function(t,x,y)  eval(driftx)
    Ay    <- function(t,x,y)  eval(drifty)
    }else{
	driftstrx <- eval(Simplify(substitute(expression(e1 - 0.5 * e2 * de2), list(e1 = driftx[[1]], e2 = diffx[[1]], de2 = Deriv(diffx,"x",cache.exp=FALSE)[[1]]))))
	driftstry <- eval(Simplify(substitute(expression(e1 - 0.5 * e2 * de2), list(e1 = drifty[[1]], e2 = diffy[[1]], de2 = Deriv(diffy,"y",cache.exp=FALSE)[[1]]))))
    Ax <- function(t,x,y) eval(driftstrx) 
    Ay <- function(t,x,y) eval(driftstry)
    # Ax <- function(t,x,y) eval(driftx) - 0.5 * eval(diffx) * eval(DSx)
    # Ay <- function(t,x,y) eval(drifty) - 0.5 * eval(diffy) * eval(DSy)
         }
    Sx    <- function(t,x,y)  eval(diffx)
    dSx   <- function(t,x,y)  eval(DSx)
    SSx   <- function(t,x,y)  Ax(t,x,y) - mu * Sx(t,x,y) * dSx(t,x,y)
    Sy    <- function(t,x,y)  eval(diffy)
    dSy   <- function(t,x,y)  eval(DSy)
    SSy   <- function(t,x,y)  Ay(t,x,y) - mu * Sy(t,x,y) * dSy(t,x,y)	
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
	XX[i + 1L,] = XX[i,] + Ax(t[i],XX[i,],Y[i,])*Dt + Sx(t[i],XX[i,],Y[i,])*Wux[i,]
	YY[i + 1L,] = YY[i,] + Ay(t[i],X[i,],YY[i,])*Dt + Sy(t[i],X[i,],YY[i,])*Wuy[i,]
	X[i + 1L,] <- X[i,] + (alpha*SSx(t[i+1],XX[i+1,],Y[i,])+(1-alpha)*SSx(t[i],X[i,],Y[i,]))*Dt+(mu*Sx(t[i+1],XX[i+1,],Y[i,])+(1-mu)*Sx(t[i],X[i,],Y[i,]))*Wdx[i,]
	Y[i + 1L,] <- Y[i,] + (alpha*SSy(t[i+1],X[i,],YY[i+1,])+(1-alpha)*SSy(t[i],X[i,],Y[i,]))*Dt+(mu*Sy(t[i+1],X[i,],YY[i+1,])+(1-mu)*Sy(t[i],X[i,],Y[i,]))*Wdy[i,]	
	}
    name <- c("X","Y")
    name <- if(M > 1) c(paste(name[1],1:M,sep=""),paste(name[2],1:M,sep=""))
    X <- ts(X, start = t0, deltat = Dt, names=name[1:M])
    Y <- ts(Y, start = t0, deltat = Dt, names=name[(M+1):(2*M)])
    return(list(X=X,Y=Y))
} 

#####
##### PredCorr3D

.PredCorr3D <- function(N =1000,M=1,x0=0,y0=0,z0=0,t0=0,T=1,Dt=NULL,alpha=0.5,mu=0.5,driftx,diffx,drifty,diffy,
                     driftz,diffz,type=c("ito","str"),...)
                       {
    DSx  <- Simplify(Deriv(diffx,"x",cache.exp=FALSE))
    DSy  <- Simplify(Deriv(diffy,"y",cache.exp=FALSE))
    DSz  <- Simplify(Deriv(diffz,"z",cache.exp=FALSE)) 
    if (type=="ito"){
    Ax    <- function(t,x,y,z)  eval(driftx)
    Ay    <- function(t,x,y,z)  eval(drifty)
    Az    <- function(t,x,y,z)  eval(driftz)
    }else{
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
    Sx    <- function(t,x,y,z)  eval(diffx)
    dSx   <- function(t,x,y,z)  eval(DSx)
    SSx   <- function(t,x,y,z)  Ax(t,x,y,z) - mu * Sx(t,x,y,z) * dSx(t,x,y,z)
    Sy    <- function(t,x,y,z)  eval(diffy)
    dSy   <- function(t,x,y,z)  eval(DSy)
    SSy   <- function(t,x,y,z)  Ay(t,x,y,z) - mu * Sy(t,x,y,z) * dSy(t,x,y,z)
    Sz    <- function(t,x,y,z)  eval(diffz)
    dSz   <- function(t,x,y,z)  eval(DSz)
    SSz   <- function(t,x,y,z)  Az(t,x,y,z) - mu * Sz(t,x,y,z) * dSz(t,x,y,z)
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
	XX[i + 1L,] = XX[i,] + Ax(t[i],XX[i,],Y[i,],Z[i,])*Dt + Sx(t[i],XX[i,],Y[i,],Z[i,])*Wux[i,]
	YY[i + 1L,] = YY[i,] + Ay(t[i],X[i,],YY[i,],Z[i,])*Dt + Sy(t[i],X[i,],YY[i,],Z[i,])*Wuy[i,]
	ZZ[i + 1L,] = ZZ[i,] + Az(t[i],X[i,],Y[i,],ZZ[i,])*Dt + Sz(t[i],X[i,],Y[i,],ZZ[i,])*Wuz[i,]
	X[i + 1L,] <- X[i,] + (alpha*SSx(t[i+1],XX[i+1,],Y[i,],Z[i,])+(1-alpha)*SSx(t[i],X[i,],Y[i,],Z[i,]))*Dt+(mu*Sx(t[i+1],XX[i+1,],Y[i,],Z[i,])+(1-mu)*Sx(t[i],X[i,],Y[i,],Z[i,]))*Wdx[i,]
	Y[i + 1L,] <- Y[i,] + (alpha*SSy(t[i+1],X[i,],YY[i+1,],Z[i,])+(1-alpha)*SSy(t[i],X[i,],Y[i,],Z[i,]))*Dt+(mu*Sy(t[i+1],X[i,],YY[i+1,],Z[i,])+(1-mu)*Sy(t[i],X[i,],Y[i,],Z[i,]))*Wdy[i,]
	Z[i + 1L,] <- Z[i,] + (alpha*SSz(t[i+1],X[i,],Y[i,],ZZ[i+1,])+(1-alpha)*SSz(t[i],X[i,],Y[i,],Z[i,]))*Dt+(mu*Sz(t[i+1],X[i,],Y[i,],ZZ[i+1,])+(1-mu)*Sz(t[i],X[i,],Y[i,],Z[i,]))*Wdz[i,]	
	}
	name <- c("X","Y","Z")
    name <- if(M > 1) c(paste(name[1],1:M,sep=""),paste(name[2],1:M,sep=""),paste(name[3],1:M,sep=""))
    X <- ts(X, start = t0, deltat = Dt, names=name[1:M])
    Y <- ts(Y, start = t0, deltat = Dt, names=name[(M+1):(2*M)])
    Z <- ts(Z, start = t0, deltat = Dt, names=name[(2*M+1):(3*M)])
    return(list(X=X,Y=Y,Z=Z))
} 

