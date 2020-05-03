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
##### Euler1D

.Euler1D <- function(N =1000,M=1,x0=0,t0=0,T=1,Dt,drift,diffusion,
                          type=c("ito","str"),...)
                       {					   
    if (type=="ito") {A    <- function(t,x)  eval(drift)}else{
	driftstr <- eval(Simplify(substitute(expression(e1 + 0.5 * e2 * de2), list(e1 = drift[[1]], e2 = diffusion[[1]], de2 = Deriv(diffusion,"x",cache.exp=FALSE)[[1]]))))
	A  <- function(t,x) eval(driftstr)
	}
    S  <- function(t,x) eval(diffusion)	
    if (missing(Dt)) {
        t <- seq(t0, T, by=Dt)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
		T <- t[N + 1]
    }
    Dt <- (T - t0)/N	
	W  <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    X <- matrix(x0, N+1, M)
    for (i in 1L:N) {X[i + 1L,] <- X[i,] + A(t[i],X[i,]) * Dt + S(t[i],X[i,]) * W[i,]  }
    name <- "X"
    name <- if(M > 1) paste("X",1:M,sep="")
    X <- ts(X, start = t0, deltat = Dt, names=name)
    return(list(X=X))
}   

#####
##### Euler2D

.Euler2D <- function(N =1000,M=1,x0=0,y0=0,t0=0,T=1,Dt,driftx,diffx,drifty,diffy,
                          type=c("ito","str"),corr=NULL,...)
   {				   
    if (type=="ito"){
          Ax <- function(t,x,y) eval(driftx)
          Ay <- function(t,x,y) eval(drifty)}else{
		     if (is.null(corr)){
		       driftstrx <- eval(Simplify(substitute(expression(e1 + 0.5 * e2 * de2), list(e1 = driftx[[1]], e2 = diffx[[1]], de2 = Deriv(diffx,"x",cache.exp=FALSE)[[1]]))))
		       driftstry <- eval(Simplify(substitute(expression(e1 + 0.5 * e2 * de2), list(e1 = drifty[[1]], e2 = diffy[[1]], de2 = Deriv(diffy,"y",cache.exp=FALSE)[[1]]))))
		     }else{
			   C <- stats::cov2cor(corr)
		       driftstrx <- eval(Simplify(substitute(expression(e1 + 0.5 * e2 * de2 + 0.5 * rho * e3 * de3), list(e1 = driftx[[1]], e2 = diffx[[1]], de2 = Deriv(diffx,"x",cache.exp=FALSE)[[1]],
			                     e3 = diffy[[1]], de3 = Deriv(diffx,"y",cache.exp=FALSE)[[1]], rho = C[1,2]))))
		       driftstry <- eval(Simplify(substitute(expression(e1 + 0.5 * e2 * de2 + 0.5 * rho * e3 * de3), list(e1 = drifty[[1]], e2 = diffy[[1]], de2 = Deriv(diffy,"y",cache.exp=FALSE)[[1]],
                                 e3 = diffx[[1]], de3 = Deriv(diffy,"x",cache.exp=FALSE)[[1]], rho = C[1,2]))))			   
		     }
          Ax <- function(t,x,y) eval(driftstrx) 
          Ay <- function(t,x,y) eval(driftstry)
         }
    Sx <- function(t,x,y) eval(diffx)
    Sy <- function(t,x,y) eval(diffy) 
    if (missing(Dt)) {
        t <- seq(t0, T, by=Dt)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
		T <- t[N + 1]
    }
    Dt <- (T - t0)/N
	Z1 <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    Z2 <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    if (!is.null(corr)){
        C <- stats::cov2cor(corr)
	    W1 <- Z1
        W2 <- C[1,2] * Z1 + sqrt(1-C[1,2]^2) * Z2	
     }else{
	    W1 <- Z1
		W2 <- Z2
	 }	
    X <- matrix(x0, N+1, M)
    Y <- matrix(y0, N+1, M)	
    for (i in 1L:N) {
       X[i + 1L,] <- X[i,] + Ax(t[i],X[i,],Y[i,]) * Dt + Sx(t[i],X[i,],Y[i,]) * W1[i,]  
       Y[i + 1L,] <- Y[i,] + Ay(t[i],X[i,],Y[i,]) * Dt + Sy(t[i],X[i,],Y[i,]) * W2[i,]}
    name <- c("X","Y")
    name <- if(M > 1) c(paste(name[1],1:M,sep=""),paste(name[2],1:M,sep=""))
    X <- ts(X, start = t0, deltat = Dt, names=name[1:M])
    Y <- ts(Y, start = t0, deltat = Dt, names=name[(M+1):(2*M)])
    return(list(X=X,Y=Y))
} 

#####
##### Euler3D


.Euler3D <- function(N =1000,M=1,x0=0,y0=0,z0=0,t0=0,T=1,Dt,driftx,diffx,drifty,diffy,
                     driftz,diffz,type=c("ito","str"),corr=NULL,...)
                       {
    if (type=="ito"){
        Ax    <- function(t,x,y,z) eval(driftx)
        Ay    <- function(t,x,y,z) eval(drifty)
        Az    <- function(t,x,y,z) eval(driftz)}else{
		   	if (is.null(corr)){
		      driftstrx <- eval(Simplify(substitute(expression(e1 + 0.5 * e2 * de2), list(e1 = driftx[[1]], e2 = diffx[[1]], de2 = Deriv(diffx,"x",cache.exp=FALSE)[[1]]))))
		      driftstry <- eval(Simplify(substitute(expression(e1 + 0.5 * e2 * de2), list(e1 = drifty[[1]], e2 = diffy[[1]], de2 = Deriv(diffy,"y",cache.exp=FALSE)[[1]]))))
	          driftstrz <- eval(Simplify(substitute(expression(e1 + 0.5 * e2 * de2), list(e1 = driftz[[1]], e2 = diffz[[1]], de2 = Deriv(diffz,"z",cache.exp=FALSE)[[1]]))))
			}else{
			  C <- stats::cov2cor(corr)
		      driftstrx <- eval(Simplify(substitute(expression(e1 + 0.5 * e2 * de2 + 0.5 * rho1 *e3 * de3 + 0.5 * rho2 *e4 * de4), list(e1 = driftx[[1]], e2 = diffx[[1]], de2 = Deriv(diffx,"x",cache.exp=FALSE)[[1]],
			   e3 = diffy[[1]], de3 = Deriv(diffx,"y",cache.exp=FALSE)[[1]],e4 = diffz[[1]], de4 = Deriv(diffx,"z",cache.exp=FALSE)[[1]],roh1=C[1,2],rho2=C[1,3]))))
		      driftstry <- eval(Simplify(substitute(expression(e1 + 0.5 * e2 * de2 + 0.5 * rho1 *e3 * de3 + 0.5 * rho3 *e4 * de4), list(e1 = drifty[[1]], e2 = diffy[[1]], de2 = Deriv(diffy,"y",cache.exp=FALSE)[[1]],
			   e3 = diffx[[1]], de3 = Deriv(diffy,"x",cache.exp=FALSE)[[1]],e4 = diffz[[1]], de4 = Deriv(diffy,"z",cache.exp=FALSE)[[1]],roh1=C[1,2],rho3=C[2,3]))))			  
	          driftstrz <- eval(Simplify(substitute(expression(e1 + 0.5 * e2 * de2 + 0.5 * rho2 *e3 * de3 + 0.5 * rho3 *e4 * de4), list(e1 = driftz[[1]], e2 = diffz[[1]], de2 = Deriv(diffz,"z",cache.exp=FALSE)[[1]],
			   e3 = diffx[[1]], de3 = Deriv(diffz,"x",cache.exp=FALSE)[[1]],e4 = diffy[[1]], de4 = Deriv(diffz,"y",cache.exp=FALSE)[[1]],roh2=C[1,3],rho3=C[2,3]))))			  
			}
        Ax <- function(t,x,y,z) eval(driftstrx) 
        Ay <- function(t,x,y,z) eval(driftstry)
        Az <- function(t,x,y,z) eval(driftstrz)
	}
    Sx <- function(t,x,y,z) eval(diffx)
    Sy <- function(t,x,y,z) eval(diffy) 
    Sz <- function(t,x,y,z) eval(diffz)
    if (missing(Dt)) {
        t <- seq(t0, T, by=Dt)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
		T <- t[N + 1]
    } 
    Dt <- (T - t0)/N
	Z1 <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    Z2 <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    Z3 <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    if (!is.null(corr)){
        C <- stats::cov2cor(corr)
	    W1 <- Z1
        W2 <- C[1,2] * Z1 + sqrt(1-C[1,2]^2) * Z2	
		W3 <- C[1,3] * Z1 + ((C[2,3]-C[1,2]*C[1,3])/(sqrt(1-C[1,2]^2))) * Z2 + ( sqrt(1-C[1,3]^2 - ((C[2,3]-C[1,2]*C[1,3])/(sqrt(1-C[1,2]^2)))^2 )  ) * Z3
     }else{
	    W1 <- Z1
		W2 <- Z2
		W3 <- Z3
	 }	
    X <- matrix(x0, N+1, M)
    Y <- matrix(y0, N+1, M)	
	Z <- matrix(z0, N+1, M)	
    for (i in 1L:N) {
       X[i + 1L,] <- X[i,] + Ax(t[i],X[i,],Y[i,],Z[i,]) * Dt + Sx(t[i],X[i,],Y[i,],Z[i,]) * W1[i,]  
       Y[i + 1L,] <- Y[i,] + Ay(t[i],X[i,],Y[i,],Z[i,]) * Dt + Sy(t[i],X[i,],Y[i,],Z[i,]) * W2[i,]
	   Z[i + 1L,] <- Z[i,] + Az(t[i],X[i,],Y[i,],Z[i,]) * Dt + Sz(t[i],X[i,],Y[i,],Z[i,]) * W3[i,]}
    name <- c("X","Y","Z")
    name <- if(M > 1) c(paste(name[1],1:M,sep=""),paste(name[2],1:M,sep=""),paste(name[3],1:M,sep=""))
    X <- ts(X, start = t0, deltat = Dt, names=name[1:M])
    Y <- ts(Y, start = t0, deltat = Dt, names=name[(M+1):(2*M)])
    Z <- ts(Z, start = t0, deltat = Dt, names=name[(2*M+1):(3*M)])
    return(list(X=X,Y=Y,Z=Z))
} 
