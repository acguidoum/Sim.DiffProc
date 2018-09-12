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
##### RK1D
         
.RK1D <- function(N =1000,M=1,x0=0,t0=0,T=1,Dt=NULL,drift,diffusion,
                          type=c("ito","str"),order=c(1,2,3),...)
                       {					   
    if (type=="ito") {A    <- function(t,x)  eval(drift)}else{
	driftstr <- eval(Simplify(substitute(expression(e1 - 0.5 * e2 * de2), list(e1 = drift[[1]], e2 = diffusion[[1]], de2 = Deriv(diffusion,"x",cache.exp=FALSE)[[1]]))))
	A  <- function(t,x) eval(driftstr)
    # A  <- function(t,x) eval(drift) - 0.5 * eval(diffusion) * eval(Deriv(diffusion,"x"))
	}
    S  <- function(t,x) eval(diffusion)
    if (is.null(Dt)) {
        Dt <- (T - t0)/N
        t <- seq(t0, T, by=Dt)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
		T <- t[N + 1]
    }
    W <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
	Wxx <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
	Wxxx <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    X=XX=XXX <- matrix(x0, N+1, M)		
    for (i in 1L:N){
    if (order==1){
      X[i + 1L,] = X[i,]+ A(t[i],X[i,])*Dt+S(t[i],X[i,])*W[i,]+(0.5/sqrt(Dt)) * (S(t[i]+Dt,X[i,]+S(t[i],X[i,])*sqrt(Dt))-S(t[i],X[i,])) * (W[i,]^2 - Dt)}
    else if (order==2){
      X[i + 1L,] = X[i,]+ 0.5 * (A(t[i],X[i,])+A(t[i]+Dt,X[i,]+A(t[i],X[i,])*Dt+S(t[i],X[i,])*W[i,])) *Dt + 0.25* (2*S(t[i],X[i,])+S(t[i]+Dt,X[i,]+A(t[i],X[i,])*Dt+
             sqrt(Dt)*S(t[i],X[i,]))+S(t[i]+Dt,X[i,]+A(t[i],X[i,])*Dt-sqrt(Dt)*S(t[i],X[i,])) ) *W[i,]+0.25 * (S(t[i]+Dt,X[i,]+A(t[i],X[i,])*Dt-sqrt(Dt)*S(t[i],X[i,])) - 
             S(t[i]+Dt,X[i,]+A(t[i],X[i,])*Dt+sqrt(Dt)*S(t[i],X[i,]))) *(sqrt(Dt) - (W[i,])^2 / sqrt(Dt))}
    else if (order==3){	          
              XX[i + 1L,]=X[i,]+0.5*Dt*A(t[i],X[i,])+S(t[i],X[i,])*Wxx[i,]
              XXX[i + 1L,]=X[i,]-A(t[i],X[i,])*Dt+2*Dt*A(t[i]+0.5*Dt,XX[i,])+(2*S(t[i]+0.5*Dt,XX[i,])-S(t[i],X[i,]))*Wxxx[i,]
              X[i + 1L,] = X[i,]+(Dt/6)*(A(t[i],X[i,])+4*A(t[i]+0.5*Dt,XX[i,])+A(t[i]+Dt,XXX[i,]))+(1/6)*(S(t[i],X[i,])+4*S(t[i,]+0.5*Dt,XX[i,])+S(t[i]+Dt,XXX[i,]))*W[i,]}
       }
    name <- "X"
    name <- if(M > 1) paste("X",1:M,sep="")
    X <- ts(X, start = t0, deltat = Dt, names=name)
    return(list(X=X))
}   

  
#####
##### RK2D

.RK2D <- function(N =1000,M=1,x0=0,y0=0,t0=0,T=1,Dt=NULL,driftx,diffx,drifty,diffy,
                          type=c("ito","str"),order=c(1,2,3),...)
                       {
    if (type=="ito"){
    Ax <- function(t,x,y) eval(driftx)
    Ay <- function(t,x,y) eval(drifty) }else{
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
    W1 <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    W2 <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
	Wxx <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
	Wxxx <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
	Wyy <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
	Wyyy <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    X=XX=XXX <- matrix(x0, N+1, M)	
    Y=YY=YYY <- matrix(y0, N+1, M)
    for (i in 1L:N){
    if (order==1){
      X[i + 1L,] = X[i,]+ Ax(t[i],X[i,],Y[i,])*Dt+Sx(t[i],X[i,],Y[i,])*W1[i,]+(0.5/sqrt(Dt)) * (Sx(t[i]+Dt,X[i,]+Sx(t[i],X[i,],Y[i,])*sqrt(Dt),Y[i,])-Sx(t[i],X[i,],Y[i,])) * (W1[i,]^2 - Dt)
	  Y[i + 1L,] = Y[i,]+ Ay(t[i],X[i,],Y[i,])*Dt+Sy(t[i],X[i,],Y[i,])*W2[i,]+(0.5/sqrt(Dt)) * (Sy(t[i]+Dt,X[i,],Y[i,]+Sy(t[i],X[i,],Y[i,])*sqrt(Dt))-Sy(t[i],X[i,],Y[i,])) * (W2[i,]^2 - Dt)}
    else if (order==2){
      X[i + 1L,] = X[i,]+ 0.5 * (Ax(t[i],X[i,],Y[i,])+Ax(t[i]+Dt,X[i,]+Ax(t[i],X[i,],Y[i,])*Dt+Sx(t[i],X[i,],Y[i,])*W1[i,],Y[i,])) *Dt + 0.25* (2*Sx(t[i],X[i,],Y[i,])+Sx(t[i]+Dt,X[i,]+Ax(t[i],X[i,],Y[i,])*Dt+
             sqrt(Dt)*Sx(t[i],X[i,],Y[i,]),Y[i,])+Sx(t[i]+Dt,X[i,]+Ax(t[i],X[i,],Y[i,])*Dt-sqrt(Dt)*Sx(t[i],X[i,],Y[i,]),Y[i,]) ) *W1[i,]+0.25 * (Sx(t[i]+Dt,X[i,]+Ax(t[i],X[i,],Y[i,])*Dt-sqrt(Dt)*Sx(t[i],X[i,],Y[i,]),Y[i,]) - 
             Sx(t[i]+Dt,X[i,]+Ax(t[i],X[i,],Y[i,])*Dt+sqrt(Dt)*Sx(t[i],X[i,],Y[i,]),Y[i,])) *(sqrt(Dt) - (W1[i,])^2 / sqrt(Dt))
	  Y[i + 1L,] = Y[i,]+ 0.5 * (Ay(t[i],X[i,],Y[i,])+Ay(t[i]+Dt,X[i,],Y[i,]+Ay(t[i],X[i,],Y[i,])*Dt+Sy(t[i],X[i,],Y[i,])*W2[i,])) *Dt + 0.25* (2*Sy(t[i],X[i,],Y[i,])+Sy(t[i]+Dt,X[i,],Y[i,]+Ay(t[i],X[i,],Y[i,])*Dt+
	         sqrt(Dt)*Sy(t[i],X[i,],Y[i,]))+Sy(t[i]+Dt,X[i,],Y[i,]+Ay(t[i],X[i,],Y[i,])*Dt-sqrt(Dt)*Sy(t[i,],X[i,],Y[i,])) ) *W2[i,]+0.25 * (Sy(t[i]+Dt,X[i,],Y[i,]+Ay(t[i],X[i,],Y[i,])*Dt-sqrt(Dt)*Sy(t[i],X[i,],Y[i,]))-
			 Sy(t[i]+Dt,X[i,],Y[i,]+Ay(t[i],X[i,],Y[i,])*Dt+sqrt(Dt)*Sy(t[i],X[i,],Y[i,]))) *(sqrt(Dt) - (W2[i,])^2 / sqrt(Dt))}
    else if (order==3){	          
      XX[i + 1L,] =X[i,]+0.5*Dt*Ax(t[i],X[i,],Y[i,])+Sx(t[i],X[i,],Y[i,])*Wxx[i,]
      YY[i + 1L,] =Y[i,]+0.5*Dt*Ay(t[i],X[i,],Y[i,])+Sy(t[i],X[i,],Y[i,])*Wyy[i,]
      XXX[i + 1L,]=X[i,]-Ax(t[i],X[i,],Y[i,])*Dt+2*Dt*Ax(t[i]+0.5*Dt,XX[i,],Y[i,])+(2*Sx(t[i]+0.5*Dt,XX[i,],Y[i,])-Sx(t[i],X[i,],Y[i,]))*Wxxx[i,]
      YYY[i + 1L,]=Y[i,]-Ay(t[i],X[i,],Y[i,])*Dt+2*Dt*Ay(t[i]+0.5*Dt,X[i,],YY[i,])+(2*Sy(t[i]+0.5*Dt,X[i,],YY[i,])-Sy(t[i],X[i,],Y[i,]))*Wyyy[i,]
      X[i + 1L,] = X[i,]+(Dt/6)*(Ax(t[i],X[i,],Y[i,])+4*Ax(t[i]+0.5*Dt,XX[i,],Y[i,])+Ax(t[i]+Dt,XXX[i,],Y[i,]))+(1/6)*(Sx(t[i],X[i,],Y[i,])+4*Sx(t[i]+0.5*Dt,XX[i,],Y[i,])+Sx(t[i]+Dt,XXX[i,],Y[i,]))*W1[i,]
      Y[i + 1L,] = Y[i,]+(Dt/6)*(Ay(t[i],X[i,],Y[i,])+4*Ay(t[i]+0.5*Dt,X[i,],YY[i,])+Ay(t[i]+Dt,X[i,],YYY[i,]))+(1/6)*(Sy(t[i],X[i,],Y[i,])+4*Sy(t[i]+0.5*Dt,X[i,],YY[i,])+Sy(t[i]+Dt,X[i,],YYY[i,]))*W2[i,]
       }}	
    name <- c("X","Y")
    name <- if(M > 1) c(paste(name[1],1:M,sep=""),paste(name[2],1:M,sep=""))
    X <- ts(X, start = t0, deltat = Dt, names=name[1:M])
    Y <- ts(Y, start = t0, deltat = Dt, names=name[(M+1):(2*M)])
    return(list(X=X,Y=Y))
}


#####
##### RK3D

.RK3D <- function(N =1000,M=1,x0=0,y0=0,z0=0,t0=0,T=1,Dt=NULL,driftx,diffx,drifty,diffy,
                     driftz,diffz,type=c("ito","str"),order=c(1,2,3),...)
                       {
    if (type=="ito"){
    Ax <- function(t,x,y,z)  eval(driftx)
    Ay <- function(t,x,y,z)  eval(drifty) 
    Az <- function(t,x,y,z)  eval(driftz)}else{
	driftstrx <- eval(Simplify(substitute(expression(e1 - 0.5 * e2 * de2), list(e1 = driftx[[1]], e2 = diffx[[1]], de2 = Deriv(diffx,"x",cache.exp=FALSE)[[1]]))))
	driftstry <- eval(Simplify(substitute(expression(e1 - 0.5 * e2 * de2), list(e1 = drifty[[1]], e2 = diffy[[1]], de2 = Deriv(diffy,"y",cache.exp=FALSE)[[1]]))))
	driftstrz <- eval(Simplify(substitute(expression(e1 - 0.5 * e2 * de2), list(e1 = driftz[[1]], e2 = diffz[[1]], de2 = Deriv(diffz,"z",cache.exp=FALSE)[[1]]))))
    Ax <- function(t,x,y,z) eval(driftstrx) 
    Ay <- function(t,x,y,z) eval(driftstry)
    Az <- function(t,x,y,z) eval(driftstrz)	
    # Ax <- function(t,x,y,z)  eval(driftx) - 0.5 * eval(diffx) * eval(Deriv(diffx,"x"))
    # Ay <- function(t,x,y,z)  eval(drifty) - 0.5 * eval(diffy) * eval(Deriv(diffy,"y"))
    # Az <- function(t,x,y,z)  eval(driftz) - 0.5 * eval(diffz) * eval(Deriv(diffz,"z"))
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
    W1 <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    W2 <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    W3 <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
	Wxx <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
	Wxxx <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
	Wyy <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
	Wyyy <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
	Wzz <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
	Wzzz <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    X=XX=XXX <- matrix(x0, N+1, M)	
    Y=YY=YYY <- matrix(y0, N+1, M)
    Z=ZZ=ZZZ <- matrix(z0, N+1, M)
    for (i in 1L:N){
    if (order==1){
      X[i + 1L,] = X[i,]+ Ax(t[i],X[i,],Y[i,],Z[i,])*Dt+Sx(t[i],X[i,],Y[i,],Z[i,])*W1[i,]+(0.5/sqrt(Dt)) * (Sx(t[i]+Dt,X[i,]+Sx(t[i],X[i,],Y[i,],Z[i,])*sqrt(Dt),Y[i,],Z[i,])-Sx(t[i],X[i,],Y[i,],Z[i,])) * (W1[i,]^2 - Dt)
	  Y[i + 1L,] = Y[i,]+ Ay(t[i],X[i,],Y[i,],Z[i,])*Dt+Sy(t[i],X[i,],Y[i,],Z[i,])*W2[i,]+(0.5/sqrt(Dt)) * (Sy(t[i]+Dt,X[i,],Y[i,]+Sy(t[i],X[i,],Y[i,],Z[i,])*sqrt(Dt),Z[i,])-Sy(t[i],X[i,],Y[i,],Z[i,])) * (W2[i,]^2 - Dt)			  
	  Z[i + 1L,] = Z[i,]+ Az(t[i],X[i,],Y[i,],Z[i,])*Dt+Sz(t[i],X[i,],Y[i,],Z[i,])*W3[i,]+(0.5/sqrt(Dt)) * (Sz(t[i]+Dt,X[i,],Y[i,],Z[i,]+Sz(t[i],X[i,],Y[i,],Z[i,])*sqrt(Dt))-Sz(t[i],X[i,],Y[i,],Z[i,])) * (W3[i,]^2 - Dt)
	  }
    else if (order==2){
      X[i + 1L,] = X[i,]+ 0.5 * (Ax(t[i],X[i,],Y[i,],Z[i,])+Ax(t[i]+Dt,X[i,]+Ax(t[i],X[i,],Y[i,],Z[i,])*Dt+Sx(t[i],X[i,],Y[i,],Z[i,])*W1[i,],Y[i,],Z[i,])) *Dt + 
			 0.25* (2*Sx(t[i],X[i,],Y[i,],Z[i,])+Sx(t[i]+Dt,X[i,]+Ax(t[i],X[i,],Y[i,],Z[i,])*Dt+sqrt(Dt)*Sx(t[i],X[i,],Y[i,],Z[i,]),Y[i,],Z[i,])+ Sx(t[i]+Dt,X[i,]+Ax(t[i],X[i,],Y[i,],Z[i,])*Dt-sqrt(Dt)*Sx(t[i],X[i,],Y[i,],Z[i,]),Y[i,],Z[i,]) ) *W1[i,]+
             0.25 * (Sx(t[i]+Dt,X[i,]+Ax(t[i],X[i,],Y[i,],Z[i,])*Dt-sqrt(Dt)*Sx(t[i],X[i,],Y[i,],Z[i,]),Y[i,],Z[i,]) - Sx(t[i]+Dt,X[i,]+Ax(t[i],X[i,],Y[i,],Z[i,])*Dt+sqrt(Dt)*Sx(t[i],X[i,],Y[i,],Z[i,]),Y[i,],Z[i,])) *(sqrt(Dt) - (W1[i,])^2 / sqrt(Dt))              
	  Y[i + 1L,] = Y[i,]+ 0.5 * (Ay(t[i],X[i,],Y[i,],Z[i,])+Ay(t[i]+Dt,X[i,],Y[i,]+Ay(t[i],X[i,],Y[i,],Z[i,])*Dt+Sy(t[i],X[i,],Y[i,],Z[i,])*W2[i,],Z[i,])) *Dt + 
			 0.25* (2*Sy(t[i],X[i,],Y[i,],Z[i,])+Sy(t[i]+Dt,X[i,],Y[i,]+Ay(t[i],X[i,],Y[i,],Z[i,])*Dt+sqrt(Dt)*Sy(t[i],X[i,],Y[i,],Z[i,]),Z[i,])+Sy(t[i]+Dt,X[i,],Y[i,]+Ay(t[i],X[i,],Y[i,],Z[i,])*Dt-sqrt(Dt)*Sy(t[i],X[i,],Y[i,],Z[i,]),Z[i,]) ) *W2[i,]+
             0.25 * (Sy(t[i]+Dt,X[i,],Y[i,]+Ay(t[i],X[i,],Y[i,],Z[i,])*Dt-sqrt(Dt)*Sy(t[i],X[i,],Y[i,],Z[i,]),Z[i,])-Sy(t[i]+Dt,X[i,],Y[i,]+Ay(t[i],X[i,],Y[i,],Z[i,])*Dt+sqrt(Dt)*Sy(t[i],X[i,],Y[i,],Z[i,]),Z[i,])) *(sqrt(Dt) - (W2[i,])^2 / sqrt(Dt))              
	  Z[i + 1L,] = Z[i,]+ 0.5 * (Az(t[i],X[i,],Y[i,],Z[i,])+Az(t[i]+Dt,X[i,],Y[i,],Z[i,]+Az(t[i],X[i,],Y[i,],Z[i,])*Dt+Sz(t[i],X[i,],Y[i,],Z[i,])*W3[i,])) *Dt + 
			 0.25* (2*Sz(t[i],X[i,],Y[i,],Z[i,])+Sz(t[i]+Dt,X[i,],Y[i,],Z[i,]+Az(t[i],X[i,],Y[i,],Z[i,])*Dt+sqrt(Dt)*Sz(t[i],X[i,],Y[i,],Z[i,]))+Sz(t[i]+Dt,X[i,],Y[i,],Z[i,]+Az(t[i,],X[i,],Y[i,],Z[i,])*Dt-sqrt(Dt)*Sz(t[i],X[i,],Y[i,],Z[i,])) ) *W3[i,]+
             0.25 * (Sz(t[i]+Dt,X[i,],Y[i,],Z[i,]+Az(t[i],X[i,],Y[i,],Z[i,])*Dt-sqrt(Dt)*Sz(t[i],X[i,],Y[i,],Z[i,]))-Sz(t[i]+Dt,X[i,],Y[i,],Z[i,]+Az(t[i],X[i,],Y[i,],Z[i,])*Dt+sqrt(Dt)*Sz(t[i],X[i,],Y[i,],Z[i,]))) *(sqrt(Dt) - (W3[i,])^2 / sqrt(Dt))  
			  }
    else if (order==3){
      XX[i + 1L,]=X[i,]+0.5*Dt*Ax(t[i],X[i,],Y[i,],Z[i,])+Sx(t[i],X[i,],Y[i,],Z[i,])*Wxx[i,]
      YY[i + 1L,]=Y[i,]+0.5*Dt*Ay(t[i],X[i,],Y[i,],Z[i,])+Sy(t[i],X[i,],Y[i,],Z[i,])*Wyy[i,]
      ZZ[i + 1L,]=Z[i,]+0.5*Dt*Az(t[i],X[i,],Y[i,],Z[i,])+Sz(t[i],X[i,],Y[i,],Z[i,])*Wzz[i,]
      XXX[i + 1L,]=X[i,]-Ax(t[i],X[i,],Y[i,],Z[i,])*Dt+2*Dt*Ax(t[i]+0.5*Dt,XX[i,],Y[i,],Z[i,])+(2*Sx(t[i]+0.5*Dt,XX[i,],Y[i,],Z[i,])-Sx(t[i],X[i,],Y[i,],Z[i,]))*Wxxx[i,]
      YYY[i + 1L,]=Y[i,]-Ay(t[i],X[i,],Y[i,],Z[i,])*Dt+2*Dt*Ay(t[i]+0.5*Dt,X[i,],YY[i,],Z[i,])+(2*Sy(t[i]+0.5*Dt,X[i,],YY[i,],Z[i,])-Sy(t[i],X[i,],Y[i,],Z[i,]))*Wyyy[i,]
      ZZZ[i + 1L,]=Z[i,]-Az(t[i],X[i,],Y[i,],Z[i,])*Dt+2*Dt*Az(t[i]+0.5*Dt,X[i,],Y[i,],ZZ[i,])+(2*Sz(t[i]+0.5*Dt,X[i,],Y[i,],ZZ[i,])-Sz(t[i],X[i,],Y[i,],Z[i,]))*Wzzz[i,]
      X[i + 1L,] = X[i,]+(Dt/6)*(Ax(t[i],X[i,],Y[i,],Z[i,])+4*Ax(t[i]+0.5*Dt,XX[i,],Y[i,],Z[i,])+Ax(t[i]+Dt,XXX[i,],Y[i,],Z[i,]))+
             (1/6)*(Sx(t[i],X[i,],Y[i,],Z[i,])+4*Sx(t[i]+0.5*Dt,XX[i,],Y[i,],Z[i,])+Sx(t[i]+Dt,XXX[i,],Y[i,],Z[i,]))*W1[i,]
      Y[i + 1L,] = Y[i,]+(Dt/6)*(Ay(t[i],X[i,],Y[i,],Z[i,])+4*Ay(t[i]+0.5*Dt,X[i,],YY[i,],Z[i,])+Ay(t[i]+Dt,X[i,],YYY[i,],Z[i,]))+
             (1/6)*(Sy(t[i],X[i,],Y[i,],Z[i,])+4*Sy(t[i]+0.5*Dt,X[i,],YY[i,],Z[i,])+Sy(t[i]+Dt,X[i,],YYY[i,],Z[i,]))*W2[i,]
      Z[i + 1L,] = Z[i,]+(Dt/6)*(Az(t[i],X[i,],Y[i,],Z[i,])+4*Az(t[i]+0.5*Dt,X[i,],Y[i,],ZZ[i,])+Az(t[i]+Dt,X[i,],Y[i,],ZZZ[i,]))+
             (1/6)*(Sz(t[i],X[i,],Y[i,],Z[i,])+4*Sz(t[i]+0.5*Dt,X[i,],Y[i,],ZZ[i,])+Sz(t[i]+Dt,X[i,],Y[i,],ZZZ[i,]))*W3[i,]
			  }
       }
    name <- c("X","Y","Z")
    name <- if(M > 1) c(paste(name[1],1:M,sep=""),paste(name[2],1:M,sep=""),paste(name[3],1:M,sep=""))
    X <- ts(X, start = t0, deltat = Dt, names=name[1:M])
    Y <- ts(Y, start = t0, deltat = Dt, names=name[(M+1):(2*M)])
    Z <- ts(Z, start = t0, deltat = Dt, names=name[(2*M+1):(3*M)])
    return(list(X=X,Y=Y,Z=Z))
} 


