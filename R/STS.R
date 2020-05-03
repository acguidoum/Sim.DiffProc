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
##### STS1D

.STS1D <- function(N =1000,M=1,x0=0,t0=0,T=1,Dt,drift,diffusion,
                          type=c("ito","str"),...)
                       {
	DSx  <- Deriv(diffusion,"x")  
    DSxx <- Deriv(DSx,"x")
    if (type=="ito"){
    A    <- function(t,x)  eval(drift)
    Ax   <- function(t,x)  eval(Deriv(drift,"x"))
    Axx  <- function(t,x)  eval(Deriv(Deriv(drift,"x"),"x"))
    }else{
    A    <- function(t,x)  eval(drift) + 0.5 * eval(diffusion) * eval(Deriv(diffusion,"x"))
    Ax   <- function(t,x)  eval(Deriv(drift,"x")) + 0.5 * (eval(Deriv(diffusion,"x")) * eval(Deriv(diffusion,"x"))+ eval(diffusion) * eval(Deriv(Deriv(diffusion,"x"),"x")))
    Axx  <- function(t,x)  eval(Deriv(Deriv(drift,"x"),"x")) + 0.5 * ( eval(Deriv(Deriv(diffusion,"x"),"x")) * eval(Deriv(diffusion,"x"))+ eval(Deriv(diffusion,"x")) * eval(Deriv(Deriv(diffusion,"x"),"x"))+
                           eval(Deriv(diffusion,"x")) * eval(Deriv(Deriv(diffusion,"x"),"x")) + eval(diffusion) * eval(Deriv(Deriv(Deriv(diffusion,"x"),"x"),"x")) )
                  }
    S    <- function(t,x)  eval(diffusion)
    Sx   <- function(t,x)  eval(DSx)
    Sxx  <- function(t,x)  eval(DSxx)
    if (missing(Dt)) {
        t <- seq(t0, T, by=Dt)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
		T <- t[N + 1]
    }
    Dt <- (T - t0)/N
	Sigma <- matrix(c(Dt, 0.5 * Dt^2, 0.5 * Dt^2, (1/3)*Dt^3), ncol=2, nrow=2)
	B <- do.call("cbind",lapply(1:M,function(i) MASS::mvrnorm(N, mu= c(0, 0), Sigma)))
	W <- B[,-seq(2,2*M,by=2)]
	Z <- B[,-seq(1,2*M-1,by=2)]
	X <- matrix(x0, N+1, M)
    for (i in 1L:N) {
	     X[i + 1L,] = X[i,]+A(t[i],X[i,])*Dt+S(t[i],X[i,])*W[i,]+ 0.5*S(t[i],X[i,])*Sx(t[i],X[i,])*((W[i,]^2)-Dt)+
            Ax(t[i],X[i,])*S(t[i],X[i,])*Z[i,]+0.5*(A(t[i],X[i,])*Ax(t[i],X[i,])+0.5*(S(t[i],X[i,])^2)*
            Axx(t[i],X[i,]))*(Dt^2)+(A(t[i],X[i,])*Sx(t[i],X[i,])+0.5*(S(t[i],X[i,])^2)*Sxx(t[i],X[i,]))*
            (W[i,]*Dt-Z[i,])+0.5*S(t[i],X[i,])*(S(t[i],X[i,])*Sxx(t[i],X[i,])+(Sx(t[i],X[i,])^2))*((1/3)*(W[i,]^2)-Dt)*W[i,]
	}
    name <- "X"
    name <- if(M > 1) paste("X",1:M,sep="")
    X <- ts(X, start = t0, deltat = Dt, names=name)
    return(list(X=X))
}   
      
             
#####
##### STS2D

.STS2D <- function(N =1000,M=1,x0=0,y0=0,t0=0,T=1,Dt,driftx,diffx,drifty,diffy,
                          type=c("ito","str"),...)
   {
    DSx  <- Deriv(diffx,"x")  
    DSxx <- Deriv(DSx,"x")
    DSy  <- Deriv(diffy,"y")  
    DSyy <- Deriv(DSy,"y")					   
    if (type=="ito"){
    Ax    <- function(t,x,y)  eval(driftx)
    dAx   <- function(t,x,y)  eval(Deriv(driftx,"x"))
    dAxx  <- function(t,x,y)  eval(Deriv(Deriv(driftx,"x"),"x"))
    Ay    <- function(t,x,y)  eval(drifty)
    dAy   <- function(t,x,y)  eval(Deriv(drifty,"y"))
    dAyy  <- function(t,x,y)  eval(Deriv(Deriv(drifty,"y"),"y"))
    }else{
    Ax    <- function(t,x,y)  eval(driftx) + 0.5 * eval(diffx) * eval(Deriv(diffx,"x"))
    dAx   <- function(t,x,y)  eval(Deriv(driftx,"x")) + 0.5 * (eval(Deriv(diffx,"x")) * eval(Deriv(diffx,"x"))+ eval(diffx) * eval(Deriv(Deriv(diffx,"x"),"x")))
    dAxx  <- function(t,x,y)  eval(Deriv(Deriv(driftx,"x"),"x")) + 0.5 * ( eval(Deriv(Deriv(diffx,"x"),"x")) * eval(Deriv(diffx,"x"))+ eval(Deriv(diffx,"x")) * eval(Deriv(Deriv(diffx,"x"),"x"))+
                              eval(Deriv(diffx,"x")) * eval(Deriv(Deriv(diffx,"x"),"x")) + eval(diffx) * eval(Deriv(Deriv(Deriv(diffx,"x"),"x"),"x")) )
    Ay    <- function(t,x,y)  eval(drifty) + 0.5 * eval(diffy) * eval(Deriv(diffy,"y"))
    dAy   <- function(t,x,y)  eval(Deriv(drifty,"y")) + 0.5 * (eval(Deriv(diffy,"y")) * eval(Deriv(diffy,"y"))+ eval(diffy) * eval(Deriv(Deriv(diffy,"y"),"y")))
    dAyy  <- function(t,x,y)  eval(Deriv(Deriv(drifty,"y"),"y")) + 0.5 * ( eval(Deriv(Deriv(diffy,"y"),"y")) * eval(Deriv(diffy,"y"))+ eval(Deriv(diffy,"y")) * eval(Deriv(Deriv(diffy,"y"),"y"))+
                              eval(Deriv(diffy,"y")) * eval(Deriv(Deriv(diffy,"y"),"y")) + eval(diffy) * eval(Deriv(Deriv(Deriv(diffy,"y"),"y"),"y")) )
                  }
    Sx    <- function(t,x,y)  eval(diffx)
    dSx   <- function(t,x,y)  eval(DSx)
    dSxx  <- function(t,x,y)  eval(DSxx)
    Sy    <- function(t,x,y)  eval(diffy)
    dSy   <- function(t,x,y)  eval(DSy)
    dSyy  <- function(t,x,y)  eval(DSyy)
    if (missing(Dt)) {
        t <- seq(t0, T, by=Dt)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
		T <- t[N + 1]
    }
    Dt <- (T - t0)/N
	Sigma <- matrix(c(Dt, 0.5 * Dt^2, 0.5 * Dt^2, (1/3)*Dt^3), ncol=2, nrow=2)
	Bx <- do.call("cbind",lapply(1:M,function(i) MASS::mvrnorm(N, mu= c(0, 0), Sigma)))
	By <- do.call("cbind",lapply(1:M,function(i) MASS::mvrnorm(N, mu= c(0, 0), Sigma)))
	Wx <- Bx[,-seq(2,2*M,by=2)]
	Wy <- By[,-seq(2,2*M,by=2)]
	Zx <- Bx[,-seq(1,2*M-1,by=2)]
	Zy <- By[,-seq(1,2*M-1,by=2)]
	X <- matrix(x0, N+1, M)
	Y <- matrix(y0, N+1, M)
    for (i in 1L:N) {
	     X[i + 1L,] = X[i,]+Ax(t[i],X[i,],Y[i,])*Dt+Sx(t[i],X[i,],Y[i,])*Wx[i,]+ 0.5*Sx(t[i],X[i,],Y[i,])*dSx(t[i],X[i,],Y[i,])*((Wx[i,]^2)-Dt)+
            dAx(t[i],X[i,],Y[i,])*Sx(t[i],X[i,],Y[i,])*Zx[i,]+0.5*(Ax(t[i],X[i,],Y[i,])*dAx(t[i],X[i,],Y[i,])+0.5*(Sx(t[i],X[i,],Y[i,])^2)*
            dAxx(t[i],X[i,],Y[i,]))*(Dt^2)+(Ax(t[i],X[i,],Y[i,])*dSx(t[i],X[i,],Y[i,])+0.5*(Sx(t[i],X[i,],Y[i,])^2)*dSxx(t[i],X[i,],Y[i,]))*
            (Wx[i,]*Dt-Zx[i,])+0.5*Sx(t[i],X[i,],Y[i,])*(Sx(t[i],X[i,],Y[i,])*dSxx(t[i],X[i,],Y[i,])+(dSx(t[i],X[i,],Y[i,])^2))*((1/3)*(Wx[i,]^2)-Dt)*Wx[i,]
		 Y[i + 1L,] = Y[i,]+Ay(t[i],X[i,],Y[i,])*Dt+Sy(t[i],X[i,],Y[i,])*Wy[i,]+ 0.5*Sy(t[i],X[i,],Y[i,])*dSy(t[i],X[i,],Y[i,])*((Wy[i,]^2)-Dt)+
            dAy(t[i],X[i,],Y[i,])*Sy(t[i],X[i,],Y[i,])*Zy[i,]+0.5*(Ay(t[i],X[i,],Y[i,])*dAy(t[i],X[i,],Y[i,])+0.5*(Sy(t[i],X[i,],Y[i,])^2)*
            dAyy(t[i],X[i,],Y[i,]))*(Dt^2)+(Ay(t[i],X[i,],Y[i,])*dSy(t[i],X[i,],Y[i,])+0.5*(Sy(t[i],X[i,],Y[i,])^2)*dSyy(t[i],X[i,],Y[i,]))*
            (Wy[i,]*Dt-Zy[i,])+0.5*Sy(t[i],X[i,],Y[i,])*(Sy(t[i],X[i,],Y[i,])*dSyy(t[i],X[i,],Y[i,])+(dSy(t[i],X[i,],Y[i,])^2))*((1/3)*(Wy[i,]^2)-Dt)*Wy[i,]
		 }
    name <- c("X","Y")
    name <- if(M > 1) c(paste(name[1],1:M,sep=""),paste(name[2],1:M,sep=""))
    X <- ts(X, start = t0, deltat = Dt, names=name[1:M])
    Y <- ts(Y, start = t0, deltat = Dt, names=name[(M+1):(2*M)])
    return(list(X=X,Y=Y))
}

#####
##### STS3D

.STS3D <- function(N =100,M=1,x0=0,y0=0,z0=0,t0=0,T=1,Dt,driftx,diffx,drifty,diffy,
                     driftz,diffz,type=c("ito","str"),...)
                       {
    if (type=="ito"){
    Ax    <- function(t,x,y,z)  eval(driftx)
    dAx   <- function(t,x,y,z)  eval(Deriv(driftx,"x"))
    dAxx  <- function(t,x,y,z)  eval(Deriv(Deriv(driftx,"x"),"x"))
    Ay    <- function(t,x,y,z)  eval(drifty)
    dAy   <- function(t,x,y,z)  eval(Deriv(drifty,"y"))
    dAyy  <- function(t,x,y,z)  eval(Deriv(Deriv(drifty,"y"),"y"))
    Az    <- function(t,x,y,z)  eval(driftz)
    dAz   <- function(t,x,y,z)  eval(Deriv(driftz,"z"))
    dAzz  <- function(t,x,y,z)  eval(Deriv(Deriv(driftz,"z"),"z"))}else{
    Ax    <- function(t,x,y,z)  eval(driftx) + 0.5 * eval(diffx) * eval(Deriv(diffx,"x"))
    dAx   <- function(t,x,y,z)  eval(Deriv(driftx,"x")) + 0.5 * (eval(Deriv(diffx,"x")) * eval(Deriv(diffx,"x"))+ eval(diffx) * eval(Deriv(Deriv(diffx,"x"),"x")))
    dAxx  <- function(t,x,y,z)  eval(Deriv(Deriv(driftx,"x"),"x")) + 0.5 * ( eval(Deriv(Deriv(diffx,"x"),"x")) * eval(Deriv(diffx,"x"))+ eval(Deriv(diffx,"x")) * eval(Deriv(Deriv(diffx,"x"),"x"))+
                                eval(Deriv(diffx,"x")) * eval(Deriv(Deriv(diffx,"x"),"x")) + eval(diffx) * eval(Deriv(Deriv(Deriv(diffx,"x"),"x"),"x")) )
    Ay    <- function(t,x,y,z)  eval(drifty) + 0.5 * eval(diffy) * eval(Deriv(diffy,"y"))
    dAy   <- function(t,x,y,z)  eval(Deriv(drifty,"y")) + 0.5 * (eval(Deriv(diffy,"y")) * eval(Deriv(diffy,"y"))+ eval(diffy) * eval(Deriv(Deriv(diffy,"y"),"y")))
    dAyy  <- function(t,x,y,z)  eval(Deriv(Deriv(drifty,"y"),"y")) + 0.5 * ( eval(Deriv(Deriv(diffy,"y"),"y")) * eval(Deriv(diffy,"y"))+ eval(Deriv(diffy,"y")) * eval(Deriv(Deriv(diffy,"y"),"y"))+
                                eval(Deriv(diffy,"y")) * eval(Deriv(Deriv(diffy,"y"),"y")) + eval(diffy) * eval(Deriv(Deriv(Deriv(diffy,"y"),"y"),"y")) )
    Az    <- function(t,x,y,z)  eval(driftz) + 0.5 * eval(diffz) * eval(Deriv(diffz,"z"))
    dAz   <- function(t,x,y,z)  eval(Deriv(driftz,"z")) + 0.5 * (eval(Deriv(diffz,"z")) * eval(Deriv(diffz,"z"))+ eval(diffz) * eval(Deriv(Deriv(diffz,"z"),"z")))
    dAzz  <- function(t,x,y,z)  eval(Deriv(Deriv(driftz,"z"),"z")) + 0.5 * ( eval(Deriv(Deriv(diffz,"z"),"z")) * eval(Deriv(diffz,"z"))+ eval(Deriv(diffz,"z")) * eval(Deriv(Deriv(diffz,"z"),"z"))+
                                eval(Deriv(diffz,"z")) * eval(Deriv(Deriv(diffz,"z"),"z")) + eval(diffz) * eval(Deriv(Deriv(Deriv(diffz,"z"),"z"),"z")) )
                  }
    DSx  <- Deriv(diffx,"x")  
    DSxx <- Deriv(DSx,"x")
    Sx    <- function(t,x,y,z)  eval(diffx)
    dSx   <- function(t,x,y,z)  eval(DSx)
    dSxx  <- function(t,x,y,z)  eval(DSxx)
    DSy  <- Deriv(diffy,"y")  
    DSyy <- Deriv(DSy,"y")
    Sy    <- function(t,x,y,z)  eval(diffy)
    dSy   <- function(t,x,y,z)  eval(DSy)
    dSyy  <- function(t,x,y,z)  eval(DSyy)
    DSz  <- Deriv(diffz,"z")  
    DSzz <- Deriv(DSz,"z")
    Sz    <- function(t,x,y,z)  eval(diffz)
    dSz   <- function(t,x,y,z)  eval(DSz)
    dSzz  <- function(t,x,y,z)  eval(DSzz)
    if (missing(Dt)) {
        t <- seq(t0, T, by=Dt)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
		T <- t[N + 1]
    } 
    Dt <- (T - t0)/N
	Sigma <- matrix(c(Dt, 0.5 * Dt^2, 0.5 * Dt^2, (1/3)*Dt^3), ncol=2, nrow=2)
	Bx <- do.call("cbind",lapply(1:M,function(i) MASS::mvrnorm(N, mu= c(0, 0), Sigma)))
	By <- do.call("cbind",lapply(1:M,function(i) MASS::mvrnorm(N, mu= c(0, 0), Sigma)))
	Bz <- do.call("cbind",lapply(1:M,function(i) MASS::mvrnorm(N, mu= c(0, 0), Sigma)))
	Wx <- Bx[,-seq(2,2*M,by=2)]
	Wy <- By[,-seq(2,2*M,by=2)]
	Wz <- Bz[,-seq(2,2*M,by=2)]
	Zx <- Bx[,-seq(1,2*M-1,by=2)]
	Zy <- By[,-seq(1,2*M-1,by=2)]
	Zz <- Bz[,-seq(1,2*M-1,by=2)]
	X <- matrix(x0, N+1, M)
	Y <- matrix(y0, N+1, M)
	Z <- matrix(z0, N+1, M)
    for (i in 1L:N) {
	     X[i + 1L,] = X[i,]+Ax(t[i],X[i,],Y[i,],Z[i,])*Dt+Sx(t[i],X[i,],Y[i,],Z[i,])*Wx[i,]+ 0.5*Sx(t[i],X[i,],Y[i,],Z[i,])*dSx(t[i],X[i,],Y[i,],Z[i,])*((Wx[i,]^2)-Dt)+
            dAx(t[i],X[i,],Y[i,],Z[i,])*Sx(t[i],X[i,],Y[i,],Z[i,])*Zx[i,]+0.5*(Ax(t[i],X[i,],Y[i,],Z[i,])*dAx(t[i],X[i,],Y[i,],Z[i,])+0.5*(Sx(t[i],X[i,],Y[i,],Z[i,])^2)*
            dAxx(t[i],X[i,],Y[i,],Z[i,]))*(Dt^2)+(Ax(t[i],X[i,],Y[i,],Z[i,])*dSx(t[i],X[i,],Y[i,],Z[i,])+0.5*(Sx(t[i],X[i,],Y[i,],Z[i,])^2)*dSxx(t[i],X[i,],Y[i,],Z[i,]))*
            (Wx[i,]*Dt-Zx[i,])+0.5*Sx(t[i],X[i,],Y[i,],Z[i,])*(Sx(t[i],X[i,],Y[i,],Z[i,])*dSxx(t[i],X[i,],Y[i,],Z[i,])+(dSx(t[i],X[i,],Y[i,],Z[i,])^2))*((1/3)*(Wx[i,]^2)-Dt)*Wx[i,]
		 Y[i + 1L,] = Y[i,]+Ay(t[i],X[i,],Y[i,],Z[i,])*Dt+Sy(t[i],X[i,],Y[i,],Z[i,])*Wy[i,]+ 0.5*Sy(t[i],X[i,],Y[i,],Z[i,])*dSy(t[i],X[i,],Y[i,],Z[i,])*((Wy[i,]^2)-Dt)+
            dAy(t[i],X[i,],Y[i,],Z[i,])*Sy(t[i],X[i,],Y[i,],Z[i,])*Zy[i,]+0.5*(Ay(t[i],X[i,],Y[i,],Z[i,])*dAy(t[i],X[i,],Y[i,],Z[i,])+0.5*(Sy(t[i],X[i,],Y[i,],Z[i,])^2)*
            dAyy(t[i],X[i,],Y[i,],Z[i,]))*(Dt^2)+(Ay(t[i],X[i,],Y[i,],Z[i,])*dSy(t[i],X[i,],Y[i,],Z[i,])+0.5*(Sy(t[i],X[i,],Y[i,],Z[i,])^2)*dSyy(t[i],X[i,],Y[i,],Z[i,]))*
            (Wy[i,]*Dt-Zy[i,])+0.5*Sy(t[i],X[i,],Y[i,],Z[i,])*(Sy(t[i],X[i,],Y[i,],Z[i,])*dSyy(t[i],X[i,],Y[i,],Z[i,])+(dSy(t[i],X[i,],Y[i,],Z[i,])^2))*((1/3)*(Wy[i,]^2)-Dt)*Wy[i,]
		 Z[i + 1L,] = Z[i,]+Az(t[i],X[i,],Y[i,],Z[i,])*Dt+Sz(t[i],X[i,],Y[i,],Z[i,])*Wz[i,]+ 0.5*Sz(t[i],X[i,],Y[i,],Z[i,])*dSz(t[i],X[i,],Y[i,],Z[i,])*((Wz[i,]^2)-Dt)+
            dAz(t[i],X[i,],Y[i,],Z[i,])*Sz(t[i],X[i,],Y[i,],Z[i,])*Zz[i,]+0.5*(Az(t[i],X[i,],Y[i,],Z[i,])*dAz(t[i],X[i,],Y[i,],Z[i,])+0.5*(Sz(t[i],X[i,],Y[i,],Z[i,])^2)*
            dAzz(t[i],X[i,],Y[i,],Z[i,]))*(Dt^2)+(Az(t[i],X[i,],Y[i,],Z[i,])*dSz(t[i],X[i,],Y[i,],Z[i,])+0.5*(Sz(t[i],X[i,],Y[i,],Z[i,])^2)*dSzz(t[i],X[i,],Y[i,],Z[i,]))*
            (Wz[i,]*Dt-Zz[i,])+0.5*Sz(t[i],X[i,],Y[i,],Z[i,])*(Sz(t[i],X[i,],Y[i,],Z[i,])*dSzz(t[i],X[i,],Y[i,],Z[i,])+(dSz(t[i],X[i,],Y[i,],Z[i,])^2))*((1/3)*(Wz[i,]^2)-Dt)*Wz[i,]
		 }
    name <- c("X","Y","Z")
    name <- if(M > 1) c(paste(name[1],1:M,sep=""),paste(name[2],1:M,sep=""),paste(name[3],1:M,sep=""))
    X <- ts(X, start = t0, deltat = Dt, names=name[1:M])
    Y <- ts(Y, start = t0, deltat = Dt, names=name[(M+1):(2*M)])
    Z <- ts(Z, start = t0, deltat = Dt, names=name[(2*M+1):(3*M)])
    return(list(X=X,Y=Y,Z=Z))
} 



